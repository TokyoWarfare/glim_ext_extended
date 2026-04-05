#include <deque>
#include <atomic>
#include <fstream>
#include <iomanip>
#include <thread>
#include <numeric>
#include <filesystem>
#include <Eigen/Core>

#define GLIM_ROS2

#include <boost/format.hpp>
#include <glim/mapping/callbacks.hpp>
#include <glim/util/logging.hpp>
#include <glim/util/concurrent_vector.hpp>

#ifdef GLIM_ROS2
#include <glim/util/extension_module_ros2.hpp>
#include <sensor_msgs/msg/nav_sat_fix.hpp>

using ExtensionModuleBase = glim::ExtensionModuleROS2;
using NavSatFix = sensor_msgs::msg::NavSatFix;
using NavSatFixConstPtr = sensor_msgs::msg::NavSatFix::ConstSharedPtr;

template <typename Stamp>
double to_sec(const Stamp& stamp) {
  return stamp.sec + stamp.nanosec / 1e9;
}
#else
#include <glim/util/extension_module_ros.hpp>
// ROS1 NavSatFix support not implemented.

using ExtensionModuleBase = glim::ExtensionModuleROS;
#endif

#include <spdlog/spdlog.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/slam/PoseTranslationPrior.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>

#include <glim/util/logging.hpp>
#include <glim/util/convert_to_string.hpp>
#include <glim_ext/util/config_ext.hpp>
#include <glim_ext/geodetic.hpp>

namespace glim {

using gtsam::symbol_shorthand::X;

/**
 * @brief GNSS global constraints from NavSatFix (lat/lon/alt → ENU).
 * @note  Subscribes to sensor_msgs/NavSatFix. The first accepted fix sets the ENU
 *        datum; all subsequent fixes are expressed as East/North/Up offsets (metres)
 *        from that origin. The ENU coordinates feed the existing SVD alignment that
 *        bootstraps T_world_enu. Fix quality is gated via min_fix_status; degraded
 *        fixes (SBAS) use an inflated noise model controlled by sbas_noise_inflation.
 */
class GNSSGlobal : public ExtensionModuleBase {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // NavSatFix status constants (sensor_msgs/NavSatStatus).
  static constexpr int STATUS_NO_FIX   = -1;  // No signal — always rejected.
  static constexpr int STATUS_FIX      =  0;  // Autonomous, no corrections.
  static constexpr int STATUS_SBAS_FIX =  1;  // SBAS/RTK-float — degraded accuracy.
  static constexpr int STATUS_GBAS_FIX =  2;  // RTK fix — nominal accuracy.

  GNSSGlobal() : logger(create_module_logger("gnss_global")) {
    logger->info("initializing GNSS global constraints (NavSatFix → ENU)");
    const std::string config_path = glim::GlobalConfigExt::get_config_path("config_gnss_global");
    logger->info("gnss_global_config_path={}", config_path);

    glim::Config config(config_path);
    gnss_topic         = config.param<std::string>("gnss", "gnss_topic", "/navsatfix");
    prior_inf_scale    = config.param<Eigen::Vector3d>("gnss", "prior_inf_scale", Eigen::Vector3d(1e3, 1e3, 1e3));
    min_baseline       = config.param<double>("gnss", "min_baseline", 5.0);
    // Minimum NavSatFix status to accept. STATUS_NO_FIX (-1) is always rejected
    // regardless of this setting. Default 1 = accept SBAS_FIX and above.
    min_fix_status     = config.param<int>("gnss", "min_fix_status", STATUS_SBAS_FIX);
    // Noise inflation factor for fixes below STATUS_GBAS_FIX (RTK float / SBAS).
    // The information scale is divided by this value squared, so larger values
    // mean softer constraints for degraded fixes. Default 10 → σ inflated ×10.
    sbas_noise_inflation = config.param<double>("gnss", "sbas_noise_inflation", 10.0);

    logger->info(
      "config: topic={} min_fix_status={} sbas_noise_inflation={:.1f} min_baseline={:.1f}",
      gnss_topic, min_fix_status, sbas_noise_inflation, min_baseline);

    transformation_initialized = false;
    T_world_enu.setIdentity();
    datum_initialized    = false;
    datum_json_written   = false;
    datum_utm_zone       = 0;
    datum_utm_easting    = 0.0;
    datum_utm_northing   = 0.0;

    kill_switch = false;
    thread = std::thread([this] { backend_task(); });

    using std::placeholders::_1;
    using std::placeholders::_2;
    using std::placeholders::_3;
    GlobalMappingCallbacks::on_insert_submap.add(std::bind(&GNSSGlobal::on_insert_submap, this, _1));
    GlobalMappingCallbacks::on_smoother_update.add(std::bind(&GNSSGlobal::on_smoother_update, this, _1, _2, _3));
  }

  ~GNSSGlobal() {
    kill_switch = true;
    thread.join();
  }

  virtual std::vector<GenericTopicSubscription::Ptr> create_subscriptions() override {
    const auto sub = std::make_shared<TopicSubscription<NavSatFix>>(
      gnss_topic, [this](const NavSatFixConstPtr msg) { navsatfix_callback(msg); });
    return {sub};
  }

  void navsatfix_callback(const NavSatFixConstPtr& msg) {
    const int status = msg->status.status;

    // STATUS_NO_FIX (-1) is hard-rejected unconditionally: no signal at all.
    if (status == STATUS_NO_FIX) {
      return;
    }
    // Reject fixes below the configured minimum quality threshold.
    if (status < min_fix_status) {
      logger->debug("GNSS fix rejected: status={} < min_fix_status={}", status, min_fix_status);
      return;
    }

    // Set ENU datum from the first accepted fix.
    if (!datum_initialized) {
      datum_lat = msg->latitude;
      datum_lon = msg->longitude;
      datum_alt = msg->altitude;
      datum_ecef = wgs84_to_ecef(datum_lat, datum_lon, datum_alt);
      R_enu_ecef = enu_rotation(datum_lat, datum_lon);

      // Precompute UTM origin at datum time.
      // UTM zone is locked to datum longitude — no zone-crossing handling.
      datum_utm_zone = ecef_to_utm_zone(datum_lat, datum_lon);
      const Eigen::Vector2d utm_xy = wgs84_to_utm_xy(datum_lat, datum_lon);
      datum_utm_easting  = utm_xy.x();
      datum_utm_northing = utm_xy.y();

      datum_initialized = true;
      logger->info(
        "ENU datum set: lat={:.9f} lon={:.9f} alt={:.3f}  UTM zone={} E={:.3f} N={:.3f}",
        datum_lat, datum_lon, datum_alt, datum_utm_zone, datum_utm_easting, datum_utm_northing);
    }

    // Convert to ENU relative to datum.
    const Eigen::Vector3d ecef = wgs84_to_ecef(msg->latitude, msg->longitude, msg->altitude);
    const Eigen::Vector3d enu  = R_enu_ecef * (ecef - datum_ecef);

    // Extract per-axis sigmas from position_covariance when available.
    // NavSatFix position_covariance is a row-major 3×3 in ENU order: [0]=σ_E², [4]=σ_N², [8]=σ_U².
    // Store -1 as a sentinel when the covariance type is unknown.
    constexpr int COVARIANCE_TYPE_UNKNOWN = 0;
    constexpr double SIGMA_FLOOR = 0.005;  // 0.5 cm
    constexpr double SIGMA_CAP   = 20.0;   // 20 m — reject measurement if exceeded

    Eigen::Vector3d sigma(-1.0, -1.0, -1.0);  // sentinel: fall back to prior_inf_scale
    if (msg->position_covariance_type != COVARIANCE_TYPE_UNKNOWN) {
      const double raw_E = std::sqrt(msg->position_covariance[0]);
      const double raw_N = std::sqrt(msg->position_covariance[4]);
      const double raw_U = std::sqrt(msg->position_covariance[8]);

      // Reject the entire measurement if any component is NaN or exceeds the cap.
      if (std::isnan(raw_E) || std::isnan(raw_N) || std::isnan(raw_U) ||
          raw_E > SIGMA_CAP || raw_N > SIGMA_CAP || raw_U > SIGMA_CAP) {
        logger->warn(
          "GNSS measurement rejected: sigma ({:.2f}, {:.2f}, {:.2f}) m exceeds {:.0f} m cap or is NaN",
          raw_E, raw_N, raw_U, SIGMA_CAP);
        return;
      }

      // Apply floor and warn if any component is clamped.
      auto floored = [&](double raw, const char* axis) -> double {
        if (raw < SIGMA_FLOOR) {
          logger->warn(
            "[WARN] GNSS covariance below 0.5cm floor — clamping {} sigma to 0.5cm (reported: {:.2f} cm)",
            axis, raw * 100.0);
          return SIGMA_FLOOR;
        }
        return raw;
      };
      sigma.x() = floored(raw_E, "E");
      sigma.y() = floored(raw_N, "N");
      sigma.z() = floored(raw_U, "U");
    }

    // Queue layout: [stamp, enu_e, enu_n, enu_u, status, σ_E, σ_N, σ_U]
    // σ < 0 = covariance unknown; factor builder falls back to prior_inf_scale.
    Eigen::Matrix<double, 8, 1> entry;
    entry << to_sec(msg->header.stamp), enu.x(), enu.y(), enu.z(),
             static_cast<double>(status), sigma.x(), sigma.y(), sigma.z();
    input_gnss_queue.push_back(entry);
  }

  void on_insert_submap(const SubMap::ConstPtr& submap) {
    input_submap_queue.push_back(submap);
  }

  void on_smoother_update(
    gtsam_points::ISAM2Ext& isam2,
    gtsam::NonlinearFactorGraph& new_factors,
    gtsam::Values& new_values)
  {
    const auto factors = output_factors.get_all_and_clear();
    if (!factors.empty()) {
      logger->debug("insert {} GNSS prior factors", factors.size());
      new_factors.add(factors);
    }
  }

  void backend_task() {
    logger->info("starting GNSS global thread");

    // enu_queue stores [stamp, enu_e, enu_n, enu_u, status, σ_E, σ_N, σ_U].
    std::deque<Eigen::Matrix<double, 8, 1>> ecef_queue;
    std::deque<SubMap::ConstPtr> submap_queue;

    while (!kill_switch) {
      // Drain incoming GNSS measurements into the local queue.
      const auto gnss_data = input_gnss_queue.get_all_and_clear();
      ecef_queue.insert(ecef_queue.end(), gnss_data.begin(), gnss_data.end());

      // Drain incoming submaps.
      const auto new_submaps = input_submap_queue.get_all_and_clear();
      if (new_submaps.empty()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        continue;
      }
      submap_queue.insert(submap_queue.end(), new_submaps.begin(), new_submaps.end());

      // Discard submaps that predate the oldest GNSS measurement — they cannot
      // be matched by interpolation and would produce stale constraints.
      while (!ecef_queue.empty() && !submap_queue.empty() &&
             submap_queue.front()->frames.front()->stamp < ecef_queue.front()[0]) {
        submap_queue.pop_front();
      }

      // Temporally interpolate ECEF coordinates to each submap's mid-frame stamp.
      while (
        !ecef_queue.empty() && !submap_queue.empty() &&
        submap_queue.front()->frames.front()->stamp > ecef_queue.front()[0] &&
        submap_queue.front()->frames.back()->stamp  < ecef_queue.back()[0])
      {
        const auto& submap = submap_queue.front();
        const double stamp = submap->frames[submap->frames.size() / 2]->stamp;

        const auto right = std::lower_bound(
          ecef_queue.begin(), ecef_queue.end(), stamp,
          [](const Eigen::Matrix<double, 8, 1>& e, double t) { return e[0] < t; });

        if (right == ecef_queue.end() || right == ecef_queue.begin()) {
          logger->warn("GNSS interpolation bounds invalid — skipping submap");
          break;
        }
        const auto left = right - 1;

        logger->debug("submap={:.6f} gnss_left={:.6f} gnss_right={:.6f}", stamp, (*left)[0], (*right)[0]);

        const double tl = (*left)[0];
        const double tr = (*right)[0];
        const double p  = (stamp - tl) / (tr - tl);

        // Linearly interpolate position and sigma; inherit the *lower* of the two statuses
        // (conservative: a mixed bracket gets the worse quality label).
        // For sigma: if either bracket has unknown covariance (sentinel < 0), propagate
        // the sentinel so the factor builder falls back to prior_inf_scale.
        Eigen::Matrix<double, 8, 1> interpolated = (1.0 - p) * (*left) + p * (*right);
        interpolated[4] = std::min((*left)[4], (*right)[4]);
        if ((*left)[5] < 0.0 || (*right)[5] < 0.0) {
          interpolated[5] = -1.0;
          interpolated[6] = -1.0;
          interpolated[7] = -1.0;
        }

        submaps.push_back(submap);
        submap_coords.push_back(interpolated);

        submap_queue.pop_front();
        ecef_queue.erase(ecef_queue.begin(), left);
      }

      // Bootstrap T_world_enu via SVD once we have enough baseline.
      if (!transformation_initialized && !submaps.empty() &&
          (submaps.front()->T_world_origin.inverse() *
           submaps.back()->T_world_origin).translation().norm() > min_baseline)
      {
        Eigen::Vector3d mean_est = Eigen::Vector3d::Zero();
        Eigen::Vector3d mean_enu = Eigen::Vector3d::Zero();
        for (size_t i = 0; i < submaps.size(); i++) {
          mean_est += submaps[i]->T_world_origin.translation();
          mean_enu += submap_coords[i].segment<3>(1);
        }
        mean_est /= submaps.size();
        mean_enu /= submaps.size();

        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        for (size_t i = 0; i < submaps.size(); i++) {
          const Eigen::Vector3d ce = submaps[i]->T_world_origin.translation() - mean_est;
          const Eigen::Vector3d cg = submap_coords[i].segment<3>(1) - mean_enu;
          cov += cg * ce.transpose();
        }
        cov /= submaps.size();

        // 2D SVD alignment (horizontal plane only; Z handled via translation).
        const Eigen::JacobiSVD<Eigen::Matrix2d> svd(
          cov.block<2, 2>(0, 0), Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Eigen::Matrix2d U = svd.matrixU();
        const Eigen::Matrix2d V = svd.matrixV();
        Eigen::Matrix2d S = Eigen::Matrix2d::Identity();
        if (U.determinant() * V.determinant() < 0.0) {
          S(1, 1) = -1.0;
        }

        Eigen::Isometry3d T_enu_world = Eigen::Isometry3d::Identity();
        T_enu_world.linear().block<2, 2>(0, 0) = U * S * V.transpose();
        T_enu_world.translation() = mean_enu - T_enu_world.linear() * mean_est;

        T_world_enu = T_enu_world.inverse();

        for (size_t i = 0; i < submaps.size(); i++) {
          const Eigen::Vector3d w = T_world_enu * submap_coords[i].segment<3>(1);
          logger->debug("submap={} enu_in_world={}",
            convert_to_string(submaps[i]->T_world_origin.translation().eval()),
            convert_to_string(w));
        }

        logger->info("T_world_enu initialized: {}", convert_to_string(T_world_enu));

        // DEBUG: print full T_enu_world rotation matrix and derived yaw relative to North.
        // T_enu_world maps world-frame vectors into ENU (East=+X, North=+Y, Up=+Z).
        // Column 0 of R is where world +X lands in ENU; column 1 is where world +Y lands.
        // Near-identity here means SVD found no meaningful heading correction.
        {
          const Eigen::Matrix3d R = T_enu_world.linear();
          const Eigen::Vector3d t = T_enu_world.translation();
          logger->info("[DEBUG] T_enu_world rotation matrix (row-major):");
          logger->info("[DEBUG]   R row0: [{:.6f}, {:.6f}, {:.6f}]", R(0,0), R(0,1), R(0,2));
          logger->info("[DEBUG]   R row1: [{:.6f}, {:.6f}, {:.6f}]", R(1,0), R(1,1), R(1,2));
          logger->info("[DEBUG]   R row2: [{:.6f}, {:.6f}, {:.6f}]", R(2,0), R(2,1), R(2,2));
          logger->info("[DEBUG]   translation: [{:.4f}, {:.4f}, {:.4f}]", t(0), t(1), t(2));
        }

        // yaw: angle world +X axis makes with geographic East (ENU +X).
        // 0 deg = world +X points East; 90 deg = world +X points North.
        const double yaw_deg = std::atan2(T_enu_world.linear()(1,0), T_enu_world.linear()(0,0)) * 180.0 / M_PI;
        logger->info("[DEBUG] yaw of world +X relative to geographic East: {:.2f} deg "
                     "(0=East, 90=North, ±180=West)", yaw_deg);

        // Write gnss_datum.json now that T_world_enu carries the real SVD-aligned heading.
        // Previously this was written on the first GPS fix when T_world_enu was still
        // identity, so the offline viewer never got the heading correction.
        if (!datum_json_written) {
          write_datum_json("/tmp/dump/config/gnss_datum.json");
          datum_json_written = true;
          logger->info("gnss_datum.json written after SVD alignment — yaw={:.2f} deg relative to geographic East",
                       yaw_deg);
        }

        transformation_initialized = true;
      }

      // Emit a PoseTranslationPrior factor for the latest matched submap.
      if (transformation_initialized) {
        const auto& coord = submap_coords.back();
        const Eigen::Vector3d xyz = T_world_enu * coord.segment<3>(1);
        const int fix_status = static_cast<int>(coord[4]);

        logger->debug("submap={} gnss_world={} fix_status={}",
          convert_to_string(submaps.back()->T_world_origin.translation().eval()),
          convert_to_string(xyz), fix_status);

        // Build precision vector from dynamic sigma if position_covariance was available,
        // otherwise fall back to the prior_inf_scale config values.
        // coord[5..7] = σ_E, σ_N, σ_U (metres); negative = sentinel (type unknown).
        Eigen::Vector3d inf_scale;
        if (coord[5] >= 0.0) {
          // Dynamic covariance path: precision = 1/σ².
          inf_scale = Eigen::Vector3d(
            1.0 / (coord[5] * coord[5]),
            1.0 / (coord[6] * coord[6]),
            1.0 / (coord[7] * coord[7]));
          logger->debug("dynamic sigma: σ_E={:.4f} σ_N={:.4f} σ_U={:.4f} m",
                        coord[5], coord[6], coord[7]);
        } else {
          // Fallback path: covariance type unknown, use configured prior_inf_scale.
          inf_scale = prior_inf_scale;
        }

        // Apply SBAS noise inflation on top of whichever source was used.
        // Information ∝ 1/σ²; inflating σ by k → divide information by k².
        if (fix_status < STATUS_GBAS_FIX) {
          const double k2 = sbas_noise_inflation * sbas_noise_inflation;
          inf_scale /= k2;
          logger->debug("degraded fix (status={}): noise inflated ×{:.1f}", fix_status, sbas_noise_inflation);
        }

        const auto model = gtsam::noiseModel::Diagonal::Precisions(inf_scale);
        gtsam::NonlinearFactor::shared_ptr factor(
          new gtsam::PoseTranslationPrior<gtsam::Pose3>(X(submaps.back()->id), xyz, model));
        output_factors.push_back(factor);
      }
    }
  }

  // Write gnss_datum.json to path. T_world_enu may be identity if SVD alignment has
  // not yet fired — the offline viewer will still get a valid UTM origin.
  void write_datum_json(const std::string& path) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream ofs(path);
    if (!ofs) {
      logger->error("failed to write gnss_datum.json to {}: check path exists", path);
      return;
    }

    const Eigen::Isometry3d T_enu_world = T_world_enu.inverse();
    const Eigen::Matrix3d R = T_enu_world.linear();
    const Eigen::Vector3d t = T_enu_world.translation();

    ofs << std::setprecision(15) << std::fixed;
    ofs << "{\n";
    ofs << "  \"latitude\": "            << datum_lat          << ",\n";
    ofs << "  \"longitude\": "           << datum_lon          << ",\n";
    ofs << "  \"altitude\": "            << datum_alt          << ",\n";
    ofs << "  \"utm_zone\": "            << datum_utm_zone     << ",\n";
    ofs << "  \"utm_easting_origin\": "  << datum_utm_easting  << ",\n";
    ofs << "  \"utm_northing_origin\": " << datum_utm_northing << ",\n";
    ofs << "  \"T_enu_world\": [\n";
    for (int r = 0; r < 3; r++) {
      ofs << "    " << R(r,0) << ", " << R(r,1) << ", " << R(r,2) << ", " << t(r);
      ofs << (r < 2 ? "," : "") << "\n";
    }
    ofs << "  ]\n";
    ofs << "}\n";

    logger->info("gnss_datum.json written to {}", path);
  }

  // Build the 3×3 rotation matrix R such that:  enu = R * (ecef_point - ecef_origin)
  // Rows are the East, North, Up unit vectors expressed in ECEF.
  static Eigen::Matrix3d enu_rotation(double lat_deg, double lon_deg) {
    const double lat = lat_deg * M_PI / 180.0;
    const double lon = lon_deg * M_PI / 180.0;
    const double slat = std::sin(lat), clat = std::cos(lat);
    const double slon = std::sin(lon), clon = std::cos(lon);
    Eigen::Matrix3d R;
    R.row(0) = Eigen::Vector3d(-slon,         clon,        0.0);   // East
    R.row(1) = Eigen::Vector3d(-slat * clon, -slat * slon, clat);  // North
    R.row(2) = Eigen::Vector3d( clat * clon,  clat * slon, slat);  // Up
    return R;
  }

private:
  std::atomic_bool kill_switch;
  std::thread thread;

  ConcurrentVector<Eigen::Matrix<double, 8, 1>> input_gnss_queue;
  ConcurrentVector<SubMap::ConstPtr> input_submap_queue;
  ConcurrentVector<gtsam::NonlinearFactor::shared_ptr> output_factors;

  std::vector<SubMap::ConstPtr> submaps;
  std::vector<Eigen::Matrix<double, 8, 1>> submap_coords;  // [stamp, enu_e, enu_n, enu_u, status, σ_E, σ_N, σ_U]
                                                            // σ < 0 = sentinel: covariance unknown, use prior_inf_scale

  std::string gnss_topic;
  Eigen::Vector3d prior_inf_scale;
  double min_baseline;
  int min_fix_status;
  double sbas_noise_inflation;

  // ENU datum — set once from the first accepted NavSatFix.
  bool datum_initialized;
  bool datum_json_written;  // true once gnss_datum.json is written post-SVD alignment; guarded by backend_task thread only
  double datum_lat, datum_lon, datum_alt;
  Eigen::Vector3d datum_ecef;
  Eigen::Matrix3d R_enu_ecef;
  // UTM origin precomputed at datum time (zone locked to datum longitude).
  int    datum_utm_zone;
  double datum_utm_easting;
  double datum_utm_northing;

  bool transformation_initialized;
  Eigen::Isometry3d T_world_enu;

  std::shared_ptr<spdlog::logger> logger;
};

}  // namespace glim

extern "C" glim::ExtensionModule* create_extension_module() {
  return new glim::GNSSGlobal();
}
