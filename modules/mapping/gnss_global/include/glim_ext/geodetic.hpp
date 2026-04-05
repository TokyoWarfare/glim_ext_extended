#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace glim {

Eigen::Vector3d ecef_to_wgs84(const Eigen::Vector3d& xyz);

Eigen::Vector3d wgs84_to_ecef(double lat, double lon, double alt);

Eigen::Isometry3d calc_T_ecef_nwz(const Eigen::Vector3d& ecef, double radius = 6378137);

double harversine(const Eigen::Vector2d& latlon1, const Eigen::Vector2d& latlon2);

/// Returns the UTM zone number (1–60) for the given longitude in degrees.
int ecef_to_utm_zone(double lat, double lon);

/// Projects (lat, lon) in degrees to UTM (easting, northing) in metres using
/// the WGS84 Transverse Mercator formulas (Snyder 1987).  The zone is the one
/// that contains the given longitude; use ecef_to_utm_zone() to obtain it.
/// False easting (500 000 m) and false northing (10 000 000 m for S hemisphere)
/// are applied so the result is a standard UTM grid coordinate.
Eigen::Vector2d wgs84_to_utm_xy(double lat, double lon);

/// Converts an ENU point (metres from datum) to UTM easting / northing / altitude.
/// The UTM zone is locked to the one that contains datum_lon.
/// utm.x = UTM easting, utm.y = UTM northing, utm.z = datum_alt + enu.z().
Eigen::Vector3d enu_to_utm(const Eigen::Vector3d& enu, double datum_lat, double datum_lon, double datum_alt);

}  // namespace gir