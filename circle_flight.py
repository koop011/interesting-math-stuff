#!/usr/bin/env python
"""
Circle point generator (WGS-84 ellipsoid).

Given a center lon/lat and a radius in km, generate evenly spaced
points around the circle. The number of points is controlled by
the --resolution argument.

Uses Vincenty's direct formula on the WGS-84 ellipsoid for
geodetically accurate destination-point computation.
"""

from __future__ import print_function

import argparse
import json
import math
import sys

# WGS-84 ellipsoid constants
a = 6378137.0              # semi-major axis (metres)
f = 1.0 / 298.257223563    # flattening
b = a * (1.0 - f)          # semi-minor axis (metres)
e2 = 2.0 * f - f ** 2      # eccentricity squared


def vincenty_direct(lat1_deg, lon1_deg, bearing_deg, distance_m):
    """
    Vincenty's direct formula.  Given a start point (lat/lon in degrees),
    a bearing (degrees clockwise from north) and a distance (metres),
    return the destination (lat, lon) in degrees on the WGS-84 ellipsoid.
    """
    phi1 = math.radians(lat1_deg)
    alpha1 = math.radians(bearing_deg)
    s = distance_m

    sin_alpha1 = math.sin(alpha1)
    cos_alpha1 = math.cos(alpha1)

    # Reduced latitude
    tan_U1 = (1.0 - f) * math.tan(phi1)
    cos_U1 = 1.0 / math.sqrt(1.0 + tan_U1 ** 2)
    sin_U1 = tan_U1 * cos_U1

    sigma1 = math.atan2(tan_U1, cos_alpha1)
    sin_alpha = cos_U1 * sin_alpha1
    cos2_alpha = 1.0 - sin_alpha ** 2

    u2 = cos2_alpha * (a ** 2 - b ** 2) / (b ** 2)
    A = 1.0 + u2 / 16384.0 * (4096.0 + u2 * (-768.0 + u2 * (320.0 - 175.0 * u2)))
    B = u2 / 1024.0 * (256.0 + u2 * (-128.0 + u2 * (74.0 - 47.0 * u2)))

    sigma = s / (b * A)
    for _iteration in range(200):
        cos_2sigma_m = math.cos(2.0 * sigma1 + sigma)
        sin_sigma = math.sin(sigma)
        cos_sigma = math.cos(sigma)

        delta_sigma = B * sin_sigma * (
            cos_2sigma_m
            + B / 4.0 * (
                cos_sigma * (-1.0 + 2.0 * cos_2sigma_m ** 2)
                - B / 6.0 * cos_2sigma_m * (-3.0 + 4.0 * sin_sigma ** 2)
                  * (-3.0 + 4.0 * cos_2sigma_m ** 2)
            )
        )
        sigma_prev = sigma
        sigma = s / (b * A) + delta_sigma
        if abs(sigma - sigma_prev) < 1e-12:
            break

    sin_sigma = math.sin(sigma)
    cos_sigma = math.cos(sigma)
    cos_2sigma_m = math.cos(2.0 * sigma1 + sigma)

    phi2 = math.atan2(
        sin_U1 * cos_sigma + cos_U1 * sin_sigma * cos_alpha1,
        (1.0 - f) * math.sqrt(
            sin_alpha ** 2
            + (sin_U1 * sin_sigma - cos_U1 * cos_sigma * cos_alpha1) ** 2
        ),
    )
    lam = math.atan2(
        sin_sigma * sin_alpha1,
        cos_U1 * cos_sigma - sin_U1 * sin_sigma * cos_alpha1,
    )
    C = f / 16.0 * cos2_alpha * (4.0 + f * (4.0 - 3.0 * cos2_alpha))
    L = lam - (1.0 - C) * f * sin_alpha * (
        sigma + C * sin_sigma * (
            cos_2sigma_m + C * cos_sigma * (-1.0 + 2.0 * cos_2sigma_m ** 2)
        )
    )
    lon2 = math.radians(lon1_deg) + L

    return math.degrees(phi2), math.degrees(lon2)


def geodetic_to_ecef(lat_deg, lon_deg, alt_m):
    """
    Convert geodetic coordinates (lat, lon in degrees; altitude in metres)
    to Earth-Centred Earth-Fixed (ECEF) cartesian coordinates (X, Y, Z)
    on the WGS-84 ellipsoid.
    """
    phi = math.radians(lat_deg)
    lam = math.radians(lon_deg)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)
    sin_lam = math.sin(lam)
    cos_lam = math.cos(lam)

    N = a / math.sqrt(1.0 - e2 * sin_phi ** 2)

    X = (N + alt_m) * cos_phi * cos_lam
    Y = (N + alt_m) * cos_phi * sin_lam
    Z = (N * (1.0 - e2) + alt_m) * sin_phi

    return X, Y, Z


def generate_circle_points(center_lat, center_lon, radius_km, resolution):
    """
    Generate `resolution` evenly spaced points around a circle.
    Returns list of dicts with point_id, lat, lon, bearing_deg.
    """
    distance_m = radius_km * 1000.0
    points = []
    for i in range(resolution):
        bearing_deg = (360.0 / resolution) * i
        lat, lon = vincenty_direct(center_lat, center_lon, bearing_deg, distance_m)
        points.append({
            "point": i + 1,
            "lat": round(lat, 7),
            "lon": round(lon, 7),
            "bearing_deg": round(bearing_deg, 2),
        })
    return points


def get_circle_coords(center_lat, center_lon, radius_km, resolution):
    """
    Returns a list of (lat, lon) tuples evenly spaced around a circle.

    Usage:
        from circle_flight import get_circle_coords

        coords = get_circle_coords(-17.367767, 116.884228, radius_km=5, resolution=8)
        for lat, lon in coords:
            print(lat, lon)
    """
    if resolution < 2:
        raise ValueError("resolution must be >= 2")
    points = generate_circle_points(center_lat, center_lon, radius_km, resolution)
    return [(p["lat"], p["lon"]) for p in points]


def draw_ascii_circle(points, center_lat, center_lon, width=60, height=30):
    """
    Draw an ASCII plot showing the circle outline and numbered point positions.
    """
    # Build a grid
    grid = [[" " for _ in range(width)] for _ in range(height)]

    # Normalize point positions to grid coordinates
    # x = lon, y = lat (invert y because terminal rows go top-to-bottom)
    lons = [p["lon"] for p in points]
    lats = [p["lat"] for p in points]
    all_lons = lons + [center_lon]
    all_lats = lats + [center_lat]

    min_lon, max_lon = min(all_lons), max(all_lons)
    min_lat, max_lat = min(all_lats), max(all_lats)

    # Add padding so points aren't on the very edge
    pad_lon = (max_lon - min_lon) * 0.15 or 0.001
    pad_lat = (max_lat - min_lat) * 0.15 or 0.001
    min_lon -= pad_lon
    max_lon += pad_lon
    min_lat -= pad_lat
    max_lat += pad_lat

    def to_grid(lat, lon):
        gx = int((lon - min_lon) / (max_lon - min_lon) * (width - 1))
        gy = int((max_lat - lat) / (max_lat - min_lat) * (height - 1))  # flip y
        gx = max(0, min(width - 1, gx))
        gy = max(0, min(height - 1, gy))
        return gx, gy

    # Draw circle outline using many intermediate angles
    steps = max(120, len(points) * 20)
    for i in range(steps):
        angle = 2 * math.pi * i / steps
        # Place relative to center using the same scale as the points
        avg_r_lat = (max(lats) - min(lats)) / 2 if len(lats) > 1 else pad_lat
        avg_r_lon = (max(lons) - min(lons)) / 2 if len(lons) > 1 else pad_lon
        clat = center_lat + avg_r_lat * math.cos(angle)
        clon = center_lon + avg_r_lon * math.sin(angle)
        gx, gy = to_grid(clat, clon)
        if grid[gy][gx] == " ":
            grid[gy][gx] = u"\u00b7"

    # Draw center
    cx, cy = to_grid(center_lat, center_lon)
    grid[cy][cx] = "+"

    # Draw points (use point number as label)
    for p in points:
        gx, gy = to_grid(p["lat"], p["lon"])
        label = str(p["point"])
        # Write label chars onto grid
        for ci, ch in enumerate(label):
            px = gx + ci
            if 0 <= px < width:
                grid[gy][px] = ch

    # Render
    border_h = "-" * (width + 2)
    print("")
    print("+%s+" % border_h)
    for row in grid:
        print("| %s |" % "".join(row))
    print("+%s+" % border_h)
    print("  + = center (%s, %s)" % (center_lat, center_lon))
    print("  . = circle outline")
    print("  # = point positions")
    print("")


def main():
    parser = argparse.ArgumentParser(
        description="Generate evenly spaced lon/lat points around a circle."
    )
    parser.add_argument(
        "--center-lat", type=float, default=-17.367767,
        help="Latitude of circle center (default: -17.367767)",
    )
    parser.add_argument(
        "--center-lon", type=float, default=116.884228,
        help="Longitude of circle center (default: 116.884228)",
    )
    parser.add_argument(
        "--radius-km", type=float, required=True,
        help="Radius of the circle in kilometres",
    )
    parser.add_argument(
        "--resolution", type=int, required=True,
        help="Number of points on the circle (must be >= 2)",
    )
    parser.add_argument(
        "--format", choices=["json", "csv"], default="json",
        help="Output format (default: json)",
    )

    args = parser.parse_args()

    if args.resolution < 2:
        print("Error: --resolution must be >= 2", file=sys.stderr)
        sys.exit(1)

    points = generate_circle_points(
        args.center_lat, args.center_lon, args.radius_km, args.resolution
    )

    if args.format == "csv":
        print("point,lat,lon,bearing_deg")
        for p in points:
            print("%s,%s,%s,%s" % (p["point"], p["lat"], p["lon"], p["bearing_deg"]))
    else:
        json.dump(points, sys.stdout, indent=2)
        print()

    draw_ascii_circle(points, args.center_lat, args.center_lon)


if __name__ == "__main__":
    main()
