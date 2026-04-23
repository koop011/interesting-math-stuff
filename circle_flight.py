#!/usr/bin/env python3
"""
Circle point generator.

Given a center lon/lat and a radius in km, generate evenly spaced
points around the circle. The number of points is controlled by
the --resolution argument.

Uses spherical Earth destination-point formula so the lon/lat
outputs are geographically accurate.
"""

import argparse
import json
import math
import sys

EARTH_RADIUS_KM = 6371.0


def destination_point(lat_deg, lon_deg, bearing_deg, distance_km):
    """
    Given a start point (lat/lon in degrees), a bearing (degrees clockwise
    from north) and a distance (km), return the destination (lat, lon) in
    degrees.  Spherical-Earth approximation — accurate enough for <100 km.
    """
    lat1 = math.radians(lat_deg)
    lon1 = math.radians(lon_deg)
    brng = math.radians(bearing_deg)
    d_r = distance_km / EARTH_RADIUS_KM

    lat2 = math.asin(
        math.sin(lat1) * math.cos(d_r)
        + math.cos(lat1) * math.sin(d_r) * math.cos(brng)
    )
    lon2 = lon1 + math.atan2(
        math.sin(brng) * math.sin(d_r) * math.cos(lat1),
        math.cos(d_r) - math.sin(lat1) * math.sin(lat2),
    )
    return math.degrees(lat2), math.degrees(lon2)


def generate_circle_points(center_lat, center_lon, radius_km, resolution):
    """
    Generate `resolution` evenly spaced points around a circle.
    Returns list of dicts with point_id, lat, lon, bearing_deg.
    """
    points = []
    for i in range(resolution):
        bearing_deg = (360.0 / resolution) * i
        lat, lon = destination_point(center_lat, center_lon, bearing_deg, radius_km)
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
        cx, cy = to_grid(center_lat, center_lon)
        # Use a point on the circle at this angle
        avg_r_lat = (max(lats) - min(lats)) / 2 if len(lats) > 1 else pad_lat
        avg_r_lon = (max(lons) - min(lons)) / 2 if len(lons) > 1 else pad_lon
        clat = center_lat + avg_r_lat * math.cos(angle)
        clon = center_lon + avg_r_lon * math.sin(angle)
        gx, gy = to_grid(clat, clon)
        if grid[gy][gx] == " ":
            grid[gy][gx] = "·"

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
    border_h = "─" * (width + 2)
    print(f"\n┌{border_h}┐")
    for row in grid:
        print(f"│ {''.join(row)} │")
    print(f"└{border_h}┘")
    print(f"  + = center ({center_lat}, {center_lon})")
    print(f"  · = circle outline")
    print(f"  # = point positions\n")


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
            print(f"{p['point']},{p['lat']},{p['lon']},{p['bearing_deg']}")
    else:
        json.dump(points, sys.stdout, indent=2)
        print()

    draw_ascii_circle(points, args.center_lat, args.center_lon)


if __name__ == "__main__":
    main()
