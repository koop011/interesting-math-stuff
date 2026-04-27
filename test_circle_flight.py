"""
Tests for circle_flight.py — Vincenty direct formula and geodetic_to_ecef.

Reference values sourced from:
  - Vincenty (1975) original paper / Geoscience Australia
  - Movable Type Scripts: https://www.movable-type.co.uk/scripts/latlong-vincenty.html
  - ECEF values verified via independent WGS-84 conversion tools
"""

from __future__ import print_function

import math
import unittest

from circle_flight import (
    geodetic_to_ecef,
    generate_circle_points,
    get_circle_coords,
    vincenty_direct,
)


class TestVincentyDirect(unittest.TestCase):
    """Verify vincenty_direct against known geodetic test cases."""

    def assertLatLonAlmostEqual(self, result, expected, places=5, msg=None):
        """Helper: compare (lat, lon) tuples to N decimal places."""
        self.assertAlmostEqual(result[0], expected[0], places=places, msg=msg)
        self.assertAlmostEqual(result[1], expected[1], places=places, msg=msg)

    def test_flinders_peak_to_buninyong(self):
        """Classic Vincenty test case from the original 1975 paper.

        Flinders Peak -> Buninyong, Australia.
        Start:   -37.95103, 144.42487
        Bearing: 306.86816 deg
        Dist:    54972.271 m
        Expect:  -37.65282, 143.92650
        """
        lat, lon = vincenty_direct(-37.95103, 144.42487, 306.86816, 54972.271)
        self.assertAlmostEqual(lat, -37.65282, places=5)
        self.assertAlmostEqual(lon, 143.92650, places=5)

    def test_due_north(self):
        """Travel due north from the equator on the prime meridian.

        Start:   0.0, 0.0
        Bearing: 0 deg (north)
        Dist:    10000 m (10 km)
        Expected latitude ~0.0904 deg (10 km / ~110.6 km per degree).
        Longitude should remain 0.
        """
        lat, lon = vincenty_direct(0.0, 0.0, 0.0, 10000.0)
        self.assertAlmostEqual(lat, 0.0904370, places=5)
        self.assertAlmostEqual(lon, 0.0, places=7)

    def test_due_east_on_equator(self):
        """Travel due east along the equator.

        Start:   0.0, 0.0
        Bearing: 90 deg (east)
        Dist:    10000 m
        Expected: lat stays ~0, lon ~0.08984 deg (10 km / ~111.32 km per degree at equator).
        """
        lat, lon = vincenty_direct(0.0, 0.0, 90.0, 10000.0)
        self.assertAlmostEqual(lat, 0.0, places=5)
        self.assertAlmostEqual(lon, 0.08984, places=4)

    def test_due_south(self):
        """Due south should decrease latitude, keep longitude constant."""
        lat, lon = vincenty_direct(0.0, 0.0, 180.0, 10000.0)
        self.assertAlmostEqual(lat, -0.0904370, places=5)
        self.assertAlmostEqual(lon, 0.0, places=7)

    def test_due_west(self):
        """Due west should decrease longitude, keep latitude ~constant."""
        lat, lon = vincenty_direct(0.0, 0.0, 270.0, 10000.0)
        self.assertAlmostEqual(lat, 0.0, places=5)
        self.assertAlmostEqual(lon, -0.08984, places=4)

    def test_symmetry_north_south(self):
        """Going north then south by the same distance should return ~start."""
        start_lat, start_lon = -27.46794, 153.02809
        mid_lat, mid_lon = vincenty_direct(start_lat, start_lon, 0.0, 5000.0)
        end_lat, end_lon = vincenty_direct(mid_lat, mid_lon, 180.0, 5000.0)
        self.assertAlmostEqual(end_lat, start_lat, places=5)
        self.assertAlmostEqual(end_lon, start_lon, places=5)

    def test_high_latitude(self):
        """Verify behaviour near the poles (lat=80 deg).

        At high latitudes, longitude degrees shrink. 10 km east at 80N
        should produce a larger lon offset than at the equator.
        """
        lat, lon = vincenty_direct(80.0, 0.0, 90.0, 10000.0)
        self.assertAlmostEqual(lat, 80.0, places=3)
        # At 80N, 1 degree lon ~ 19.4 km, so 10 km ~ 0.516 deg
        self.assertAlmostEqual(lon, 0.5156, places=2)

    def test_brisbane_5km_circle_north(self):
        """Brisbane center, 5 km due north — sanity check."""
        lat, lon = vincenty_direct(-27.46794, 153.02809, 0.0, 5000.0)
        # Should be ~0.045 deg north of start
        self.assertTrue(lat > -27.46794)
        self.assertAlmostEqual(lat, -27.42285, places=4)
        self.assertAlmostEqual(lon, 153.02809, places=5)

    def test_zero_distance(self):
        """Zero distance should return the start point."""
        lat, lon = vincenty_direct(-27.46794, 153.02809, 42.0, 0.0)
        self.assertAlmostEqual(lat, -27.46794, places=7)
        self.assertAlmostEqual(lon, 153.02809, places=7)


class TestGeodeticToECEF(unittest.TestCase):
    """Verify geodetic_to_ecef against known ECEF coordinates."""

    def test_equator_prime_meridian(self):
        """(0, 0, 0) should be on the x-axis at distance = semi-major axis a."""
        X, Y, Z = geodetic_to_ecef(0.0, 0.0, 0.0)
        self.assertAlmostEqual(X, 6378137.0, places=1)
        self.assertAlmostEqual(Y, 0.0, places=1)
        self.assertAlmostEqual(Z, 0.0, places=1)

    def test_north_pole(self):
        """North pole should be on the z-axis at distance = semi-minor axis b."""
        X, Y, Z = geodetic_to_ecef(90.0, 0.0, 0.0)
        self.assertAlmostEqual(X, 0.0, places=1)
        self.assertAlmostEqual(Y, 0.0, places=1)
        self.assertAlmostEqual(Z, 6356752.3142, places=1)

    def test_south_pole(self):
        """South pole: negative z, x and y near zero."""
        X, Y, Z = geodetic_to_ecef(-90.0, 0.0, 0.0)
        self.assertAlmostEqual(X, 0.0, places=1)
        self.assertAlmostEqual(Y, 0.0, places=1)
        self.assertAlmostEqual(Z, -6356752.3142, places=1)

    def test_equator_90E(self):
        """(0, 90, 0) should be on the y-axis."""
        X, Y, Z = geodetic_to_ecef(0.0, 90.0, 0.0)
        self.assertAlmostEqual(X, 0.0, places=1)
        self.assertAlmostEqual(Y, 6378137.0, places=1)
        self.assertAlmostEqual(Z, 0.0, places=1)

    def test_brisbane(self):
        """Brisbane, Australia at sea level."""
        X, Y, Z = geodetic_to_ecef(-27.46794, 153.02809, 0.0)
        # Verified via external WGS-84 ECEF converter
        self.assertAlmostEqual(X, -5047170.58, places=0)
        self.assertAlmostEqual(Y, 2568545.79, places=0)
        self.assertAlmostEqual(Z, -2924318.63, places=0)

    def test_altitude(self):
        """Adding altitude should increase distance from Earth centre."""
        X0, Y0, Z0 = geodetic_to_ecef(0.0, 0.0, 0.0)
        X1, Y1, Z1 = geodetic_to_ecef(0.0, 0.0, 1000.0)
        dist0 = math.sqrt(X0 ** 2 + Y0 ** 2 + Z0 ** 2)
        dist1 = math.sqrt(X1 ** 2 + Y1 ** 2 + Z1 ** 2)
        self.assertAlmostEqual(dist1 - dist0, 1000.0, places=1)


class TestGenerateCirclePoints(unittest.TestCase):
    """Integration tests for the circle point generator."""

    def test_point_count(self):
        """Should generate exactly `resolution` points."""
        points = generate_circle_points(-17.367767, 116.884228, 5.0, 8)
        self.assertEqual(len(points), 8)

    def test_bearings(self):
        """Bearings should be evenly spaced starting at 0."""
        points = generate_circle_points(0.0, 0.0, 10.0, 4)
        bearings = [p["bearing_deg"] for p in points]
        self.assertEqual(bearings, [0.0, 90.0, 180.0, 270.0])

    def test_point_ids_sequential(self):
        """Point IDs should be 1-based sequential."""
        points = generate_circle_points(0.0, 0.0, 5.0, 6)
        ids = [p["point"] for p in points]
        self.assertEqual(ids, [1, 2, 3, 4, 5, 6])

    def test_north_point_higher_latitude(self):
        """Point at bearing 0 (north) should have higher latitude than centre."""
        points = generate_circle_points(-27.0, 153.0, 5.0, 4)
        north_point = points[0]
        self.assertGreater(north_point["lat"], -27.0)

    def test_south_point_lower_latitude(self):
        """Point at bearing 180 (south) should have lower latitude than centre."""
        points = generate_circle_points(-27.0, 153.0, 5.0, 4)
        south_point = points[2]  # bearing 180
        self.assertLess(south_point["lat"], -27.0)


class TestGetCircleCoords(unittest.TestCase):
    """Tests for the public get_circle_coords API."""

    def test_returns_tuples(self):
        coords = get_circle_coords(0.0, 0.0, 5.0, 4)
        self.assertEqual(len(coords), 4)
        for item in coords:
            self.assertIsInstance(item, tuple)
            self.assertEqual(len(item), 2)

    def test_resolution_too_low(self):
        with self.assertRaises(ValueError):
            get_circle_coords(0.0, 0.0, 5.0, 1)


if __name__ == "__main__":
    unittest.main()
