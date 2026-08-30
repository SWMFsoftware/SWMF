#!/usr/bin/env python3
"""Focused regression checks for the RIM-to-IPE export path."""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src" / "RIM_grid_comp.f90"


def subroutine_source(name):
    source = SOURCE.read_text(encoding="ascii")
    start = source.index(f"  subroutine {name}")
    end = source.index(f"  end subroutine {name}", start)
    return source[start:end]


def refresh_snapshot(data, north, south, nlat=181, nlon=361):
    """Refresh a combined grid with update_export_state's inverse mapping."""
    for ilat in range(1, 2 * nlat):
        is_north = ilat >= nlat
        itheta = 2 * nlat - ilat if is_north else nlat - ilat + 1
        source = north if is_north else south
        for ilon in range(1, nlon + 1):
            ipsi = (ilon + nlon // 2 - 1) % (nlon - 1) + 1
            data[ilon, ilat] = source(itheta, ipsi)


class RimExportMappingTest(unittest.TestCase):
    def test_each_run_refreshes_after_coordinate_setup(self):
        source = subroutine_source("my_run")
        source = source[source.index("deallocate(Data_VII") :]
        match = re.search(
            r"if\(DoShiftDataCoupling\) then(?P<shifted>.*?)"
            r"else(?P<unshifted>.*?)end if(?P<after>.*?)"
            r"call write_log\(\"RIM_grid_comp:run routine returned\"\)",
            source,
            re.DOTALL,
        )
        self.assertIsNotNone(match)
        self.assertNotIn("update_export_state", match.group("shifted"))
        self.assertNotIn("update_export_state", match.group("unshifted"))
        self.assertEqual(
            match.group("after").count("call update_export_state(gComp, iError)"),
            1,
        )
        self.assertEqual(source.count("call update_export_state(gComp, iError)"), 1)

    def test_global_mapping_and_two_successive_refreshes(self):
        source = subroutine_source("update_export_state")
        self.assertIn("iTheta = 2*nLat - j", source)
        self.assertIn("iTheta = nLat - j + 1", source)
        self.assertIn(
            "iPsi = modulo(i + nLon/2 - 1, nLon - 1) + 1", source
        )

        data = {}
        refresh_snapshot(
            data,
            lambda itheta, ipsi: 1_000_000 + 1_000 * itheta + ipsi,
            lambda itheta, ipsi: -1_000_000 - 1_000 * itheta - ipsi,
        )
        first = data.copy()
        refresh_snapshot(
            data,
            lambda itheta, ipsi: 2_000_000 + 1_000 * itheta + ipsi,
            lambda itheta, ipsi: -2_000_000 - 1_000 * itheta - ipsi,
        )
        second = data.copy()

        self.assertEqual(first[181, 181], 1_181_001)
        self.assertEqual(first[181, 361], 1_001_001)
        self.assertEqual(first[181, 180], -1_002_001)
        self.assertEqual(first[181, 1], -1_181_001)
        self.assertEqual(first[1, 181], first[361, 181])
        self.assertEqual(first[1, 180], first[361, 180])
        self.assertEqual(second[181, 181], 2_181_001)
        self.assertEqual(second[181, 180], -2_002_001)
        self.assertNotEqual(first[181, 181], second[181, 181])
        self.assertNotEqual(first[181, 180], second[181, 180])

    def test_real_branch_does_not_use_analytic_values(self):
        source = subroutine_source("update_export_state")
        formula = "Data_VII(iVar,i,j) = abs(Lon_I(i))*abs(Lat_I(j))*Coef"
        test_start = source.index("if(DoTest) then")
        formula_at = source.index(formula)
        real_start = source.index("else", formula_at)
        self.assertLess(test_start, formula_at)
        self.assertLess(formula_at, real_start)
        self.assertNotIn(formula, source[real_start:])


if __name__ == "__main__":
    unittest.main()
