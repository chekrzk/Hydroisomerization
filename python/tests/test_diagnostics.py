from __future__ import annotations

import sys
import unittest
from pathlib import Path

from scipy.special import ive

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from hydroiso.plots import _first_order_cylinder_eta


class DiagnosticTests(unittest.TestCase):
    def test_cylindrical_effectiveness_matches_analytical_solution(self) -> None:
        phi = 2.0
        analytical = 2.0*ive(1, phi)/(phi*ive(0, phi))
        numerical = _first_order_cylinder_eta(phi, 200)
        self.assertLess(abs(numerical-analytical)/analytical, 2.0e-5)


if __name__ == "__main__":
    unittest.main()
