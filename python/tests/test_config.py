from __future__ import annotations

import copy
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from hydroiso.config import ConfigurationError, apply_overrides, load_config, write_fortran_namelist


ROOT = Path(__file__).resolve().parents[2]


class ConfigTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = load_config(ROOT / "config.yaml")

    def test_reaction_15_is_irreversible(self) -> None:
        reaction = self.config["kinetics"]["reactions"]["R15"]
        self.assertFalse(reaction["reversible"])
        self.assertNotIn("reverse", reaction)

    def test_scientific_rejects_demo_values(self) -> None:
        changed = copy.deepcopy(self.config)
        changed["run_mode"] = "scientific"
        with self.assertRaises(ConfigurationError):
            apply_overrides(changed, [])

    def test_multiple_overrides(self) -> None:
        changed = apply_overrides(self.config, [
            "conditions.temperature.value=613.15",
            "numerics.spatial.radial_cells=25",
        ])
        self.assertEqual(changed["conditions"]["temperature"]["value"], 613.15)
        self.assertEqual(changed["numerics"]["spatial"]["radial_cells"], 25)
        expected_total = 3.0e6/(8.314462618*613.15)
        self.assertAlmostEqual(
            changed["kinetics"]["demo_embedding"]["normalization_concentration"]["value"],
            expected_total,
        )
        self.assertAlmostEqual(changed["bulk"]["concentrations"]["C1"]["value"], 0.07*expected_total)

    def test_pore_properties_are_recalculated(self) -> None:
        changed = apply_overrides(self.config, ["pellet.mesopore_volume.value=3.0e-4"])
        vpore = 3.7e-4
        self.assertAlmostEqual(changed["pellet"]["total_pore_volume"]["value"], vpore)
        self.assertAlmostEqual(
            changed["pellet"]["pore_radius"]["value"],
            2.0*vpore/changed["pellet"]["bet_surface_area"]["value"],
        )

    def test_hydrogen_pretreated_initial_state(self) -> None:
        changed = apply_overrides(self.config, ["initial_conditions.model=hydrogen_pretreated"])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.nml"
            write_fortran_namelist(changed, path)
            text = path.read_text(encoding="utf-8")
        initial_line = next(line for line in text.splitlines() if line.strip().startswith("c_initial"))
        self.assertIn("5.657", initial_line)

    def test_namelist_converts_activation_energy_to_si(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.nml"
            write_fortran_namelist(self.config, path)
            text = path.read_text(encoding="utf-8")
        self.assertIn("activation_forward", text)
        self.assertIn("2.18404800000000017d+05", text)


if __name__ == "__main__":
    unittest.main()
