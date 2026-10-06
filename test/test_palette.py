"""Cohort colours must survive changes in order and visible subsets."""
import json
import tempfile
import unittest
from pathlib import Path

from distinctipy.colorblind import colorblind_filter
from cami_amber.palette import categorical_colors, colors_for_labels, method_palette, save_palette


class PaletteTests(unittest.TestCase):
    def test_complete_cohort_and_filtered_order(self):
        labels = [f'tool_{i:02}' for i in range(64)]
        palette = method_palette(labels)
        self.assertEqual(palette, method_palette(list(reversed(labels))))
        self.assertEqual(len(set(palette.values())), 64)
        subset = [labels[40], labels[3]]
        self.assertEqual(colors_for_labels(subset, palette), [palette[label] for label in subset])
        with self.assertRaises(KeyError):
            colors_for_labels(['missing'], palette)

    def test_seed_ranges_simulations_and_receipt(self):
        colors = categorical_colors(32)
        self.assertEqual(colors, categorical_colors(32))
        for color in colors:
            self.assertTrue(all(0 <= channel <= 1 for channel in color))
            self.assertNotIn(color, [(0, 0, 0), (1, 1, 1)])
            for mode in ('Deuteranopia', 'Protanopia', 'Tritanopia'):
                self.assertTrue(all(0 <= channel <= 1 for channel in colorblind_filter(color, mode)))
        palette = method_palette(['a', 'b'])
        with tempfile.TemporaryDirectory() as output:
            save_palette(output, palette)
            receipt = json.loads(Path(output, 'method_palette.json').read_text())
            self.assertEqual(receipt['seed'], 42)
            self.assertEqual(receipt['version'], '1.3.4')
            self.assertEqual(receipt['colors'], palette)


if __name__ == '__main__':
    unittest.main()
