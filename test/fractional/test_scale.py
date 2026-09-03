import time
import unittest

from cami_amber.fractional_metrics import compute_sample_metrics
from test.fractional.helpers import genome_ids_json, load_assignments, load_truth_from_rows


class TestFractionalScale(unittest.TestCase):
    def test_thousands_of_components_without_per_base_expansion(self):
        n = 4000
        rows = []
        pred = []
        for i in range(n):
            seq = 'c{}'.format(i)
            rows.append({
                'SEQUENCEID': seq,
                '_LENGTH': 1000,
                'COMPONENT_BP': 700,
                'COMPONENT_TYPE': 'unique',
                'GENOME_IDS': genome_ids_json(['G{}'.format(i % 50)]),
            })
            rows.append({
                'SEQUENCEID': seq,
                '_LENGTH': 1000,
                'COMPONENT_BP': 300,
                'COMPONENT_TYPE': 'unique',
                'GENOME_IDS': genome_ids_json(['G{}'.format((i + 1) % 50)]),
            })
            pred.append({'SEQUENCEID': seq, 'BINID': 'b{}'.format(i % 80)})
        truth = load_truth_from_rows('scale', rows)
        assignments = load_assignments('scale', pred, truth)
        started = time.perf_counter()
        result = compute_sample_metrics(truth, assignments, [0.5], [0.1])
        elapsed = time.perf_counter() - started
        self.assertLess(elapsed, 30)
        self.assertEqual(result['truth_summary']['n_sequences'], n)
        self.assertGreater(result['metrics']['precision_weighted_bp'], 0)


if __name__ == '__main__':
    unittest.main()
