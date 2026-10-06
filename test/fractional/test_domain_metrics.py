import unittest

import pandas as pd

from cami_amber.fractional_metrics import build_support_matrices, compute_sample_metrics
from test.fractional.helpers import load_assignments, load_truth_from_rows
from test.fractional.test_fractional_metrics import _row


class DomainMetricsTest(unittest.TestCase):
    def setUp(self):
        self.truth = load_truth_from_rows('s1', [
            _row('mixed', 100, 70, 'unique', ['bacteria']),
            _row('mixed', 100, 30, 'unique', ['eukaryote']),
            _row('unbinned', 100, 100, 'unique', ['missing_bacteria']),
        ])
        self.assignments = load_assignments('s1', [{'SEQUENCEID': 'mixed', 'BINID': 'bin1'}], self.truth)

    def evaluate(self, selected=None, cache=None):
        return compute_sample_metrics(self.truth, self.assignments, [.9], [.05],
                                      selected_genomes=selected, support_cache=cache)

    def test_disabling_eukaryotes_preserves_foreign_bp(self):
        result = self.evaluate({'bacteria', 'missing_bacteria'})
        bin_row = result['precision_df'].iloc[0]
        self.assertEqual(bin_row['total_length'], 100)
        self.assertEqual(bin_row['foreign_bp'], 30)
        self.assertEqual(bin_row['precision_bp'], .7)
        self.assertEqual(result['recovered'][0]['count'], 0)
        self.assertEqual(result['metrics']['recall_avg_bp'], .5)
        self.assertAlmostEqual(result['metrics']['recall_weighted_bp'], 70 / 170)
        self.assertAlmostEqual(result['metrics']['percentage_of_assigned_bps'], 70 / 170)

    def test_bin_is_not_rematched_when_its_domain_is_disabled(self):
        result = self.evaluate({'eukaryote'})
        self.assertTrue(result['precision_df'].empty)
        self.assertEqual(result['genome_df'].iloc[0]['best_recall_bp'], 1.0)

    def test_all_selected_preserves_every_metric(self):
        expected = self.evaluate()
        actual = self.evaluate(set(self.truth.genomes))
        pd.testing.assert_frame_equal(pd.DataFrame([actual['metrics']]), pd.DataFrame([expected['metrics']]))
        pd.testing.assert_frame_equal(actual['precision_df'], expected['precision_df'])

    def test_empty_selection_is_undefined_not_a_perfect_score(self):
        result = self.evaluate(set())
        self.assertTrue(result['precision_df'].empty)
        self.assertTrue(pd.isna(result['metrics']['precision_weighted_bp']))
        self.assertTrue(pd.isna(result['metrics']['recall_weighted_bp']))

    def test_cached_support_is_unchanged(self):
        cache = build_support_matrices(self.truth, self.assignments)
        before = [frame.copy(deep=True) for frame in cache]
        actual = self.evaluate({'bacteria'}, cache)
        expected = self.evaluate({'bacteria'})
        pd.testing.assert_frame_equal(pd.DataFrame([actual['metrics']]), pd.DataFrame([expected['metrics']]))
        for frame, original in zip(cache, before):
            pd.testing.assert_frame_equal(frame, original)


if __name__ == '__main__':
    unittest.main()
