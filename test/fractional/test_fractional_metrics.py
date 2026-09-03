import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from cami_amber.binning_classes import Metrics
from cami_amber.fractional_metrics import (
    MATCH_COMPATIBLE_TIE,
    FractionalPredictionError,
    compute_sample_fari,
    compute_sample_metrics,
)
from test.fractional.helpers import genome_ids_json, load_assignments, load_truth_from_rows


def _row(seq, length, bp, kind, genomes):
    return {
        'SEQUENCEID': seq,
        '_LENGTH': length,
        'COMPONENT_BP': bp,
        'COMPONENT_TYPE': kind,
        'GENOME_IDS': genome_ids_json(genomes),
    }


class TestFractionalMetrics(unittest.TestCase):
    def _run(self, truth_rows, pred_rows, sample_id='s1'):
        truth = load_truth_from_rows(sample_id, truth_rows)
        assignments = load_assignments(sample_id, pred_rows, truth)
        return compute_sample_metrics(truth, assignments, min_completeness=[0.5], max_contamination=[0.1])

    def test_pure_hard(self):
        result = self._run(
            [_row('c1', 100, 100, 'unique', ['A'])],
            [{'SEQUENCEID': 'c1', 'BINID': 'binA'}],
        )
        bin_row = result['precision_df'].iloc[0]
        self.assertEqual(bin_row['precision_bp'], 1.0)
        self.assertEqual(bin_row['precision_seq'], 1.0)
        self.assertEqual(bin_row['recall_bp'], 1.0)

    def test_chimera_70_30_in_a(self):
        result = self._run(
            [_row('c1', 100, 70, 'unique', ['A']), _row('c1', 100, 30, 'unique', ['B'])],
            [{'SEQUENCEID': 'c1', 'BINID': 'binA'}],
        )
        bin_row = result['precision_df'].iloc[0]
        self.assertAlmostEqual(bin_row['precision_bp'], 0.70)
        self.assertAlmostEqual(bin_row['precision_seq'], 0.70)
        self.assertEqual(bin_row['matched_genome_id'], 'A')
        self.assertEqual(bin_row['unique_tp_length'], 70)
        self.assertEqual(bin_row['foreign_bp'], 30)

    def test_chimera_70_30_in_b(self):
        """A lone 70/30 contig is matched by max compatible support (A), so purity is 0.70.

        To match B, the bin needs strictly more unique B support than A.
        """
        result = self._run(
            [
                _row('c1', 100, 70, 'unique', ['A']),
                _row('c1', 100, 30, 'unique', ['B']),
                _row('c2', 80, 80, 'unique', ['B']),
            ],
            [{'SEQUENCEID': 'c1', 'BINID': 'binB'}, {'SEQUENCEID': 'c2', 'BINID': 'binB'}],
        )
        bin_row = result['precision_df'].iloc[0]
        self.assertEqual(bin_row['matched_genome_id'], 'B')
        self.assertAlmostEqual(bin_row['precision_bp'], 110 / 180)

    def test_fully_compatible(self):
        result = self._run(
            [_row('c1', 100, 100, 'compatible', ['A', 'B'])],
            [{'SEQUENCEID': 'c1', 'BINID': 'binX'}],
        )
        bin_row = result['precision_df'].iloc[0]
        self.assertEqual(bin_row['precision_bp'], 1.0)
        self.assertEqual(bin_row['match_status'], MATCH_COMPATIBLE_TIE)
        self.assertTrue(np.isnan(bin_row['recall_bp']))
        self.assertEqual(result['recovered'][0]['count'], 0)
        self.assertEqual(result['genome_df'].loc[result['genome_df']['genome_id'] == 'A', 'identifiable_truth_bp'].iloc[0], 0)

    def test_compatible_in_foreign_bin_with_unique_c(self):
        result = self._run(
            [
                _row('c1', 100, 100, 'compatible', ['A', 'B']),
                _row('c2', 200, 200, 'unique', ['C']),
            ],
            [{'SEQUENCEID': 'c1', 'BINID': 'binC'}, {'SEQUENCEID': 'c2', 'BINID': 'binC'}],
        )
        bin_row = result['precision_df'].iloc[0]
        self.assertEqual(str(bin_row['matched_genome_id']), 'C')
        self.assertAlmostEqual(bin_row['precision_bp'], 200 / 300)

    def test_mixed_unique_compatible(self):
        result = self._run(
            [
                _row('c1', 100, 60, 'unique', ['A']),
                _row('c1', 100, 30, 'compatible', ['A', 'B']),
                _row('c1', 100, 10, 'unique', ['B']),
            ],
            [{'SEQUENCEID': 'c1', 'BINID': 'binA'}],
        )
        bin_row = result['precision_df'].iloc[0]
        self.assertAlmostEqual(bin_row['precision_bp'], 0.90)
        self.assertEqual(bin_row['unique_tp_length'], 60)
        self.assertAlmostEqual(bin_row['recall_bp'], 1.0)

    def test_unresolved(self):
        result = self._run(
            [_row('c1', 100, 70, 'unique', ['A']), _row('c1', 100, 30, 'unresolved', [])],
            [{'SEQUENCEID': 'c1', 'BINID': 'binA'}],
        )
        bin_row = result['precision_df'].iloc[0]
        self.assertAlmostEqual(bin_row['precision_bp'], 0.70)
        self.assertAlmostEqual(bin_row['precision_bp_resolved'], 1.0)

    def test_tie_broken_by_unique(self):
        result = self._run(
            [
                _row('c1', 100, 40, 'compatible', ['A', 'B']),
                _row('c1', 100, 40, 'unique', ['A']),
                _row('c1', 100, 20, 'unique', ['B']),
            ],
            [{'SEQUENCEID': 'c1', 'BINID': 'binX'}],
        )
        self.assertEqual(result['precision_df'].iloc[0]['matched_genome_id'], 'A')
        self.assertEqual(result['precision_df'].iloc[0]['match_status'], 'resolved')

    def test_fractional_sequence_units(self):
        result = self._run(
            [
                _row('c1', 100, 70, 'unique', ['A']),
                _row('c1', 100, 30, 'unique', ['B']),
                _row('c2', 100, 100, 'unique', ['A']),
            ],
            [{'SEQUENCEID': 'c1', 'BINID': 'binA'}, {'SEQUENCEID': 'c2', 'BINID': 'binA'}],
        )
        self.assertAlmostEqual(result['precision_df'].iloc[0]['precision_seq'], 0.85)

    def test_partially_binned(self):
        result = self._run(
            [
                _row('c1', 100, 100, 'unique', ['A']),
                _row('c2', 100, 100, 'unique', ['B']),
            ],
            [{'SEQUENCEID': 'c1', 'BINID': 'binA'}],
        )
        self.assertAlmostEqual(result['metrics']['percentage_of_assigned_bps'], 0.5)
        self.assertAlmostEqual(result['metrics']['accuracy_bp'], 0.5)
        self.assertAlmostEqual(result['metrics']['precision_weighted_bp'], 1.0)
        self.assertLess(result['metrics']['recall_weighted_bp'], 1.0)

    def test_zero_identifiable_genome(self):
        result = self._run(
            [_row('c1', 100, 100, 'compatible', ['A', 'B'])],
            [{'SEQUENCEID': 'c1', 'BINID': 'binX'}],
        )
        a = result['genome_df'][result['genome_df']['genome_id'] == 'A'].iloc[0]
        self.assertEqual(a['identifiable_truth_bp'], 0)
        self.assertTrue(np.isnan(a['best_recall_bp']))

    def test_unknown_prediction_sequence_is_fatal(self):
        truth = load_truth_from_rows('s1', [_row('c1', 10, 10, 'unique', ['A'])])
        with self.assertRaises(Exception):
            load_assignments('s1', [{'SEQUENCEID': 'missing', 'BINID': 'b'}], truth)

    def test_repeated_unique_components_accumulate_membership(self):
        result = self._run(
            [
                _row('c1', 100, 40, 'unique', ['A']),
                _row('c1', 100, 30, 'unique', ['A']),
                _row('c1', 100, 30, 'unique', ['B']),
            ],
            [{'SEQUENCEID': 'c1', 'BINID': 'binA'}],
        )
        self.assertAlmostEqual(result['precision_df'].iloc[0]['precision_bp'], 0.70)
        self.assertEqual(result['precision_df'].iloc[0]['unique_tp_length'], 70)

    def test_identifiable_fraction_not_assembly_length(self):
        result = self._run(
            [
                _row('tiny', 10, 10, 'unique', ['A']),
                _row('big', 1000, 1000, 'unique', ['B']),
            ],
            [{'SEQUENCEID': 'tiny', 'BINID': 'binA'}, {'SEQUENCEID': 'big', 'BINID': 'binB'}],
        )
        a_bin = result['precision_df'][result['precision_df']['matched_genome_id'] == 'A'].iloc[0]
        self.assertAlmostEqual(a_bin['matched_genome_identifiable_fraction'], 1.0)
        self.assertNotAlmostEqual(a_bin['matched_genome_identifiable_fraction'], 10 / 1010)

    def test_min_length_drops_short_prediction(self):
        from cami_amber.fractional_truth import load_fractional_truth_file, write_fractional_truth
        fd, path = tempfile.mkstemp(suffix='.tsv')
        os.close(fd)
        try:
            write_fractional_truth(path, 's1', [
                _row('short', 10, 10, 'unique', ['A']),
                _row('long', 100, 100, 'unique', ['B']),
            ])
            truth = load_fractional_truth_file(path)['s1'].filter_min_length(50)
            assignments = load_assignments('s1', [
                {'SEQUENCEID': 'short', 'BINID': 'b'},
                {'SEQUENCEID': 'long', 'BINID': 'b'},
            ], truth)
            self.assertEqual(list(assignments['SEQUENCEID']), ['long'])
        finally:
            os.remove(path)

    def test_unknown_prediction_still_fatal_after_min_length(self):
        from cami_amber.fractional_truth import load_fractional_truth_file, write_fractional_truth
        fd, path = tempfile.mkstemp(suffix='.tsv')
        os.close(fd)
        try:
            write_fractional_truth(path, 's1', [_row('long', 100, 100, 'unique', ['A'])])
            truth = load_fractional_truth_file(path)['s1'].filter_min_length(50)
            with self.assertRaises(FractionalPredictionError):
                load_assignments('s1', [{'SEQUENCEID': 'ghost', 'BINID': 'b'}], truth)
        finally:
            os.remove(path)

    def test_compatible_not_in_completeness(self):
        result = self._run(
            [
                _row('c1', 100, 50, 'unique', ['A']),
                _row('c1', 100, 50, 'compatible', ['A', 'B']),
            ],
            [{'SEQUENCEID': 'c1', 'BINID': 'binA'}],
        )
        a = result['genome_df'][result['genome_df']['genome_id'] == 'A'].iloc[0]
        self.assertEqual(a['identifiable_truth_bp'], 50)
        self.assertAlmostEqual(a['best_recall_bp'], 1.0)
        self.assertEqual(a['compatible_truth_bp_involving_genome'], 50)

    def test_unresolved_not_correct(self):
        result = self._run(
            [_row('c1', 100, 70, 'unique', ['A']), _row('c1', 100, 30, 'unresolved', [])],
            [{'SEQUENCEID': 'c1', 'BINID': 'binA'}],
        )
        self.assertAlmostEqual(result['precision_df'].iloc[0]['precision_bp'], 0.70)
        self.assertLess(result['precision_df'].iloc[0]['precision_bp'], 1.0)

    def test_identifiable_bp_ari_matches_explicit_base_expansion(self):
        result = self._run(
            [
                _row('c1', 100, 70, 'unique', ['A']),
                _row('c1', 100, 30, 'unique', ['B']),
                _row('c2', 100, 100, 'unique', ['A']),
                _row('c3', 100, 100, 'unique', ['B']),
            ],
            [
                {'SEQUENCEID': 'c1', 'BINID': 'bin1'},
                {'SEQUENCEID': 'c2', 'BINID': 'bin1'},
                {'SEQUENCEID': 'c3', 'BINID': 'bin2'},
            ],
        )
        true_lab = ['A'] * 70 + ['B'] * 30 + ['A'] * 100 + ['B'] * 100
        pred_lab = ['bin1'] * 100 + ['bin1'] * 100 + ['bin2'] * 100
        confusion = pd.DataFrame({'BINID': pred_lab, 'genome_id': true_lab, 'bp': 1}).groupby(
            ['BINID', 'genome_id'], as_index=False
        ).sum()
        ri, ari = Metrics.compute_rand_index(confusion, 'BINID', 'genome_id', 'bp')
        self.assertAlmostEqual(result['metrics']['adjusted_rand_index_bp_identifiable'], ari, places=12)
        self.assertAlmostEqual(result['metrics']['rand_index_bp_identifiable'], ri, places=12)
        self.assertAlmostEqual(result['metrics']['ari_bp_identifiable_fraction'], 1.0)

    def test_length_weighted_fari_is_not_identifiable_bp_ari(self):
        truth_rows = [
            _row('c1', 100, 70, 'unique', ['A']),
            _row('c1', 100, 30, 'unique', ['B']),
            _row('c2', 100, 100, 'unique', ['A']),
            _row('c3', 100, 100, 'unique', ['B']),
        ]
        pred_rows = [
            {'SEQUENCEID': 'c1', 'BINID': 'bin1'},
            {'SEQUENCEID': 'c2', 'BINID': 'bin1'},
            {'SEQUENCEID': 'c3', 'BINID': 'bin2'},
        ]
        result = self._run(truth_rows, pred_rows)
        truth = load_truth_from_rows('s1', truth_rows)
        assignments = load_assignments('s1', pred_rows, truth)
        weighted_fari, _, _, _ = compute_sample_fari(truth, assignments, weighted=True)
        self.assertTrue(np.isfinite(weighted_fari))
        self.assertGreater(
            abs(weighted_fari - result['metrics']['adjusted_rand_index_bp_identifiable']),
            1e-6,
        )

    def test_chimera_not_whole_contig_from_a(self):
        result = self._run(
            [_row('c1', 100, 70, 'unique', ['A']), _row('c1', 100, 30, 'unique', ['B'])],
            [{'SEQUENCEID': 'c1', 'BINID': 'anything'}],
        )
        self.assertAlmostEqual(result['precision_df'].iloc[0]['precision_bp'], 0.70)
        self.assertNotAlmostEqual(result['precision_df'].iloc[0]['precision_bp'], 1.0)


if __name__ == '__main__':
    unittest.main()
