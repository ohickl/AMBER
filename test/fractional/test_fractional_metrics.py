import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from cami_amber.binning_classes import Metrics
from cami_amber.fractional_metrics import MATCH_COMPATIBLE_TIE, compute_fari, compute_sample_metrics
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

    def test_fari_hard_equals_sklearn(self):
        labels_true = np.array([0, 0, 1, 1, 1, 2])
        labels_pred = np.array([1, 1, 0, 0, 0, 2])
        k = 3
        n = len(labels_true)
        U = np.zeros((n, k))
        V = np.zeros((n, k))
        U[np.arange(n), labels_true] = 1
        V[np.arange(n), labels_pred] = 1
        fari = compute_fari(U, V)
        confusion = pd.DataFrame({
            'BINID': labels_pred,
            'genome_id': labels_true,
            'n': 1,
        }).groupby(['BINID', 'genome_id'], as_index=False).sum()
        _ri, ari = Metrics.compute_rand_index(confusion, 'BINID', 'genome_id', 'n')
        self.assertAlmostEqual(fari, ari, places=12)

    def test_fari_reflexive(self):
        U = np.array([[0.7, 0.3], [1.0, 0.0], [0.2, 0.8]])
        self.assertAlmostEqual(compute_fari(U, U), 1.0, places=12)

    def test_fari_label_permutation(self):
        U = np.array([[0.7, 0.3], [1.0, 0.0], [0.2, 0.8]])
        V = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        Vperm = V[:, [1, 0]]
        self.assertAlmostEqual(compute_fari(U, V), compute_fari(U, Vperm), places=12)

    def test_weighted_fari_matches_replication(self):
        U = np.array([[0.7, 0.3], [1.0, 0.0]])
        V = np.array([[1.0, 0.0], [0.0, 1.0]])
        weights = np.array([3, 2], dtype=float)
        closed = compute_fari(U, V, weights=weights)
        Urep = np.repeat(U, weights.astype(int), axis=0)
        Vrep = np.repeat(V, weights.astype(int), axis=0)
        replicated = compute_fari(Urep, Vrep)
        self.assertAlmostEqual(closed, replicated, places=12)


if __name__ == '__main__':
    unittest.main()
