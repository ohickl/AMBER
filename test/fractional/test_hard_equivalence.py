import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from cami_amber.binning_classes import GenomeQuery, Metrics, Options
from cami_amber.fractional_metrics import compute_fari, compute_sample_metrics
from cami_amber.fractional_truth import write_fractional_truth
from cami_amber.utils import load_data
from test.fractional.helpers import genome_ids_json, write_prediction


class TestHardEquivalence(unittest.TestCase):
    def test_one_hot_truth_matches_hard_scorer(self):
        sample_id = 'tiny'
        gs_rows = [
            ('c1', 'A', 100),
            ('c2', 'A', 50),
            ('c3', 'B', 80),
            ('c4', 'B', 20),
        ]
        pred_rows = [
            {'SEQUENCEID': 'c1', 'BINID': 'b1'},
            {'SEQUENCEID': 'c2', 'BINID': 'b1'},
            {'SEQUENCEID': 'c3', 'BINID': 'b2'},
        ]
        fd_gs, gs_path = tempfile.mkstemp(suffix='.binning')
        os.close(fd_gs)
        fd_pred, pred_path = tempfile.mkstemp(suffix='.binning')
        os.close(fd_pred)
        fd_frac, frac_path = tempfile.mkstemp(suffix='.tsv')
        os.close(fd_frac)
        try:
            with open(gs_path, 'w', encoding='utf-8') as handle:
                handle.write('@Version:0.9.1\n@SampleID:{}\n@@SEQUENCEID\tBINID\t_LENGTH\n'.format(sample_id))
                for seq, genome, length in gs_rows:
                    handle.write('{}\t{}\t{}\n'.format(seq, genome, length))
            write_prediction(pred_path, sample_id, pred_rows)
            write_fractional_truth(frac_path, sample_id, [
                {
                    'SEQUENCEID': seq,
                    '_LENGTH': length,
                    'COMPONENT_BP': length,
                    'COMPONENT_TYPE': 'unique',
                    'GENOME_IDS': genome_ids_json([genome]),
                }
                for seq, genome, length in gs_rows
            ])
            options = Options(skip_gs=True, skip_heatmap=True)
            gs_meta = load_data.read_metadata((gs_path, 'gs'))[0]
            pred_meta = load_data.read_metadata((pred_path, 'tool'))[0]
            gs_query = GenomeQuery('Gold standard', sample_id, options, gs_meta, True)
            gs_query.gold_standard = gs_query
            pred_query = GenomeQuery('tool', sample_id, options, pred_meta, False)
            pred_query.gold_standard = gs_query
            pred_query.compute_metrics()

            from cami_amber.fractional_truth import load_fractional_truth_file
            from cami_amber.fractional_metrics import load_prediction_assignments
            truth = load_fractional_truth_file(frac_path)[sample_id]
            assignments = load_prediction_assignments(pred_meta, truth)
            frac = compute_sample_metrics(truth, assignments, [0.5], [0.1])

            self.assertAlmostEqual(pred_query.metrics.percentage_of_assigned_bps, frac['metrics']['percentage_of_assigned_bps'])
            self.assertAlmostEqual(pred_query.metrics.percentage_of_assigned_seqs, frac['metrics']['percentage_of_assigned_seqs'])
            self.assertAlmostEqual(pred_query.metrics.precision_weighted_bp, frac['metrics']['precision_weighted_bp'])
            self.assertAlmostEqual(pred_query.metrics.precision_weighted_seq, frac['metrics']['precision_weighted_seq'])
            self.assertAlmostEqual(pred_query.metrics.recall_weighted_bp, frac['metrics']['recall_weighted_bp'])
            self.assertAlmostEqual(pred_query.metrics.accuracy_bp, frac['metrics']['accuracy_bp'])
            self.assertAlmostEqual(pred_query.metrics.accuracy_seq, frac['metrics']['accuracy_seq'])
            self.assertTrue(np.isnan(frac['metrics']['adjusted_rand_index_seq']))
            self.assertAlmostEqual(
                pred_query.metrics.adjusted_rand_index_bp,
                frac['metrics']['adjusted_rand_index_bp_identifiable'],
                places=12,
            )
            self.assertAlmostEqual(
                pred_query.metrics.rand_index_bp,
                frac['metrics']['rand_index_bp_identifiable'],
                places=12,
            )
        finally:
            os.remove(gs_path)
            os.remove(pred_path)
            os.remove(frac_path)


if __name__ == '__main__':
    unittest.main()
