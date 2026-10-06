"""Domain controls must retain foreign contamination and missed truth genomes."""
from pathlib import Path
import tempfile
import unittest

import amber


class OfficialDomainsTests(unittest.TestCase):
    def test_full_bin_contamination_and_missing_recall(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'truth.tsv').write_text('@SampleID:s0\n@@SEQUENCEID\tBINID\t_LENGTH\nx\tB\t70\ny\tE\t30\nz\tM\t100\n')
            (root / 'pred.tsv').write_text('@SampleID:s0\n@@SEQUENCEID\tBINID\nx\tmixed\ny\tmixed\n')
            (root / 'domains.tsv').write_text('SampleID\tGenomeID\tDomain\ns0\tB\tBacteria\ns0\tE\tEukaryotes\ns0\tM\tBacteria\n')
            result = amber.main(['-g', str(root / 'truth.tsv'), '-o', str(root / 'out'), '--skip_gs',
                                 '--genome-domains', str(root / 'domains.tsv'), '-l', 'tool', str(root / 'pred.tsv')], render_html=False)
            summary, bins = result['options'].domain_profiles[2]
            self.assertEqual(len(bins), 1)
            self.assertAlmostEqual(bins['precision_bp'].iloc[0], .7)
            self.assertAlmostEqual(summary['recall_avg_bp'].iloc[0], .5)
            self.assertAlmostEqual(summary['recall_weighted_bp'].iloc[0], 70 / 170)
            self.assertAlmostEqual(summary['percentage_of_assigned_bps'].iloc[0], 70 / 170)
            selected_euk_summary, selected_euk_bins = result['options'].domain_profiles[8]
            self.assertTrue(selected_euk_bins.empty)  # Original mixed-bin matching stays B.
            self.assertEqual(selected_euk_summary['recall_avg_bp'].iloc[0], 1)  # AMBER best-bin recall is independent of bin matching.
            full, _ = result['options'].domain_profiles[10]
            for name in ('precision_avg_bp', 'recall_avg_bp', 'recall_weighted_bp', 'percentage_of_assigned_bps'):
                self.assertEqual(full[name].iloc[0], result['summary'][name].iloc[0])


if __name__ == '__main__':
    unittest.main()
