import os
import tempfile
import unittest

from cami_amber.fractional_truth import FractionalTruthError, load_fractional_truth_file, write_fractional_truth
from test.fractional.helpers import genome_ids_json


class TestFractionalTruthParser(unittest.TestCase):
    def _write(self, rows, sample_id='s1', extra_header=None):
        fd, path = tempfile.mkstemp(suffix='.tsv')
        os.close(fd)
        with open(path, 'w', encoding='utf-8') as handle:
            handle.write('@Version:0.1.0\n')
            handle.write('@SampleID:{}\n'.format(sample_id))
            if extra_header:
                handle.write(extra_header)
            else:
                handle.write('@TruthModel:fractional-origin-v1\n')
            handle.write('@@SEQUENCEID\t_LENGTH\tCOMPONENT_BP\tCOMPONENT_TYPE\tGENOME_IDS\n')
            for row in rows:
                handle.write(row + '\n')
        return path

    def test_valid_mixed(self):
        path = self._write([
            'c1\t100\t70\tunique\t{}'.format(genome_ids_json(['A'])),
            'c1\t100\t25\tunique\t{}'.format(genome_ids_json(['B'])),
            'c1\t100\t5\tcompatible\t{}'.format(genome_ids_json(['A', 'B'])),
        ])
        try:
            sample = load_fractional_truth_file(path)['s1']
            self.assertEqual(sample.assembly_bp(), 100)
            self.assertEqual(sample.unique_truth_bp(), 95)
            self.assertEqual(sample.compatible_truth_bp(), 5)
        finally:
            os.remove(path)

    def test_wrong_truth_model(self):
        path = self._write(['c1\t10\t10\tunique\t["A"]'], extra_header='@TruthModel:nope\n')
        try:
            with self.assertRaises(FractionalTruthError):
                load_fractional_truth_file(path)
        finally:
            os.remove(path)

    def test_sum_mismatch(self):
        path = self._write(['c1\t10\t9\tunique\t["A"]'])
        try:
            with self.assertRaises(FractionalTruthError):
                load_fractional_truth_file(path)
        finally:
            os.remove(path)

    def test_duplicate_same_genome_components_are_merged(self):
        path = self._write([
            'c1\t100\t40\tunique\t["A"]',
            'c1\t100\t30\tunique\t["A"]',
            'c1\t100\t30\tunique\t["B"]',
        ])
        try:
            sample = load_fractional_truth_file(path)['s1']
            unique_a = sample.sequences['c1'].unique_bp_for('A')
            self.assertEqual(unique_a, 70)
        finally:
            os.remove(path)

    def test_wrong_schema_version(self):
        fd, path = tempfile.mkstemp(suffix='.tsv')
        os.close(fd)
        with open(path, 'w', encoding='utf-8') as handle:
            handle.write('@Version:9.9.9\n@SampleID:s1\n@TruthModel:fractional-origin-v1\n')
            handle.write('@@SEQUENCEID\t_LENGTH\tCOMPONENT_BP\tCOMPONENT_TYPE\tGENOME_IDS\n')
            handle.write('c1\t10\t10\tunique\t["A"]\n')
        try:
            with self.assertRaises(FractionalTruthError):
                load_fractional_truth_file(path)
        finally:
            os.remove(path)

    def test_unique_cardinality(self):
        path = self._write(['c1\t10\t10\tunique\t["A","B"]'])
        try:
            with self.assertRaises(FractionalTruthError):
                load_fractional_truth_file(path)
        finally:
            os.remove(path)

    def test_compatible_cardinality(self):
        path = self._write(['c1\t10\t10\tcompatible\t["A"]'])
        try:
            with self.assertRaises(FractionalTruthError):
                load_fractional_truth_file(path)
        finally:
            os.remove(path)

    def test_unresolved_must_be_empty(self):
        path = self._write(['c1\t10\t10\tunresolved\t["A"]'])
        try:
            with self.assertRaises(FractionalTruthError):
                load_fractional_truth_file(path)
        finally:
            os.remove(path)

    def test_malformed_json(self):
        path = self._write(['c1\t10\t10\tunique\t[A]'])
        try:
            with self.assertRaises(FractionalTruthError):
                load_fractional_truth_file(path)
        finally:
            os.remove(path)

    def test_zero_bp(self):
        path = self._write(['c1\t10\t0\tunique\t["A"]'])
        try:
            with self.assertRaises(FractionalTruthError):
                load_fractional_truth_file(path)
        finally:
            os.remove(path)

    def test_inconsistent_length(self):
        path = self._write([
            'c1\t10\t5\tunique\t["A"]',
            'c1\t11\t5\tunique\t["B"]',
        ])
        try:
            with self.assertRaises(FractionalTruthError):
                load_fractional_truth_file(path)
        finally:
            os.remove(path)

    def test_remove_genomes_compatible_collapses(self):
        path = self._write([
            'c1\t100\t70\tunique\t["A"]',
            'c1\t100\t30\tcompatible\t["A","B"]',
        ])
        try:
            sample = load_fractional_truth_file(path)['s1'].remove_genomes(['A'])
            kinds = {(c.kind, tuple(sorted(c.genome_ids)), c.bp) for c in sample.sequences['c1'].components}
            self.assertIn(('unresolved', (), 70), kinds)
            self.assertIn(('unique', ('B',), 30), kinds)
        finally:
            os.remove(path)


if __name__ == '__main__':
    unittest.main()
