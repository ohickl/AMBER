"""Complete sample coverage, valid empty panels and sealed model restoration."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import amber_cohort
from cami_amber.utils import load_data


class CohortTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        blocks = [('s0__control', 'c1', 'A'), ('s0__metacarvel', 'c2', 'A')]
        gold, fractional, predicted = [], [], []
        for sample, seq, genome in blocks:
            header = f'@Version:0.10.0\n@SampleID:{sample}\n'
            gold.append(header + f'@@SEQUENCEID\tBINID\t_LENGTH\n{seq}\t{genome}\t100\n')
            fractional.append('@Version:0.1.0\n@SampleID:' + sample + '\n@TruthModel:fractional-origin-v1\n@SchemaVersion:0.1.0\n' +
                              '@@SEQUENCEID\t_LENGTH\tCOMPONENT_BP\tCOMPONENT_TYPE\tGENOME_IDS\n' +
                              f'{seq}\t100\t100\tunique\t["A"]\n')
            predicted.append(header + f'@@SEQUENCEID\tBINID\n{seq}\tbinA\n')
        (self.root / 'truth.tsv').write_text(''.join(gold))
        (self.root / 'fractional.tsv').write_text(''.join(fractional))
        (self.root / 'pred.tsv').write_text(''.join(predicted))
        (self.root / 'empty.tsv').write_text(''.join(f'@Version:0.10.0\n@SampleID:{sample}\n@@SEQUENCEID\tBINID\n' for sample, _, _ in blocks))
        (self.root / 'domains.tsv').write_text('SampleID\tGenomeID\tDomain\ns0\tA\tBacteria\n')
        (self.root / 'observations.tsv').write_text('SampleID\tBiologicalSampleID\tAssemblyVariant\n' +
                                                 's0__control\ts0\tcontrol\ns0__metacarvel\ts0\tmetacarvel\n')
        manifest = dict(schema='amber-cohort-v1', domains='domains.tsv', observations='observations.tsv', models={})
        for model, truth in [('official', 'truth.tsv'), ('fractional', 'fractional.tsv')]:
            manifest['models'][model] = dict(truth=truth, min_length=0 if model == 'fractional' else None, predictions=[dict(label='one', path='pred.tsv'), dict(label='empty', path='empty.tsv')])
        self.manifest = self.root / 'manifest.json'
        self.manifest.write_text(json.dumps(manifest))

    def test_empty_sections_and_coverage(self):
        sections = load_data.read_metadata((str(self.root / 'empty.tsv'), 'empty'))
        self.assertEqual(len(sections), 2)
        for section in sections:
            self.assertTrue(load_data.load_sample(section).empty)
        amber_cohort.load_manifest(self.manifest)
        (self.root / 'empty.tsv').write_text('@SampleID:s0__control\n@@SEQUENCEID\tBINID\n')
        with self.assertRaisesRegex(ValueError, 'every observation'):
            amber_cohort.load_manifest(self.manifest)

    def test_cached_scores_and_corruption_detection(self):
        cache = self.root / 'cache'
        with patch('amber_cohort.create_cohort_report'):
            first = amber_cohort.run(self.manifest, self.root / 'first', cache)
            self.assertTrue(all(not row['restored'] for row in first['models'].values()))
            with patch('amber_cohort.amber.main', side_effect=AssertionError('scores recomputed')):
                second = amber_cohort.run(self.manifest, self.root / 'second', cache)
            self.assertTrue(all(row['restored'] for row in second['models'].values()))
            for model in ('official', 'fractional'):
                a = pd.read_csv(self.root / 'first' / model / 'results.tsv', sep='\t')
                b = pd.read_csv(self.root / 'second' / model / 'results.tsv', sep='\t')
                pd.testing.assert_frame_equal(a, b)
            key = first['models']['official']['cache_key']
            (cache / key / 'summary.parquet').write_bytes(b'corrupted')
            with self.assertRaisesRegex(ValueError, 'digest mismatch'):
                amber_cohort.run(self.manifest, self.root / 'third', cache)


if __name__ == '__main__':
    unittest.main()
