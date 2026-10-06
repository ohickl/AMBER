import json
import subprocess
import unittest

import pandas as pd
from bokeh.embed import file_html
from bokeh.resources import INLINE

from cami_amber.recovered_genomes import (
    RECOVERY_JS, count_recovered, create_recovery_controls, recovery_records,
)


class RecoveryControlsTest(unittest.TestCase):
    def setUp(self):
        self.bins = pd.DataFrame([
            ('s1', 'a', .91, .96), ('s1', 'a', .90, 1.0),
            ('s1', 'a', 1.0, .95), ('s1', 'b', float('nan'), 1.0),
            ('s2', 'a', .99, .99), ('s2', 'b', .8, .99),
        ], columns=['sample_id', 'Tool', 'recall_bp', 'precision_bp'])

    def test_strict_boundaries_and_undefined_scores(self):
        counts = count_recovered(recovery_records(self.bins), 90, 5, ['s1'], ['a', 'b'])
        self.assertEqual(counts, {'Sample': ['s1', 's1'], 'Tool': ['a', 'b'], 'Count': [1, 0]})

    def test_filter_combinations_and_empty_selection(self):
        records = recovery_records(self.bins)
        self.assertEqual(count_recovered(records, 70, 5, ['s2'], ['b'])['Count'], [1])
        self.assertEqual(count_recovered(records, 90, 5, [], ['a'])['Count'], [])

    def test_javascript_matches_backend_after_control_changes(self):
        records = recovery_records(self.bins)
        for samples, tools, completeness, contamination in [
            (['s1', 's2'], ['a', 'b'], 90, 5), (['s2'], ['b'], 70, 5),
            ([], ['a'], 90, 5), (['s1'], ['a'], 100, 0),
        ]:
            script = 'const records = ' + json.dumps(records) + ';\n'
            for name, value in [('samples', samples), ('tools', tools),
                                ('completeness', completeness), ('contamination', contamination)]:
                script += 'const ' + name + ' = {value: ' + json.dumps(value) + '};\n'
            script += 'const source = {data: {}, change: {emit() {}}};\n'
            script += RECOVERY_JS + '\nconsole.log(JSON.stringify(source.data));'
            actual = json.loads(subprocess.check_output(['node', '-e', script], text=True))
            self.assertEqual(actual, count_recovered(records, completeness, contamination, samples, tools))

    def test_standalone_html_contains_connected_controls(self):
        controls = create_recovery_controls(self.bins, [.5, .7, .9], [.1, .05])
        html = file_html(controls, INLINE, 'Recovery controls')
        self.assertIn('Completeness greater than (%)', html)
        self.assertIn('Contamination less than (%)', html)
        self.assertEqual(html.count('new Map()'), 1)


if __name__ == '__main__':
    unittest.main()
