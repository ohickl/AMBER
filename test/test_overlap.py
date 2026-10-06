import json
import subprocess
import unittest

from cami_amber.overlap import ALL_SAMPLES, OVERLAP_FUNCTION_JS, overlap_counts


class OverlapTest(unittest.TestCase):
    def setUp(self):
        self.records = [dict(sample=sample, genome=genome, tool=tool, completeness=.99, purity=.99)
                        for sample, genome, tool in [('s1','g','a'),('s1','g','a'),('s1','g','b'),
                                                     ('s1','unique','a'),('s2','g','a')]]

    def test_duplicate_bins_and_sample_identity(self):
        result = overlap_counts(self.records, ['a','b'], ALL_SAMPLES, 90, 5)
        self.assertEqual(result['tools'], dict(Tool=['a','b'], Recovered=[3,1], Unique=[2,0]))
        self.assertEqual(result['intersections']['Pattern'], ['10','11'])
        self.assertEqual(result['intersections']['Count'], [2,1])
        self.assertEqual(len(result['genomes']['Genome']), 3)

    def test_individual_sample_and_selected_tool_unique(self):
        result = overlap_counts(self.records, ['b'], 's1', 90, 5)
        self.assertEqual(result['tools']['Unique'], [1])
        self.assertEqual(overlap_counts(self.records, [], ALL_SAMPLES, 90, 5)['genomes']['Genome'], [])

    def test_javascript_matches_python_with_more_than_32_tools(self):
        for tools, sample in [(['a','b'], ALL_SAMPLES), (['a'], 's1'),
                              (['a','b'] + ['tool'+str(i) for i in range(40)], ALL_SAMPLES)]:
            script = OVERLAP_FUNCTION_JS + '\nconsole.log(JSON.stringify(overlapCounts(' + ','.join(
                json.dumps(value) for value in [self.records, tools, sample, 90, 5]) + ')));'
            actual = json.loads(subprocess.check_output(['node','-e',script], text=True))
            self.assertEqual(actual, overlap_counts(self.records, tools, sample, 90, 5))

    def test_strict_thresholds(self):
        self.records.append(dict(sample='s1',genome='boundary',tool='a',completeness=.90,purity=1.0))
        self.assertEqual(overlap_counts(self.records, ['a'], ALL_SAMPLES,90,5)['tools']['Recovered'], [3])


if __name__ == '__main__':
    unittest.main()
