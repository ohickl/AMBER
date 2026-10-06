"""Standalone entity decoding must match Bokeh rather than display HTML text."""
import json
import unittest
from refresh_cohort_report import decode_document_json, packed


class DocumentDecodeTests(unittest.TestCase):
    def test_decode_table_and_callback_exactly_once(self):
        original = {'table':'<table><tr><td>literal &lt; &amp; &#123;</td></tr></table>',
                    'callback':'if (a > 1 && b < 2) { value = "ok"; }'}
        escaped = json.dumps(original).replace('&','&amp;').replace('<','&lt;').replace('>','&gt;').replace("'",'&#x27;')
        self.assertEqual(decode_document_json(escaped),original)

    def test_serialized_maps_use_traversable_json_lists(self):
        value = packed({'callback': {'id': 'action'}})
        self.assertIsInstance(value['entries'][0], list)
        self.assertEqual(value, json.loads(json.dumps(value)))

    def test_unescaped_json_item_stays_plain_except_encoded_transport_entities(self):
        original = {'value': 'scores unchanged', 'score': .75}
        self.assertEqual(decode_document_json(json.dumps(original)),original)


if __name__ == '__main__':
    unittest.main()
