from __future__ import annotations

import json
import os
import tempfile
from typing import Iterable, List, Mapping, Sequence

from cami_amber.fractional_truth import write_fractional_truth
from cami_amber.utils import load_data


def genome_ids_json(ids: Sequence[str]) -> str:
    return json.dumps(list(ids))


def write_prediction(path: str, sample_id: str, rows: Iterable[Mapping[str, str]]) -> None:
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('@Version:0.10.0\n')
        handle.write('@SampleID:{}\n'.format(sample_id))
        handle.write('@@SEQUENCEID\tBINID\n')
        for row in rows:
            handle.write('{SEQUENCEID}\t{BINID}\n'.format(**row))


def load_truth_from_rows(sample_id: str, rows: List[dict]):
    fd, path = tempfile.mkstemp(suffix='.fractional.tsv')
    os.close(fd)
    try:
        write_fractional_truth(path, sample_id, rows)
        samples = __import__('cami_amber.fractional_truth', fromlist=['load_fractional_truth_file']).load_fractional_truth_file(path)
        return samples[sample_id]
    finally:
        os.remove(path)


def load_assignments(sample_id: str, pred_rows, truth):
    fd, path = tempfile.mkstemp(suffix='.binning')
    os.close(fd)
    try:
        write_prediction(path, sample_id, pred_rows)
        metadata = load_data.read_metadata((path, 'pred'))[0]
        from cami_amber.fractional_metrics import load_prediction_assignments
        return load_prediction_assignments(metadata, truth)
    finally:
        os.remove(path)
