#!/usr/bin/env python3
"""Evaluate a complete multisample comparison once and render one cohort report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import amber
from cami_amber.cohort_cache import cache_key, restore, seal
from cami_amber.cohort_report import create_cohort_report
from cami_amber.domain_metadata import read_observation_metadata
from cami_amber.utils.load_data import read_metadata


SCHEMA = 'amber-cohort-v1'


def sample_ids(path: Path) -> list[str]:
    samples = [row[2]['SAMPLEID'] for row in read_metadata((str(path), str(path)))]
    if len(samples) != len(set(samples)):
        raise ValueError('Duplicate SampleID sections: {}'.format(path))
    return samples


def load_manifest(path: Path) -> dict:
    manifest = json.loads(path.read_text())
    if manifest.get('schema') != SCHEMA:
        raise ValueError('Unsupported cohort manifest schema')
    def resolve(value):
        candidate = Path(value)
        return str((path.parent / candidate).resolve())
    manifest['domains'] = resolve(manifest['domains'])
    manifest['observations'] = resolve(manifest['observations'])
    observations = read_observation_metadata(manifest['observations'])
    if not manifest.get('models') or set(manifest['models']) != {'official', 'fractional'}:
        raise ValueError('A cohort requires official and fractional scoring models')
    labels = None
    for model, definition in manifest['models'].items():
        definition['truth'] = resolve(definition['truth'])
        if set(sample_ids(Path(definition['truth']))) != set(observations):
            raise ValueError('Truth and observation sample coverage differ: ' + model)
        predictions = definition['predictions']
        current = [row['label'] for row in predictions]
        if not current or len(current) != len(set(current)) or any(not isinstance(label, str) or not label or ',' in label for label in current):
            raise ValueError('Method labels must be nonempty, unique and comma-free')
        if labels is not None and current != labels:
            raise ValueError('Both models must use the same ordered methods')
        labels = current
        for row in predictions:
            row['path'] = resolve(row['path'])
            if set(sample_ids(Path(row['path']))) != set(observations):
                raise ValueError('Method does not cover every observation: ' + row['label'])
    manifest['labels'] = labels
    return manifest


def run(manifest_path: Path, output: Path, model_cache: Path) -> dict:
    manifest = load_manifest(manifest_path)
    output.mkdir(parents=True, exist_ok=True)
    models, receipt = {}, dict(schema=SCHEMA, models={})
    for model in ('official', 'fractional'):
        definition = manifest['models'][model]
        minimum_length = definition.get('min_length', 1500 if model == 'fractional' else None)
        if minimum_length is not None and (not isinstance(minimum_length, int) or minimum_length < 0):
            raise ValueError('Minimum length must be a nonnegative integer')
        inputs = dict(truth=definition['truth'], domains=manifest['domains'], observations=manifest['observations'])
        inputs.update({'prediction/' + row['label']: row['path'] for row in definition['predictions']})
        key = cache_key(inputs, dict(model=model, labels=manifest['labels'], min_completeness='90', max_contamination='5', min_length=minimum_length))
        result = restore(model_cache, key, output / model)
        restored = result is not None
        if result is None:
            truth_flag = '-g' if model == 'official' else '--fractional-gold-standard'
            args = [truth_flag, definition['truth'], '-o', str(output / model), '--skip_gs',
                    '--genome-domains', manifest['domains'], '--sample-metadata', manifest['observations'],
                    '-x', '90', '-y', '5', '-l', ','.join(manifest['labels'])]
            if minimum_length is not None:
                args.extend(['--min_length', str(minimum_length)])
            args.extend(row['path'] for row in definition['predictions'])
            result = amber.main(args, render_html=False)
            seal(model_cache, key, result)
        models[model] = result
        receipt['models'][model] = dict(cache_key=key, restored=restored, observations=len(sample_ids(Path(definition['truth']))),
                                        methods=len(manifest['labels']))
    create_cohort_report(models, output)
    (output / 'cohort.receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--model-cache', type=Path, required=True)
    args = parser.parse_args()
    run(args.manifest, args.output, args.model_cache)


if __name__ == '__main__':
    main()
