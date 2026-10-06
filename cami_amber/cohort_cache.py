"""Sealed, content-addressed model results survive interrupted cohort reports."""

from __future__ import annotations

import hashlib
import json
import shutil
import uuid
from pathlib import Path
from types import SimpleNamespace

import pandas as pd


SCHEMA = 'amber-cohort-cache-v1'


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def source_digest() -> str:
    from importlib import metadata, util
    root = Path(__file__).resolve().parent
    records = [(str(path.relative_to(root)), sha256(path)) for path in sorted(root.rglob('*.py'))]
    # Installed command scripts live in bin/, while the package lives in
    # site-packages. Never assume the checkout's requirements.txt is present.
    for entrypoint in ('amber', 'amber_cohort'):
        spec = util.find_spec(entrypoint)
        if spec is None or spec.origin is None:
            raise ValueError('Cannot locate AMBER entrypoint: ' + entrypoint)
        records.append(('entrypoint/' + entrypoint, sha256(Path(spec.origin))))
    dependencies = {name: metadata.version(name) for name in
                    ('numpy', 'pandas', 'bokeh', 'matplotlib', 'seaborn', 'pyarrow', 'distinctipy')}
    payload = dict(source=records, dependencies=dependencies)
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def cache_key(inputs: dict, parameters: dict) -> str:
    payload = dict(schema=SCHEMA, source=source_digest(), parameters=parameters,
                   inputs={label: sha256(Path(path)) for label, path in sorted(inputs.items())})
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def restore(cache_root: Path, key: str, output: Path):
    leaf = cache_root / key
    if not leaf.exists():
        return None
    receipt = json.loads((leaf / 'receipt.json').read_text())
    if receipt.get('schema') != SCHEMA or receipt.get('key') != key:
        raise ValueError('Invalid sealed model cache identity: {}'.format(leaf))
    for relative, digest in receipt['files'].items():
        if sha256(leaf / relative) != digest:
            raise ValueError('Model cache digest mismatch: {}'.format(leaf / relative))
    state = json.loads((leaf / 'state.json').read_text())
    output.mkdir(parents=True, exist_ok=True)
    shutil.copytree(leaf / 'native', output, dirs_exist_ok=True)
    options = SimpleNamespace(output_dir=str(output), min_completeness=state['min_completeness'],
                              max_contamination=state['max_contamination'])
    options.domain_profiles = {int(mask): (pd.read_parquet(leaf / f'profiles/{mask}.summary.parquet'),
                                          pd.read_parquet(leaf / f'profiles/{mask}.bins.parquet'))
                               for mask in state['masks']}
    return dict(summary=pd.read_parquet(leaf / 'summary.parquet'),
                bins=pd.read_parquet(leaf / 'bins.parquet'), labels=state['labels'],
                samples=state['samples'], options=options)


def seal(cache_root: Path, key: str, result: dict) -> None:
    cache_root.mkdir(parents=True, exist_ok=True)
    if (cache_root / key).exists():
        raise ValueError('Refusing to overwrite an existing sealed model cache')
    pending = cache_root / ('.pending-' + uuid.uuid4().hex)
    pending.mkdir()
    shutil.copytree(result['options'].output_dir, pending / 'native')
    result['summary'].to_parquet(pending / 'summary.parquet', index=False)
    result['bins'].to_parquet(pending / 'bins.parquet', index=False)
    (pending / 'profiles').mkdir()
    options = result['options']
    for mask, (summary, bins) in options.domain_profiles.items():
        summary.to_parquet(pending / f'profiles/{mask}.summary.parquet', index=False)
        bins.to_parquet(pending / f'profiles/{mask}.bins.parquet', index=False)
    state = dict(labels=result['labels'], samples=result['samples'], masks=list(options.domain_profiles),
                 min_completeness=options.min_completeness, max_contamination=options.max_contamination)
    (pending / 'state.json').write_text(json.dumps(state, sort_keys=True) + '\n')
    files = {str(path.relative_to(pending)): sha256(path) for path in sorted(pending.rglob('*')) if path.is_file()}
    receipt = dict(schema=SCHEMA, key=key, files=files)
    (pending / 'receipt.json').write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    pending.rename(cache_root / key)
