# Copyright 2026 Department of Computational Biology for Infection Research - Helmholtz Centre for Infection Research
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

import json
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import FrozenSet, Iterable, List, Mapping, Optional, Sequence, Tuple

from cami_amber.utils import load_data
from cami_amber.version import FRACTIONAL_SCHEMA_VERSION, TRUTH_MODEL_FRACTIONAL

KIND_UNIQUE = 'unique'
KIND_COMPATIBLE = 'compatible'
KIND_UNRESOLVED = 'unresolved'
VALID_KINDS = (KIND_UNIQUE, KIND_COMPATIBLE, KIND_UNRESOLVED)
REQUIRED_COLUMNS = ('SEQUENCEID', '_LENGTH', 'COMPONENT_BP', 'COMPONENT_TYPE', 'GENOME_IDS')


class FractionalTruthError(ValueError):
    """Fail-closed fractional truth validation error."""


@dataclass(frozen=True)
class TruthComponent:
    bp: int
    kind: str
    genome_ids: FrozenSet[str]


@dataclass
class SequenceTruth:
    sequence_id: str
    length: int
    components: List[TruthComponent] = field(default_factory=list)

    @property
    def unique_bp(self) -> int:
        return sum(c.bp for c in self.components if c.kind == KIND_UNIQUE)

    @property
    def compatible_bp(self) -> int:
        return sum(c.bp for c in self.components if c.kind == KIND_COMPATIBLE)

    @property
    def unresolved_bp(self) -> int:
        return sum(c.bp for c in self.components if c.kind == KIND_UNRESOLVED)

    def unique_bp_for(self, genome_id: str) -> int:
        return sum(c.bp for c in self.components if c.kind == KIND_UNIQUE and genome_id in c.genome_ids)

    def compatible_bp_for(self, genome_id: str) -> int:
        return sum(c.bp for c in self.components if genome_id in c.genome_ids)

    def is_fari_eligible(self) -> bool:
        return self.compatible_bp == 0 and self.unresolved_bp == 0 and self.unique_bp == self.length


@dataclass
class FractionalTruthSample:
    sample_id: str
    sequences: OrderedDict[str, SequenceTruth]
    genomes: List[str]
    truth_model: str = TRUTH_MODEL_FRACTIONAL
    schema_version: str = FRACTIONAL_SCHEMA_VERSION
    component_count: int = 0
    original_sequence_ids: Optional[FrozenSet[str]] = None
    _genome_agg: Optional[dict] = None

    def sequence_ids(self) -> List[str]:
        return list(self.sequences.keys())

    def assembly_bp(self) -> int:
        return sum(seq.length for seq in self.sequences.values())

    def unique_truth_bp(self) -> int:
        return sum(seq.unique_bp for seq in self.sequences.values())

    def compatible_truth_bp(self) -> int:
        return sum(seq.compatible_bp for seq in self.sequences.values())

    def unresolved_truth_bp(self) -> int:
        return sum(seq.unresolved_bp for seq in self.sequences.values())

    def unique_truth_seq_units(self, genome_id: str) -> float:
        return float(self.genome_aggregates()['unique_seq'].get(genome_id, 0.0))

    def identifiable_truth_bp(self, genome_id: str) -> int:
        return int(self.genome_aggregates()['identifiable_bp'].get(genome_id, 0))

    def compatible_truth_bp_involving(self, genome_id: str) -> int:
        return int(self.genome_aggregates()['compatible_bp'].get(genome_id, 0))

    def genome_aggregates(self) -> dict:
        if self._genome_agg is not None:
            return self._genome_agg
        identifiable_bp = {}
        unique_seq = {}
        compatible_bp = {}
        for seq in self.sequences.values():
            inv = (1.0 / seq.length) if seq.length else 0.0
            for component in seq.components:
                if component.kind == KIND_UNIQUE:
                    genome_id = next(iter(component.genome_ids))
                    identifiable_bp[genome_id] = identifiable_bp.get(genome_id, 0) + component.bp
                    unique_seq[genome_id] = unique_seq.get(genome_id, 0.0) + component.bp * inv
                elif component.kind == KIND_COMPATIBLE:
                    for genome_id in component.genome_ids:
                        compatible_bp[genome_id] = compatible_bp.get(genome_id, 0) + component.bp
        self._genome_agg = {
            'identifiable_bp': identifiable_bp,
            'unique_seq': unique_seq,
            'compatible_bp': compatible_bp,
        }
        return self._genome_agg

    def identifiable_fraction(self, genome_id: str) -> float:
        unique_bp = self.identifiable_truth_bp(genome_id)
        compatible_bp = self.compatible_truth_bp_involving(genome_id)
        den = unique_bp + compatible_bp
        if den == 0:
            return float('nan')
        return unique_bp / den

    def filter_min_length(self, min_length: int) -> 'FractionalTruthSample':
        original = self.original_sequence_ids or frozenset(self.sequences)
        if not min_length:
            self.original_sequence_ids = original
            return self
        kept = OrderedDict((sid, seq) for sid, seq in self.sequences.items() if seq.length >= min_length)
        rebuilt = _rebuild_sample(self.sample_id, kept, self.truth_model, self.schema_version)
        rebuilt.original_sequence_ids = original
        return rebuilt

    def remove_genomes(self, genome_ids: Optional[Sequence[str]]) -> 'FractionalTruthSample':
        if not genome_ids:
            return self
        removed = set(genome_ids)
        kept: OrderedDict[str, SequenceTruth] = OrderedDict()
        for sid, seq in self.sequences.items():
            new_components: List[TruthComponent] = []
            for component in seq.components:
                remaining = frozenset(g for g in component.genome_ids if g not in removed)
                if component.kind == KIND_UNIQUE:
                    if remaining:
                        new_components.append(component)
                    else:
                        new_components.append(TruthComponent(bp=component.bp, kind=KIND_UNRESOLVED, genome_ids=frozenset()))
                elif component.kind == KIND_COMPATIBLE:
                    if len(remaining) >= 2:
                        new_components.append(TruthComponent(bp=component.bp, kind=KIND_COMPATIBLE, genome_ids=remaining))
                    elif len(remaining) == 1:
                        new_components.append(TruthComponent(bp=component.bp, kind=KIND_UNIQUE, genome_ids=remaining))
                    else:
                        new_components.append(TruthComponent(bp=component.bp, kind=KIND_UNRESOLVED, genome_ids=frozenset()))
                else:
                    new_components.append(component)
            kept[sid] = SequenceTruth(sequence_id=sid, length=seq.length, components=_merge_like_components(new_components))
        return _rebuild_sample(self.sample_id, kept, self.truth_model)

    def summary_row(self) -> dict:
        assembly = self.assembly_bp()
        unique_bp = self.unique_truth_bp()
        compatible_bp = self.compatible_truth_bp()
        unresolved_bp = self.unresolved_truth_bp()
        return {
            'sample_id': self.sample_id,
            'n_sequences': len(self.sequences),
            'assembly_bp': assembly,
            'n_genomes': len(self.genomes),
            'unique_bp': unique_bp,
            'compatible_bp': compatible_bp,
            'unresolved_bp': unresolved_bp,
            'unique_truth_fraction': _safe_div(unique_bp, assembly),
            'compatible_truth_fraction': _safe_div(compatible_bp, assembly),
            'unresolved_truth_fraction': _safe_div(unresolved_bp, assembly),
            'component_count': self.component_count,
            'validation_status': 'ok',
            'truth_model': self.truth_model,
        }


def _safe_div(num: float, den: float) -> float:
    if den == 0:
        return float('nan')
    return num / den


def _merge_like_components(components: Sequence[TruthComponent]) -> List[TruthComponent]:
    buckets: OrderedDict[Tuple[str, FrozenSet[str]], int] = OrderedDict()
    for component in components:
        key = (component.kind, component.genome_ids)
        buckets[key] = buckets.get(key, 0) + component.bp
    return [TruthComponent(bp=bp, kind=kind, genome_ids=gids) for (kind, gids), bp in buckets.items()]


def _rebuild_sample(sample_id: str, sequences: OrderedDict[str, SequenceTruth], truth_model: str, schema_version: str = FRACTIONAL_SCHEMA_VERSION) -> FractionalTruthSample:
    genomes: List[str] = []
    seen = set()
    component_count = 0
    for seq in sequences.values():
        component_count += len(seq.components)
        for component in seq.components:
            for genome_id in sorted(component.genome_ids):
                if genome_id not in seen:
                    seen.add(genome_id)
                    genomes.append(genome_id)
    return FractionalTruthSample(
        sample_id=sample_id,
        sequences=sequences,
        genomes=genomes,
        truth_model=truth_model,
        schema_version=schema_version,
        component_count=component_count,
    )


def _parse_genome_ids(raw: str, line_no: int) -> FrozenSet[str]:
    text = (raw or '').strip()
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError as exc:
        raise FractionalTruthError('Malformed GENOME_IDS JSON at data line {}: {}'.format(line_no, exc)) from exc
    if not isinstance(parsed, list):
        raise FractionalTruthError('GENOME_IDS must be a JSON array at data line {}'.format(line_no))
    genomes = []
    seen = set()
    for item in parsed:
        if not isinstance(item, str):
            raise FractionalTruthError('GENOME_IDS entries must be strings at data line {}'.format(line_no))
        if item.strip() == '' or any(ch in item for ch in ('\t', '\n', '\r')):
            raise FractionalTruthError('Invalid genome ID at data line {}'.format(line_no))
        if item in seen:
            raise FractionalTruthError('Duplicate genome ID in GENOME_IDS at data line {}'.format(line_no))
        seen.add(item)
        genomes.append(item)
    return frozenset(genomes)


def _validate_component(kind: str, genome_ids: FrozenSet[str], bp: int, line_no: int) -> None:
    if bp <= 0:
        raise FractionalTruthError('COMPONENT_BP must be a positive integer at data line {}'.format(line_no))
    if kind == KIND_UNIQUE and len(genome_ids) != 1:
        raise FractionalTruthError('unique components require exactly one genome ID at data line {}'.format(line_no))
    if kind == KIND_COMPATIBLE and len(genome_ids) < 2:
        raise FractionalTruthError('compatible components require >=2 genome IDs at data line {}'.format(line_no))
    if kind == KIND_UNRESOLVED and len(genome_ids) != 0:
        raise FractionalTruthError('unresolved components require empty GENOME_IDS at data line {}'.format(line_no))
    if kind not in VALID_KINDS:
        raise FractionalTruthError('Unknown COMPONENT_TYPE {!r} at data line {}'.format(kind, line_no))


def load_fractional_truth_file(path: str) -> OrderedDict[str, FractionalTruthSample]:
    samples_metadata = load_data.read_metadata((path, 'fractional gold standard'))
    samples: OrderedDict[str, FractionalTruthSample] = OrderedDict()
    for metadata in samples_metadata:
        sample = load_fractional_truth_sample(metadata)
        if sample.sample_id in samples:
            raise FractionalTruthError('Duplicate @SampleID {}'.format(sample.sample_id))
        samples[sample.sample_id] = sample
    if not samples:
        raise FractionalTruthError('No fractional truth samples found in {}'.format(path))
    return samples


def load_fractional_truth_sample(metadata) -> FractionalTruthSample:
    data_start, data_end, header, columns_list, file_path, _label = metadata
    sample_id = header.get('SAMPLEID')
    if not sample_id:
        raise FractionalTruthError('Fractional truth is missing @SampleID')
    truth_model = header.get('TRUTHMODEL') or header.get('TruthModel')
    if truth_model != TRUTH_MODEL_FRACTIONAL:
        raise FractionalTruthError(
            'Fractional truth @TruthModel must be exactly {} (got {!r})'.format(TRUTH_MODEL_FRACTIONAL, truth_model)
        )
    schema_version = header.get('VERSION')
    if schema_version != FRACTIONAL_SCHEMA_VERSION:
        raise FractionalTruthError(
            'Fractional truth @Version must be exactly {} (got {!r})'.format(FRACTIONAL_SCHEMA_VERSION, schema_version)
        )
    missing = [col for col in REQUIRED_COLUMNS if col not in columns_list]
    if missing:
        raise FractionalTruthError('Fractional truth missing columns: {}'.format(', '.join(missing)))
    col_index = {name: i for i, name in enumerate(columns_list)}

    sequences: OrderedDict[str, SequenceTruth] = OrderedDict()
    genomes: List[str] = []
    seen_genomes = set()
    component_count = 0
    with load_data.open_generic(file_path) as handle:
        for offset, line in enumerate(handle):
            line_no = offset
            if line_no < data_start or line_no > data_end:
                continue
            stripped = line.strip('\n')
            if not stripped or stripped.startswith('#'):
                continue
            fields = stripped.split('\t')
            if len(fields) < len(columns_list):
                raise FractionalTruthError('Truncated fractional truth row at file line {}'.format(line_no + 1))
            sequence_id = fields[col_index['SEQUENCEID']]
            if sequence_id.strip() == '':
                raise FractionalTruthError('Blank SEQUENCEID at file line {}'.format(line_no + 1))
            try:
                length = int(fields[col_index['_LENGTH']])
                bp = int(fields[col_index['COMPONENT_BP']])
            except ValueError as exc:
                raise FractionalTruthError('Non-integer length or COMPONENT_BP at file line {}'.format(line_no + 1)) from exc
            if length <= 0:
                raise FractionalTruthError('Sequence _LENGTH must be positive at file line {}'.format(line_no + 1))
            kind = fields[col_index['COMPONENT_TYPE']].strip()
            genome_ids = _parse_genome_ids(fields[col_index['GENOME_IDS']], line_no + 1)
            _validate_component(kind, genome_ids, bp, line_no + 1)
            if sequence_id not in sequences:
                sequences[sequence_id] = SequenceTruth(sequence_id=sequence_id, length=length, components=[])
            elif sequences[sequence_id].length != length:
                raise FractionalTruthError('Inconsistent _LENGTH for {} at file line {}'.format(sequence_id, line_no + 1))
            existing = None
            for component in sequences[sequence_id].components:
                if component.kind == kind and component.genome_ids == genome_ids:
                    existing = component
                    break
            if existing is not None:
                sequences[sequence_id].components = [
                    TruthComponent(bp=c.bp + bp, kind=c.kind, genome_ids=c.genome_ids) if c is existing else c
                    for c in sequences[sequence_id].components
                ]
            else:
                sequences[sequence_id].components.append(TruthComponent(bp=bp, kind=kind, genome_ids=genome_ids))
            component_count += 1
            for genome_id in sorted(genome_ids):
                if genome_id not in seen_genomes:
                    seen_genomes.add(genome_id)
                    genomes.append(genome_id)

    if not sequences:
        raise FractionalTruthError('Fractional truth sample {} has no component rows'.format(sample_id))
    for seq in sequences.values():
        total = sum(c.bp for c in seq.components)
        if total != seq.length:
            raise FractionalTruthError(
                'Component bp sum {} != _LENGTH {} for sequence {}'.format(total, seq.length, seq.sequence_id)
            )
    return FractionalTruthSample(
        sample_id=sample_id,
        sequences=sequences,
        genomes=genomes,
        truth_model=truth_model,
        schema_version=schema_version,
        component_count=component_count,
        original_sequence_ids=frozenset(sequences),
    )


def write_fractional_truth(path: str, sample_id: str, rows: Iterable[Mapping[str, object]]) -> None:
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('@Version:0.1.0\n')
        handle.write('@SampleID:{}\n'.format(sample_id))
        handle.write('@TruthModel:{}\n'.format(TRUTH_MODEL_FRACTIONAL))
        handle.write('@@SEQUENCEID\t_LENGTH\tCOMPONENT_BP\tCOMPONENT_TYPE\tGENOME_IDS\n')
        for row in rows:
            handle.write('{SEQUENCEID}\t{_LENGTH}\t{COMPONENT_BP}\t{COMPONENT_TYPE}\t{GENOME_IDS}\n'.format(**row))
