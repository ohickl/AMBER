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

from collections import defaultdict
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from cami_amber.fractional_truth import KIND_COMPATIBLE, KIND_UNIQUE, KIND_UNRESOLVED, FractionalTruthSample
from cami_amber.version import FARI_IMPLEMENTATION_VERSION, SCORING_MODEL_VERSION, TRUTH_MODEL_FRACTIONAL

MATCH_RESOLVED = 'resolved'
MATCH_COMPATIBLE_TIE = 'compatible_tie'
MATCH_UNMATCHED = 'unmatched'


class FractionalPredictionError(ValueError):
    """Fatal fractional-mode prediction contract error."""


def _nan() -> float:
    return float('nan')


def _finite_or_nan(value: float) -> float:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return _nan()
    return float(value)


def _div(num: float, den: float) -> float:
    if den == 0:
        return _nan()
    return _finite_or_nan(num / den)


def load_prediction_assignments(metadata, truth: FractionalTruthSample) -> pd.DataFrame:
    df = pd.read_csv(
        metadata[4],
        sep='\t',
        comment='#',
        skiprows=metadata[0],
        nrows=metadata[1] - metadata[0] + 1,
        header=None,
        names=metadata[3],
        dtype=str,
    )
    if 'SEQUENCEID' not in df.columns or 'BINID' not in df.columns:
        raise FractionalPredictionError('Prediction file must contain SEQUENCEID and BINID')
    df = df[['SEQUENCEID', 'BINID']].copy()
    df['SEQUENCEID'] = df['SEQUENCEID'].astype(str)
    df['BINID'] = df['BINID'].astype(str)
    if df['SEQUENCEID'].duplicated().any():
        dupes = df.loc[df['SEQUENCEID'].duplicated(), 'SEQUENCEID'].unique().tolist()
        raise FractionalPredictionError('Duplicate prediction assignment for sequences: {}'.format(dupes[:10]))
    blank = df['BINID'].isna() | (df['BINID'].str.strip() == '') | (df['BINID'] == 'nan')
    if blank.any():
        raise FractionalPredictionError('Blank prediction BINID is fatal in fractional mode')
    unknown = set(df['SEQUENCEID']) - set(truth.sequences)
    if unknown:
        raise FractionalPredictionError(
            'Prediction sequences outside the truth universe: {}'.format(sorted(unknown)[:10])
        )
    return df


def build_support_matrices(
    truth: FractionalTruthSample,
    assignments: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    unique_rows = []
    compatible_rows = []
    bin_rows = []
    for _, row in assignments.iterrows():
        seq = truth.sequences[row['SEQUENCEID']]
        bin_id = row['BINID']
        unresolved = seq.unresolved_bp
        compatible_physical = seq.compatible_bp
        bin_rows.append({
            'BINID': bin_id,
            'SEQUENCEID': seq.sequence_id,
            'length': seq.length,
            'unresolved_bp': unresolved,
            'compatible_physical_bp': compatible_physical,
        })
        for component in seq.components:
            if component.kind == KIND_UNIQUE:
                genome_id = next(iter(component.genome_ids))
                unique_rows.append({'BINID': bin_id, 'genome_id': genome_id, 'bp': component.bp, 'SEQUENCEID': seq.sequence_id})
                compatible_rows.append({'BINID': bin_id, 'genome_id': genome_id, 'bp': component.bp, 'SEQUENCEID': seq.sequence_id})
            elif component.kind == KIND_COMPATIBLE:
                for genome_id in component.genome_ids:
                    compatible_rows.append({'BINID': bin_id, 'genome_id': genome_id, 'bp': component.bp, 'SEQUENCEID': seq.sequence_id})
    unique_df = pd.DataFrame(unique_rows)
    compatible_df = pd.DataFrame(compatible_rows)
    bin_seq_df = pd.DataFrame(bin_rows)

    def _aggregate(df, bp_name, seq_name):
        if df.empty:
            return pd.DataFrame(columns=['BINID', 'genome_id', bp_name, seq_name])
        rows = []
        for (bin_id, genome_id), grp in df.groupby(['BINID', 'genome_id']):
            rows.append({
                'BINID': bin_id,
                'genome_id': genome_id,
                bp_name: int(grp['bp'].sum()),
                seq_name: sum(r.bp / truth.sequences[r.SEQUENCEID].length for r in grp.itertuples()),
            })
        return pd.DataFrame(rows)

    unique_support = _aggregate(unique_df, 'unique_support_bp', 'unique_support_seq')
    compatible_support = _aggregate(compatible_df, 'compatible_support_bp', 'compatible_support_seq')
    return unique_support, compatible_support, bin_seq_df


def match_bins(unique_support: pd.DataFrame, compatible_support: pd.DataFrame, bin_ids: Sequence[str]) -> pd.DataFrame:
    records = []
    for bin_id in bin_ids:
        cand = compatible_support[compatible_support['BINID'] == bin_id]
        if cand.empty:
            records.append({
                'BINID': bin_id,
                'match_status': MATCH_UNMATCHED,
                'matched_genome_id': pd.NA,
                'matched_genome_ids': '',
                'correct_bp': 0,
                'correct_seq_units': 0.0,
            })
            continue
        max_bp = cand['compatible_support_bp'].max()
        tied = cand[cand['compatible_support_bp'] == max_bp]
        if len(tied) > 1:
            uniq = unique_support[(unique_support['BINID'] == bin_id) & (unique_support['genome_id'].isin(tied['genome_id']))]
            if not uniq.empty:
                max_u = uniq['unique_support_bp'].max()
                uniq_tied = uniq[uniq['unique_support_bp'] == max_u]
                if len(uniq_tied) == 1:
                    genome_id = uniq_tied.iloc[0]['genome_id']
                    seq_units = float(tied.loc[tied['genome_id'] == genome_id, 'compatible_support_seq'].iloc[0])
                    records.append({
                        'BINID': bin_id,
                        'match_status': MATCH_RESOLVED,
                        'matched_genome_id': genome_id,
                        'matched_genome_ids': genome_id,
                        'correct_bp': int(max_bp),
                        'correct_seq_units': seq_units,
                    })
                    continue
            genomes = sorted(tied['genome_id'].tolist())
            seq_units = float(tied['compatible_support_seq'].max())
            records.append({
                'BINID': bin_id,
                'match_status': MATCH_COMPATIBLE_TIE,
                'matched_genome_id': pd.NA,
                'matched_genome_ids': ','.join(genomes),
                'correct_bp': int(max_bp),
                'correct_seq_units': seq_units,
            })
            continue
        genome_id = tied.iloc[0]['genome_id']
        records.append({
            'BINID': bin_id,
            'match_status': MATCH_RESOLVED,
            'matched_genome_id': genome_id,
            'matched_genome_ids': genome_id,
            'correct_bp': int(max_bp),
            'correct_seq_units': float(tied.iloc[0]['compatible_support_seq']),
        })
    return pd.DataFrame(records)


def compute_fari(truth_membership: np.ndarray, pred_membership: np.ndarray, weights: Optional[np.ndarray] = None) -> float:
    """Frobenius ARI reducing to Hubert-Arabie ARI on hard partitions.

    Andrews et al., Journal of Classification 2022, DOI 10.1007/s00357-021-09407-3.
    Implementation version: {}.
    """.format(FARI_IMPLEMENTATION_VERSION)
    if truth_membership.size == 0:
        return _nan()
    if weights is None:
        weights = np.ones(truth_membership.shape[0], dtype=float)
    n = float(weights.sum())
    if n <= 1:
        return 1.0
    w = weights[:, None]
    contingency = truth_membership.T @ (pred_membership * w)
    truth_gram = truth_membership.T @ (truth_membership * w)
    pred_gram = pred_membership.T @ (pred_membership * w)
    pair_agree = 0.5 * (float(np.square(contingency).sum()) - n)
    pred_pairs = 0.5 * (float(np.square(pred_gram).sum()) - n)
    truth_pairs = 0.5 * (float(np.square(truth_gram).sum()) - n)
    total_pairs = 0.5 * (n * n - n)
    if total_pairs == 0:
        return 1.0
    expected = pred_pairs * (truth_pairs / total_pairs)
    denominator = ((pred_pairs + truth_pairs) / 2.0) - expected
    if denominator == 0:
        return 1.0
    return _finite_or_nan((pair_agree - expected) / denominator)


def _fari_memberships(truth: FractionalTruthSample, assignments: pd.DataFrame, weighted: bool) -> Tuple[float, int, float, float]:
    assigned_map = dict(zip(assignments['SEQUENCEID'], assignments['BINID']))
    eligible = []
    for sid, seq in truth.sequences.items():
        if seq.is_fari_eligible() and seq.sequence_id in assigned_map:
            eligible.append(seq)
    n_all = len(truth.sequences)
    bp_all = truth.assembly_bp()
    if not eligible:
        return _nan(), 0, 0.0, 0.0
    genomes = sorted({g for seq in eligible for c in seq.components for g in c.genome_ids})
    bins = sorted(set(assigned_map.get(seq.sequence_id, '__unbinned__') for seq in eligible))
    bin_index = {b: i for i, b in enumerate(bins)}
    genome_index = {g: i for i, g in enumerate(genomes)}
    n = len(eligible)
    k = len(genomes)
    l = len(bins)
    U = np.zeros((n, k), dtype=float)
    V = np.zeros((n, l), dtype=float)
    weights = np.zeros(n, dtype=float)
    for i, seq in enumerate(eligible):
        for component in seq.components:
            if component.kind == KIND_UNIQUE:
                U[i, genome_index[next(iter(component.genome_ids))]] = component.bp / seq.length
        V[i, bin_index[assigned_map.get(seq.sequence_id, '__unbinned__')]] = 1.0
        weights[i] = seq.length if weighted else 1.0
    value = compute_fari(U, V, weights=weights)
    seq_frac = n / n_all if n_all else _nan()
    bp_frac = sum(seq.length for seq in eligible) / bp_all if bp_all else _nan()
    return value, n, seq_frac, bp_frac


def compute_sample_metrics(
    truth: FractionalTruthSample,
    assignments: pd.DataFrame,
    min_completeness: Sequence[float],
    max_contamination: Sequence[float],
    filter_tail_percentage: float = 0.0,
) -> dict:
    unique_support, compatible_support, bin_seq_df = build_support_matrices(truth, assignments)
    if bin_seq_df.empty:
        bin_ids: List[str] = []
    else:
        bin_ids = sorted(bin_seq_df['BINID'].unique())
    matches = match_bins(unique_support, compatible_support, bin_ids)

    bin_totals = bin_seq_df.groupby('BINID', as_index=False).agg(
        total_length=('length', 'sum'),
        total_seq_counts=('SEQUENCEID', 'count'),
        unresolved_bp=('unresolved_bp', 'sum'),
        compatible_bp=('compatible_physical_bp', 'sum'),
    ) if not bin_seq_df.empty else pd.DataFrame(columns=['BINID', 'total_length', 'total_seq_counts', 'unresolved_bp', 'compatible_bp'])

    unique_lookup = {}
    if not unique_support.empty:
        for row in unique_support.itertuples(index=False):
            unique_lookup[(row.BINID, row.genome_id)] = (row.unique_support_bp, row.unique_support_seq)

    precision_rows = []
    for row in matches.itertuples(index=False):
        totals = bin_totals[bin_totals['BINID'] == row.BINID].iloc[0]
        total_bp = int(totals['total_length'])
        total_seq = int(totals['total_seq_counts'])
        unresolved_bp = int(totals['unresolved_bp'])
        compatible_physical = int(totals['compatible_bp'])
        genome_id = row.matched_genome_id if pd.notna(row.matched_genome_id) else None
        unique_tp_bp, unique_tp_seq = unique_lookup.get((row.BINID, genome_id), (0, 0.0)) if genome_id else (0, 0.0)
        identifiable = truth.identifiable_truth_bp(genome_id) if genome_id else 0
        unique_seq_den = truth.unique_truth_seq_units(genome_id) if genome_id else 0.0
        precision_bp = _div(row.correct_bp, total_bp)
        resolved_den = total_bp - unresolved_bp
        precision_bp_resolved = _div(row.correct_bp, resolved_den)
        precision_seq = _div(row.correct_seq_units, total_seq)
        if row.match_status != MATCH_RESOLVED:
            recall_bp = _nan()
            recall_seq = _nan()
        else:
            recall_bp = _div(unique_tp_bp, identifiable)
            recall_seq = _div(unique_tp_seq, unique_seq_den)
        foreign_bp = total_bp - int(row.correct_bp)
        precision_rows.append({
            'BINID': row.BINID,
            'match_status': row.match_status,
            'matched_genome_id': row.matched_genome_id,
            'matched_genome_ids': row.matched_genome_ids,
            'genome_id': row.matched_genome_id,
            'total_length': total_bp,
            'total_seq_counts': total_seq,
            'compatible_tp_length': int(row.correct_bp),
            'tp_length': int(row.correct_bp),
            'unique_tp_length': int(unique_tp_bp),
            'tp_seq_counts': row.correct_seq_units,
            'foreign_bp': foreign_bp,
            'compatible_bp': compatible_physical,
            'unresolved_bp': unresolved_bp,
            'precision_bp': precision_bp,
            'precision_bp_resolved': precision_bp_resolved,
            'precision_seq': precision_seq,
            'recall_bp': recall_bp,
            'recall_seq': recall_seq,
            'length_gs': identifiable if genome_id else _nan(),
            'seq_counts_gs': unique_seq_den if genome_id else _nan(),
            'matched_genome_identifiable_truth_bp': identifiable if genome_id else _nan(),
            'matched_genome_identifiable_fraction': _div(identifiable, truth.assembly_bp()) if genome_id else _nan(),
            'rank': 'NA',
        })
    precision_df = pd.DataFrame(precision_rows)

    if filter_tail_percentage and not precision_df.empty:
        precision_df = precision_df.copy()
        precision_df['total_length_pct'] = precision_df['total_length'] / precision_df['total_length'].sum()
        precision_df = precision_df.sort_values(by='total_length')
        precision_df['cumsum_length_pct'] = precision_df['total_length_pct'].cumsum()
        mask = precision_df['cumsum_length_pct'] <= filter_tail_percentage / 100.0
        precision_df.loc[mask, ['precision_bp', 'precision_seq']] = np.nan
        precision_df.drop(columns=['cumsum_length_pct', 'total_length_pct'], inplace=True)

    genome_rows = []
    for genome_id in truth.genomes:
        identifiable = truth.identifiable_truth_bp(genome_id)
        unique_seq_units = truth.unique_truth_seq_units(genome_id)
        compatible_involving = truth.compatible_truth_bp_involving(genome_id)
        ident_frac = _div(identifiable, identifiable + compatible_involving + 0)
        # identifiable_fraction vs assembly of that genome's unique+compatible involvement
        genome_mass = identifiable + compatible_involving
        ident_frac = _div(identifiable, genome_mass) if genome_mass else _nan()
        support = unique_support[unique_support['genome_id'] == genome_id] if not unique_support.empty else pd.DataFrame()
        if support.empty or identifiable == 0:
            best_bin = pd.NA
            best_unique_tp = 0
            best_unique_seq = 0.0
            best_recall_bp = _nan() if identifiable == 0 else 0.0
            best_recall_seq = _nan() if unique_seq_units == 0 else 0.0
        else:
            best = support.loc[support['unique_support_bp'].idxmax()]
            best_bin = best['BINID']
            best_unique_tp = int(best['unique_support_bp'])
            best_unique_seq = float(best['unique_support_seq'])
            best_recall_bp = _div(best_unique_tp, identifiable)
            best_recall_seq = _div(best_unique_seq, unique_seq_units)
        genome_rows.append({
            'genome_id': genome_id,
            'identifiable_truth_bp': identifiable,
            'unique_truth_seq_units': unique_seq_units,
            'compatible_truth_bp_involving_genome': compatible_involving,
            'identifiable_fraction': ident_frac,
            'best_bin_id': best_bin,
            'best_unique_tp_bp': best_unique_tp,
            'best_recall_bp': best_recall_bp,
            'best_recall_seq': best_recall_seq,
        })
    genome_df = pd.DataFrame(genome_rows)

    assigned_bp = int(bin_seq_df['length'].sum()) if not bin_seq_df.empty else 0
    assigned_seq = int(bin_seq_df['SEQUENCEID'].nunique()) if not bin_seq_df.empty else 0
    assembly_bp = truth.assembly_bp()
    n_seq = len(truth.sequences)

    correct_bp = float(precision_df['tp_length'].sum()) if not precision_df.empty else 0.0
    correct_seq = float(precision_df['tp_seq_counts'].sum()) if not precision_df.empty else 0.0
    total_binned_bp = float(precision_df['total_length'].sum()) if not precision_df.empty else 0.0
    total_binned_seq = float(precision_df['total_seq_counts'].sum()) if not precision_df.empty else 0.0

    best_unique_tp_bp = float(genome_df['best_unique_tp_bp'].fillna(0).sum()) if not genome_df.empty else 0.0
    identifiable_all = float(genome_df['identifiable_truth_bp'].sum()) if not genome_df.empty else 0.0
    best_unique_seq = 0.0
    unique_seq_all = 0.0
    if not genome_df.empty:
        for row in genome_df.itertuples(index=False):
            unique_seq_all += row.unique_truth_seq_units
            if pd.notna(row.best_bin_id):
                best_unique_seq += (unique_support.loc[
                    (unique_support['BINID'] == row.best_bin_id) & (unique_support['genome_id'] == row.genome_id),
                    'unique_support_seq'
                ].sum() if not unique_support.empty else 0.0)

    fari_seq, fari_n, fari_seq_frac, fari_bp_frac_eligible = _fari_memberships(truth, assignments, weighted=False)
    fari_bp, _, _, _ = _fari_memberships(truth, assignments, weighted=True)

    recovered = []
    for min_c, max_cont in ((a, b) for a in min_completeness for b in max_contamination):
        if precision_df.empty:
            count = 0
        else:
            eligible_bins = precision_df[
                (precision_df['match_status'] == MATCH_RESOLVED)
                & (precision_df['recall_bp'] > min_c)
                & (precision_df['precision_bp'] > (1 - max_cont))
            ]
            count = int(eligible_bins.shape[0])
        recovered.append({
            'min_completeness': min_c,
            'max_contamination': max_cont,
            'count': count,
        })

    summary = truth.summary_row()
    metrics = {
        'truth_model': TRUTH_MODEL_FRACTIONAL,
        'scoring_model_version': SCORING_MODEL_VERSION,
        'fari_implementation_version': FARI_IMPLEMENTATION_VERSION,
        'percentage_of_assigned_bps': _div(assigned_bp, assembly_bp),
        'percentage_of_assigned_seqs': _div(assigned_seq, n_seq),
        'precision_avg_bp': float(precision_df['precision_bp'].mean()) if not precision_df.empty else _nan(),
        'precision_avg_seq': float(precision_df['precision_seq'].mean()) if not precision_df.empty else _nan(),
        'precision_avg_bp_sem': 0.0 if precision_df.empty or pd.isna(precision_df['precision_bp'].sem()) else float(precision_df['precision_bp'].sem()),
        'precision_avg_seq_sem': 0.0 if precision_df.empty or pd.isna(precision_df['precision_seq'].sem()) else float(precision_df['precision_seq'].sem()),
        'precision_weighted_bp': _div(correct_bp, total_binned_bp),
        'precision_weighted_seq': _div(correct_seq, total_binned_seq),
        'recall_avg_bp': float(genome_df['best_recall_bp'].mean()) if not genome_df.empty else _nan(),
        'recall_avg_seq': float(genome_df['best_recall_seq'].mean()) if not genome_df.empty else _nan(),
        'recall_avg_bp_sem': 0.0 if genome_df.empty or pd.isna(genome_df['best_recall_bp'].sem()) else float(genome_df['best_recall_bp'].sem()),
        'recall_avg_seq_sem': 0.0 if genome_df.empty or pd.isna(genome_df['best_recall_seq'].sem()) else float(genome_df['best_recall_seq'].sem()),
        'recall_weighted_bp': _div(best_unique_tp_bp, identifiable_all),
        'recall_weighted_seq': _div(best_unique_seq, unique_seq_all),
        'accuracy_bp': _div(correct_bp, assembly_bp),
        'accuracy_seq': _div(correct_seq, n_seq),
        'misclassification_bp': 1 - _div(correct_bp, total_binned_bp) if total_binned_bp else _nan(),
        'misclassification_seq': 1 - _div(correct_seq, total_binned_seq) if total_binned_seq else _nan(),
        'rand_index_bp': _nan(),
        'rand_index_seq': _nan(),
        'adjusted_rand_index_bp': _nan(),
        'adjusted_rand_index_seq': _nan(),
        'fari_seq': fari_seq,
        'fari_bp': fari_bp,
        'fari_seq_n_sequences': fari_n,
        'fari_seq_sequence_fraction': fari_seq_frac,
        'fari_seq_bp_fraction': fari_bp_frac_eligible,
        'fari_sequence_fraction': fari_seq_frac,
        'fari_bp_fraction': fari_bp_frac_eligible,
        'truth_unique_bp': summary['unique_bp'],
        'truth_compatible_bp': summary['compatible_bp'],
        'truth_unresolved_bp': summary['unresolved_bp'],
        'truth_unique_fraction': summary['unique_truth_fraction'],
        'truth_compatible_fraction': summary['compatible_truth_fraction'],
        'truth_unresolved_fraction': summary['unresolved_truth_fraction'],
    }
    for key, value in list(metrics.items()):
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            metrics[key] = _finite_or_nan(value) if not isinstance(value, int) else value

    return {
        'metrics': metrics,
        'precision_df': precision_df,
        'genome_df': genome_df,
        'unique_support': unique_support,
        'compatible_support': compatible_support,
        'recovered': recovered,
        'truth_summary': summary,
    }
