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
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from cami_amber.fractional_truth import KIND_COMPATIBLE, KIND_UNIQUE, FractionalTruthSample
from cami_amber.version import FARI_IMPLEMENTATION_VERSION, SCORING_MODEL_VERSION, TRUTH_MODEL_FRACTIONAL

MATCH_RESOLVED = 'resolved'
MATCH_COMPATIBLE_TIE = 'compatible_tie'
MATCH_UNMATCHED = 'unmatched'
HEATMAP_MAX_CELLS = 250000


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
    pred_sample_id = metadata[2].get('SAMPLEID')
    if pred_sample_id != truth.sample_id:
        raise FractionalPredictionError(
            'Prediction @SampleID {!r} does not match fractional truth sample {!r}'.format(
                pred_sample_id, truth.sample_id
            )
        )
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
    original = set(truth.original_sequence_ids) if truth.original_sequence_ids is not None else set(truth.sequences)
    unknown = set(df['SEQUENCEID']) - original
    if unknown:
        raise FractionalPredictionError(
            'Prediction sequences outside the original truth universe: {}'.format(sorted(unknown)[:10])
        )
    df = df[df['SEQUENCEID'].isin(truth.sequences)].copy()
    return df


def build_support_matrices(
    truth: FractionalTruthSample,
    assignments: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    unique_acc = defaultdict(lambda: [0, 0.0])
    compat_acc = defaultdict(lambda: [0, 0.0])
    bin_len = defaultdict(int)
    bin_n = defaultdict(int)
    bin_unres = defaultdict(int)
    bin_compat_phys = defaultdict(int)
    for rec in assignments.itertuples(index=False):
        seq = truth.sequences[rec.SEQUENCEID]
        bin_id = rec.BINID
        bin_len[bin_id] += seq.length
        bin_n[bin_id] += 1
        bin_unres[bin_id] += seq.unresolved_bp
        bin_compat_phys[bin_id] += seq.compatible_bp
        inv = 1.0 / seq.length if seq.length else 0.0
        for component in seq.components:
            if component.kind == KIND_UNIQUE:
                genome_id = next(iter(component.genome_ids))
                key = (bin_id, genome_id)
                unique_acc[key][0] += component.bp
                unique_acc[key][1] += component.bp * inv
                compat_acc[key][0] += component.bp
                compat_acc[key][1] += component.bp * inv
            elif component.kind == KIND_COMPATIBLE:
                for genome_id in component.genome_ids:
                    key = (bin_id, genome_id)
                    compat_acc[key][0] += component.bp
                    compat_acc[key][1] += component.bp * inv
    unique_support = pd.DataFrame(
        [{'BINID': b, 'genome_id': g, 'unique_support_bp': v[0], 'unique_support_seq': v[1]}
         for (b, g), v in unique_acc.items()]
    )
    compatible_support = pd.DataFrame(
        [{'BINID': b, 'genome_id': g, 'compatible_support_bp': v[0], 'compatible_support_seq': v[1]}
         for (b, g), v in compat_acc.items()]
    )
    bin_seq_df = pd.DataFrame(
        [{'BINID': b, 'total_length': bin_len[b], 'total_seq_counts': bin_n[b],
          'unresolved_bp': bin_unres[b], 'compatible_bp': bin_compat_phys[b]}
         for b in bin_len]
    )
    return unique_support, compatible_support, bin_seq_df


def match_bins(unique_support: pd.DataFrame, compatible_support: pd.DataFrame, bin_ids: Sequence[str]) -> pd.DataFrame:
    records = []
    if compatible_support.empty:
        return pd.DataFrame([{
            'BINID': bin_id,
            'match_status': MATCH_UNMATCHED,
            'matched_genome_id': pd.NA,
            'matched_genome_ids': '',
            'correct_bp': 0,
            'correct_seq_units': 0.0,
        } for bin_id in bin_ids])
    compat = compatible_support.set_index('BINID')
    uniq = unique_support.set_index('BINID') if not unique_support.empty else pd.DataFrame()
    for bin_id, cand in compat.groupby(level=0):
        max_bp = cand['compatible_support_bp'].max()
        tied = cand[cand['compatible_support_bp'] == max_bp]
        if len(tied) == 1:
            genome_id = tied.iloc[0]['genome_id']
            records.append({
                'BINID': bin_id,
                'match_status': MATCH_RESOLVED,
                'matched_genome_id': genome_id,
                'matched_genome_ids': genome_id,
                'correct_bp': int(max_bp),
                'correct_seq_units': float(tied.iloc[0]['compatible_support_seq']),
            })
            continue
        if not uniq.empty and bin_id in uniq.index:
            uniq_cand = uniq.loc[[bin_id]]
            uniq_tied = uniq_cand[uniq_cand['genome_id'].isin(set(tied['genome_id']))]
            if not uniq_tied.empty:
                max_u = uniq_tied['unique_support_bp'].max()
                uniq_best = uniq_tied[uniq_tied['unique_support_bp'] == max_u]
                if len(uniq_best) == 1:
                    genome_id = uniq_best.iloc[0]['genome_id']
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
        records.append({
            'BINID': bin_id,
            'match_status': MATCH_COMPATIBLE_TIE,
            'matched_genome_id': pd.NA,
            'matched_genome_ids': ','.join(genomes),
            'correct_bp': int(max_bp),
            'correct_seq_units': float(tied['compatible_support_seq'].max()),
        })
    seen = {r['BINID'] for r in records}
    for bin_id in bin_ids:
        if bin_id not in seen:
            records.append({
                'BINID': bin_id,
                'match_status': MATCH_UNMATCHED,
                'matched_genome_id': pd.NA,
                'matched_genome_ids': '',
                'correct_bp': 0,
                'correct_seq_units': 0.0,
            })
    return pd.DataFrame(records)


def compute_fari_from_stats(n: float, sA: float, sB: float, qA: float, qB: float, qAB: float, tA: float, tB: float) -> float:
    """Andrews et al. 2022 FARI from bonding-matrix sufficient statistics.

    Line-for-line equivalent of its-likeli-jeff/FARI ``R/fari.R`` (commit
    9f2e7e8769120a26d8456242ba5ed77bf539f2c2) without materializing n×n
    bonding matrices A=UU^T or B=VV^T.
    """
    if n <= 1:
        return 1.0
    if qA == 0 or qB == 0:
        return _nan()
    pair_den = n * (n - 1.0)
    sum_NaNb = (sA / qA) * (sB / qB) * qAB
    sum_j_minus = (n * n) - (sA * sA / qA) - (sB * sB / qB) + sum_NaNb
    ri = (sum_NaNb + sum_j_minus - n) / pair_den
    sum_MA = sA / n
    sum_MB = sB / n
    sum_RA = tA - sA / n
    sum_RB = tB - sB / n
    eri = (
        (2.0 * sA * sB / (qA * qB)) * (sum_MA * sum_MB + (1.0 / (n - 1.0)) * sum_RA * sum_RB)
        - (sA * sA) / qA
        - (sB * sB) / qB
        + n * n
        - n
    ) / pair_den
    if eri == 1.0:
        return 1.0
    return _finite_or_nan((ri - eri) / (1.0 - eri))


def fari_stats_from_memberships(U: np.ndarray, V: np.ndarray, weights: Optional[np.ndarray] = None) -> Tuple[float, ...]:
    """Sufficient statistics for FARI; used by tests. Does not build n×n bonding matrices."""
    if weights is None:
        weights = np.ones(U.shape[0], dtype=float)
    n = float(weights.sum())
    w = weights[:, None]
    gram_uu = U.T @ (U * w)
    gram_vv = V.T @ (V * w)
    gram_uv = U.T @ (V * w)
    col_u = U.T @ weights
    col_v = V.T @ weights
    sA = float(np.dot(col_u, col_u))
    sB = float(np.dot(col_v, col_v))
    qA = float(np.square(gram_uu).sum())
    qB = float(np.square(gram_vv).sum())
    qAB = float(np.square(gram_uv).sum())
    tA = float((np.square(U).sum(axis=1) * weights).sum())
    tB = float((np.square(V).sum(axis=1) * weights).sum())
    return n, sA, sB, qA, qB, qAB, tA, tB


def compute_fari(U: np.ndarray, V: np.ndarray, weights: Optional[np.ndarray] = None) -> float:
    return compute_fari_from_stats(*fari_stats_from_memberships(U, V, weights=weights))


def compute_fari_nxn_reference(U: np.ndarray, V: np.ndarray) -> float:
    """Literal n×n translation of its-likeli-jeff/FARI R/fari.R. Tests only; never use in scoring."""
    n = U.shape[0]
    A = U @ U.T
    B = V @ V.T
    j = np.ones((n, n))
    qA = float(np.sum(A * A))
    qB = float(np.sum(B * B))
    sA = float(np.sum(A))
    sB = float(np.sum(B))
    Na = (sA / qA) * A
    Nb = (sB / qB) * B
    ri = (np.sum(Na * Nb) + np.sum((j - Na) * (j - Nb)) - n) / (n * (n - 1))
    M = j / n
    R = np.eye(n) - M
    eri = (
        (2 * sA * sB / (qA * qB)) * (np.sum(M * A) * np.sum(M * B) + (1.0 / (n - 1)) * np.sum(R * A) * np.sum(R * B))
        - (sA ** 2) / qA
        - (sB ** 2) / qB
        + n * n
        - n
    ) / (n * (n - 1))
    return float((ri - eri) / (1.0 - eri))


def compute_identifiable_bp_rand(unique_support: pd.DataFrame, assembly_bp: int) -> Tuple[float, float, float]:
    """Hubert-Arabie RI/ARI on unique-origin bp contingency; same formula as AMBER Metrics.compute_rand_index."""
    from cami_amber.binning_classes import Metrics

    if unique_support is None or unique_support.empty:
        return _nan(), _nan(), 0.0
    confusion = unique_support[['BINID', 'genome_id', 'unique_support_bp']].rename(
        columns={'unique_support_bp': 'bp'}
    )
    confusion['bp'] = confusion['bp'].astype(int)
    participating = int(confusion['bp'].sum())
    ri, ari = Metrics.compute_rand_index(confusion, 'BINID', 'genome_id', 'bp')
    return float(ri), float(ari), _div(participating, assembly_bp)


def compute_sample_fari(truth: FractionalTruthSample, assignments: pd.DataFrame, weighted: bool) -> Tuple[float, int, float, float]:
    assigned_map = dict(zip(assignments['SEQUENCEID'], assignments['BINID']))
    eligible = [seq for seq in truth.sequences.values() if seq.is_fari_eligible() and seq.sequence_id in assigned_map]
    n_all = len(truth.sequences)
    bp_all = truth.assembly_bp()
    if not eligible:
        return _nan(), 0, 0.0, 0.0
    genomes = []
    genome_index = {}
    for seq in eligible:
        for component in seq.components:
            if component.kind == KIND_UNIQUE:
                g = next(iter(component.genome_ids))
                if g not in genome_index:
                    genome_index[g] = len(genomes)
                    genomes.append(g)
    bins = []
    bin_index = {}
    for seq in eligible:
        b = assigned_map[seq.sequence_id]
        if b not in bin_index:
            bin_index[b] = len(bins)
            bins.append(b)
    k = len(genomes)
    gram_uu = np.zeros((k, k), dtype=float)
    col_u = np.zeros(k, dtype=float)
    bin_w = np.zeros(len(bins), dtype=float)
    gram_uv = np.zeros((k, len(bins)), dtype=float)
    n = 0.0
    tA = 0.0
    tB = 0.0
    for seq in eligible:
        w = float(seq.length if weighted else 1.0)
        n += w
        contrib = {}
        for component in seq.components:
            if component.kind == KIND_UNIQUE:
                gi = genome_index[next(iter(component.genome_ids))]
                contrib[gi] = contrib.get(gi, 0.0) + component.bp / seq.length
        b = bin_index[assigned_map[seq.sequence_id]]
        norm2 = 0.0
        items = list(contrib.items())
        for gi, ug in items:
            col_u[gi] += w * ug
            gram_uv[gi, b] += w * ug
            norm2 += ug * ug
            for hj, uh in items:
                gram_uu[gi, hj] += w * ug * uh
        bin_w[b] += w
        tA += w * norm2
        tB += w
    sA = float(np.dot(col_u, col_u))
    sB = float(np.dot(bin_w, bin_w))
    qA = float(np.square(gram_uu).sum())
    qB = float(np.square(bin_w).sum())
    qAB = float(np.square(gram_uv).sum())
    value = compute_fari_from_stats(n, sA, sB, qA, qB, qAB, tA, tB)
    seq_frac = len(eligible) / n_all if n_all else _nan()
    bp_frac = sum(seq.length for seq in eligible) / bp_all if bp_all else _nan()
    return value, len(eligible), seq_frac, bp_frac


def compute_sample_metrics(
    truth: FractionalTruthSample,
    assignments: pd.DataFrame,
    min_completeness: Sequence[float],
    max_contamination: Sequence[float],
    filter_tail_percentage: float = 0.0,
) -> dict:
    unique_support, compatible_support, bin_totals = build_support_matrices(truth, assignments)
    bin_ids = sorted(bin_totals['BINID'].tolist()) if not bin_totals.empty else []
    matches = match_bins(unique_support, compatible_support, bin_ids)
    agg = truth.genome_aggregates()

    unique_lookup = {}
    if not unique_support.empty:
        for row in unique_support.itertuples(index=False):
            unique_lookup[(row.BINID, row.genome_id)] = (row.unique_support_bp, row.unique_support_seq)

    totals_idx = bin_totals.set_index('BINID') if not bin_totals.empty else pd.DataFrame()
    precision_rows = []
    for row in matches.itertuples(index=False):
        totals = totals_idx.loc[row.BINID]
        total_bp = int(totals['total_length'])
        total_seq = int(totals['total_seq_counts'])
        unresolved_bp = int(totals['unresolved_bp'])
        compatible_physical = int(totals['compatible_bp'])
        genome_id = row.matched_genome_id if pd.notna(row.matched_genome_id) else None
        unique_tp_bp, unique_tp_seq = unique_lookup.get((row.BINID, genome_id), (0, 0.0)) if genome_id else (0, 0.0)
        identifiable = truth.identifiable_truth_bp(genome_id) if genome_id else 0
        unique_seq_den = truth.unique_truth_seq_units(genome_id) if genome_id else 0.0
        ident_frac = truth.identifiable_fraction(genome_id) if genome_id else _nan()
        precision_bp = _div(row.correct_bp, total_bp)
        precision_bp_resolved = _div(row.correct_bp, total_bp - unresolved_bp)
        precision_seq = _div(row.correct_seq_units, total_seq)
        if row.match_status != MATCH_RESOLVED:
            recall_bp = _nan()
            recall_seq = _nan()
        else:
            recall_bp = _div(unique_tp_bp, identifiable)
            recall_seq = _div(unique_tp_seq, unique_seq_den)
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
            'foreign_bp': total_bp - int(row.correct_bp),
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
            'matched_genome_identifiable_fraction': ident_frac,
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

    best_by_genome = pd.DataFrame()
    if not unique_support.empty:
        idx = unique_support.groupby('genome_id')['unique_support_bp'].idxmax()
        best_by_genome = unique_support.loc[idx].set_index('genome_id')

    genome_rows = []
    for genome_id in truth.genomes:
        identifiable = int(agg['identifiable_bp'].get(genome_id, 0))
        unique_seq_units = float(agg['unique_seq'].get(genome_id, 0.0))
        compatible_involving = int(agg['compatible_bp'].get(genome_id, 0))
        ident_frac = _div(identifiable, identifiable + compatible_involving)
        has_best = (not best_by_genome.empty) and (genome_id in best_by_genome.index)
        if has_best and identifiable > 0:
            best = best_by_genome.loc[genome_id]
            best_bin = best['BINID']
            best_unique_tp = int(best['unique_support_bp'])
            best_unique_seq = float(best['unique_support_seq'])
            best_recall_bp = _div(best_unique_tp, identifiable)
            best_recall_seq = _div(best_unique_seq, unique_seq_units)
        else:
            best_bin = pd.NA
            best_unique_tp = 0
            best_unique_seq = 0.0
            best_recall_bp = _nan() if identifiable == 0 else 0.0
            best_recall_seq = _nan() if unique_seq_units == 0 else 0.0
        genome_rows.append({
            'genome_id': genome_id,
            'identifiable_truth_bp': identifiable,
            'unique_truth_seq_units': unique_seq_units,
            'compatible_truth_bp_involving_genome': compatible_involving,
            'identifiable_fraction': ident_frac,
            'best_bin_id': best_bin,
            'best_unique_tp_bp': best_unique_tp,
            'best_unique_tp_seq': best_unique_seq,
            'best_recall_bp': best_recall_bp,
            'best_recall_seq': best_recall_seq,
        })
    genome_df = pd.DataFrame(genome_rows)

    assigned_bp = int(bin_totals['total_length'].sum()) if not bin_totals.empty else 0
    assigned_seq = int(bin_totals['total_seq_counts'].sum()) if not bin_totals.empty else 0
    assembly_bp = truth.assembly_bp()
    n_seq = len(truth.sequences)
    correct_bp = float(precision_df['tp_length'].sum()) if not precision_df.empty else 0.0
    correct_seq = float(precision_df['tp_seq_counts'].sum()) if not precision_df.empty else 0.0
    total_binned_bp = float(precision_df['total_length'].sum()) if not precision_df.empty else 0.0
    total_binned_seq = float(precision_df['total_seq_counts'].sum()) if not precision_df.empty else 0.0
    best_unique_tp_bp = float(genome_df['best_unique_tp_bp'].fillna(0).sum()) if not genome_df.empty else 0.0
    identifiable_all = float(genome_df['identifiable_truth_bp'].sum()) if not genome_df.empty else 0.0
    unique_seq_all = float(genome_df['unique_truth_seq_units'].sum()) if not genome_df.empty else 0.0
    best_unique_seq = float(genome_df['best_unique_tp_seq'].fillna(0).sum()) if not genome_df.empty else 0.0

    fari_seq, fari_n, fari_seq_frac, fari_seq_bp_frac = compute_sample_fari(truth, assignments, weighted=False)
    ri_bp_id, ari_bp_id, ari_bp_id_frac = compute_identifiable_bp_rand(unique_support, assembly_bp)

    recovered = []
    for min_c in min_completeness:
        for max_cont in max_contamination:
            if precision_df.empty:
                count = 0
            else:
                count = int(precision_df[
                    (precision_df['match_status'] == MATCH_RESOLVED)
                    & (precision_df['recall_bp'] > min_c)
                    & (precision_df['precision_bp'] > (1 - max_cont))
                ].shape[0])
            recovered.append({'min_completeness': min_c, 'max_contamination': max_cont, 'count': count})

    summary = truth.summary_row()
    nan = _nan()
    metrics = {
        'truth_model': TRUTH_MODEL_FRACTIONAL,
        'scoring_model_version': SCORING_MODEL_VERSION,
        'fari_implementation_version': FARI_IMPLEMENTATION_VERSION,
        'percentage_of_assigned_bps': _div(assigned_bp, assembly_bp),
        'percentage_of_assigned_seqs': _div(assigned_seq, n_seq),
        'precision_avg_bp': float(precision_df['precision_bp'].mean()) if not precision_df.empty else nan,
        'precision_avg_seq': float(precision_df['precision_seq'].mean()) if not precision_df.empty else nan,
        'precision_avg_bp_sem': 0.0 if precision_df.empty or pd.isna(precision_df['precision_bp'].sem()) else float(precision_df['precision_bp'].sem()),
        'precision_avg_seq_sem': 0.0 if precision_df.empty or pd.isna(precision_df['precision_seq'].sem()) else float(precision_df['precision_seq'].sem()),
        'precision_weighted_bp': _div(correct_bp, total_binned_bp),
        'precision_weighted_seq': _div(correct_seq, total_binned_seq),
        'recall_avg_bp': float(genome_df['best_recall_bp'].mean()) if not genome_df.empty else nan,
        'recall_avg_seq': float(genome_df['best_recall_seq'].mean()) if not genome_df.empty else nan,
        'recall_avg_bp_sem': 0.0 if genome_df.empty or pd.isna(genome_df['best_recall_bp'].sem()) else float(genome_df['best_recall_bp'].sem()),
        'recall_avg_seq_sem': 0.0 if genome_df.empty or pd.isna(genome_df['best_recall_seq'].sem()) else float(genome_df['best_recall_seq'].sem()),
        'recall_weighted_bp': _div(best_unique_tp_bp, identifiable_all),
        'recall_weighted_seq': _div(best_unique_seq, unique_seq_all),
        'accuracy_bp': _div(correct_bp, assembly_bp),
        'accuracy_seq': _div(correct_seq, n_seq),
        'misclassification_bp': 1 - _div(correct_bp, total_binned_bp) if total_binned_bp else nan,
        'misclassification_seq': 1 - _div(correct_seq, total_binned_seq) if total_binned_seq else nan,
        'rand_index_bp': nan,
        'rand_index_seq': nan,
        'adjusted_rand_index_bp': nan,
        'adjusted_rand_index_seq': nan,
        'recall_avg_bp_cami1': nan,
        'recall_avg_seq_cami1': nan,
        'recall_avg_bp_sem_cami1': nan,
        'recall_avg_seq_sem_cami1': nan,
        'recall_avg_bp_var_cami1': nan,
        'f1_score_bp_cami1': nan,
        'f1_score_seq_cami1': nan,
        'fari_seq': fari_seq,
        'fari_seq_n_sequences': fari_n,
        'fari_seq_sequence_fraction': fari_seq_frac,
        'fari_sequence_fraction': fari_seq_frac,
        'fari_seq_bp_fraction': fari_seq_bp_frac,
        'rand_index_bp_identifiable': ri_bp_id,
        'adjusted_rand_index_bp_identifiable': ari_bp_id,
        'ari_bp_identifiable_fraction': ari_bp_id_frac,
        'truth_unique_bp': summary['unique_bp'],
        'truth_compatible_bp': summary['compatible_bp'],
        'truth_unresolved_bp': summary['unresolved_bp'],
        'truth_unique_fraction': summary['unique_truth_fraction'],
        'truth_compatible_fraction': summary['compatible_truth_fraction'],
        'truth_unresolved_fraction': summary['unresolved_truth_fraction'],
    }
    for key, value in list(metrics.items()):
        if isinstance(value, float):
            metrics[key] = _finite_or_nan(value)
    return {
        'metrics': metrics,
        'precision_df': precision_df,
        'genome_df': genome_df,
        'unique_support': unique_support,
        'compatible_support': compatible_support,
        'recovered': recovered,
        'truth_summary': summary,
    }
