"""Domain views retain full-bin matching and contamination denominators."""

import numpy as np
import pandas as pd

from cami_amber.binning_classes import Metrics
from cami_amber.utils import labels
from cami_amber.domain_metadata import DOMAINS


def official_domain_metrics(precision, recall, confusion, truth_sizes, assigned, selected):
    """Scope target genomes without rematching or cleaning any predicted bin.

    Recall and assignment include all selected truth genomes, including genomes
    missing from predictions. Rand metrics use selected-origin confusion cells.
    """
    precision = precision[precision['genome_id'].isin(selected)].copy()
    recall = recall[recall['genome_id'].isin(selected)].copy()
    sizes = truth_sizes[truth_sizes.index.isin(selected)]
    cells = confusion[confusion.index.get_level_values('genome_id').isin(selected)]
    assigned = assigned[assigned.index.isin(selected)]
    metrics = Metrics()
    for attribute in list(vars(metrics)):
        setattr(metrics, attribute, np.nan)

    def divide(numerator, denominator):
        return numerator / denominator if denominator else np.nan

    for unit, true_length, tp, bin_total, cell_weight in [
        ('bp', 'length_gs', 'tp_length', 'total_length', 'genome_length'),
        ('seq', 'seq_counts_gs', 'tp_seq_counts', 'total_seq_counts', 'genome_seq_counts'),
    ]:
        setattr(metrics, 'precision_avg_' + unit, precision['precision_' + unit].mean())
        setattr(metrics, 'precision_avg_' + unit + '_sem', precision['precision_' + unit].sem())
        setattr(metrics, 'precision_weighted_' + unit, divide(precision[tp].sum(), precision[bin_total].sum()))
        setattr(metrics, 'recall_avg_' + unit, recall['recall_' + unit].mean())
        setattr(metrics, 'recall_avg_' + unit + '_sem', recall['recall_' + unit].sem())
        setattr(metrics, 'recall_weighted_' + unit, divide(recall[cell_weight].sum(), sizes[true_length].sum()))
        setattr(metrics, 'accuracy_' + unit, divide(precision[tp].sum(), sizes[true_length].sum()))
        missing = set(selected) - set(precision['genome_id'])
        cami1 = pd.Series(precision['recall_' + unit].tolist() + [0.0] * len(missing), dtype=float)
        setattr(metrics, 'recall_avg_' + unit + '_cami1', cami1.mean())
        setattr(metrics, 'recall_avg_' + unit + '_sem_cami1', cami1.sem())
        if unit == 'bp':
            metrics.recall_avg_bp_var_cami1 = cami1.var()
        if not cells.empty:
            ri, ari = Metrics.compute_rand_index(cells, 'BINID', 'genome_id', cell_weight)
            setattr(metrics, 'rand_index_' + unit, ri)
            setattr(metrics, 'adjusted_rand_index_' + unit, ari)
    metrics.precision_avg_bp_var = precision['precision_bp'].var()
    metrics.recall_avg_bp_var = recall['recall_bp'].var()
    metrics.percentage_of_assigned_bps = divide(assigned['seq_length'].sum(), sizes['length_gs'].sum())
    metrics.percentage_of_assigned_seqs = divide(assigned['SEQUENCEID'].sum(), sizes['seq_counts_gs'].sum())
    return metrics.get_ordered_dict(), precision


def with_query_identity(metrics, query):
    return dict(metrics, **{labels.TOOL: query.label, labels.SAMPLE: query.sample_id,
                            labels.BINNING_TYPE: 'genome', labels.RANK: 'NA'})


def selected_genomes(query, genomes, mask):
    sample = getattr(query.options, 'biological_samples', {}).get(query.sample_id, query.sample_id)
    domains = query.options.genome_domains
    return {str(genome) for genome in genomes
            if mask & (1 << DOMAINS.index(domains[sample, str(genome)]))}


def capture_official_profiles(query, confusion, truth_sizes, assigned):
    if not hasattr(query.options, 'genome_domains'):
        return
    query.domain_profiles = {}
    genomes = set(truth_sizes.index)
    for mask in query.options.domain_masks:
        selected = selected_genomes(query, genomes, mask)
        if selected == genomes:
            query.domain_profiles[mask] = (query.get_metrics_df(), query.precision_df.reset_index())
            continue
        metrics, bins = official_domain_metrics(query.precision_df, query.recall_df, confusion,
                                                truth_sizes, assigned, selected)
        query.domain_profiles[mask] = (pd.DataFrame([with_query_identity(metrics, query)]), bins.reset_index())


def capture_fractional_profiles(query, assignments, result):
    if not hasattr(query.options, 'genome_domains'):
        return
    from cami_amber.fractional_metrics import compute_sample_metrics
    query.domain_profiles = {}
    genomes = set(query.fractional_truth.genomes)
    for mask in query.options.domain_masks:
        selected = selected_genomes(query, genomes, mask)
        if selected == genomes:
            query.domain_profiles[mask] = (query.get_metrics_df(), query.precision_df.reset_index())
            continue
        view = compute_sample_metrics(query.fractional_truth, assignments,
                                      query.options.min_completeness, query.options.max_contamination,
                                      query.options.filter_tail_percentage,
                                      selected_genomes=selected, support_cache=result['support_cache'])
        standard = Metrics()
        for key, value in view['metrics'].items():
            if hasattr(standard, key):
                setattr(standard, key, value)
        metrics = standard.get_ordered_dict()
        metrics.update({key: value for key, value in view['metrics'].items() if key not in metrics})
        bins = view['precision_df'].copy()
        bins[labels.TOOL], bins['sample_id'] = query.label, query.sample_id
        query.domain_profiles[mask] = (pd.DataFrame([with_query_identity(metrics, query)]), bins)


def combine_profiles(queries_by_sample, masks):
    profiles = {}
    for mask in masks:
        frames = [query.domain_profiles[mask] for queries in queries_by_sample.values()
                  for query in queries if query.eval_success]
        profiles[mask] = (pd.concat([frame[0] for frame in frames], ignore_index=True, sort=True),
                          pd.concat([frame[1] for frame in frames], ignore_index=True, sort=True))
    return profiles
