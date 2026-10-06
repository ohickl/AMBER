"""Validated truth-genome domains; taxonomy never comes from predicted bins."""

import csv
import itertools


DOMAINS = ('Archaea', 'Bacteria', 'Viruses', 'Eukaryotes', 'Plasmids', 'Unknown')


def read_genome_domains(path):
    domains = {}
    with open(path, encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != ['SampleID', 'GenomeID', 'Domain']:
            raise ValueError('Domain manifest must have SampleID, GenomeID, Domain columns')
        for row in reader:
            sample, genome, domain = (row[key] for key in reader.fieldnames)
            if not sample or not genome or domain not in DOMAINS:
                raise ValueError('Invalid genome-domain manifest row: {!r}'.format(row))
            key = sample, genome
            if key in domains:
                raise ValueError('Duplicate genome-domain identity: {!r}'.format(key))
            domains[key] = domain
    if not domains:
        raise ValueError('Empty genome-domain manifest')
    return domains


def validate_domain_coverage(domains, sample_genomes):
    missing = sorted((sample, str(genome)) for sample, genomes in sample_genomes.items()
                     for genome in genomes if (sample, str(genome)) not in domains)
    if missing:
        raise ValueError('Truth genomes lack domain labels: {!r} ({} total)'.format(missing[:10], len(missing)))


def read_observation_metadata(path):
    metadata = {}
    with open(path, encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != ['SampleID', 'BiologicalSampleID', 'AssemblyVariant']:
            raise ValueError('Observation metadata requires SampleID, BiologicalSampleID, AssemblyVariant')
        for row in reader:
            sample = row['SampleID']
            if sample in metadata or any(not value for value in row.values()):
                raise ValueError('Duplicate or empty observation metadata row')
            metadata[sample] = row
    if not metadata:
        raise ValueError('Empty observation metadata')
    return metadata


def configure_domain_views(queries, options, options_gs, path, observations=None):
    domains = read_genome_domains(path)
    biological = {sample: row['BiologicalSampleID'] for sample, row in (observations or {}).items()}
    sample_genomes = {}
    for sample, sample_queries in queries.items():
        if observations is not None and sample not in observations:
            raise ValueError('Observation metadata missing sample {!r}'.format(sample))
        for query in sample_queries:
            genomes = (query.fractional_truth.genomes if hasattr(query, 'fractional_truth')
                       else query.gold_standard.df['BINID'].unique())
            sample_genomes.setdefault(biological.get(sample, sample), set()).update(map(str, genomes))
    validate_domain_coverage(domains, sample_genomes)
    present = sorted({DOMAINS.index(domains[sample, genome])
                      for sample, genomes in sample_genomes.items() for genome in genomes})
    masks = [sum(1 << bit for bit, enabled in zip(present, flags) if enabled)
             for flags in itertools.product((False, True), repeat=len(present))]
    for opt in (options, options_gs):
        opt.genome_domains, opt.biological_samples = domains, biological
        opt.domain_masks = masks


def report_observation_views(summary, bins, observations):
    """Compare variants as distinct series while preserving each truth denominator."""
    from cami_amber.utils import labels
    summary, bins = summary.copy(), bins.copy()
    for frame, sample_column in ((summary, labels.SAMPLE), (bins, 'sample_id')):
        if frame.empty:
            continue
        native_samples = frame[sample_column].copy()
        frame['Observation'] = native_samples
        frame['SourceTool'] = frame[labels.TOOL]
        frame['AssemblyVariant'] = native_samples.map(lambda sample: observations[sample]['AssemblyVariant'])
        frame[labels.TOOL] = frame[labels.TOOL] + ' / ' + frame['AssemblyVariant']
        frame[sample_column] = native_samples.map(lambda sample: observations[sample]['BiologicalSampleID'])
    return summary, bins
