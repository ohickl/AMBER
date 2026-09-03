# Changelog

## 2.0.8+fractional.1

- Standard AMBER 2.0.8 CLI, CAMI hard gold-standard files, prediction files, taxonomic binning, and hard-truth metrics are unchanged when `-g/--gold_standard_file` is used.
- Opt-in `--fractional-gold-standard` adds fractional / multi-origin genome-truth scoring (`truth_model=fractional-origin-v1`, `scoring_model-v1`).
- Compatible `{A,B}` sequence may count toward bin purity for either genome and is excluded from primary identifiable completeness denominators.
- Sequence-count metrics use fractional contig units. Recovered-genome counts require a single resolved matched genome.
- Frobenius Adjusted Rand Index (FARI) is reported for fully unique-origin sequences only. Legacy Rand/ARI columns are NA in fractional mode and are not reused as FARI labels.
- This is a fork extension, not an official CAMI/AMBER metric set.

## 2.0.8

Upstream CAMI AMBER 2.0.8.
