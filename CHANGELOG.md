# Changelog

## 2.0.8+fractional.1

- Standard AMBER 2.0.8 CLI, CAMI hard gold-standard files, prediction files, taxonomic binning, and hard-truth metrics are unchanged when `-g/--gold_standard_file` is used.
- Opt-in `--fractional-gold-standard` is mutually exclusive with `-g`.
- Compatible `{A,B}` sequence may count toward bin purity for either genome and is excluded from primary identifiable completeness denominators. Predicted BINID strings have no matching semantics.
- Sequence clustering agreement is FARI (Andrews et al. 2022 / `its-likeli-jeff/FARI` `R/fari.R`) on identifiable contig memberships.
- Base-pair clustering agreement is ordinary Hubert–Arabie RI/ARI on unique-origin bp (`adjusted_rand_index_bp_identifiable`), not length-weighted fuzzy FARI.
- ARI coverage is split into participating-of-assembly vs assignment-of-identifiable-truth. SEM is NA when undefined. Heatmap includes `__UNASSIGNED_IDENTIFIABLE__`. `matched_genome_ids` is JSON. `component_count` is post-merge.
- `--min_length` drops matching prediction rows; IDs absent from the original truth remain fatal.
- `--remove_genomes` is unsupported in fractional-v1 (fail-fast).
- Repeated `(sequence, kind, genome-set)` components are merged by summing bp.
- Legacy Rand/ARI and CAMI1 completeness fields are NA in fractional mode; HTML shows FARI.
- This is a fork extension, not an official CAMI/AMBER metric set.

## 2.0.8

Upstream CAMI AMBER 2.0.8.
