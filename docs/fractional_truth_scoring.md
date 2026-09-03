# Fractional / multi-origin truth scoring

**Truth model:** `fractional-origin-v1`  
**Scoring model:** `scoring-model-v1`  
**FARI implementation:** `fari-frobenius-v1`  
**Upstream base:** AMBER 2.0.8 (GPL-3.0-or-later)

This document specifies the opt-in evaluator in this fork. Standard AMBER with `-g` is unchanged.

## Schema

```
@Version:0.1.0
@SampleID:<assembly-artifact-id>
@TruthModel:fractional-origin-v1
@@SEQUENCEID	_LENGTH	COMPONENT_BP	COMPONENT_TYPE	GENOME_IDS
```

`GENOME_IDS` is a JSON array of strings.

| COMPONENT_TYPE | GENOME_IDS | Role |
| --- | --- | --- |
| unique | exactly one | Identifiable origin; completeness denominator |
| compatible | >= 2 distinct | Purity-correct if the bin match is in the set; not in primary completeness |
| unresolved | `[]` | Conservatively non-correct in primary purity |

Invariants: positive integer `COMPONENT_BP`; per-sequence `sum(COMPONENT_BP) == _LENGTH`; no duplicate component rows; prediction sequences must be a subset of the truth universe; duplicate prediction assignments and blank BINID are fatal.

## Matching

`compatible_support_bp[b,g]` sums component bp of sequences in bin `b` for which `g` is allowed (unique or compatible).

`unique_support_bp[b,g]` sums unique-origin bp only.

Match: max compatible support, then unique-support tie-break. Remaining ties are `match_status=compatible_tie` with `matched_genome_id=NA`. Purity remains defined; the bin cannot recover a genome.

## Purity and completeness

- `precision_bp = max_g compatible_support_bp / bin total bp`
- `precision_bp_resolved` omits unresolved bp from the denominator
- `precision_seq` uses `compatible_bp / length` summed over contigs in the bin
- Completeness uses uniquely attributable truth only
- Weighted recall chooses, per genome, the bin with greatest unique support

Compatible `{A,B}` is set-valued uncertainty, not 50/50 ancestry.

## FARI

Sequence FARI uses the Frobenius / Hubert–Arabie formula of Andrews et al. (Journal of Classification 2022, DOI 10.1007/s00357-021-09407-3). Only sequences whose entire length is unique-origin are included. Compatible/unresolved sequences are excluded (not converted to 0.5/0.5).

bp-weighted FARI is the same closed form with sequence-length weights and is tested against explicit row replication.

Hard one-hot memberships reduce to AMBER `Metrics.compute_rand_index` ARI.

Independent oracle used in tests: AMBER's own hard ARI plus replication equivalence. R `MoEClust::FARI` was not executed in this environment.

## `min_length` and `remove_genomes`

`min_length` drops entire sequences. Removing genome `g` turns unique `{g}` into unresolved for scoring, collapses compatible `{g,h}` to unique `{h}`, and turns compatible sets whose candidates are all removed into unresolved.

## Limitations

See the implementation plan: resolver quality, identifiable-only completeness, FARI eligibility, and the requirement to label these metrics as fork extensions in publications.
