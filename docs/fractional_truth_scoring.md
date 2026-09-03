# Fractional / multi-origin truth scoring

**Truth model:** `fractional-origin-v1`  
**Schema version:** `@Version:0.1.0` (required; distinct from the truth-model tag)  
**Scoring model:** `scoring-model-v1`  
**FARI implementation:** `fari-andrews-2022-v1`  
**FARI reference:** `its-likeli-jeff/FARI` `R/fari.R` commit `9f2e7e8769120a26d8456242ba5ed77bf539f2c2`; Andrews, Browne & Hvingelby, *Journal of Classification* 2022, DOI 10.1007/s00357-021-09407-3  
**Upstream base:** AMBER 2.0.8 (GPL-3.0-or-later)

These metrics are a **fork extension**, not official CAMI/AMBER metrics.

Standard AMBER with `-g` is unchanged.

---

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

Repeated rows with the same `(SEQUENCEID, COMPONENT_TYPE, genome set)` are **merged by summing `COMPONENT_BP`** (composition table, not genomic coordinates). After merge, `sum(COMPONENT_BP)` must equal `_LENGTH`.

---

## Matching (predicted BINID has no truth meaning)

A predicted bin name such as `binB` does **not** force a match to genome B.

`correct_bp[b] = max_g compatible_support_bp[b,g]`.

A lone contig `70 unique A + 30 unique B` therefore matches **A** and has purity **0.70** even if the predicted BINID is `binB`.

A fully compatible `{A,B}` contig may have purity 1 with `match_status=compatible_tie` and **must not** recover A or B.

---

## Purity, completeness, recovered genomes

- Primary purity uses max compatible support over physical bin bp.
- `precision_bp_resolved` omits unresolved bp from the denominator.
- Completeness uses uniquely attributable truth only.
- Compatible-tie bins have NA recall and cannot count as recovered genomes.
- Zero-identifiable genomes have completeness NA.
- Thresholds use upstream AMBER's strict `>` comparisons.
- Per-genome and per-bin `identifiable_fraction` is `unique(g) / (unique(g) + compatible involving g)`, **not** unique/assembly.

---

## `min_length`

Filtering removes entire sequences from the **evaluation** universe. Prediction rows for those IDs are dropped. A prediction ID that was **never** in the original truth file remains fatal.

## `remove_genomes`

**Unsupported** with `--fractional-gold-standard` in fractional-v1 (fail-fast). Use standard `-g` mode for genome exclusion.

---

## Clustering agreement (v1 taxonomy)

**Purity/completeness** use fractional bp and fractional sequence-unit metrics (compatibility-aware purity; unique-origin completeness).

**Sequence clustering agreement** is FARI (`fari_seq`) on **assigned** sequences whose entire length is unique-origin (FARI-eligible). Each contig is one observation. A 70/30 A/B contig is the fuzzy membership row `[A=0.7, B=0.3]`. Compatible/unresolved sequences are excluded. Coverage:

- `fari_seq_participating_fraction_of_assembly` / `fari_seq_bp_fraction`: assigned FARI-eligible bp / assembly bp
- `fari_seq_assignment_fraction_of_identifiable_truth`: assigned FARI-eligible sequences / all FARI-eligible sequences in truth
- `fari_seq_sequence_fraction`: assigned FARI-eligible sequences / all truth sequences

These are not “fraction of the assembly that is identifiable” alone.

Bonding matrices `A = U Uᵀ` and `B = V Vᵀ` are not materialized. Sufficient statistics:

```
sA = ||Uᵀ 1||²    qA = ||Uᵀ U||_F²    qAB = ||Uᵀ V||_F²    tA = Σ_i ||u_i||²
```

plus the `Na`, `Nb`, FRI, expected-FRI, FARI equations of `R/fari.R`.

**Base-pair clustering agreement** is ordinary Hubert–Arabie RI/ARI on uniquely attributable physical bp (`rand_index_bp_identifiable`, `adjusted_rand_index_bp_identifiable`), using the same combinatorics as standard AMBER `adjusted_rand_index_bp`. The contingency is

```
C[bin, genome] = sum unique-origin COMPONENT_BP for sequences assigned to that bin
```

A 70/30 contig therefore contributes **70 A observations and 30 B observations**, not 100 fuzzy `[.7,.3]` bases. Compatible/unresolved bp are excluded. Unbinned identifiable bp are **not** placed in the ARI contingency.

Coverage:

- `ari_bp_participating_fraction_of_assembly` = assigned unique-origin bp / assembly bp
- `ari_bp_assignment_fraction_of_identifiable_truth` = assigned unique-origin bp / total unique-origin truth bp

The unique-origin heatmap has predicted-bin rows plus `__UNASSIGNED_IDENTIFIABLE__`. Each genome column sums to `identifiable_truth_bp[genome]`. Compatible/unresolved bp are omitted from the heatmap; sample-level unique/compatible/unresolved fractions remain in the summary.

`matched_genome_ids` is a JSON array (safe if genome IDs contain commas).

`component_count` is the **canonical post-merge** component count. `input_component_rows` is the raw input row count.

Average SEM values are NA when fewer than two observations exist (not 0).

Length-weighted fuzzy FARI (replicating `[.7,.3]` once per base) is **not** a v1 headline metric and is not named `fari_bp`.

Legacy hard Rand/ARI and CAMI1 completeness fields remain **NA** in fractional mode. HTML shows FARI (seq) and identifiable-bp ARI.

Production never expands one row per base.

---

## Limitations

- Completeness measures identifiable unique-origin recovery only.
- Compatible sequence is set-valued, not 50/50 ancestry.
- Sequence FARI and identifiable-bp ARI both exclude compatible/unresolved sequence; coverage is reported separately.
- Fractional mode does not synthesize gold-standard-vs-self.
- Keep a standard AMBER `-g` run for external comparability.
