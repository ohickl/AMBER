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

## FARI

Eligible sequences: assigned, 100% unique-origin (no compatible, no unresolved). Coverage is reported as `fari_seq_n_sequences`, `fari_sequence_fraction`, `fari_bp_fraction`.

Bonding matrices `A = U Uᵀ`, `B = V Vᵀ` are **not** materialized. Sufficient statistics:

```
sA = ||Uᵀ 1||²    qA = ||Uᵀ U||_F²    qAB = ||Uᵀ V||_F²    tA = Σ_i ||u_i||²
```

and the `Na`, `Nb`, FRI, expected-FRI, FARI equations of `R/fari.R`.

bp-weighted FARI uses the same equations with conceptual row replication by integer length `w_i`:

```
n = Σ w_i
Uᵀ w,  Uᵀ diag(w) U,  trace terms Σ w_i ||u_i||²
```

Legacy Rand/ARI and CAMI1 completeness fields are **NA** in fractional mode. HTML ranks FARI, not ARI.

---

## Limitations

- Completeness measures identifiable unique-origin recovery only.
- Compatible sequence is set-valued, not 50/50 ancestry.
- FARI v1 excludes compatible/unresolved sequences.
- Fractional mode does not synthesize gold-standard-vs-self.
- Keep a standard AMBER `-g` run for external comparability.
