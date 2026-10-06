# Cohort comparison report

`amber_cohort.py --manifest cohort.json --output report --model-cache cache`
evaluates native multi-SampleID inputs for official and fractional scoring,
then writes one self-contained `report/index.html`. Each assembly observation
has a distinct SampleID. The metadata maps those observations to biological
samples and assembly variants for presentation; truth denominators are never
pooled across different assemblies.

The manifest has schema `amber-cohort-v1`, `domains` and `observations` TSV paths,
and `models.official` / `models.fractional` entries with a `truth` path and an
ordered `predictions` list of `{ "label": "method", "path": "method.tsv" }`.
Paths resolve relative to the manifest. Both models must cover exactly the same
observations and ordered methods, including header-only zero-bin sections.
Domain TSV columns are `SampleID`, `GenomeID`, `Domain`; SampleID refers to the
biological sample. Observation columns are `SampleID`, `BiologicalSampleID`,
`AssemblyVariant`.

All report pages share domain checkboxes for Archaea, Bacteria, Viruses,
Eukaryotes, Plasmids and Unknown. Authoritative truth supplies the labels.
Filtering retains each bin's original matching and complete contamination
denominator. Recall includes missing selected truth genomes. Fractional
compatible components touching selected genomes count once; FARI excludes
sequences spanning selected and excluded domains. The all-domain scores
preserve the original scoring models. Undefined denominators remain N/A.

Recovered-genome and Overlap pages share strict completeness and contamination
threshold sliders: completeness greater than the selected value and
contamination less than the selected value, matching AMBER's recovery policy.
Recovery allows sample/tool filtering; overlap selects a sample or sums over
sample/genome identities. Multiple predicted bins recovering the same truth
genome count once for overlap. Venn diagrams support two or three selected
methods, with schematic areas and exact region counts. Larger selections use
UpSet graphics for the 30 largest exclusive intersections plus complete
intersection, per-tool unique-recovery and sample/genome membership tables.

Categorical colours use distinctipy 1.3.4 with seed 42 and sorted full-cohort
labels. `method_palette.json` records the fixed mapping, reused by filtered
views and both scoring models. Colour-vision simulations are validation aids;
method labels and tables remain available independently of colour.

Completed model scores and domain profiles are sealed atomically in a
content-addressed cache bound to input hashes, scoring parameters and source
hashes. Restore verifies every sealed file. Corruption fails explicitly rather
than silently recomputing. An interrupted report can reuse a completed model.
