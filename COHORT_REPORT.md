# Cohort presentation

Cohort reports retain separate official and fractional models and all six tabs.
The selected domain combination is a compressed, self-contained Bokeh document.
Only that document is decompressed and rendered; switching domains removes the
old views and document and restores the active model/tab and every tool/sample
selection, including Metrics and Metrics per bin. Empty selections are retained.
Choices unavailable in a domain view remain remembered until they return; a new
explicit selection replaces that remembered choice. Restoration waits for binner
callbacks to rebuild sample options before applying the saved sample. Modern browsers with `DecompressionStream` can open the HTML directly
from disk without a web server or network resources.

Recovery and overlap filters use responsive columns. Selected tool/sample tags
wrap into columns inside a 180-pixel scroll area. Narrow containers stack the
selectors. Tool and recovery-sample selectors provide Select all and Clear buttons.
All six global domain choices remain visible, including domains absent from this
truth set. Recovery defaults to `By sample`; `Total across selected samples`
produces one `Sample=Total` row per selected tool, summing qualifying bins across
selected observations. Strict AMBER completeness/purity thresholds and zero
counts remain unchanged. Empty sample/tool selections produce no rows. These
are bin counts; the overlap tab counts distinct sample/genome identities.

## Upgrade an accepted eager report

With the pinned requirements installed:

```sh
python refresh_cohort_report.py accepted/index.html accepted/index-interactive.html
```

The upgrader refuses an existing output path, preserves the original report,
and does not run scoring or QC. It accepts the Bokeh 3.8.2 cohort schema,
decodes the transport entities exactly once, matching Bokeh's standalone loader,
then checks every scientific `ColumnDataSource` against the accepted document, and
requires complete source coverage across domain views. Only the old toolbar's
single domain-state source is replaced. The adjacent receipt records input and
output SHA256, data equivalence, sizes and model counts. Keep the original
export receipt separate from the derivative presentation receipt.

## Focused validation

```sh
python -m unittest test.test_report_refresh test.test_recovered_genomes test.test_overlap test.test_cohort
python test/validate_cohort_browser.py accepted/index-interactive.html
```

The browser check needs development-only Selenium, Firefox and geckodriver.
It checks rendered metric tables and their dropdown callbacks, overlap selections
against independently counted records, plot legends, ranking sorting, Select all/Clear
buttons in both models, wide/narrow layout bounds, bounded scrolling, selected-sample totals,
empty selections, every tool/sample selection across populated and empty domain
views in both scoring models, model/tab/filter preservation, and release of old domain
views. It writes an adjacent browser receipt. Report-only validation does not
qualify a rebuilt production container or replace the paired pipeline microgate.
