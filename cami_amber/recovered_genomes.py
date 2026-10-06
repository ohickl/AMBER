"""Interactive recovery counts using AMBER's strict threshold convention."""

import math

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, DataTable, Div, MultiChoice, Slider, TableColumn


def recovery_records(pd_bins):
    """Keep empty groups and serialize undefined scores safely for JavaScript."""
    if pd_bins.empty:
        return []
    return [
        dict(sample=str(sample), tool=str(tool),
             completeness=float(recall) if math.isfinite(recall) else None,
             purity=float(precision) if math.isfinite(precision) else None)
        for sample, tool, recall, precision in pd_bins[
            ['sample_id', 'Tool', 'recall_bp', 'precision_bp']
        ].itertuples(index=False, name=None)
    ]


def count_recovered(records, completeness, contamination, samples, tools):
    """Count bins (including zero results), as in AMBER's existing table."""
    counts = {}
    for record in records:
        if record['sample'] not in samples or record['tool'] not in tools:
            continue
        key = record['sample'], record['tool']
        counts.setdefault(key, 0)
        if (record['completeness'] is not None and record['purity'] is not None
                and record['completeness'] > completeness / 100
                and record['purity'] > 1 - contamination / 100):
            counts[key] += 1
    return dict(Sample=[key[0] for key in counts], Tool=[key[1] for key in counts],
                Count=list(counts.values()))


RECOVERY_JS = """
const counts = new Map();
for (const record of records) {
    if (!samples.value.includes(record.sample) || !tools.value.includes(record.tool)) continue;
    const key = JSON.stringify([record.sample, record.tool]);
    if (!counts.has(key)) counts.set(key, {sample: record.sample, tool: record.tool, count: 0});
    if (record.completeness !== null && record.purity !== null &&
        record.completeness > completeness.value / 100 &&
        record.purity > 1 - contamination.value / 100) counts.get(key).count++;
}
source.data = {Sample: [], Tool: [], Count: []};
for (const value of counts.values()) {
    source.data.Sample.push(value.sample);
    source.data.Tool.push(value.tool);
    source.data.Count.push(value.count);
}
source.change.emit();
"""


def create_recovery_controls(pd_bins, min_completeness, max_contamination, groups=None):
    records = recovery_records(pd_bins)
    known = {(record['sample'], record['tool']) for record in records}
    for sample, tool in groups or []:
        if (str(sample), str(tool)) not in known:
            records.append(dict(sample=str(sample), tool=str(tool), completeness=None, purity=None))
    sample_ids = list(dict.fromkeys(record['sample'] for record in records))
    tool_ids = list(dict.fromkeys(record['tool'] for record in records))
    completeness = Slider(title='Completeness greater than (%)', start=0, end=100,
                          step=1, value=max(min_completeness) * 100, name='recovery_completeness')
    contamination = Slider(title='Contamination less than (%)', start=0, end=100,
                           step=1, value=min(max_contamination) * 100, name='recovery_contamination')
    samples = MultiChoice(title='Samples', options=sample_ids, value=sample_ids, name='recovery_samples')
    tools = MultiChoice(title='Tools', options=tool_ids, value=tool_ids, name='recovery_tools')
    source = ColumnDataSource(count_recovered(records, completeness.value,
                                             contamination.value, sample_ids, tool_ids))
    callback = CustomJS(args=dict(records=records, source=source, samples=samples,
                                  tools=tools, completeness=completeness,
                                  contamination=contamination), code=RECOVERY_JS)
    for control in (samples, tools, completeness, contamination):
        control.js_on_change('value', callback)
    table = DataTable(source=source, columns=[TableColumn(field=field, title=field)
                                             for field in ('Sample', 'Tool', 'Count')],
                      sizing_mode='stretch_width', height=500, index_position=None)
    return column(Div(text='Counts use strict AMBER thresholds. Undefined scores do not qualify. '
                           'Clear a filter to select no rows. Counts refer to qualifying bins.'),
                  row(completeness, contamination), row(samples, tools), table,
                  sizing_mode='stretch_width')
