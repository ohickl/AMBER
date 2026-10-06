"""Recovered truth-genome intersections, keyed by sample rather than bin name."""

import math
from collections import defaultdict

import pandas as pd
from bokeh.layouts import column, row
from bokeh.models import (ColumnDataSource, CustomJS, DataTable, Div, FactorRange,
                          FixedTicker, MultiChoice, Select, Slider, TableColumn)
from bokeh.plotting import figure
from bokeh.events import DocumentReady
from cami_amber.report_controls import bounded_choices, responsive_filters


ALL_SAMPLES = '[sum over samples]'


def overlap_records(bins):
    result = []
    if bins.empty:
        return result
    for sample, tool, genome, recall, purity in bins[
            ['sample_id', 'Tool', 'genome_id', 'recall_bp', 'precision_bp']].itertuples(index=False, name=None):
        if pd.isna(genome) or not math.isfinite(recall) or not math.isfinite(purity):
            continue
        result.append(dict(sample=str(sample), tool=str(tool), genome=str(genome),
                           completeness=float(recall), purity=float(purity)))
    return result


def overlap_counts(records, tools, sample, completeness, contamination):
    memberships = defaultdict(set)
    for record in records:
        if (record['tool'] in tools and (sample == ALL_SAMPLES or record['sample'] == sample)
                and record['completeness'] > completeness / 100
                and record['purity'] > 1 - contamination / 100):
            memberships[record['sample'], record['genome']].add(record['tool'])
    intersections = defaultdict(int)
    counts, unique = dict.fromkeys(tools, 0), dict.fromkeys(tools, 0)
    genomes = []
    for (sample_id, genome), members in sorted(memberships.items()):
        pattern = ''.join('1' if tool in members else '0' for tool in tools)
        intersections[pattern] += 1
        for tool in members:
            counts[tool] += 1
            unique[tool] += len(members) == 1
        genomes.append(dict(Sample=sample_id, Genome=genome,
                            Tools='; '.join(tool for tool in tools if tool in members)))
    patterns = sorted(intersections, key=lambda pattern: (-intersections[pattern], pattern))
    return dict(tools=dict(Tool=tools, Recovered=[counts[t] for t in tools],
                           Unique=[unique[t] for t in tools]),
                intersections=dict(Pattern=patterns, Count=[intersections[p] for p in patterns],
                                   Tools=['; '.join(t for t, bit in zip(tools, p) if bit == '1') for p in patterns]),
                genomes={key: [r[key] for r in genomes] for key in ('Sample', 'Genome', 'Tools')})


OVERLAP_FUNCTION_JS = """
function overlapCounts(records, tools, sample, completeness, contamination) {
    const memberships = new Map();
    for (const record of records) {
        if (!tools.includes(record.tool) || (sample !== '[sum over samples]' && record.sample !== sample) ||
            record.completeness <= completeness / 100 || record.purity <= 1 - contamination / 100) continue;
        const key = JSON.stringify([record.sample, record.genome]);
        if (!memberships.has(key)) memberships.set(key, new Set());
        memberships.get(key).add(record.tool);
    }
    const intersections = new Map();
    const counts = new Map(tools.map(tool => [tool, 0]));
    const unique = new Map(tools.map(tool => [tool, 0]));
    const genomes = {Sample: [], Genome: [], Tools: []};
    const keys = [...memberships.keys()].sort((a, b) => {
        const aa = JSON.parse(a), bb = JSON.parse(b);
        return aa[0] < bb[0] ? -1 : aa[0] > bb[0] ? 1 : aa[1] < bb[1] ? -1 : aa[1] > bb[1] ? 1 : 0;
    });
    for (const key of keys) {
        const members = memberships.get(key);
        const [sampleId, genome] = JSON.parse(key);
        const pattern = tools.map(tool => members.has(tool) ? '1' : '0').join('');
        intersections.set(pattern, (intersections.get(pattern) || 0) + 1);
        for (const tool of members) {
            counts.set(tool, counts.get(tool) + 1);
            if (members.size === 1) unique.set(tool, unique.get(tool) + 1);
        }
        genomes.Sample.push(sampleId); genomes.Genome.push(genome);
        genomes.Tools.push(tools.filter(tool => members.has(tool)).join('; '));
    }
    const patterns = [...intersections.keys()].sort((a, b) => intersections.get(b) - intersections.get(a) || (a < b ? -1 : a > b ? 1 : 0));
    return {tools: {Tool: tools, Recovered: tools.map(tool => counts.get(tool)), Unique: tools.map(tool => unique.get(tool))},
            intersections: {Pattern: patterns, Count: patterns.map(p => intersections.get(p)),
                            Tools: patterns.map(p => tools.filter((tool, i) => p[i] === '1').join('; '))}, genomes};
}
"""


OVERLAP_UPDATE_JS = OVERLAP_FUNCTION_JS + """
const selected = tools.value;
const result = overlapCounts(records, selected, sample.value, completeness.value, contamination.value);
totals.data = result.tools; totals.change.emit();
intersections.data = result.intersections; intersections.change.emit();
genomes.data = result.genomes; genomes.change.emit();
const patterns = result.intersections.Pattern.slice(0, 30);
bars.data = {Pattern: patterns, Count: result.intersections.Count.slice(0, 30)}; bars.change.emit();
bar_plot.x_range.factors = patterns;
const dots = {Pattern: [], y: [], color: []};
const lines = {Pattern: [], ymin: [], ymax: []};
for (const pattern of patterns) {
    const enabled = [];
    for (let i = 0; i < selected.length; i++) {
        dots.Pattern.push(pattern); dots.y.push(i); dots.color.push(pattern[i] === '1' ? '#222222' : '#dddddd');
        if (pattern[i] === '1') enabled.push(i);
    }
    lines.Pattern.push(pattern); lines.ymin.push(Math.min(...enabled)); lines.ymax.push(Math.max(...enabled));
}
matrix.data = dots; matrix.change.emit(); connectors.data = lines; connectors.change.emit();
matrix_plot.y_range.start = -0.5; matrix_plot.y_range.end = Math.max(selected.length - 0.5, 0.5);
matrix_axis.ticker.ticks = selected.map((tool, i) => i);
matrix_axis.major_label_overrides = Object.fromEntries(selected.map((tool, i) => [i, tool]));
matrix_plot.height = Math.max(160, selected.length * 24 + 40);
const useVenn = selected.length === 2 || selected.length === 3;
venn_plot.visible = useVenn; bar_plot.visible = !useVenn; matrix_plot.visible = !useVenn;
if (useVenn) {
    const centers = selected.length === 2 ? [[-0.4, 0], [0.4, 0]] : [[-0.4, 0.3], [0.4, 0.3], [0, -0.4]];
    circles.data = {x: centers.map(p => p[0]), y: centers.map(p => p[1]), color: selected.map(tool => palette[tool])};
    circles.change.emit();
    const positions = selected.length === 2 ? {'10': [-0.75,0], '01': [0.75,0], '11': [0,0]} :
        {'100': [-0.75,0.6], '010': [0.75,0.6], '001': [0,-0.95], '110': [0,0.6],
         '101': [-0.45,-0.3], '011': [0.45,-0.3], '111': [0,0]};
    const regionCounts = new Map(result.intersections.Pattern.map((p, i) => [p, result.intersections.Count[i]]));
    const keys = Object.keys(positions);
    venn_counts.data = {x: keys.map(p => positions[p][0]), y: keys.map(p => positions[p][1]), text: keys.map(p => String(regionCounts.get(p) || 0))};
    venn_counts.change.emit();
    venn_names.data = {x: [-1.15,1.15,0].slice(0,selected.length), y: [1.35,1.35,-1.6].slice(0,selected.length), text: selected};
    venn_names.change.emit();
}
note.text = 'Recovered truth genomes: ' + result.genomes.Genome.length +
    '. Unique means recovered by one selected tool. ' + (useVenn ? 'Venn areas are schematic; region counts are exact.' :
    'UpSet shows the 30 largest exclusive intersections; the tables contain every intersection.');
"""


def create_overlap_panel(bins, summary, recovery_completeness, recovery_contamination, palette=None):
    records = overlap_records(bins)
    tool_ids = list(dict.fromkeys(summary['Tool']))
    from cami_amber.palette import method_palette
    palette = method_palette(tool_ids) if palette is None else palette
    for tool in tool_ids:
        palette[tool]
    samples = list(dict.fromkeys(summary['Sample']))
    tools = bounded_choices(MultiChoice(title='Tools', options=tool_ids, value=tool_ids, name='overlap_tools'))
    sample = Select(title='Sample', options=[ALL_SAMPLES] + samples, value=ALL_SAMPLES, name='overlap_sample')
    completeness = Slider(title='Completeness greater than (%)', start=0, end=100, step=1,
                          value=recovery_completeness.value, name='overlap_completeness')
    contamination = Slider(title='Contamination less than (%)', start=0, end=100, step=1,
                           value=recovery_contamination.value, name='overlap_contamination')
    for control, peer in [(completeness, recovery_completeness), (contamination, recovery_contamination)]:
        control.js_on_change('value', CustomJS(args=dict(peer=peer), code='if (peer.value !== cb_obj.value) peer.value = cb_obj.value;'))
        peer.js_on_change('value', CustomJS(args=dict(peer=control), code='if (peer.value !== cb_obj.value) peer.value = cb_obj.value;'))
    initial = overlap_counts(records, tool_ids, ALL_SAMPLES, completeness.value, contamination.value)
    totals = ColumnDataSource(initial['tools']); intersections = ColumnDataSource(initial['intersections'])
    genomes = ColumnDataSource(initial['genomes'])
    patterns = initial['intersections']['Pattern'][:30]
    bars = ColumnDataSource(dict(Pattern=patterns, Count=initial['intersections']['Count'][:30]))
    bar_plot = figure(x_range=FactorRange(factors=patterns), width=1000, height=300,
                      title='Exclusive intersection size', min_border_left=300)
    bar_plot.vbar(x='Pattern', top='Count', width=.8, source=bars)
    bar_plot.xaxis.visible = False
    matrix = ColumnDataSource(dict(Pattern=[p for p in patterns for _ in tool_ids],
                                   y=[i for _ in patterns for i in range(len(tool_ids))],
                                   color=['#222222' if bit == '1' else '#dddddd' for p in patterns for bit in p]))
    connectors = ColumnDataSource(dict(Pattern=patterns,
                                       ymin=[min(i for i, bit in enumerate(p) if bit == '1') for p in patterns],
                                       ymax=[max(i for i, bit in enumerate(p) if bit == '1') for p in patterns]))
    matrix_plot = figure(x_range=bar_plot.x_range, y_range=(-.5, max(len(tool_ids)-.5, .5)),
                         width=1000, height=max(160, len(tool_ids)*24+40), min_border_left=300)
    matrix_plot.segment(x0='Pattern', x1='Pattern', y0='ymin', y1='ymax', source=connectors, line_color='#222222')
    matrix_plot.scatter(x='Pattern', y='y', color='color', size=8, source=matrix)
    matrix_plot.yaxis.ticker = FixedTicker(ticks=list(range(len(tool_ids))))
    matrix_plot.yaxis.major_label_overrides = dict(enumerate(tool_ids))
    matrix_plot.xaxis.visible = False
    venn_plot = figure(x_range=(-2,2), y_range=(-2,2), width=1000, height=500, visible=False,
                       title='Shared and unique recovered genomes (schematic areas)')
    venn_plot.axis.visible = False; venn_plot.grid.visible = False
    circles = ColumnDataSource(dict(x=[], y=[], color=[]))
    venn_counts = ColumnDataSource(dict(x=[], y=[], text=[])); venn_names = ColumnDataSource(dict(x=[], y=[], text=[]))
    venn_plot.ellipse(x='x', y='y', width=1.9, height=1.9, fill_color='color', fill_alpha=.35, source=circles)
    venn_plot.text(x='x', y='y', text='text', source=venn_counts, text_align='center')
    venn_plot.text(x='x', y='y', text='text', source=venn_names, text_align='center', text_font_size='10pt')
    note = Div(text='Unique means recovered by one selected tool. Aggregate counts preserve sample identity. '
                    'UpSet shows the 30 largest exclusive intersections; tables contain all intersections.')
    callback = CustomJS(args=dict(palette=palette, records=records, tools=tools, sample=sample, completeness=completeness,
                                  contamination=contamination, totals=totals, intersections=intersections,
                                  genomes=genomes, bars=bars, matrix=matrix, connectors=connectors,
                                  bar_plot=bar_plot, matrix_plot=matrix_plot, matrix_axis=matrix_plot.yaxis[0], venn_plot=venn_plot,
                                  circles=circles, venn_counts=venn_counts, venn_names=venn_names, note=note), code=OVERLAP_UPDATE_JS)
    for control in (tools, sample, completeness, contamination):
        control.js_on_change('value', callback)
    tables = [DataTable(source=source, columns=[TableColumn(field=field, title=field) for field in fields],
                        sizing_mode='stretch_width', height=300, index_position=None)
              for source, fields in [(totals, ['Tool','Recovered','Unique']),
                                      (intersections, ['Tools','Count']), (genomes, ['Sample','Genome','Tools'])]]
    for plot in (venn_plot, bar_plot, matrix_plot):
        plot.sizing_mode = 'stretch_width'
        plot.min_width = 0
    layout = column(responsive_filters(tools, sample), responsive_filters(completeness, contamination), note, venn_plot, bar_plot, matrix_plot,
                  Div(text='Per-tool recovery and unique genomes'), tables[0],
                  Div(text='All exclusive intersections'), tables[1],
                  Div(text='Recovered sample–genome pairs'), tables[2], sizing_mode='stretch_width')
    layout.js_on_event(DocumentReady, callback)
    return layout
