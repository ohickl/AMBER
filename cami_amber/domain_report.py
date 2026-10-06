"""One domain selector controls every genome-binning report tab."""

from bokeh.layouts import column
from bokeh.models import CheckboxGroup, ColumnDataSource, CustomJS, Div

from cami_amber.domain_metadata import DOMAINS


DOMAIN_JS = """
let mask = 0;
for (const index of domains.active) mask |= 1 << index;
mask &= present_mask;
const key = String(mask);
const previous = String(state.data.key[0]);
if (key === previous) return;
const old_panel = panels[previous];
const new_panel = panels[key];
new_panel.active = old_panel.active;
for (let i = 0; i < old_panel.tabs.length; i++) {
    const old_child = old_panel.tabs[i].child;
    const new_child = new_panel.tabs[i].child;
    if ('active' in old_child && 'active' in new_child) new_child.active = old_child.active;
}
for (const name of Object.keys(controls[previous])) {
    const old_control = controls[previous][name];
    const new_control = controls[key][name];
    const value = old_control.value;
    new_control.value = Array.isArray(value) ? [...value] : value;
}
view.children = [new_panel];
state.data = {key: [key]};
state.change.emit();
"""


def create_domain_report(profiles, labels, samples, options, make_panel):
    panels, controls = {}, {}
    for mask, (summary, bins) in profiles.items():
        key = str(mask)
        panels[key] = make_panel(summary, bins, labels, samples, options)
        controls[key] = {name: panels[key].select_one({'name': name})
                         for name in ('recovery_completeness', 'recovery_contamination',
                                      'recovery_samples', 'recovery_tools', 'overlap_tools', 'overlap_sample',
                                      'overlap_completeness', 'overlap_contamination')}
    return create_domain_selector(panels, controls, max(profiles))


def create_domain_selector(panels, controls, present_mask):
    initial = str(present_mask)
    view = column(panels[initial], sizing_mode='stretch_width')
    domains = CheckboxGroup(labels=list(DOMAINS), active=list(range(len(DOMAINS))),
                            name='global_domains', inline=True)
    state = ColumnDataSource(dict(key=[initial]))
    domains.js_on_change('active', CustomJS(args=dict(domains=domains, state=state, view=view,
                                                     panels=panels, controls=controls,
                                                     present_mask=present_mask), code=DOMAIN_JS))
    note = Div(text='Domain selection applies to all tabs. Bins keep their complete contamination '
                    'denominators and original matched genomes. Recall includes missing selected truth genomes. '
                    'Plasmids are separate; Unknown means no authoritative domain label. '
                    'In fractional views, compatible components touching a selected genome count once; '
                    'FARI excludes sequences spanning selected and excluded domains.', sizing_mode='stretch_width')
    return column(domains, note, view, sizing_mode='stretch_width')
