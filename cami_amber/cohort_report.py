"""One HTML report, one domain selector, scientifically separate scoring models."""

from __future__ import annotations

from pathlib import Path
from bokeh.embed import file_html
from bokeh.layouts import column
from bokeh.models import Div, TabPanel, Tabs
from bokeh.resources import INLINE

from cami_amber import amber_html
from cami_amber.domain_report import create_domain_selector
from cami_amber.palette import method_palette, save_palette


CONTROLS = ('recovery_completeness', 'recovery_contamination', 'recovery_samples', 'recovery_tools',
            'overlap_tools', 'overlap_sample', 'overlap_completeness', 'overlap_contamination')


def create_cohort_report(models: dict, output: Path) -> None:
    palette = method_palette(label for result in models.values() for label in result['labels'])
    save_palette(output, palette)
    amber_html.create_heatmap_bar(str(output))
    present = 0
    for result in models.values():
        present |= max(result['options'].domain_profiles)
    masks = [mask for mask in range(present + 1) if mask & ~present == 0]
    panels, controls, model_panels = {}, {}, {}
    for model, result in models.items():
        options = result['options']
        options.tool_palette = palette
        options.report_native_prefix = model
        for mask, (summary, bins) in options.domain_profiles.items():
            model_panels[model, mask] = amber_html.create_genome_binning_html(
                summary, bins, result['labels'], result['samples'], options)
    for mask in masks:
        pages = []
        controls[str(mask)] = {}
        for model, result in models.items():
            local_mask = mask & max(result['options'].domain_profiles)
            panel = model_panels[model, local_mask]
            pages.append(TabPanel(child=panel, title=model.capitalize()))
            for name in CONTROLS:
                controls[str(mask)][model + '/' + name] = panel.select_one({'name': name})
        panels[str(mask)] = Tabs(tabs=pages)
    selector = create_domain_selector(panels, controls, present)
    title = Div(text='<h1>AMBER cohort comparison</h1><p>Every biological sample and assembly variant is '
                     'included. Official and fractional scoring remain separate. Method colours stay fixed '
                     'across samples and domain selections.</p>')
    html = file_html(column(title, selector, sizing_mode='stretch_width'), resources=INLINE,
                     title='AMBER cohort comparison', template=amber_html.TEMPLATE)
    (output / 'index.html').write_text(html)
