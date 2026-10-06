"""One HTML report, one domain selector, scientifically separate scoring models."""

from __future__ import annotations

from pathlib import Path
from bokeh.embed import json_item
from bokeh.models import TabPanel, Tabs

from cami_amber import amber_html
from cami_amber.lazy_report import encode_view, lazy_html
from cami_amber.report_controls import TAB_STYLE
from cami_amber.palette import method_palette, save_palette


CONTROLS = ('recovery_completeness', 'recovery_contamination', 'recovery_samples', 'recovery_tools', 'recovery_mode',
            'overlap_tools', 'overlap_sample', 'overlap_completeness', 'overlap_contamination')


def create_cohort_report(models: dict, output: Path) -> None:
    palette = method_palette(label for result in models.values() for label in result['labels'])
    save_palette(output, palette)
    amber_html.create_heatmap_bar(str(output))
    present = 0
    for result in models.values():
        present |= max(result['options'].domain_profiles)
    masks = [mask for mask in range(present + 1) if mask & ~present == 0]
    views = {}
    for mask in masks:
        pages, controls = [], {}
        for model, result in models.items():
            options = result['options']
            options.tool_palette = palette
            options.report_native_prefix = model
            # Each serialized view owns its models; retain no attached Bokeh cache.
            options.bin_view_cache = {}
            local_mask = mask & max(options.domain_profiles)
            summary, bins = options.domain_profiles[local_mask]
            panel = amber_html.create_genome_binning_html(
                summary, bins, result['labels'], result['samples'], options)
            pages.append(TabPanel(child=panel, title=model.capitalize()))
            for name in CONTROLS:
                control = panel.select_one({'name': name})
                if control is None:
                    raise ValueError(f'Missing cohort control: {model}/{name}')
                controls[model + '/' + name] = control.id
        root = Tabs(tabs=pages, sizing_mode='stretch_width', min_width=0, stylesheets=[TAB_STYLE])
        item = json_item(root)
        views[str(mask)] = encode_view(item['doc'], item['root_id'], controls)
        for result in models.values():
            result['options'].bin_view_cache = {}
    (output / 'index.html').write_text(lazy_html(views, present))
