#!/usr/bin/env python3
"""Upgrade an accepted cohort HTML to lazy views without rescoring any data."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

from cami_amber.recovered_genomes import RECOVERY_JS
from cami_amber.report_controls import CHOICE_STYLE, FILTER_STYLES, TAB_STYLE
from cami_amber.lazy_report import encode_view, lazy_html


def model_definitions(value: object) -> dict:
    models = {}
    def walk(item):
        if isinstance(item, dict):
            if item.get('type') == 'object' and 'id' in item and 'name' in item:
                models[item['id']] = item
            for child in item.values():
                walk(child)
        elif isinstance(item, list):
            for child in item:
                walk(child)
    walk(value)
    return models


def mapping(value: dict) -> dict:
    if value.get('type') != 'map':
        raise ValueError('Expected a serialized Bokeh map')
    return dict(value['entries'])


def packed(value: dict) -> dict:
    return dict(type='map', entries=list(value.items()))


def data_digest(models: dict) -> str:
    data = {key: model.get('attributes', {}).get('data') for key, model in models.items()
            if model['name'] == 'ColumnDataSource'}
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def refresh(source: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    html = source.read_text()
    matches = list(re.finditer(r'<script type="application/json"[^>]*>(.*?)</script>', html, re.S))
    if len(matches) != 1:
        raise ValueError('Expected one standalone Bokeh document')
    match = matches[0]
    payload = json.loads(match[1])
    if len(payload) != 1 or next(iter(payload.values())).get('version') != '3.8.2':
        raise ValueError('Expected the pinned Bokeh 3.8.2 report schema')
    models = model_definitions(payload)
    before = data_digest(models)
    originals = {key: json.dumps(model.get('attributes', {}).get('data'),sort_keys=True) for key, model in models.items() if model['name']=='ColumnDataSource'}
    named = {key: model for key, model in models.items() if model.get('attributes', {}).get('name')}
    modes = {}
    for key, model in named.items():
        attrs = model['attributes']
        if attrs['name'] != 'recovery_samples':
            continue
        callbacks = mapping(attrs['js_property_callbacks'])['change:value']
        callback = models[callbacks[0]['id']]
        mode_id = key + '_count_mode'
        mode = dict(type='object', name='Select', id=mode_id, attributes=dict(
            name='recovery_mode', title='Show counts', width=300,
            value='By sample', options=['By sample', 'Total across selected samples'],
            js_property_callbacks=packed({'change:value': [dict(id=callback['id'])]})))
        modes[key] = mode
        args = mapping(callback['attributes']['args'])
        args['mode'] = mode
        callback['attributes']['args'] = packed(args)
        callback['attributes']['code'] = RECOVERY_JS
    if not modes:
        raise ValueError('No recovery controls found')
    for model in list(models.values()):
        attrs = model.setdefault('attributes', {})
        children = attrs.get('children', [])
        if model['name'] == 'Tabs':
            attrs.update(sizing_mode='stretch_width', min_width=0)
            attrs['stylesheets'] = [*attrs.get('stylesheets', []), TAB_STYLE]
        if model['name'] == 'Row' and len(children) == 2:
            controls = [models[child['id']] for child in children]
            names = {control.get('attributes', {}).get('name') for control in controls}
            if names in ({'recovery_tools','recovery_samples'}, {'overlap_tools','overlap_sample'},
                         {'recovery_completeness','recovery_contamination'},
                         {'overlap_completeness','overlap_contamination'}):
                attrs.update(sizing_mode='stretch_width', min_width=0, styles=packed(FILTER_STYLES))
                for control in controls:
                    control['attributes'].update(sizing_mode='stretch_width', min_width=0)
                    control['attributes']['stylesheets'] = [*control['attributes'].get('stylesheets', []), ':host {min-width:0; width:100% !important; max-width:100%;}']
                if 'recovery_tools' in names or 'overlap_tools' in names:
                    attrs['children'] = sorted(children, key=lambda child:
                        0 if models[child['id']]['attributes']['name'].endswith('_tools') else 1)
        if model['name'] == 'Column':
            for index, child in enumerate(list(children)):
                row = models.get(child['id'], {})
                for item in row.get('attributes', {}).get('children', []):
                    if item['id'] in modes:
                        children.insert(index+1, dict(id=modes[item['id']]['id']))
                        break
        if attrs.get('name') in ('recovery_samples','recovery_tools','overlap_tools'):
            attrs['stylesheets'] = [*attrs.get('stylesheets', []), CHOICE_STYLE]
        if model['name'] == 'CustomJS' and 'args' in attrs:
            args = mapping(attrs['args'])
            if 'controls' in args and 'present_mask' in args:
                for panel_controls in mapping(args['controls']).values():
                    controls = mapping(panel_controls)
                    for key, value in list(controls.items()):
                        if key.endswith('recovery_samples'):
                            controls[key.replace('recovery_samples','recovery_mode')] = dict(id=modes[value['id']]['id'])
                    panel_controls.update(packed(controls))
            for key in ('bar_plot','matrix_plot','venn_plot'):
                if key in args:
                    models[args[key]['id']]['attributes'].update(sizing_mode='stretch_width', min_width=0)
        if model['name'] == 'TabPanel' and attrs.get('title') in ('#Recovered genomes','Overlap'):
            pending=[attrs['child']['id']]
            visited=set()
            while pending:
                key=pending.pop()
                if key in visited: continue
                visited.add(key)
                layout=models[key]
                if layout['name'] not in ('Column','Row'): continue
                properties=layout['attributes']
                properties.update(sizing_mode='stretch_width',min_width=0)
                properties['stylesheets']=[*properties.get('stylesheets',[]),
                    ':host {min-width:0; width:100% !important; max-width:100%;}']
                pending.extend(child['id'] for child in properties.get('children',[]))
    models.update({mode['id']: mode for mode in modes.values()})
    if before != data_digest(models):
        raise ValueError('Embedded source data changed during control refresh')
    domain = next(model for model in models.values() if model['name']=='CustomJS'
                  and 'present_mask' in mapping(model.get('attributes',{}).get('args',packed({}))))
    domain_args=mapping(domain['attributes']['args'])
    panels=mapping(domain_args['panels'])
    controls=mapping(domain_args['controls'])
    views, counts, covered = {}, {}, set()
    def document_for(root_id):
        seen=set()
        def expand(value):
            if isinstance(value, dict):
                if value.get('id') in models:
                    key=value['id']
                    if key in seen: return dict(id=key)
                    seen.add(key)
                    definition=models[key]
                    return {k:expand(v) for k,v in definition.items()}
                return {k:expand(v) for k,v in value.items()}
            if isinstance(value,list): return [expand(v) for v in value]
            return value
        document=dict(version='3.8.2',title='AMBER cohort comparison',roots=[expand(dict(id=root_id))])
        return document,seen
    for mask,panel in panels.items():
        models[panel['id']]['attributes'].update(sizing_mode='stretch_width',min_width=0)
        document,seen=document_for(panel['id'])
        scoped=model_definitions(document)
        for key,model in scoped.items():
            if model['name']=='ColumnDataSource':
                if json.dumps(model.get('attributes',{}).get('data'),sort_keys=True)!=originals[key]:
                    raise ValueError('View source data changed: '+key)
                covered.add(key)
        ids={name:control['id'] for name,control in mapping(controls[mask]).items()}
        views[mask]=encode_view(document,panel['id'],ids)
        counts[mask]=len(seen)
    excluded=set(originals)-covered
    # The global domain-state source is replaced by the HTML toolbar.
    if excluded != {domain_args['state']['id']}:
        raise ValueError('Scientific sources omitted: '+str(excluded))
    output.write_text(lazy_html(views,domain_args['present_mask']))
    receipt = dict(schema='amber-cohort-presentation-refresh-v2', source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                   output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(), embedded_data_sha256=before,
                   embedded_data_unchanged=True, recovery_views=len(modes), models=len(models),
                   live_models_by_domain=counts, scientific_sources_verified=len(covered),
                   original_bytes=source.stat().st_size, output_bytes=output.stat().st_size, rescored=0)
    output.with_suffix('.receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(f'cohort-presentation-refresh-v2 domains={len(views)} bytes={output.stat().st_size} embedded_data_unchanged=true rescored=0')
    return receipt


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args=parser.parse_args()
    refresh(args.source, args.output)
