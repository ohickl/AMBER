"""Standalone, compressed domain views with one live Bokeh document at a time."""
from __future__ import annotations

import base64
import gzip
import json
from bokeh.resources import INLINE
from cami_amber.domain_metadata import DOMAINS


def encode_view(document, root_id, controls):
    raw = json.dumps(dict(doc=document, root_id=root_id, controls=controls),
                     separators=(',', ':'), ensure_ascii=False).encode()
    return base64.b64encode(gzip.compress(raw, compresslevel=6, mtime=0)).decode()


def lazy_html(views, present_mask):
    controls = ''.join(f'<label><input type="checkbox" value="{i}" checked>{name}</label>'
                       for i, name in enumerate(DOMAINS))
    return '''<!DOCTYPE html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>AMBER cohort comparison</title>''' + INLINE.render() + '''
<style>body{font:15px system-ui,sans-serif;margin:16px;color:#263238}
#domains{display:flex;flex-wrap:wrap;gap:18px;margin:18px 0}label{cursor:pointer}
#status{min-height:24px;color:#52616b}#cohort-view{width:100%;min-width:0;overflow-x:auto}
h1{font-size:26px}p{max-width:1100px;line-height:1.5}</style></head><body>
<h1>AMBER cohort comparison</h1><p>Every biological sample and assembly variant is included.
Official and fractional scoring remain separate. Method colours stay fixed across selections.</p>
<div id="domains">''' + controls + '''</div>
<p>Domain selection applies to all tabs. Bins keep their complete contamination denominators and
original matched genomes. Recall includes missing selected truth genomes. Plasmids are separate.
In fractional views, compatible components touching a selected genome count once; FARI excludes
sequences spanning selected and excluded domains.</p><div id="status" role="status">Loading selected view…</div>
<div id="cohort-view"></div><script>
const bundles = ''' + json.dumps(views, separators=(',', ':')) + ''';
const report = window.cohortReport = {ready:false, mask:null, root:null, controls:{}, timings:[], error:null};
let activeViews = [], saved = null, busy = false, desired = ''' + str(present_mask) + ''';
function saveState(){
 if (!report.root) return null;
 const values = {};
 for (const [name, control] of Object.entries(report.controls)) values[name] = Array.isArray(control.value) ? [...control.value] : control.value;
 return {outer:report.root.active, inner:report.root.tabs.map(t=>t.child.active), values};
}
async function decode(encoded){
 const bytes=Uint8Array.from(atob(encoded), c=>c.charCodeAt(0));
 const stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
 return JSON.parse(await new Response(stream).text());
}
async function load(){
 if (busy) return;
 busy=true;
 try {
  while (report.mask !== desired) {
   const mask=desired, started=performance.now();
   report.ready=false;
   document.getElementById('status').textContent='Loading selected view…';
   saved=saveState() || saved;
   for (const view of activeViews) view.remove();
   activeViews=[];
   if (report.root) {
    const doc=report.root.document;
    doc.clear();
    const index=Bokeh.documents.indexOf(doc);
    if(index>=0) Bokeh.documents.splice(index,1);
   }
   report.root=null; report.controls={};
   document.getElementById('cohort-view').replaceChildren();
   const item=await decode(bundles[String(mask)]);
   const views=await Bokeh.embed.embed_item(item,'cohort-view');
   activeViews=Array.from(views);
   const root=activeViews[0].model, doc=root.document;
   report.root=root;
   report.controls=Object.fromEntries(Object.entries(item.controls).map(([k,id])=>[k,doc.get_model_by_id(id)]));
   if(saved){
    root.active=saved.outer;
    root.tabs.forEach((tab,i)=>tab.child.active=saved.inner[i]);
    for(const [name,value] of Object.entries(saved.values)) if(report.controls[name]) report.controls[name].value=value;
   }
   report.mask=mask;
   report.timings.push({mask,milliseconds:performance.now()-started,models:doc._all_models.size});
  }
  report.ready=true;
  document.getElementById('status').textContent='Selected view ready';
 } catch(error){report.error=String(error);document.getElementById('status').textContent='Unable to load view: '+error;console.error(error);}
 finally {busy=false;}
}
document.getElementById('domains').addEventListener('change',()=>{
 desired=Array.from(document.querySelectorAll('#domains input:checked')).reduce((mask,input)=>mask | (1<<Number(input.value)),0) & ''' + str(present_mask) + ''';
 load();
});
load();
</script></body></html>'''
