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
let activeViews = [], saved = null, busy = false, restoring = false, deferred = new Set(), desired = ''' + str(present_mask) + ''';
function saveState(){
 if (!report.root) return null;
 // Keep remembered values for controls absent from an empty domain view.
 const values = {...(saved?.values || {})};
 for (const [name, control] of Object.entries(report.controls)) {
  if (deferred.has(name)) continue;
  values[name] = Array.isArray(control.value) ? [...control.value] : control.value;
 }
 return {outer:report.root.active, inner:report.root.tabs.map(t=>t.child.active), values};
}
function registerControls(root){
 const slug = text=>text.toLowerCase().replace(/[^a-z0-9]+/g,'_').replace(/^_|_$/g,'');
 const registered = new Set(Object.values(report.controls));
 for(const model of root.tabs){
  for(const tab of model.child.tabs){
   for(const control of tab.child.references()){
    if(control.type !== 'Select' || registered.has(control)) continue;
    const key=slug(model.title)+'/'+slug(tab.title)+'/'+slug(control.title);
    if(report.controls[key]) throw Error('Ambiguous report control: '+key);
    report.controls[key]=control;
    registered.add(control);
   }
  }
 }
 for(const [name,control] of Object.entries(report.controls)){
  control.properties.value.change.connect(()=>{if(!restoring)deferred.delete(name);});
 }
}
async function restoreState(root){
 if(!saved) return;
 restoring=true;
 deferred=new Set();
 try{
  root.active=saved.outer;
  root.tabs.forEach((tab,i)=>tab.child.active=saved.inner[i]);
  // Bokeh CustomJS callbacks run asynchronously. Let binner callbacks rebuild
  // sample options before restoring samples, and suppress deferred-choice edits
  // until every restoration callback has settled.
  const entries=Object.entries(saved.values).sort(([a],[b])=>Number(b.endsWith('/binner'))-Number(a.endsWith('/binner')));
  for(const [name,value] of entries){
   const control=report.controls[name];
   if(!control) continue;
   if(control.type==='Select'){
    const options=control.options.map(option=>Array.isArray(option)?option[0]:option);
    if(!options.includes(value)){deferred.add(name);continue;}
   }
   control.value=Array.isArray(value)?[...value]:value;
   if(name.endsWith('/binner')) await new Promise(resolve=>setTimeout(resolve,0));
  }
  await new Promise(resolve=>setTimeout(resolve,0));
 }finally{restoring=false;}
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
   registerControls(root);
   await restoreState(root);
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
