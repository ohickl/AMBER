"""Firefox checks for lazy views, responsive selectors and recovery totals.

Usage: python test/validate_cohort_browser.py PATH/TO/index.html
Selenium and Firefox are development-only dependencies.
"""
from pathlib import Path
import json
import sys
import time
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.ui import WebDriverWait

options = Options()
options.add_argument('-headless')
with webdriver.Firefox(options=options) as driver:
    driver.set_window_size(1440, 1000)
    started = time.monotonic()
    driver.get(Path(sys.argv[1]).resolve().as_uri())
    def ready(mask=None):
        WebDriverWait(driver, 90).until(lambda d: d.execute_script(
            'return window.cohortReport && (cohortReport.error || (cohortReport.ready && (arguments[0]===null || cohortReport.mask===arguments[0])))', mask))
        assert driver.execute_script('return cohortReport.error') is None, driver.execute_script('return cohortReport.error')
    ready()
    initial_seconds = time.monotonic() - started
    present_mask = driver.execute_script('return cohortReport.mask')
    driver.execute_script('''
        window.control = (name,model='official')=>cohortReport.controls[model+'/'+name];
        window.findView = model => {
            const stack=[Bokeh.index[cohortReport.root.id]], seen=new Set();
            while(stack.length){const view=stack.pop();if(seen.has(view))continue;seen.add(view);
                if(view.model===model)return view;stack.push(...view.children());}
            throw Error('Control view not found: '+model.name);
        };
        window.checkTotals = model => {
            const samples=control('recovery_samples',model), tools=control('recovery_tools',model);
            const args=samples.js_property_callbacks['change:value'][0].args;
            const expected=new Map();
            for(const row of args.records){
                if(!samples.value.includes(row.sample)||!tools.value.includes(row.tool))continue;
                if(!expected.has(row.tool))expected.set(row.tool,0);
                if(row.completeness!==null && row.purity!==null && row.completeness>args.completeness.value/100 && row.purity>1-args.contamination.value/100)
                    expected.set(row.tool,expected.get(row.tool)+1);
            }
            const data=args.source.data;
            if(data.Count.length!==expected.size || data.Sample.some(s=>s!=='Total')) throw Error('Incorrect total rows');
            data.Tool.forEach((tool,i)=>{if(data.Count[i]!==expected.get(tool))throw Error('Incorrect total count');});
            return data.Count.reduce((a,b)=>a+b,0);
        };
        cohortReport.root.tabs[0].child.active=4;
        control('recovery_samples').value=control('recovery_samples').options.slice(0,3);
        control('recovery_tools').value=control('recovery_tools').options.slice(0,2);
        control('recovery_mode').value='Total across selected samples';
    ''')
    labels = driver.execute_script("return [...document.querySelectorAll('#domains label')].map(label=>label.textContent)")
    assert labels == ['Archaea','Bacteria','Viruses','Eukaryotes','Plasmids','Unknown'], labels
    for model in (0,1):
        driver.execute_script('cohortReport.root.active=arguments[0]',model)
        for tab in (0,2):
            driver.execute_script('cohortReport.root.tabs[arguments[0]].child.active=arguments[1]',model,tab)
            time.sleep(.25)
            driver.execute_script("""
                const panel=cohortReport.root.tabs[arguments[0]].child.tabs[arguments[1]].child;
                window.tableModel=[...panel.references()].find(item=>item.type==='Div' && item.text.includes('<table'));
                if(!tableModel)throw Error('Expected actual table markup');
                if(!findView(tableModel).shadow_el.querySelector('td'))throw Error('Table rendered as text');
                window.tableSelectors=[...panel.references()].filter(item=>item.type==='Select' &&
                    Object.values(item.js_property_callbacks).flat().some(callback=>callback.args.mytable===tableModel));
            """,model,tab)
            count = driver.execute_script('return tableSelectors.length')
            for index in range(count):
                driver.execute_script("""const select=tableSelectors[arguments[0]];
                    if(select.options.length>1){const option=select.options[1];select.value=Array.isArray(option)?option[0]:option;}""",index)
                time.sleep(.25)
                assert driver.execute_script("return !!findView(tableModel).shadow_el.querySelector('td')")
    driver.execute_script('cohortReport.root.active=0;cohortReport.root.tabs[0].child.active=4')
    time.sleep(.4)
    assert driver.execute_script("return checkTotals('official')") >= 0
    driver.execute_script("control('recovery_tools').value=[]")
    time.sleep(.2)
    assert driver.execute_script("return checkTotals('official')") == 0
    driver.execute_script("control('recovery_tools').value=control('recovery_tools').options;control('recovery_samples').value=control('recovery_samples').options")
    for model in ('official','fractional'):
        driver.execute_script("cohortReport.root.active=arguments[0];cohortReport.root.tabs[arguments[0]].child.active=5",0 if model=='official' else 1)
        driver.execute_script("""
            window.overlapArgs=control('overlap_tools',arguments[0]).js_property_callbacks['change:value'][0].args;
            window.checkOverlap=()=>{
                const a=overlapArgs, expected=new Map(a.tools.value.map(tool=>[tool,new Set()]));
                for(const row of a.records){
                    if(!expected.has(row.tool)||(a.sample.value!=='[sum over samples]'&&row.sample!==a.sample.value))continue;
                    if(row.completeness>a.completeness.value/100&&row.purity>1-a.contamination.value/100)
                        expected.get(row.tool).add(JSON.stringify([row.sample,row.genome]));
                }
                if(JSON.stringify(a.totals.data.Tool)!==JSON.stringify(a.tools.value))throw Error('Overlap tool selection did not update');
                a.totals.data.Tool.forEach((tool,i)=>{if(a.totals.data.Recovered[i]!==expected.get(tool).size)throw Error('Overlap recovery differs');});
                if(a.sample.value!=='[sum over samples]'&&a.genomes.data.Sample.some(sample=>sample!==a.sample.value))throw Error('Overlap sample selection did not update');
                return {sample:a.sample.value,tools:a.totals.data.Tool,rows:a.genomes.data.Sample.length,
                        venn:a.venn_plot.visible,upset:a.bar_plot.visible};
            };
            control('overlap_tools',arguments[0]).value=control('overlap_tools',arguments[0]).options.slice(0,2);
            control('overlap_completeness',arguments[0]).value=70;
        """,model)
        time.sleep(.3)
        result=driver.execute_script('return checkOverlap()')
        assert len(result['tools'])==2 and result['venn'] and not result['upset'],result
        driver.execute_script("control('overlap_sample',arguments[0]).value=control('overlap_sample',arguments[0]).options[1]",model)
        time.sleep(.3)
        result=driver.execute_script('return checkOverlap()')
        assert result['sample']!='[sum over samples]',result
        driver.execute_script("control('overlap_tools',arguments[0]).value=control('overlap_tools',arguments[0]).options",model)
        time.sleep(.3)
        result=driver.execute_script('return checkOverlap()')
        if len(result['tools'])>3:assert result['upset'] and not result['venn'],result
        driver.execute_script("control('overlap_tools',arguments[0]).value=[]",model)
        time.sleep(.3)
        assert driver.execute_script('return checkOverlap().rows')==0
        driver.execute_script("control('overlap_tools',arguments[0]).value=control('overlap_tools',arguments[0]).options;control('overlap_sample',arguments[0]).value=control('overlap_sample',arguments[0]).options[0]",model)
        # Plot legends use Bokeh's own hide interaction; exercise the real view.
        driver.execute_script("cohortReport.root.tabs[arguments[0]].child.active=1",0 if model=='official' else 1)
        time.sleep(.25)
        driver.execute_script("""
            const panel=cohortReport.root.tabs[arguments[0]].child.tabs[1].child;
            const legend=[...panel.references()].find(item=>item.type==='Legend'&&item.click_policy==='hide'&&item.items.length);
            if(!legend)throw Error('Missing interactive plot legend');
            window.legendView=findView(legend); window.legendRenderer=legend.items[0].renderers[0];
            window.previousVisible=legendRenderer.visible;
            legendView.entries[0].el.dispatchEvent(new PointerEvent('pointerdown',{bubbles:true}));
        """,0 if model=='official' else 1)
        time.sleep(.15)
        assert driver.execute_script('return legendRenderer.visible!==previousVisible')
        driver.execute_script("legendView.entries[0].el.dispatchEvent(new PointerEvent('pointerdown',{bubbles:true}))")
    for model_index in (0,1):
        driver.execute_script('cohortReport.root.active=arguments[0];cohortReport.root.tabs[arguments[0]].child.active=3',model_index)
        time.sleep(.2)
        driver.execute_script("""
            const panel=cohortReport.root.tabs[arguments[0]].child.tabs[3].child;
            const table=[...panel.references()].find(item=>item.type==='DataTable');
            window.rankingView=findView(table);
            window.rankingHeaders=rankingView.shadow_el.querySelectorAll('.slick-header-column');
            if(rankingHeaders.length<2)throw Error('Missing ranking headers');
            rankingHeaders[1].dispatchEvent(new MouseEvent('click',{bubbles:true}));
        """,model_index)
        time.sleep(.2)
        assert driver.execute_script("return !!rankingView.shadow_el.querySelector('.slick-header-column-sorted')")
    # Click actual Bokeh buttons, then verify the corresponding selection and data.
    for model_index, model in enumerate(('official','fractional')):
        for tab, names in ((4,('recovery_tools','recovery_samples')),(5,('overlap_tools',))):
            driver.execute_script('cohortReport.root.active=arguments[0];cohortReport.root.tabs[arguments[0]].child.active=arguments[1]',model_index,tab)
            time.sleep(.2)
            for name in names:
                for suffix in ('clear','all'):
                    driver.execute_script("""
                        const choice=control(arguments[0],arguments[1]);
                        const panel=cohortReport.root.tabs[arguments[1]==='official'?0:1].child;
                        const button=[...panel.references()].find(item=>item.name===arguments[0]+'_'+arguments[2]);
                        if(!button)throw Error('Missing selection button');
                        findView(button).shadow_el.querySelector('button').click();
                    """,name,model,suffix)
                    time.sleep(.2)
                    assert driver.execute_script("const choice=control(arguments[0],arguments[1]);return arguments[2]==='clear'?choice.value.length===0:choice.value.length===choice.options.length",name,model,suffix), (name,model,suffix)
    driver.execute_script('cohortReport.root.active=0')
    widths = []
    for width in (1440, 650):
        driver.set_window_size(width, 1000)
        for tab, names in ((4, ('recovery_tools','recovery_samples')), (5, ('overlap_tools','overlap_sample'))):
            driver.execute_script('cohortReport.root.tabs[0].child.active=arguments[0]', tab)
            time.sleep(.4)
            measured = driver.execute_script('''return arguments[0].map(name=>{
                const view=findView(control(name)),rect=view.el.getBoundingClientRect();
                const list=view.shadow_el.querySelector('.choices__list--multiple');
                return {name,x:rect.x,y:rect.y,width:rect.width,right:rect.right,viewport:innerWidth,
                        listHeight:list?.clientHeight,scrollHeight:list?.scrollHeight};
            })''', list(names))
            assert all(row['width'] > 0 and row['right'] <= row['viewport'] + 1 for row in measured), measured
            if width == 650:
                assert measured[1]['y'] > measured[0]['y'], measured
            assert all(row['listHeight'] is None or row['listHeight'] <= 180 for row in measured), measured
            widths.append(dict(width=width, tab=tab, controls=measured))
    # Give both models different subsets, including empty selections, and test
    # every sample/binner selector across populated and empty domain views.
    driver.execute_script('''
        window.expectedSelections={};
        for(const model of ['official','fractional']){
            control('recovery_tools',model).value=control('recovery_tools',model).options.slice(1,3);
            control('recovery_samples',model).value=model==='official'?control('recovery_samples',model).options.slice(1,3):[];
            control('overlap_tools',model).value=model==='official'?[]:control('overlap_tools',model).options.slice(0,2);
            const overlap=control('overlap_sample',model);overlap.value=overlap.options[overlap.options.length-1];
            const metrics=control('metrics/sample',model);metrics.value=metrics.options[metrics.options.length-1];
            const binner=control('metrics_per_bin/binner',model);
            if(binner){
                binner.value=binner.options[binner.options.length-1];
                const sample=control('metrics_per_bin/sample',model);
                const option=sample.options[sample.options.length-1];
                sample.value=Array.isArray(option)?option[0]:option;
            }
        }
        for(const [key,item] of Object.entries(cohortReport.controls)){
            if(item.type==='Select'||Array.isArray(item.value))expectedSelections[key]=item.value;
        }
        window.checkSelections=()=>{
            for(const [key,expected] of Object.entries(expectedSelections)){
                const item=cohortReport.controls[key];
                if(!item)continue;
                if(item.type==='Select'){
                    const options=item.options.map(option=>Array.isArray(option)?option[0]:option);
                    if(!options.includes(expected))continue;
                }
                if(JSON.stringify(item.value)!==JSON.stringify(expected))throw Error('Selection reset: '+key);
            }
        };
    ''')
    driver.execute_script('''
        cohortReport.root.active=1;
        cohortReport.root.tabs[1].child.active=4;
        control('recovery_mode','fractional').value='Total across selected samples';
        expectedSelections['fractional/recovery_mode']='Total across selected samples';
        control('overlap_contamination').value=35;
        window.selectedSamples=control('recovery_samples').value;
        const boxes=document.querySelectorAll('#domains input');
        boxes.forEach(box=>box.checked=box.value==='1');
        boxes[0].dispatchEvent(new Event('change',{bubbles:true}));
    ''')
    ready(2)
    state = driver.execute_script('''return {model:cohortReport.root.active,
        official:cohortReport.root.tabs[0].child.active, fractional:cohortReport.root.tabs[1].child.active,
        mode:control('recovery_mode').value,fractionalMode:control('recovery_mode','fractional').value,
        samples:control('recovery_samples').value,expected:window.selectedSamples,
        contamination:control('overlap_contamination').value,documents:Bokeh.documents.length}''')
    assert state == dict(model=1, official=5, fractional=4, mode='Total across selected samples',
                         fractionalMode='Total across selected samples', samples=state['expected'],
                         expected=state['expected'], contamination=35, documents=1), state
    time.sleep(.4)
    driver.execute_script("checkSelections()")
    driver.execute_script("checkTotals('official');checkTotals('fractional')")
    for mask in (0,present_mask,4 & present_mask,present_mask):
        driver.execute_script('''const boxes=document.querySelectorAll('#domains input');
            boxes.forEach(box=>box.checked=!!(arguments[0]&(1<<Number(box.value))));
            boxes[0].dispatchEvent(new Event('change',{bubbles:true}));''', mask)
        ready(mask)
        time.sleep(.2)
        assert driver.execute_script('return Bokeh.documents.length') == 1
        assert driver.execute_script('return Object.keys(Bokeh.index).length') == 2
        driver.execute_script("checkSelections();checkTotals('official');checkTotals('fractional')")
    # Every formerly absent selector must return with its original selection.
    assert driver.execute_script("return Object.keys(expectedSelections).every(key=>!!cohortReport.controls[key])")
    receipt = dict(schema='cohort-browser-v4',initial_seconds=initial_seconds, metric_tables='pass',domain_choices=6,overlap_data='pass',plot_legends='pass',ranking_sort='pass',selection_buttons='pass',
                   timings=driver.execute_script('return cohortReport.timings'),
                   responsive=widths, totals='pass', empty_filters='pass', state_preserved='all_selectors_both_models',live_documents=1)
    Path(sys.argv[1]).with_suffix('.browser.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'cohort-browser-v4 metric_tables=pass domain_choices=6 initial_seconds={initial_seconds:.2f} models=2 tabs=6 widths=1440,650 totals=pass domain_switch=pass live_documents=1')
