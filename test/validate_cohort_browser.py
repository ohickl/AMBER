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
        assert driver.execute_script('return cohortReport.error') is None
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
    time.sleep(.4)
    assert driver.execute_script("return checkTotals('official')") >= 0
    driver.execute_script("control('recovery_tools').value=[]")
    time.sleep(.2)
    assert driver.execute_script("return checkTotals('official')") == 0
    driver.execute_script("control('recovery_tools').value=control('recovery_tools').options;control('recovery_samples').value=control('recovery_samples').options")
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
    driver.execute_script('''
        cohortReport.root.active=1;
        cohortReport.root.tabs[1].child.active=4;
        control('recovery_mode','fractional').value='Total across selected samples';
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
    driver.execute_script("checkTotals('official');checkTotals('fractional')")
    for mask in (0,present_mask,4 & present_mask,present_mask):
        driver.execute_script('''const boxes=document.querySelectorAll('#domains input');
            boxes.forEach(box=>box.checked=!!(arguments[0]&(1<<Number(box.value))));
            boxes[0].dispatchEvent(new Event('change',{bubbles:true}));''', mask)
        ready(mask)
        time.sleep(.2)
        assert driver.execute_script('return Bokeh.documents.length') == 1
        assert driver.execute_script('return Object.keys(Bokeh.index).length') == 2
        driver.execute_script("checkTotals('official');checkTotals('fractional')")
    receipt = dict(schema='cohort-browser-v2',initial_seconds=initial_seconds,
                   timings=driver.execute_script('return cohortReport.timings'),
                   responsive=widths, totals='pass', empty_filters='pass', state_preserved='pass',live_documents=1)
    Path(sys.argv[1]).with_suffix('.browser.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'cohort-browser-v2 initial_seconds={initial_seconds:.2f} models=2 tabs=6 widths=1440,650 totals=pass domain_switch=pass live_documents=1')
