"""Real Firefox regression for the self-contained cohort fixture report.

Usage: python test/validate_cohort_browser.py PATH/TO/index.html
Requires development-only selenium and Firefox; production has neither dependency.
"""
from pathlib import Path
import sys

from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.ui import WebDriverWait


options = Options()
options.add_argument('-headless')
with webdriver.Firefox(options=options) as driver:
    driver.get(Path(sys.argv[1]).resolve().as_uri())
    WebDriverWait(driver, 30).until(lambda d: d.execute_script('return typeof Bokeh !== "undefined" && Bokeh.documents.length > 0'))
    driver.execute_script('''
        const root = Bokeh.documents[0].roots()[0];
        window.selector = root.children[1];
        window.domains = selector.children[0];
        window.models = () => selector.children[2].children[0];
        window.page = (model=0) => models().tabs[model].child;
        window.control = (name, model=0) => [...page(model).references()].find(item => item.name === name);
        models().active = 1;
        page(0).active = 5;
        page(1).active = 5;
        window.expectedTools = control('overlap_tools').options.slice(0, 2);
        window.expectedSample = control('overlap_sample').options[1];
        control('overlap_tools').value = expectedTools;
        control('overlap_contamination').value = 35;
        control('overlap_sample').value = expectedSample;
        domains.active = [1];
    ''')
    actual = driver.execute_script('''
        return {model: models().active, officialTab: page(0).active, fractionalTab: page(1).active,
                sample: control('overlap_sample').value,
                contamination: control('overlap_contamination').value,
                tools: control('overlap_tools').value};
    ''')
    expected = driver.execute_script('return {tools: expectedTools, sample: expectedSample}')
    assert actual == dict(model=1, officialTab=5, fractionalTab=5, sample=expected['sample'], contamination=35,
                          tools=expected['tools']), actual
    driver.execute_script("domains.active = [];")
    assert driver.execute_script("return control('overlap_contamination').value") == 35
    print('cohort-browser-v1 models=2 tabs=6 domain_switch=pass controls_preserved=pass')
