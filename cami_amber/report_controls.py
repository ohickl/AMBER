"""Responsive, bounded filter controls shared by recovery and overlap."""

from bokeh.layouts import row


CHOICE_STYLE = """
:host { min-width: 0; width: 100% !important; max-width: 100%; }
.choices { min-width: 0; width: 100%; }
.choices__inner { box-sizing: border-box; min-width: 0; width: 100%; }
.choices__list--multiple {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(min(100%, 220px), 1fr));
  gap: 5px;
  max-height: 180px;
  overflow-y: auto;
  overflow-x: hidden;
  scrollbar-gutter: stable;
}
.choices__list--multiple .choices__item {
  box-sizing: border-box; min-width: 0; max-width: 100%; margin: 0;
  white-space: normal; overflow-wrap: anywhere;
}
.choices__input { box-sizing: border-box; max-width: 100%; }
.choices__list--dropdown .choices__list { max-height: 240px; overflow-y: auto; }
"""

TAB_STYLE = """
:host { min-width:0; width:100% !important; max-width:100%;
        grid-template-columns:minmax(0,1fr) !important; }
.bk-header { flex-wrap:wrap !important; min-width:0; }
"""

FILTER_STYLES = {'display': 'grid', 'gap': '16px', 'align-items': 'start',
                 'grid-template-columns': 'repeat(auto-fit, minmax(min(100%, 340px), 1fr))'}


def responsive_filters(*controls):
    """Fit two columns when space permits; stack within narrow containers."""
    for control in controls:
        control.sizing_mode = 'stretch_width'
        control.min_width = 0
        control.styles = {**control.styles, 'min-width': '0', 'width': '100%'}
        control.stylesheets = [*control.stylesheets, ':host {min-width:0; width:100% !important; max-width:100%;}']
    return row(*controls, sizing_mode='stretch_width', min_width=0,
               styles=FILTER_STYLES)


def bounded_choices(control):
    control.stylesheets = [*control.stylesheets, CHOICE_STYLE]
    return control


def responsive_panel(panel):
    """Constrain nested layout hosts as well as their individual controls."""
    from bokeh.models import Column, Row
    for model in panel.references():
        if isinstance(model, (Column, Row)):
            model.sizing_mode = 'stretch_width'
            model.min_width = 0
            model.stylesheets = [*model.stylesheets,
                                 ':host {min-width:0; width:100% !important; max-width:100%;}']
    return panel
