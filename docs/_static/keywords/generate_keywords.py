from koopmans.settings import (ConvergenceSettingsDict, MLSettingsDict,
                               PlotSettingsDict,
                               UnfoldAndInterpolateSettingsDict,
                               WorkflowSettingsDict)


def print_settings_as_html(settings, table_id, fd):
    header = [f'<table id="{table_id}" style="width:100%; text-align:left">',
              '   <tr>',
              '      <th>Keyword</th>',
              '      <th>Description</th>',
              '      <th>Type</th>',
              '      <th>Default</th>',
              '   </tr>']

    fd.write('\n'.join(header))
    for s in settings.settings:

        # Parsing type
        if isinstance(s.kind, (tuple, list)):
            stype = '/'.join([f'<code>{t.__name__}</code>' for t in s.kind])
        else:
            stype = f'<code>{s.kind.__name__}</code>'

        # Parsing default
        if s.default is None:
            default = ''
        else:
            default = f'<code>{s.default}</code>'

        # Adding options if necessary
        if isinstance(s.options, (tuple, list)) and s.kind is not bool:
            default += ' (must be ' + '/'.join([f'<code>{o}</code>' for o in s.options]) + ')'

        fd.write(f'\n   <tr>')
        fd.write(f'\n      <td><code>{s.name}</code></td>')
        fd.write(f'\n      <td>{s.description.replace("<","&lt;").replace(">","&gt;")}</td>')
        fd.write(f'\n      <td>{stype}</td>')
        fd.write(f'\n      <td>{default}</td>')
        fd.write(f'\n   </tr>')
    fd.write('\n</table>')


if __name__ == '__main__':

    with open('workflow_keywords.html', 'w') as fd:
        print_settings_as_html(WorkflowSettingsDict(), 'workflowTable', fd)

    with open('ui_keywords.html', 'w') as fd:
        print_settings_as_html(UnfoldAndInterpolateSettingsDict(), 'uiTable', fd)

    with open('plot_keywords.html', 'w') as fd:
        print_settings_as_html(PlotSettingsDict(), 'plotTable', fd)

    with open('ml_keywords.html', 'w') as fd:
        print_settings_as_html(MLSettingsDict(), 'mlTable', fd)

    with open('convergence_keywords.html', 'w') as fd:
        print_settings_as_html(ConvergenceSettingsDict(), 'convTable', fd)
