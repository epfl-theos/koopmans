'''
Script for running the bulk Si example

Written by Edward Linscott, Mar 2021
'''

from koopmans import utils, io


if __name__ == '__main__':
    names = ['dscf', 'dfpt', 'dfpt_with_dscf_screening']

    wfs = {}
    for name in names:
        with utils.chdir(name):
            wf = io.read_json(f'si.json')
            wf.from_scratch = True
            if 'screening' in name:
                wf.calculate_alpha = False
                dscf_alphas = wfs['dscf'].bands.alphas
                wf.alpha_guess = dscf_alphas[len(dscf_alphas) // 2 - 4: len(dscf_alphas) // 2 + 4]
            wf.run()
        wfs[name] = wf
