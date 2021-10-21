#! /usr/bin/env python3

import glob
import os
import copy
from pathlib import Path
from koopmans import utils
from koopmans.io import read, write
from koopmans.io import write_kwf as write_encoded_json

readin_exceptions = {''}


if __name__ == '__main__':

    # Run tests
    for test_json in sorted(glob.glob('test_??/test*.json')):
        test_json = Path(test_json)
        print(f'{test_json}...', end='', flush=True)
        wf = read(test_json)
        with utils.chdir(test_json.parent):
            wf.run()
        kwf_name = (test_json.parent / test_json.stem).with_suffix('.kwf')
        write(wf, kwf_name)
        print(' done')

    # Construct json file (missing input files)
    kwf_names = [Path(f) for f in sorted(glob.glob('test_??/*.kwf'))]

    data = {}
    for fname in kwf_names:
        fname_without_suffix = fname.parent / fname.stem

        wf = read(fname)

        for calc in wf.calculations:
            cname = str(calc.directory.relative_to(Path.cwd()) / calc.prefix)

            # Convert the SettingsDict to a plain dictionary to store in the json, making sure we
            # are going to store all paths as relative paths
            calc.parameters.use_relative_paths = True
            parameters = dict(calc.parameters)

            if cname in data:
                raise ValueError(f'Encountered a duplicate for {cname}')

            # Don't include walltime info
            calc.results.pop('walltime', None)

            # Here we will create a minimal representation of the calculator to store in the benchmarks
            data[cname] = {'parameters': parameters, 'results': calc.results}

            # Load alphas if do_orbdep = True and the alpha file is older than the input file
            if hasattr(calc, 'alphas'):
                data[cname]['alphas'] = calc.alphas

    with open('benchmarks.json', 'w') as f:
        write_encoded_json(data, f)
    os.system('cp benchmarks.json benchmarks_backup.json')

    # Add input files
    utils.system_call("sed -i 's/construct_exceptions = False/construct_exceptions = True/g' conftest.py")
    os.chdir('..')
    utils.system_call('pytest -m "mock" tests/test_?? --pdb')
    os.chdir('tests')
    utils.system_call("sed -i 's/construct_exceptions = True/construct_exceptions = False/g' conftest.py")
