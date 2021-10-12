#! /usr/bin/env python3

import glob
import os
import copy
from pathlib import Path
from koopmans import utils
from koopmans.io import read
from koopmans.io import write_kwf as write_encoded_json

readin_exceptions = {''}


if __name__ == '__main__':

    # Run tests
    for test_json in sorted(glob.glob('test_??/test*.json')):
        print(test_json + '...', end='', flush=True)
        test_directory, test_json = test_json.rsplit('/', 1)
        with utils.chdir(test_directory):
            utils.system_call(f'koopmans {test_json} > {test_json.replace("json", "stdout")} 2>&1')
        print(' done')

    # Construct json file (missing input files)
    fnames = [Path(f) for ext in ['cpi', 'pwi', 'win', 'p2wi', 'uii', 'w2ki', 'ksi', 'khi'] for f
              in glob.glob(f'test_??/**/*.{ext}', recursive=True)]

    fnames.sort(key=os.path.getmtime)

    data = {}
    for fname in fnames:
        fname_without_suffix = fname.parent / fname.stem
        calc = read(fname_without_suffix)
        if fname_without_suffix in data:
            raise ValueError(f'Encountered a duplicate for {fname}')

        # Convert the SettingsDict to a plain dictionary to store in the json, making sure we
        # are going to store all paths as relative paths
        calc.parameters.use_relative_paths = True
        parameters = dict(calc.parameters)

        results = {}
        for k, v in calc.results.items():
            if k == 'walltime':
                continue
            results[k] = v

        data[str(fname_without_suffix)] = {'parameters': parameters, 'results': results}

        # Load alphas if do_orbdep = True and the alpha file is older than the input file
        if calc.parameters.get('do_orbdep', False):
            if os.path.getctime(fname.parent / 'file_alpharef.txt') < os.path.getctime(fname):
                data[str(fname_without_suffix)]['alphas'] = calc.read_alphas()[0]

    # Create dummy workflow object in order to write to file
    with open('benchmarks.json', 'w') as f:
        write_encoded_json(data, f)
    os.system('cp benchmarks.json benchmarks_backup.json')

    # Add input files
    utils.system_call("sed -i 's/construct_exceptions = False/construct_exceptions = True/g' conftest.py")
    os.chdir('..')
    utils.system_call('pytest -m "mock" tests/ --pdb')
    os.chdir('tests')
    utils.system_call("sed -i 's/construct_exceptions = True/construct_exceptions = False/g' conftest.py")
