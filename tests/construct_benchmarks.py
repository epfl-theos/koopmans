#! /usr/bin/env python3

import glob
import os
from koopmans import utils
from koopmans.io import read
from koopmans.io import write_json as write_encoded_json
from koopmans.io import read_json as read_encoded_json

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
    fnames = [f for ext in ['cpi', 'pwi', 'win', 'p2wi', 'uii', 'w2ki', 'ksi', 'khi'] for f
              in glob.glob(f'test_??/**/*.{ext}', recursive=True)]

    fnames.sort(key=os.path.getmtime)

    data = {}
    for fname in fnames:
        print(fname)
        fname_without_ext, _ = fname.lstrip('./').rsplit('.', 1)
        calc = read(fname_without_ext)
        if fname_without_ext in data:
            raise ValueError(f'Encountered a duplicate for {fname}')

        settings = calc._settings
        for key in ['pseudo_dir', 'outdir', 'dft_ham_file', 'kc_ham_file', 'dft_smooth_ham_file', 'w90_seedname']:
            if key in settings:
                test_dir = fname.split('/')[0]
                if settings[key][0] == '/':
                    # Convert all absolute paths to be relative to the input file
                    fname_folder = fname.rsplit('/', 1)[0]
                    settings[key] = os.path.relpath(settings[key], os.getcwd() + '/' + fname_folder)

        results = {}
        for k, v in calc.results.items():
            if k == 'walltime':
                continue
            results[k] = v

        data[fname_without_ext] = {'settings': settings, 'results': results}

        # Load alphas if do_orbdep = True and the alpha file is older than the input file
        if getattr(calc, 'do_orbdep', False):
            directory, _ = fname_without_ext.rsplit('/', 1)
            if os.path.getctime(directory + '/' + 'file_alpharef.txt') < os.path.getctime(fname):
                data[fname_without_ext]['alphas'] = calc.read_alphas()[0]

    with open('benchmarks.json', 'w') as f:
        write_encoded_json(f, data)
    os.system('cp benchmarks.json benchmarks_backup.json')

    # Add input files
    utils.system_call("sed -i 's/construct_exceptions = False/construct_exceptions = True/g' conftest.py")
    os.chdir('..')
    utils.system_call('pytest -m "mock" tests/')
    os.chdir('tests')
    utils.system_call("sed -i 's/construct_exceptions = True/construct_exceptions = False/g' conftest.py")
