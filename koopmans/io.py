"""

I/O module for python_KI

Written by Edward Linscott Jan 2020

"""


def cpi_diff(calcs):

    # If calcs is a dict, convert it to a list (we only need the values)
    if isinstance(calcs, dict):
        calcs = calcs.values()

    diffs = []

    blocks = set([b for c in calcs for b in c.parameters['input_data'].keys()])
    for block in sorted(blocks):
        keys = set(
            [k for c in calcs for k in c.parameters['input_data'].get(block, {}).keys()])
        for key in sorted(keys):
            vals = [c.parameters['input_data']
                    [block].get(key, None) for c in calcs]
            if len(set(vals)) > 1:
                print(f'{block}.{key}: ' + ', '.join(map(str, vals)))
                diffs.append(key)

    return diffs


def print_summary(alpha_df, error_df):
    # Printing out a progress summary
    print('\nalpha')
    print(alpha_df)
    print('\ndelta E - lambda^alpha_ii (eV)')
    print(error_df)
    print('')


def write_alpharef(alphas, filling, directory='.'):
    '''
    Generates file_alpharef.txt and file_alpharef_empty.txt from a list of alpha values

    Arguments:
       alphas    -- a list of alpha values (floats)
       filling   -- a list of booleans; true if the corresponding orbital is filled
       directory -- the directory within which to write the files
    '''

    a_filled = [a for a, f in zip(alphas, filling) if f]
    a_empty = [a for a, f in zip(alphas, filling) if not f]
    for alphas, suffix in zip([a_filled, a_empty], ['', '_empty']):
        with open(f'{directory}/file_alpharef{suffix}.txt', 'w') as fd:
            fd.write('{}\n'.format(2*len(alphas)))
            fd.writelines(['{} {} 1.0\n'.format(i+1, a)
                           for i, a in enumerate(alphas)])
            fd.writelines(['{} {} 1.0\n'.format(i+1+len(alphas), a)
                           for i, a in enumerate(alphas)])


def read_alpharef(calc=None, directory=None):
    '''
    Reads in file_alpharef.txt and file_alpharef_empty.txt from a calculation's directory

    Arguments:
       calc      -- an ASE calculator
       directory -- a directory

    Output:
       alphas -- a list of alpha values (1 per orbital)
    '''

    if calc is not None:
        directory = calc.label.rsplit('/', 1)[0]
    elif directory is None:
        raise ValueError(
            'read_alpharef called without a calculator or a directory. Please provide at least one.')

    alphas = []
    for suffix in ['', '_empty']:
        with open(f'{directory}/file_alpharef{suffix}.txt', 'r') as fd:
            flines = fd.readlines()
            n_orbs = int(flines[0]) // 2
            alphas += [float(line.split()[1]) for line in flines[1:n_orbs + 1]]
    return alphas
