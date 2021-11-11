from koopmans import utils
calc_directory = 'screening'
utils.symlink(f'wannier/occ/wann_u.mat', f'{calc_directory}/', exist_ok=True)
