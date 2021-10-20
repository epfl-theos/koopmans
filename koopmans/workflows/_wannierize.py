"""

Workflow module for koopmans, containing the workflow for generating maximally localized
Wannier functions (MLWFs), using the Wannier90 code

Written by Riccardo De Gennaro Nov 2020

"""

from pathlib import Path
from typing import List
from koopmans.pseudopotentials import nelec_from_pseudos
from ._generic import Workflow
from koopmans import utils
from ase.io.wannier90 import num_wann_from_projections
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


def list_to_formatted_str(values: List[int]):
    # Converts a list of integers into the format expected by Wannier90
    # e.g. list_to_formatted_str([1, 2, 3, 4, 5, 7]) = "1-5,7"
    assert all(a > b for a, b in zip(values[1:], values[:-1])), 'values must be monotonically increasing'
    indices = [None] + [i + 1 for i in range(len(values) - 1) if values[i + 1] != values[i] + 1] + [None]
    sectors = [values[slice(a, b)] for a, b in zip(indices[:-1], indices[1:])]
    out = []
    for sector in sectors:
        if len(sector) == 1:
            out.append(str(sector[0]))
        else:
            out.append(f'{sector[0]}-{sector[-1]}')
    return ','.join(out)


class WannierizeWorkflow(Workflow):

    def __init__(self, *args, nspin=1, check_wannierisation=None, **kwargs):
        super().__init__(*args, **kwargs)

        if check_wannierisation is not None:
            # In certain cases we never want to check the wannierisation, even if this is requested by the
            # JSON input file (e.g. with the smooth interpolation)
            self.parameters.check_wannierisation = check_wannierisation

        if 'pw' not in self.master_calc_params:
            raise ValueError(
                'You need to provide a pw block in your input when init_orbitals = "mlwfs" or "projwfs"')
        if 'w90_occ' not in self.master_calc_params and 'w90_emp' not in self.master_calc_params:
            raise ValueError(
                'You need to provide a w90 block in your input when init_orbitals = "mlwfs" or "projwfs"')

        # Make sure num_wann (occ/empty), num_bands (occ/empty), and nbnd are present and consistent
        pw_params = self.master_calc_params['pw']
        pw2w_params = self.master_calc_params['pw2wannier']
        w90_occ_params = self.master_calc_params['w90_occ']
        w90_emp_params = self.master_calc_params['w90_emp']

        # Filling out missing fields
        extra_core_bands = 0
        extra_conduction_bands = 0
        n_filled_bands = nelec_from_pseudos(self.atoms, self.pseudopotentials, pw_params.pseudo_dir) // 2
        for params in [w90_occ_params, w90_emp_params]:
            # Populate num_wann based on the projections provided
            if params.num_wann is None:
                assert 'projections' in params, 'Missing projections block for the Wannier90 calculation'
                params.num_wann = num_wann_from_projections(params.projections, self.atoms)
        if w90_occ_params.num_bands is None:
            # If num_bands has not been defined, this should just match the number of occupied Wannier functions
            # we want to define
            w90_occ_params.num_bands = w90_occ_params.num_wann
        if w90_emp_params.num_bands is None:
            # If num_bands has not been defined, this should just match the number of empty bands from the pw calculation
            w90_emp_params.num_bands = pw_params.nbnd - n_filled_bands
        if w90_occ_params.exclude_bands is None:
            # If exclude_bands hasn't been defined for the occupied calculation, this should exclude...
            exclude_bands = []
            extra_core_bands = n_filled_bands - w90_occ_params.num_bands
            if extra_core_bands > 0:
                # (a) the core bands if there are more filled bands than w90 num_bands, and
                exclude_bands += list(range(1, extra_core_bands + 1))
            # (b) all empty bands
            if n_filled_bands != pw_params.nbnd:
                exclude_bands += list(range(n_filled_bands + 1, pw_params.nbnd + 1))
            w90_occ_params.exclude_bands = list_to_formatted_str(exclude_bands)
        if w90_emp_params.exclude_bands is None:
            extra_conduction_bands = pw_params.nbnd - extra_core_bands - w90_occ_params.num_bands - w90_emp_params.num_bands
            # If exclude bands hasn't been defined for the empty calculation, this should exclude all filled bands
            w90_emp_params.exclude_bands = f'1-{n_filled_bands}'

        # Sanity checking
        w90_nbnd = w90_occ_params.num_bands + w90_emp_params.num_bands + extra_core_bands + extra_conduction_bands
        if pw_params.nbnd != w90_nbnd:
            raise ValueError(f'Number of bands disagrees between pw ({pw_params.nbnd}) and wannier90 ({w90_nbnd})')
        if w90_emp_params.num_wann == 0:
            raise ValueError('Cannot run a wannier90 calculation with num_wann = 0. Please set empty_states_nbnd > 0 '
                             'in the setup block, or num_wann > 0 in the wannier90 empty subblock')
        if nspin == 1:
            pw_params.nspin = 1
        else:
            pw_params.nspin = 2
            pw_params.tot_magnetization = 0.0
            pw2w_params.spin_component = 'up'

    def run(self):
        '''

        Wrapper for the calculation of (maximally localized) Wannier functions
        using PW and Wannier90

        '''

        self.print('Wannierisation', style='heading')

        if self.parameters.from_scratch:
            utils.system_call("rm -rf wannier", False)

        # Run PW scf and nscf calculations
        # PWscf needs only the valence bands
        calc_pw = self.new_calculator('pw')
        calc_pw.parameters.pop('nbnd', None)
        calc_pw.directory = 'wannier'
        calc_pw.prefix = 'scf'
        self.run_calculator(calc_pw)

        calc_pw = self.new_calculator('pw', calculation='nscf', nosym=True, noinv=True)
        calc_pw.directory = 'wannier'
        calc_pw.prefix = 'nscf'
        self.run_calculator(calc_pw)

        for typ in ['occ', 'emp']:
            # 1) pre-processing Wannier90 calculation
            calc_w90 = self.new_calculator('w90_' + typ, directory='wannier/' + typ)
            calc_w90.prefix = 'wann_preproc'
            calc_w90.command.flags = '-pp'
            self.run_calculator(calc_w90)
            utils.system_call(f'rsync -a {calc_w90.directory}/wann_preproc.nnkp {calc_w90.directory}/wann.nnkp')

            # 2) standard pw2wannier90 calculation
            calc_p2w = self.new_calculator('pw2wannier', directory=calc_w90.directory,
                                           outdir=calc_pw.parameters.outdir)
            calc_p2w.prefix = 'pw2wan'
            self.run_calculator(calc_p2w)

            # 3) Wannier90 calculation
            calc_w90 = self.new_calculator('w90_' + typ, directory='wannier/' + typ,
                                           bands_plot=self.parameters.check_wannierisation)
            calc_w90.prefix = 'wann'
            self.run_calculator(calc_w90)

        if self.parameters.check_wannierisation:
            # Run a "bands" calculation
            calc_pw = self.new_calculator('pw', calculation='bands', kpts=self.kpath)
            calc_pw.directory = 'wannier'
            calc_pw.prefix = 'bands'
            self.run_calculator(calc_pw)

            # Plot the bandstructures on top of one another
            ax = None
            labels = ['interpolation (occ)', 'interpolation (emp)', 'explicit']
            colour_cycle = plt.rcParams["axes.prop_cycle"]()
            selected_calcs = [c for c in self.calculations if 'band structure' in c.results]
            emin = np.min(selected_calcs[0].results['band structure'].energies) - 1
            emax = np.max(selected_calcs[1].results['band structure'].energies) + 1
            for calc, label in zip(selected_calcs, labels):
                if 'band structure' in calc.results:
                    # Load the bandstructure
                    bs = calc.results['band structure']

                    # Tweaking the plot aesthetics
                    colours = [next(colour_cycle)['color'] for _ in range(bs.energies.shape[0])]
                    kwargs = {}
                    if 'explicit' in label:
                        kwargs['ls'] = 'none'
                        kwargs['marker'] = 'x'

                    # Plot
                    ax = bs.plot(ax=ax, emin=emin, emax=emax, colors=colours, label=label, **kwargs)

            # Move the legend
            plt.legend(bbox_to_anchor=(1, 1), loc="lower right", ncol=2)

            # Save the comparison to file
            plt.savefig('interpolated_bandstructure_{}x{}x{}.png'.format(*self.kgrid))

        return

    def new_calculator(self, calc_type, *args, **kwargs):
        calc = super().new_calculator(calc_type, *args, **kwargs)

        # Extra tweaks for Wannier90 calculations
        if calc_type.startswith('w90'):
            if calc.parameters.num_bands != calc.parameters.num_wann and calc.parameters.dis_num_iter is None:
                calc.parameters.dis_num_iter = 5000
            if self.parameters.init_orbitals == 'projwfs':
                calc.parameters.num_iter = 0

        # Use a unified tmp directory
        if 'outdir' in calc.parameters.valid:
            calc.parameters.outdir = Path('wannier/TMP').resolve()

        return calc
