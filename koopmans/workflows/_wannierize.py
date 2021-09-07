"""

Workflow module for koopmans, containing the workflow for generating maximally localized
Wannier functions (MLWFs), using the Wannier90 code

Written by Riccardo De Gennaro Nov 2020

"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from koopmans import utils
from ._generic import Workflow


class WannierizeWorkflow(Workflow):

    def __init__(self, workflow_settings, calcs_dct, nspin=1, check_wannierisation=None):
        super().__init__(workflow_settings, calcs_dct)

        if check_wannierisation is not None:
            # In certain cases we never want to check the wannierisation, even if this is requested by the
            # JSON input file (e.g. with the smooth interpolation)
            self.check_wannierisation = check_wannierisation

        if 'pw' not in self.master_calcs:
            raise ValueError(
                'You need to provide a pw block in your input when init_orbitals = "mlwfs" or "projwfs"')
        if 'w90_occ' not in self.master_calcs and 'w90_emp' not in self.master_calcs:
            raise ValueError(
                'You need to provide a w90 block in your input when init_orbitals = "mlwfs" or "projwfs"')

        # Make sure num_wann (occ/empty), num_bands (occ/empty), and nbnd are present and consistent
        pw_calc = self.master_calcs['pw']
        pw2w_calc = self.master_calcs['pw2wannier']
        w90_occ_calc = self.master_calcs['w90_occ']
        w90_emp_calc = self.master_calcs['w90_emp']

        # Filling out missing fields
        if w90_emp_calc.num_bands is None:
            w90_emp_calc.num_bands = pw_calc.nbnd - w90_occ_calc.num_bands
        if w90_occ_calc.exclude_bands is None:
            w90_occ_calc.exclude_bands = f'{w90_occ_calc.num_bands + 1}-{pw_calc.nbnd}'

        # Sanity checking
        w90_nbnd = w90_occ_calc.num_bands + w90_emp_calc.num_bands
        if pw_calc.nbnd != w90_nbnd:
            raise ValueError('Number of bands disagrees between pw ({pw_calc.nbnd}) and wannier90 ({w90_nbnd})')
        if w90_emp_calc.num_wann == 0:
            raise ValueError('Cannot run a wannier90 calculation with num_wann = 0. Please set empty_states_nbnd > 0 '
                             'in the setup block, or num_wann > 0 in the wannier90 empty subblock')
        if nspin == 1:
            pw_calc.nspin = 1
        else:
            pw_calc.nspin = 2
            pw_calc.tot_magnetization = 0.0
            pw2w_calc.spin_component = 'up'

    def run(self):
        '''

        Wrapper for the calculation of (maximally localized) Wannier functions
        using PW and Wannier90

        '''

        self.print('Wannierisation', style='heading')

        if self.from_scratch:
            utils.system_call("rm -rf wannier", False)

        # Run PW scf and nscf calculations
        # PWscf needs only the valence bands
        calc_pw = self.new_calculator('pw', nbnd=None)
        calc_pw.directory = 'wannier'
        calc_pw.name = 'scf'
        self.run_calculator(calc_pw)

        calc_pw = self.new_calculator('pw', calculation='nscf', nosym=True, noinv=True)
        calc_pw.directory = 'wannier'
        calc_pw.name = 'nscf'
        self.run_calculator(calc_pw)

        for typ in ['occ', 'emp']:
            # 1) pre-processing Wannier90 calculation
            calc_w90 = self.new_calculator('w90_' + typ, directory='wannier/' + typ, name='wann_preproc')
            calc_w90.calc.command.flags = '-pp'
            self.run_calculator(calc_w90)
            utils.system_call(f'rsync -a {calc_w90.directory}/wann_preproc.nnkp {calc_w90.directory}/wann.nnkp')

            # 2) standard pw2wannier90 calculation
            calc_p2w = self.new_calculator('pw2wannier', directory=calc_w90.directory,
                                           outdir=calc_pw.outdir, name='pw2wan', write_unk=self.check_wannierisation)
            self.run_calculator(calc_p2w)

            # 3) Wannier90 calculation
            calc_w90 = self.new_calculator('w90_' + typ, directory='wannier/' + typ, name='wann',
                                           wannier_plot=self.check_wannierisation, bands_plot=self.check_wannierisation)
            self.run_calculator(calc_w90)

        if self.check_wannierisation:
            # Run a "bands" calculation
            calc_pw = self.new_calculator('pw', calculation='bands')
            calc_pw.directory = 'wannier'
            calc_pw.name = 'bands'
            self.run_calculator(calc_pw)

            # Plot the bandstructures on top of one another
            ax = None
            labels = ['interpolation (occ)', 'interpolation (emp)', 'explicit']
            colour_cycle = plt.rcParams["axes.prop_cycle"]()
            selected_calcs = [c for c in self.all_calcs if 'band structure' in c.results]
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
            plt.savefig('interpolated_bandstructure_{}x{}x{}.png'.format(*calc_w90.mp_grid))

        return

    def new_calculator(self, calc_type, *args, **kwargs):
        calc = super().new_calculator(calc_type, *args, **kwargs)

        # Extra tweaks for Wannier90 calculations
        if calc_type.startswith('w90'):
            if calc.num_bands != calc.num_wann and calc.dis_num_iter is None:
                calc.dis_num_iter = 5000
            if self.init_orbitals == 'projwfs':
                calc.num_iter = 0

        # Use a unified tmp directory
        if hasattr(calc, 'outdir'):
            calc.outdir = os.path.abspath('wannier/TMP')

        return calc
