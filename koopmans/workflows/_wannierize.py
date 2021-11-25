"""

Workflow module for koopmans, containing the workflow for generating maximally localized
Wannier functions (MLWFs), using the Wannier90 code

Written by Riccardo De Gennaro Nov 2020

"""

import itertools
import pickle
from ._generic import Workflow
from koopmans import utils
from koopmans.pseudopotentials import nelec_from_pseudos
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import List, Union
import matplotlib
matplotlib.use('Agg')


class WannierizeWorkflow(Workflow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

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
        if w90_occ_params.num_bands is None:
            # If num_bands has not been defined, this should just match the number of occupied Wannier functions
            # we want to define
            w90_occ_params.num_bands = w90_occ_params.num_wann
        if w90_emp_params.num_bands is None:
            # If num_bands has not been defined, this should just match the number of empty bands from the pw
            # calculation
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
            if len(exclude_bands) > 0:
                w90_occ_params.exclude_bands = utils.list_to_formatted_str(exclude_bands)
        if w90_emp_params.exclude_bands is None:
            extra_conduction_bands = pw_params.nbnd - extra_core_bands - w90_occ_params.num_bands \
                - w90_emp_params.num_bands
            # If exclude bands hasn't been defined for the empty calculation, this should exclude all filled bands
            w90_emp_params.exclude_bands = f'1-{n_filled_bands}'

        # Update the projections_blocks to account for additional bands
        self.parameters.w90_projections_blocks.add_excluded_bands(
            w90_occ_params.num_bands + extra_core_bands - w90_occ_params.num_wann, above=False)
        self.parameters.w90_projections_blocks.add_excluded_bands(
            w90_emp_params.num_bands - w90_emp_params.num_wann, above=True)

        # Sanity checking
        w90_nbnd = w90_occ_params.num_bands + w90_emp_params.num_bands + extra_core_bands + extra_conduction_bands
        if pw_params.nbnd != w90_nbnd:
            raise ValueError(f'Number of bands disagrees between pw ({pw_params.nbnd}) and wannier90 ({w90_nbnd})')

        # Spin-polarisation
        if self.parameters.spin_polarised:
            pw_params.nspin = 2
        else:
            pw_params.nspin = 1

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

        fillings = []
        for filling in ['occ', 'emp']:
            if self.master_calc_params['w90_' + filling].num_wann > 0:
                fillings.append(filling)

        if self.parameters.spin_polarised:
            spins = ['up', 'down']
        else:
            spins = [None]

        for filling, spin in itertools.product(fillings, spins):
            # First, construct any extra settings in case we are Wannierising block-by-block
            filled = (filling == 'occ')
            block_projs = [b for b in self.parameters.w90_projections_blocks if b.filled == filled]
            if len(block_projs) == 1:
                # We are not Wannierising block-by-block; don't provide any extra arguments
                blocks_kwargs = [{}]
            else:
                # We are Wannierising block-by-block; provide extra arguments that will overwrite projections,
                # num_wann, and exclude_bands appropriately
                blocks_kwargs = [b.kwargs for b in block_projs]

            # Loop over each block that will be Wannierised separately
            for i_block, block_kwargs in enumerate(blocks_kwargs):
                # Construct the subdirectory label
                if spin is None:
                    typ = filling
                else:
                    typ = f'{filling}_{spin}'
                if block_kwargs != {}:
                    typ += f'_block{i_block + 1}'

                # 1) pre-processing Wannier90 calculation
                calc_w90 = self.new_calculator('w90_' + filling, directory='wannier/'
                                               + typ, spin_component=spin, **block_kwargs)
                calc_w90.prefix = 'wann_preproc'
                calc_w90.command.flags = '-pp'
                self.run_calculator(calc_w90)
                utils.system_call(f'rsync -a {calc_w90.directory}/wann_preproc.nnkp {calc_w90.directory}/wann.nnkp')

                # 2) standard pw2wannier90 calculation
                calc_p2w = self.new_calculator('pw2wannier', directory=calc_w90.directory,
                                               outdir=calc_pw.parameters.outdir)
                calc_p2w.prefix = 'pw2wan'
                if spin is not None:
                    calc_p2w.spin_component = spin
                self.run_calculator(calc_p2w)

                # 3) Wannier90 calculation
                calc_w90 = self.new_calculator('w90_' + filling, directory='wannier/' + typ,
                                               bands_plot=self.parameters.check_wannierisation, spin_component=spin, **block_kwargs)
                calc_w90.prefix = 'wann'
                self.run_calculator(calc_w90)

        if self.parameters.check_wannierisation:
            # Run a "bands" calculation, making sure we don't overwrite the scf/nscf tmp files by setting a different prefix
            calc_pw_bands = self.new_calculator('pw', calculation='bands', kpts=self.kpath)
            calc_pw_bands.directory = 'wannier'
            calc_pw_bands.prefix = 'bands'
            calc_pw_bands.parameters.prefix += '_bands'
            # Link the save directory so that the bands calculation can use the old density
            if self.parameters.from_scratch:
                [src, dest] = [(c.parameters.outdir / c.parameters.prefix).with_suffix('.save')
                               for c in [calc_pw, calc_pw_bands]]
                utils.symlink(src, dest)
            self.run_calculator(calc_pw_bands)

            # Select those calculations that generated a band structure
            selected_calcs = [c for c in self.calculations[:-1] if 'band structure' in c.results]

            # Work out the vertical shift to set the valence band edge to zero
            w90_emp_num_bands = self.master_calc_params['w90_emp'].num_bands
            if w90_emp_num_bands > 0:
                vbe = np.max(calc_pw_bands.results['band structure'].energies[:, :, :-w90_emp_num_bands])
            else:
                vbe = np.max(calc_pw_bands.results['band structure'].energies)

            # Work out the energy ranges for plotting
            emin = np.min(selected_calcs[0].results['band structure'].energies) - 1 - vbe
            emax = np.max(selected_calcs[-1].results['band structure'].energies) + 1 - vbe

            # Plot the bandstructures on top of one another
            ax = None
            labels = ['explicit'] \
                + [f'interpolation ({c.directory.name.replace("_",", ").replace("block", "block ")})'
                   for c in selected_calcs]
            colour_cycle = plt.rcParams["axes.prop_cycle"]()
            for calc, label in zip([calc_pw_bands] + selected_calcs, labels):
                if 'band structure' in calc.results:
                    # Load the bandstructure
                    bs = calc.results['band structure']

                    # Unfortunately once a bandstructure object is created you cannot tweak it, so we must alter
                    # this private variable
                    bs._energies -= vbe

                    # Tweaking the plot aesthetics
                    colours = [next(colour_cycle)['color'] for _ in range(bs.energies.shape[0])]
                    kwargs = {}
                    if 'explicit' in label:
                        kwargs['ls'] = 'none'
                        kwargs['marker'] = 'x'

                    # Plot
                    ax = bs.plot(ax=ax, emin=emin, emax=emax, colors=colours, label=label, **kwargs)

            # Move the legend
            lgd = ax.legend(bbox_to_anchor=(1, 1), loc="lower right", ncol=2)

            # Save the comparison to file (as png and also in editable form)
            with open('interpolated_bandstructure_{}x{}x{}.fig.pkl'.format(*self.kgrid), 'wb') as fd:
                pickle.dump(plt.gcf(), fd)
            # The "bbox_extra_artists" and "bbox_inches" mean that the legend is not cropped out
            plt.savefig('interpolated_bandstructure_{}x{}x{}.png'.format(*self.kgrid),
                        bbox_extra_artists=(lgd,), bbox_inches='tight')

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
