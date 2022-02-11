"""

Workflow module for koopmans, containing the workflow for generating maximally localized
Wannier functions (MLWFs), using the Wannier90 code

Written by Riccardo De Gennaro Nov 2020

"""

import numpy as np
import math
from pathlib import Path
import pickle
import shutil
from typing import List, TypeVar
import koopmans.mpl_config
import matplotlib.pyplot as plt
from koopmans import utils, projections, calculators
from koopmans.pseudopotentials import nelec_from_pseudos
from ._workflow import Workflow


CalcExtType = TypeVar('CalcExtType', bound='calculators.CalculatorExt')


class WannierizeWorkflow(Workflow):

    def __init__(self, *args, force_nspin2=False, scf_kgrid=None, **kwargs):
        super().__init__(*args, **kwargs)

        if 'pw' not in self.master_calc_params:
            raise ValueError(
                'You need to provide a pw block in your input')

        pw_params = self.master_calc_params['pw']

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:

            if self.parameters.spin_polarised:
                spins = [0, 1]
            else:
                spins = [None]

            for spin in spins:
                # Update the projections_blocks to account for additional occupied bands
                num_wann_occ = self.projections.num_bands(occ=True, spin=spin)
                nelec = nelec_from_pseudos(self.atoms, self.pseudopotentials, pw_params.pseudo_dir)
                if self.parameters.spin_polarised:
                    num_bands_occ = nelec
                    if spin == 0:
                        num_bands_occ += pw_params.tot_magnetization
                    else:
                        num_bands_occ -= pw_params.tot_magnetization
                    num_bands_occ //= 2
                else:
                    num_bands_occ = nelec // 2
                self.projections.add_bands(num_bands_occ - num_wann_occ, above=False, spin=spin)

                # Update the projections_blocks to account for additional empty bands
                num_wann_emp = self.projections.num_bands(occ=False, spin=spin)
                num_bands_emp = pw_params.nbnd - num_bands_occ
                self.projections.add_bands(num_bands_emp - num_wann_emp, above=True, spin=spin)

                # Sanity checking
                w90_nbnd = self.projections.num_bands(spin=spin)
                if pw_params.nbnd != w90_nbnd:
                    raise ValueError(
                        f'Number of bands disagrees between pw ({pw_params.nbnd}) and wannier90 ({w90_nbnd})')

        elif self.parameters.init_orbitals == 'kohn-sham':
            pass

        else:
            raise NotImplementedError('WannierizeWorkflow supports only init_orbitals = '
                                      '"mlwfs", "projwfs" or "kohn-sham"')

        # Spin-polarisation
        self._force_nspin2 = force_nspin2
        if force_nspin2:
            pw_params.nspin = 2
        elif self.parameters.spin_polarised:
            pw_params.nspin = 2
        else:
            pw_params.nspin = 1

        # When running a smooth PW calculation the k-grid for the scf calculation
        # must match the original k-grid
        self._scf_kgrid = scf_kgrid

    def run(self):
        '''

        Wrapper for the calculation of (maximally localized) Wannier functions
        using PW and Wannier90

        '''
        if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
            self.print('Wannierisation', style='heading')
        else:
            self.print('Kohn-Sham orbitals', style='heading')

        if self.parameters.from_scratch:
            utils.system_call("rm -rf wannier", False)

        # Run PW scf and nscf calculations
        # PWscf needs only the valence bands
        calc_pw = self.new_calculator('pw')
        calc_pw.parameters.pop('nbnd', None)
        calc_pw.directory = 'wannier'
        calc_pw.prefix = 'scf'
        if self._scf_kgrid:
            calc_pw.parameters.kpts = self._scf_kgrid
        self.run_calculator(calc_pw)

        calc_pw = self.new_calculator('pw', calculation='nscf', nosym=True, noinv=True)
        calc_pw.directory = 'wannier'
        calc_pw.prefix = 'nscf'
        self.run_calculator(calc_pw)

        if self.parameters.init_orbitals in ['mlwfs', 'projwfs']:
            # Loop over the various subblocks that we must wannierise separately
            for block in self.projections:
                # Construct the subdirectory label
                w90_dir = Path('wannier') / block.directory

                # 1) pre-processing Wannier90 calculation
                calc_w90 = self.new_calculator(block.calc_type, directory=w90_dir,
                                               **block.w90_kwargs)
                calc_w90.prefix = 'wann_preproc'
                calc_w90.command.flags = '-pp'
                self.run_calculator(calc_w90)
                utils.system_call(f'rsync -a {calc_w90.directory}/wann_preproc.nnkp {calc_w90.directory}/wann.nnkp')

                # 2) standard pw2wannier90 calculation
                calc_p2w = self.new_calculator('pw2wannier', directory=w90_dir,
                                               spin_component=block.spin,
                                               outdir=calc_pw.parameters.outdir)
                calc_p2w.prefix = 'pw2wan'
                self.run_calculator(calc_p2w)

                # 3) Wannier90 calculation
                calc_w90 = self.new_calculator(block.calc_type, directory=w90_dir,
                                               bands_plot=self.parameters.calculate_bands,
                                               **block.w90_kwargs)
                calc_w90.prefix = 'wann'
                self.run_calculator(calc_w90)

            # Merging Hamiltonian files, if necessary
            for block in self.projections.to_merge():
                self.merge_hr_files(block, prefix=calc_w90.prefix)

        if self.parameters.calculate_bands:
            # Run a "bands" calculation, making sure we don't overwrite
            # the scf/nscf tmp files by setting a different prefix
            calc_pw_bands = self.new_calculator('pw', calculation='bands', kpts=self.kpath)
            calc_pw_bands.directory = 'wannier'
            calc_pw_bands.prefix = 'bands'
            calc_pw_bands.parameters.prefix += '_bands'
            # Link the save directory so that the bands calculation can use the old density
            if self.parameters.from_scratch:
                [src, dest] = [(c.parameters.outdir / c.parameters.prefix).with_suffix('.save')
                               for c in [calc_pw, calc_pw_bands]]

                shutil.copytree(src, dest)
            self.run_calculator(calc_pw_bands)

            # Select those calculations that generated a band structure
            selected_calcs = [c for c in self.calculations[:-1] if 'band structure' in c.results]

            # Work out the vertical shift to set the valence band edge to zero
            w90_emp_num_bands = self.projections.num_bands(occ=False)
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

                    # Undo the vertical shift (we don't want the stored band structure to be vertically shifted
                    # depending on the value of self.parameters.calculate_bands!)
                    bs._energies += vbe

            # Move the legend
            lgd = ax.legend(bbox_to_anchor=(1, 1), loc="lower right", ncol=2)

            # Save the comparison to file (as png and also in editable form)
            with open(self.name + '_interpolated_bandstructure_{}x{}x{}.fig.pkl'.format(*self.kgrid), 'wb') as fd:
                pickle.dump(plt.gcf(), fd)
            # The "bbox_extra_artists" and "bbox_inches" mean that the legend is not cropped out
            plt.savefig('interpolated_bandstructure_{}x{}x{}.png'.format(*self.kgrid),
                        bbox_extra_artists=(lgd,), bbox_inches='tight')

        return

    def new_calculator(self, calc_type, *args, **kwargs) -> CalcExtType:
        calc: CalcExtType = super().new_calculator(calc_type, *args, **kwargs)

        # Extra tweaks for Wannier90 calculations
        if calc_type.startswith('w90'):
            if calc.parameters.num_bands != calc.parameters.num_wann and calc.parameters.dis_num_iter is None:
                calc.parameters.dis_num_iter = 5000
            if self.parameters.init_orbitals == 'projwfs':
                calc.parameters.num_iter = 0
        if calc_type == 'pw2wannier':
            if self._force_nspin2 and not self.parameters.spin_polarised:
                calc.parameters.spin_component = 'up'

        # Use a unified tmp directory
        if 'outdir' in calc.parameters.valid:
            calc.parameters.outdir = Path('wannier/TMP').resolve()

        return calc

    def merge_hr_files(self, block: List[projections.ProjectionBlock], prefix: str = 'wann'):
        """
        Merges the hr (Hamiltonian) files of a collection of blocks that share the same filling and spin
        """

        # Working out the files to read in and where to write out to
        fnames_in: List[Path] = []
        for b in block:
            assert b.directory is not None, 'The block which you are trying to merge is missing a directory; this ' \
                'should not happen'
            fnames_in.append(Path('wannier') / b.directory / (prefix + '_hr.dat'))
        assert b.merge_directory is not None, 'The block which you are trying to merge is missing a ' \
            'merge_directory; this should not happen'
        fname_out = Path('wannier') / b.merge_directory / (prefix + '_hr.dat')

        # Reading in each hr file in turn
        hr_list = []
        weights_out = None
        rvect_out = None
        for fname_in in fnames_in:
            # Reading the hr file
            hr, rvect, weights, nrpts = utils.read_hr_file(fname_in)

            # Sanity checking
            if weights_out is None:
                weights_out = weights
            elif weights != weights_out:
                raise ValueError(f'{fname_in} contains weights that differ from the other blocks. This should not '
                                 'happen.')
            if rvect_out is None:
                rvect_out = rvect
            elif np.all(rvect != rvect_out):
                raise ValueError(f'{fname_in} contains a set of R-vectors that differ from the other blocks. This '
                                 'should not happen.')

            # Reshaping this block of the Hamiltonian in preparation for constructing the block matrix, and storing it
            num_wann2 = hr.size // nrpts
            num_wann = int(math.sqrt(num_wann2))
            hr_list.append(hr.reshape(nrpts, num_wann, num_wann))

        # Constructing the block matrix hr_out which is dimensions (nrpts, num_wann_tot, num_wann_tot)
        num_wann_tot = sum([hr.shape[-1] for hr in hr_list])
        hr_out = np.zeros((nrpts, num_wann_tot, num_wann_tot), dtype=complex)
        start = 0
        for hr in hr_list:
            end = start + hr.shape[1]
            for irpt in range(nrpts):
                hr_out[irpt, start:end, start:end] = hr[irpt, :, :]
            start = end

        assert rvect_out is not None
        assert weights_out is not None

        utils.write_hr_file(fname_out, hr_out, rvect_out.tolist(), weights_out)
