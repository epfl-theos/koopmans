"""

Workflow module for performing a single DFT calculation with either kcp.x or pw.x

Written by Edward Linscott Oct 2020

"""

from koopmans import utils, pseudopotentials, mpl_config
from ._generic import Workflow
import matplotlib.pyplot as plt


class DFTCPWorkflow(Workflow):

    def run(self):

        calc = self.new_calculator('kcp')

        # Removing old directories
        if self.parameters.from_scratch:
            utils.system_call(f'rm -r {calc.parameters.outdir} 2>/dev/null', False)

        calc.prefix = 'dft'
        calc.directory = '.'
        calc.parameters.ndr = 50
        calc.parameters.ndw = 51
        calc.parameters.restart_mode = 'from_scratch'
        calc.parameters.do_orbdep = False
        calc.parameters.fixed_state = False
        calc.parameters.do_outerloop = True
        calc.parameters.which_compensation = 'tcc'

        if calc.parameters.maxiter is None:
            calc.parameters.maxiter = 300
        if calc.parameters.empty_states_nbnd is not None and calc.parameters.empty_states_nbnd > 0:
            calc.parameters.do_outerloop_empty = True
            if calc.parameters.empty_states_maxstep is None:
                calc.parameters.empty_states_maxstep = 300

        self.run_calculator(calc, enforce_ss=self.parameters.fix_spin_contamination)

        return calc


class DFTPWWorkflow(Workflow):

    def run(self):

        # Create the calculator
        calc = self.new_calculator('pw')

        # Update keywords
        calc.prefix = 'dft'
        calc.parameters.ndr = 50
        calc.parameters.ndw = 51
        calc.parameters.restart_mode = 'from_scratch'

        # Remove old directories
        if self.parameters.from_scratch:
            utils.system_call(f'rm -r {calc.parameters.outdir} 2>/dev/null', False)

        # Run the calculator
        self.run_calculator(calc)

        return


class PWBandStructureWorkflow(Workflow):

    def run(self):

        # First, a scf calculation
        calc_scf = self.new_calculator('pw', nbnd=None)
        calc_scf.prefix = 'scf'
        self.run_calculator(calc_scf)

        # Second, a bands calculation
        calc_bands = self.new_calculator('pw', calculation='bands', kpts=self.kpath)
        calc_bands.prefix = 'bands'
        self.run_calculator(calc_bands)

        # Third, a PDOS calculation
        calc_dos = self.new_calculator('projwfc', filpdos=self.name)
        pseudo_dir = calc_scf.parameters.pseudo_dir
        calc_dos.expected_orbitals = {}
        z_core_to_first_orbital = {0: '1s', 2: '2s', 4: '2p', 10: '3s', 12: '3p', 18: '4s', 20: '3d', 30: '4p'}
        for atom in self.atoms:
            if atom.symbol in calc_dos.expected_orbitals:
                continue
            pseudo_file = self.pseudopotentials[atom.symbol]
            upf = pseudopotentials.read_pseudo_file(pseudo_dir / pseudo_file)
            z_valence = float(upf.find('PP_HEADER').get('z_valence'))
            z_core = atom.number - z_valence
            first_orbital = z_core_to_first_orbital[z_core]
            all_orbitals = list(z_core_to_first_orbital.values())
            expected_orbitals = all_orbitals[all_orbitals.index(first_orbital):]
            calc_dos.expected_orbitals[atom.symbol] = expected_orbitals
        self.run_calculator(calc_dos)

        # Prepare the band structure for plotting
        bs = calc_bands.results['band structure']
        n_filled = pseudopotentials.nelec_from_pseudos(
            self.atoms, self.pseudopotentials, calc_scf.parameters.pseudo_dir) // 2
        vbe = bs._energies[:, :, :n_filled].max()
        bs._energies -= vbe

        # Prepare the projected density of states for plotting
        dc = calc_dos.results['dos']
        dc_summed = dc.sum_by('symbol', 'n', 'l', 'spin')
        dc_summed._energies -= vbe
        dc_up = dc_summed.select(spin='up')
        dc_down = dc_summed.select(spin='down')
        dc_down._weights *= -1

        # Plot the band structure and DOS
        _, axes = plt.subplots(1, 2, sharey=True, gridspec_kw={'width_ratios': [3, 1]})
        ax_bs = axes[0]
        ax_dos = axes[1]

        # Plot the band structure
        bs.plot(ax=ax_bs, colors=plt.rcParams['axes.prop_cycle'].by_key()['color'])
        [xmin, xmax] = ax_bs.get_ylim()

        # Plot the spin-up DOS
        spdf_order = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
        for d in sorted(dc_up, key=lambda x: (x.info['symbol'], x.info['n'], spdf_order[x.info['l']])):
            label = f'{d.info["symbol"]} {d.info["n"]}{d.info["l"]}'
            d.plot_dos(ax=ax_dos, xmin=xmin, xmax=xmax, orientation='vertical', mplargs={'label': label})

        # Reset color cycle
        ax_dos.set_prop_cycle(None)

        # Plot the spin-down DOS
        for d in sorted(dc_down, key=lambda x: (x.info['symbol'], x.info['n'], spdf_order[x.info['l']])):
            d.plot_dos(ax=ax_dos, xmin=xmin, xmax=xmax, orientation='vertical', mplargs={'label': None})

        # Tweaking the DOS figure aesthetics
        ax_dos.text(0.25, 0.10, 'up', ha='center', va='top', transform=ax_dos.transAxes)
        ax_dos.text(0.75, 0.10, 'down', ha='center', va='top', transform=ax_dos.transAxes)
        maxval = 1.1 * dc_summed._weights[:, [e >= xmin and e <= xmax for e in dc_summed._energies]].max()
        ax_dos.set_xlim([maxval, -maxval])
        ax_dos.set_xticks([])
        ax_dos.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.subplots_adjust(right=0.85, wspace=0.05)

        # Saving the figure
        plt.savefig(fname=f'{self.name}_bands.png')
