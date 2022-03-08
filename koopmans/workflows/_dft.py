"""

Workflow module for performing a single DFT calculation with either kcp.x or pw.x

Written by Edward Linscott Oct 2020

"""

from koopmans import utils, pseudopotentials, mpl_config
from ._workflow import Workflow
import matplotlib.pyplot as plt


class DFTCPWorkflow(Workflow):

    def _run(self):

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

    def _run(self):

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

    def _run(self):

        # First, a scf calculation
        calc_scf = self.new_calculator('pw', nbnd=None)
        calc_scf.prefix = 'scf'
        self.run_calculator(calc_scf)

        # Second, a bands calculation
        calc_bands = self.new_calculator('pw', calculation='bands', kpts=self.kpath)
        calc_bands.prefix = 'bands'
        self.run_calculator(calc_bands)

        # Third, a PDOS calculation
        calc_dos = self.new_calculator('projwfc', filpdos=self.name, pseudopotentials=self.pseudopotentials,
                                       spin_polarised=self.parameters.spin_polarised,
                                       pseudo_dir=calc_scf.parameters.pseudo_dir)
        self.run_calculator(calc_dos)

        # Prepare the band structure for plotting
        bs = calc_bands.results['band structure']
        n_filled = pseudopotentials.nelec_from_pseudos(
            self.atoms, self.pseudopotentials, calc_scf.parameters.pseudo_dir) // 2
        vbe = bs._energies[:, :, :n_filled].max()
        bs._energies -= vbe

        # Prepare the projected density of states for plotting
        dos = calc_dos.results['dos']
        dos_summed = dos.sum_by('symbol', 'n', 'l', 'spin')
        dos_summed._energies -= vbe

        if self.parameters.spin_polarised:
            dos_up = dos_summed.select(spin='up')
            dos_down = dos_summed.select(spin='down')
            dos_down._weights *= -1
            doss = [dos_up, dos_down]
        else:
            doss = [dos_summed]

        # Plot the band structure and DOS
        _, axes = plt.subplots(1, 2, sharey=True, gridspec_kw={'width_ratios': [3, 1]})
        ax_bs = axes[0]
        ax_dos = axes[1]

        # Plot the band structure
        bs.plot(ax=ax_bs, colors=plt.rcParams['axes.prop_cycle'].by_key()['color'])
        [xmin, xmax] = ax_bs.get_ylim()

        # Plot the DOSs
        spdf_order = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
        for dos in doss:
            for d in sorted(dos, key=lambda x: (x.info['symbol'], x.info['n'], spdf_order[x.info['l']])):
                if not self.parameters.spin_polarised or d.info.get('spin') == 'up':
                    label = f'{d.info["symbol"]} {d.info["n"]}{d.info["l"]}'
                else:
                    label = None
                d.plot_dos(ax=ax_dos, xmin=xmin, xmax=xmax, orientation='vertical', mplargs={'label': label})

            # Reset color cycle
            ax_dos.set_prop_cycle(None)

        # Tweaking the DOS figure aesthetics
        maxval = 1.1 * dos_summed._weights[:, [e >= xmin and e <= xmax for e in dos_summed._energies]].max()
        if self.parameters.spin_polarised:
            ax_dos.set_xlim([maxval, -maxval])
            ax_dos.text(0.25, 0.10, 'up', ha='center', va='top', transform=ax_dos.transAxes)
            ax_dos.text(0.75, 0.10, 'down', ha='center', va='top', transform=ax_dos.transAxes)
        else:
            ax_dos.set_xlim([0, maxval])
        ax_dos.set_xticks([])
        ax_dos.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
        plt.subplots_adjust(right=0.85, wspace=0.05)

        # Saving the figure
        plt.savefig(fname=f'{self.name}_bands.png')
