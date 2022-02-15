"""

Workflow module for performing a single DFT calculation with either kcp.x or pw.x

Written by Edward Linscott Oct 2020

"""

import numpy as np
from koopmans import utils, pseudopotentials, mpl_config
from ._generic import Workflow
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
from matplotlib import transforms


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
        self.run_calculator(calc_dos)

        # Prepare the band structure for plotting
        bs = calc_bands.results['band structure']
        n_filled = pseudopotentials.nelec_from_pseudos(self.atoms, self.pseudopotentials, calc_scf.parameters.pseudo_dir) // 2
        vbe = bs._energies[:, :, :n_filled].max()
        bs._energies -= vbe

        # Prepare the projected density of states for plotting
        dc = calc_dos.results['dos']
        dc_summed = dc.sum_by('symbol', 'l', 'spin')
        dc_summed._energies -= vbe
        dc_up = dc_summed.select(spin = 'up')
        dc_down = dc_summed.select(spin = 'down')
        dc_down._weights *=-1

        # Plot the band structure and DOS
        fig,axes= plt.subplots(1,2,sharey=True,gridspec_kw={'width_ratios':[3,1]})
        ax_bs=axes[0]
        ax_dos= axes[1]
        bs.plot(ax=ax_bs)
        [xmin,xmax] = ax_bs.get_ylim()

        dc_up.plot(ax=ax_dos, xmin=xmin , xmax=xmax, orientation = 'vertical')
        
        ax_dos.set_prop_cycle(None)
        dc_down.plot(ax=ax_dos, xmin = xmin, orientation = 'vertical')
        
        #ax_dos.axes.get_xaxis().set_ticklabels([]) #No Tick Labels
        #ax_dos.axes.xaxis.set_visible(False) #No Tick Labels
        #ax_dos.axes.yaxis.set_visible(False) 
        ax_dos.xaxis.set_major_locator(plt.NullLocator())
        ax_dos.xaxis.set_major_formatter(plt.NullFormatter())
        ax_dos.xticklabels=('down', 'up')
        labels = [item.get_text() for item in ax_dos.get_xticklabels()]
        labels[1] = 'down'
        labels[2] = 'up'
        ax_dos.set_xticklabels(labels)

        ax_dos.legend(loc= 'upper left', bbox_to_anchor=(1,1))
        plt.savefig(fname=f'{self.name}_bands.png')
