"""
Interpolation module for the UI calculator

Originally written by Riccardo De Gennaro as part of the standalone 'unfolding and interpolate' code
Integrated within koopmans by Edward Linscott Jan 2021

"""

import numpy as np
from time import time
from ._utils import crys_to_cart
from ase.dft.dos import DOS
from ase.spectrum.band_structure import BandStructure


def interpolate(self, start_time):
    """
    interpolate is the main program in this module and it calls consecutively
                the three independent functions:
                - map_wannier
                - calc_bands
                - calc_dos

    """

    # Step 1: map the WFs
    if self.do_map:
        self.map_wannier()
        self.f_out.write(f'\tBuilding the map |i> --> |Rn> in:\t{time()-start_time:.3f} sec\n')
    reset = time()

    # Step 2: calculate the electronic bands along kpath
    self.calc_bands()
    self.f_out.write(f'\tCalculating bands in: {time()-reset:22.3f} sec\n')
    reset = time()

    # Step 3 (optional) : calculate the density-of-states
    if self.do_dos:
        self.calc_dos()
        self.f_out.write(f'\tCalculating DOS in: {time()-reset:24.3f} sec\n')

    return


def map_wannier(self):
    """
    map_wannier builds the map |i> --> |Rn> between the WFs in the SC and in the PC.
    """

    centers = []
    spreads = []
    index = []
    num_wann = 0

    # here we identify the WFs within the R=0 cell
    for n in range(self.num_wann_sc):
        # shift the WFs within the SC
        self.centers[n] = self.centers[n] / self.sc_dim
        self.centers[n] = self.centers[n] - np.floor(self.centers[n])

        # converting centers from crystal units of SC to crystal units of PC
        self.centers[n] = self.centers[n] * self.sc_dim

        if (self.centers[n][0] - 1 < 1.e-3) and (self.centers[n][1] - 1 < 1.e-3) and (self.centers[n][2] - 1 < 1.e-3):
            centers.append(self.centers[n])
            spreads.append(self.spreads[n])
            index.append(n)
            num_wann += 1

    # check on the WFs found in the R=0 cell
    assert num_wann == self.num_wann, 'Did not find the right number of WFs in the R=0 cell'

    # here we identify with |Rn> the WFs in the rest of the SC
    # the WFs are now ordered as (R0,1),(R0,2),...,(R0,n),(R1,1),...
    for rvect in self.Rvec[1:]:
        count = 0
        for m in range(self.num_wann):
            for n in range(self.num_wann_sc):
                wf_dist = self.centers[n] - centers[m]
                if (np.linalg.norm(wf_dist - rvect) < 1.e-3) and (abs(self.spreads[n] - spreads[m] < 1.e-3)):
                    centers.append(self.centers[n])
                    spreads.append(self.spreads[n])
                    index.append(n)
                    count += 1

        assert count == self.num_wann, f'Found {count} WFs in the {rvect} cell'

    # redefine phases and hr in order to follow the new order of WFs
    phases = []
    hr = np.zeros((self.num_wann_sc, self.num_wann_sc), dtype=complex)
    for n in range(self.num_wann_sc):
        phases.append(self.phases[index[n]])
        for m in range(self.num_wann_sc):
            hr[m, n] = self.hr[index[m], index[n]]

    self.centers = centers
    self.spreads = spreads
    self.hr = hr
    self.phases = phases

    return


def calc_bands(self):
    """
    calc_bands interpolates the k-space hamiltonian along the input path, by Fourier
               transforming the Wannier hamiltonian H(R). The function generates two
               new attributes:
               - self.hk containing H(k) for any k-vector in the input path
               - self.results['band structure'] containing the interpolated electronic energies

    """

    # when smooth interpolation is on, we remove the DFT part from hr
    hr = np.array(self.hr[:, :self.num_wann], dtype=complex)
    if self.parameters.do_smooth_interpolation:
        hr = hr - self.hr_coarse

    # renormalize H(R) on the WF phases
    for m in range(self.num_wann_sc):
        for n in range(self.num_wann):
            hr[m, n] = self.phases[m].conjugate() * hr[m, n] * self.phases[n]

    # here we build the interpolated H(k)
    hk = np.zeros((len(self.kpath.kpts), self.num_wann, self.num_wann), dtype=complex)
    bands = []
    self.f_out.write(f"\n\t\tTotal number of k-points: {len(self.kpath.kpts):6d}\n\n")
    for ik, kpt in enumerate(self.kpath.kpts):
        self.f_out.write(f"\t\t      calculating point # {ik+1}\n")
        for m in range(self.num_wann):
            for n in range(self.num_wann):
                for ir in range(len(self.Rvec)):

                    mm = ir * self.num_wann + m
                    phase = self.correct_phase(self.centers[n], self.centers[mm], self.Rvec[ir], kpt)

                    hk[ik, m, n] = hk[ik, m, n] + np.exp(1j * 2 * np.pi * np.dot(kpt, self.Rvec[ir])) * \
                        phase * hr[mm, n]

                if self.parameters.do_smooth_interpolation:
                    hk_smooth = 0
                    for jr in range(len(self.Rsmooth)):
                        Rs = self.Rsmooth[jr]
                        hk_smooth = hk_smooth + np.exp(1j * 2 * np.pi * np.dot(kpt, Rs)) * \
                            self.hr_smooth[jr, m, n] / self.wRs[jr]
                    hk[ik, m, n] = hk[ik, m, n] + hk_smooth

        bands.append(np.linalg.eigvalsh(hk[ik]))

    self.f_out.write('\n')

    self.hk = hk
    self.results['band structure'] = BandStructure(self.kpath, [bands])

    return


def calc_dos(self):
    """
    calc_dos calculates the density of states using a gaussian smearing. The DOS is saved
             as a list [ [E1, DOS(E1)], [E2, DOS[E2]], ... , [En, DOS(En)] ]
    """

    self.results['dos'] = DOS(self, width=self.degauss, window=(self.Emin, self.Emax), npts=self.nstep + 1)

    return


def correct_phase(self, center_ref, center, rvect, kvect):
    """
    correct_phase calculate the correct phase factor to put in the Fourier transform
                  to get the interpolated k-space hamiltonian. The correction consists
                  of finding the right distance, i.e. the right R-vector, considering
                  also the BVK boundary conditions.
                  if use_ws_distance=True, the function accounts also for the intracell
                  distance between Wannier functions, otherwise only the intercell
                  distances are considered.

       IMPORTANT: the input vectors (center_ref, center, Rvec and kvec) must all
                  be in crystal units otherwise the distances are not properly evaluated.
    """

    if self.use_ws_distance:
        wf_dist = crys_to_cart(center - center_ref, self.at, +1)
    else:
        wf_dist = crys_to_cart(rvect, self.at, +1)

    dist_min = np.linalg.norm(wf_dist)
    Tlist = []

    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                tvect = np.array([i, j, k]) * self.sc_dim
                Tvec = crys_to_cart(tvect, self.at, +1)
                dist = np.linalg.norm(wf_dist + Tvec)

                if (abs(dist - dist_min) < 1.e-3):
                    #
                    # an equidistant replica is found
                    #
                    Tlist.append(tvect)

                elif (dist < dist_min):
                    #
                    # a new nearest replica is found
                    # reset dist_min and reinitialize Tlist
                    #
                    dist_min = dist
                    Tlist = [tvect]

                else:
                    #
                    # this replica is rejected
                    #
                    continue

    phase = sum(np.exp(1j * 2 * np.pi * np.dot(Tlist, kvect))) / len(Tlist)

    return phase
