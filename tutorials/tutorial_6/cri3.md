
  koopmans
  ========

  *Koopmans spectral functional calculations with `Quantum ESPRESSO`*

  📦 **Version:** 1.1.0  
  🧑 **Authors:** Edward Linscott, Nicola Colonna, Riccardo De Gennaro, Ngoc Linh Nguyen, Giovanni Borghi, Andrea 
  Ferretti, Ismaila Dabo, and Nicola Marzari  
  📚 **Documentation:** https://koopmans-functionals.org  
  ❓ **Support:** https://groups.google.com/g/koopmans-users  
  🐛 **Report a bug:** https://github.com/epfl-theos/koopmans/issues/new

  > [!NOTE]  
  > Please cite the papers listed in `cri3.bib` in work involving this calculation

  cri3
  ----
  - **Koopmans DSCF**
    - **Initialization**
      - **Wannierize**
        - ✅ `01-scf` completed  
        - ✅ `02-nscf` completed  
        - **Wannierize spin up Block 1**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize spin up Block 2**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize spin up Block 3**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize spin up Block 4**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize spin up Block 5**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize spin down Block 1**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize spin down Block 2**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize spin down Block 3**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize spin down Block 4**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize spin down Block 5**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize spin down Block 6**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - ✅ `14-merge_occ_up_wannier_hamiltonian` completed  
        - ✅ `15-merge_occ_down_wannier_hamiltonian` completed  
        - ✅ `16-merge_emp_down_wannier_hamiltonian` completed  
        - ✅ `17-bands` completed  
        - ✅ `18-projwfc` completed  
      - **Fold To Supercell**
        - ✅ `01-convert_spin_up_block_1_to_supercell` completed  
        - ✅ `02-convert_spin_up_block_2_to_supercell` completed  
        - ✅ `03-convert_spin_up_block_3_to_supercell` completed  
        - ✅ `04-convert_spin_up_block_4_to_supercell` completed  
        - ✅ `05-convert_spin_up_block_5_to_supercell` completed  
        - ✅ `06-convert_spin_down_block_1_to_supercell` completed  
        - ✅ `07-convert_spin_down_block_2_to_supercell` completed  
        - ✅ `08-convert_spin_down_block_3_to_supercell` completed  
        - ✅ `09-convert_spin_down_block_4_to_supercell` completed  
        - ✅ `10-convert_spin_down_block_5_to_supercell` completed  
        - ✅ `11-convert_spin_down_block_6_to_supercell` completed  
        - ✅ `12-merging_wavefunctions_for_occupied_spin_up` completed  
        - ✅ `13-merging_wavefunctions_for_occupied_spin_down` completed  
        - ✅ `14-merging_wavefunctions_for_empty_spin_down` completed  
      - ✅ `03-dft_dummy` completed  
      - ✅ `04-dft_init` completed  
    Skipping calculation of screening parameters
**α**

    **spin up**
    |    |   267 |   268 |   269 |   270 |   271 |   272 |   273 |   274 |   275 |   276 |   277 |   278 |   279 |   280 |   281 |   282 |   283 |   284 |   285 |   286 |   287 |   288 |   289 |   290 |   291 |   292 |   293 |   294 |   295 |   296 |   297 |   298 |   299 |   300 |   301 |   302 |   303 |   304 |   305 |   306 |   307 |   308 |
    |---:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
    |  0 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 |

    **spin down**
    |    |   225 |   226 |   227 |   228 |   229 |   230 |   231 |   232 |   233 |   234 |   235 |   236 |   237 |   238 |   239 |   240 |   241 |   242 |   243 |   244 |   245 |   246 |   247 |   248 |   249 |   250 |   251 |   252 |   253 |   254 |   255 |   256 |   257 |   258 |   259 |   260 |   261 |   262 |   263 |   264 |   265 |   266 |
    |---:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
    |  0 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 | 0.122 |

    - ✅ `02-ki_final` completed  
    - **Unfold And Interpolate**
      - ✅ `01-unfold_and_interpolate_occ_up` completed  
      - ✅ `02-unfold_and_interpolate_emp_up` completed  
      - ✅ `03-unfold_and_interpolate_occ_down` completed  
      - ✅ `04-unfold_and_interpolate_emp_down` completed  

  **Workflow complete** 🎉
