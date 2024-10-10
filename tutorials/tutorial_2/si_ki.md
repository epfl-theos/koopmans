
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
  > Please cite the papers listed in `si.bib` in work involving this calculation

  si
  --

  > [!WARNING]  
  > Makov-Payne corrections are not being used; do this with caution for periodic systems


  > [!WARNING]  
  > `eps_inf` missing in input; it will default to 1.0. Proceed with caution for periodic systems; consider setting
  > `eps_inf == "auto"` to calculate it automatically.

  - **Koopmans DSCF**
    - **Initialization**
      - **Wannierize**
        - ✅ `01-scf` completed  
        - ✅ `02-nscf` completed  
        - **Wannierize Block 1**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize Block 2**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
      - **Fold To Supercell**
        - ✅ `01-convert_block_1_to_supercell` completed  
        - ✅ `02-convert_block_2_to_supercell` completed  

        > [!WARNING]  
        > Small box parameters `nrb` not provided in input: these will be automatically set to safe default values.
        > These values can probably be decreased, but this would require convergence tests. Estimated real mesh
        > dimension `(nr1, nr2, nr3) = 72 72 72`. Small box mesh dimension `(nr1b, nr2b, nr3b) = 30 30 30`.

      - ✅ `03-dft_dummy` completed  
      - ✅ `04-dft_init` completed  
    - **Calculate Screening Via DSCF**
      - **Iteration 1**
        - ✅ `01-ki` completed  
        - **Orbital 32**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 33**
          - ✅ `01-dft_n+1_dummy` completed  
          - ✅ `02-pz_print` completed  
          - ✅ `03-dft_n+1` completed  

        **α**
        |    |       32 |        33 |
        |---:|---------:|----------:|
        |  0 | 0.077    | 0.077     |
        |  1 | 0.133007 | 0.0404255 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        32 |        33 |
        |---:|----------:|----------:|
        |  0 | -0.311871 | -0.113775 |


        > [!WARNING]  
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final` completed  
    - **Unfold And Interpolate**
      - **Wannierize**
        - ✅ `01-scf` completed  
        - ✅ `02-nscf` completed  
        - **Wannierize Block 1**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - **Wannierize Block 2**
          - ✅ `01-wannier90_preproc` completed  
          - ✅ `02-pw2wannier90` completed  
          - ✅ `03-wannier90` completed  
        - ✅ `05-bands` completed  
        - ✅ `06-projwfc` completed  
      - ✅ `02-unfold_and_interpolate_occ` completed  
      - ✅ `03-unfold_and_interpolate_emp` completed  

  **Workflow complete** 🎉
