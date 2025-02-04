
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
  > Please cite the papers listed in `ozone.bib` in work involving this calculation

  ozone
  -----
  - **Koopmans DSCF**
    - **Initialization**
      - ✅ `01-dft_init_nspin1` completed  
      - ✅ `02-dft_init_nspin2_dummy` completed  
      - ✅ `03-convert_files_from_spin1to2` completed  
      - ✅ `04-dft_init_nspin2` completed  
      - ✅ `05-pz_init` completed  
    - **Calculate Screening Via DSCF**
      - **Iteration 1**
        - ✅ `01-kipz` completed  
        - **Orbital 1**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 2**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 3**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 4**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 5**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 6**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 7**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 8**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 9**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 10**
          - ✅ `01-dft_n+1_dummy` completed  
          - ✅ `02-kipz_print` completed  
          - ✅ `03-kipz_n+1` completed  

        **α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |        7 |        8 |        9 |       10 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      |
        |  1 | 0.467158 | 0.467126 | 0.470452 | 0.567696 | 0.567582 | 0.522102 | 0.522062 | 0.418538 | 0.470476 | 0.715769 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |       1 |       2 |       3 |        4 |        5 |        6 |        7 |       8 |      9 |       10 |
        |---:|--------:|--------:|--------:|---------:|---------:|---------:|---------:|--------:|-------:|---------:|
        |  0 | 1.30604 | 1.30636 | 1.25483 | 0.230445 | 0.231262 | 0.677356 | 0.677703 | 1.87798 | 1.2546 | 0.844239 |

      - **Iteration 2**
        - ✅ `01-kipz_nspin1_dummy` completed  
        - ✅ `02-convert_files_from_spin2to1` completed  
        - ✅ `03-kipz_nspin1` completed  
        - ✅ `04-kipz_nspin2_dummy` completed  
        - ✅ `05-convert_files_from_spin1to2` completed  
        - ✅ `06-kipz_nspin2` completed  
        - **Orbital 1**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 2**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 3**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 4**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 5**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 6**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 7**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 8**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 9**
          - ✅ `01-kipz_n-1` completed  
        - **Orbital 10**
          - ✅ `01-kipz_print` completed  
          - ✅ `02-kipz_n+1` completed  

        **α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |        7 |        8 |        9 |       10 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      |
        |  1 | 0.467158 | 0.467126 | 0.470452 | 0.567696 | 0.567582 | 0.522102 | 0.522062 | 0.418538 | 0.470476 | 0.715769 |
        |  2 | 0.470115 | 0.470178 | 0.475533 | 0.567498 | 0.56759  | 0.529347 | 0.529279 | 0.425569 | 0.475533 | 0.790019 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |          1 |          2 |          3 |          4 |            5 |          6 |          7 |          8 |          9 |       10 |
        |---:|-----------:|-----------:|-----------:|-----------:|-------------:|-----------:|-----------:|-----------:|-----------:|---------:|
        |  0 |  1.30604   |  1.30636   |  1.25483   | 0.230445   |  0.231262    |  0.677356  |  0.677703  |  1.87798   |  1.2546    | 0.844239 |
        |  1 | -0.0288726 | -0.0297991 | -0.0486399 | 0.00144054 | -6.20404e-05 | -0.0621426 | -0.0619079 | -0.0722085 | -0.0484106 | 0.240338 |


        > [!WARNING]  
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-kipz_final` completed  

  **Workflow complete** 🎉
