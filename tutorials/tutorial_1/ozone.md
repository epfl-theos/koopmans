
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
    - **Calculate Screening Via DSCF**
      - **Iteration 1**
        - ✅ `01-ki` completed  
        - **Orbital 1**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 2**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 3**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 4**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 5**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 6**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 7**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 8**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 9**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 10**
          - ✅ `01-dft_n+1_dummy` completed  
          - ✅ `02-pz_print` completed  
          - ✅ `03-dft_n+1` completed  

        **α**
        |    |        1 |        2 |       3 |        4 |        5 |        6 |        7 |        8 |        9 |       10 |
        |---:|---------:|---------:|--------:|---------:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.6      | 0.6      | 0.6     | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      |
        |  1 | 0.655691 | 0.727571 | 0.78386 | 0.663859 | 0.772354 | 0.726848 | 0.729968 | 0.741899 | 0.779264 | 0.717389 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |         2 |        3 |         4 |        5 |         6 |         7 |         8 |        9 |       10 |
        |---:|---------:|----------:|---------:|----------:|---------:|----------:|----------:|----------:|---------:|---------:|
        |  0 | -0.47736 | -0.920499 | -1.12409 | -0.452149 | -1.05055 | -0.830893 | -0.785139 | -0.888597 | -1.05016 | 0.711209 |


        > [!WARNING]  
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final` completed  

  **Workflow complete** 🎉
