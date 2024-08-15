
  koopmans
  ========

  *Koopmans spectral functional calculations with `Quantum ESPRESSO`*

  📦 **Version:** 1.1.0  
  🧑 **Authors:** Edward Linscott, Nicola Colonna, Riccardo De Gennaro, Ngoc Linh Nguyen, Giovanni Borghi, Andrea 
  Ferretti, Ismaila Dabo, and Nicola Marzari  
  📝 **Documentation:** https://koopmans-functionals.org  
  ❓ **Support:** https://groups.google.com/g/koopmans-users  
  🐛 **Report a bug:** https://github.com/epfl-theos/koopmans/issues/new

  > [!NOTE] 
  > Please cite the papers listed in `h2o_train.bib` in work involving this calculation

  h2o_train
  ---------

  > [!WARNING] 
  > This system is not cubic and will therefore not have a uniform dielectric tensor. However, the image-correction
  > schemes that are currently implemented assume a uniform dielectric. Proceed with caution


  > [!WARNING] 
  > Predicting screening parameters with machine-learning is an experimental feature; proceed with caution

  - **Koopmans DSCF Snapshot 1 of 5**
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
      - ✅ `03-dft_dummy` completed  
      - ✅ `04-dft_init` completed  
    - **Calculate Screening Via DSCF**
      - **Iteration 1**
        - ✅ `01-ki` completed  
        - **Power Spectrum Decomposition**
          - **Convert Orbital Files To XML**
            - ✅ `01-bin2xml_total_density` completed  
            - ✅ `02-bin2xml_occ_spin_0_orb_1_density` completed  
            - ✅ `03-bin2xml_occ_spin_0_orb_2_density` completed  
            - ✅ `04-bin2xml_occ_spin_0_orb_3_density` completed  
            - ✅ `05-bin2xml_occ_spin_0_orb_4_density` completed  
            - ✅ `06-bin2xml_emp_spin_0_orb_5_density` completed  
            - ✅ `07-bin2xml_emp_spin_0_orb_6_density` completed  
          - ✅ `02-extract_coefficients_from_xml` completed  
          - ✅ `03-compute_power_spectrum_orbital_1` completed  
          - ✅ `04-compute_power_spectrum_orbital_2` completed  
          - ✅ `05-compute_power_spectrum_orbital_3` completed  
          - ✅ `06-compute_power_spectrum_orbital_4` completed  
          - ✅ `07-compute_power_spectrum_orbital_5` completed  
          - ✅ `08-compute_power_spectrum_orbital_6` completed  
        - **Orbital 1**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 2**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 3**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 4**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 5**
          - ✅ `01-pz_print` completed  
          - ✅ `02-dft_n+1_dummy` completed  
          - ✅ `03-dft_n+1` completed  
        - **Orbital 6**
          - ✅ `01-pz_print` completed  
          - ✅ `02-dft_n+1_dummy` completed  
          - ✅ `03-dft_n+1` completed  

        **α**
        |    |        1 |       2 |       3 |        4 |        5 |       6 |
        |---:|---------:|--------:|--------:|---------:|---------:|--------:|
        |  0 | 0.6      | 0.6     | 0.6     | 0.6      | 0.6      | 0.6     |
        |  1 | 0.574716 | 0.56594 | 0.49966 | 0.499231 | 0.566399 | 0.54828 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |        3 |        4 |          5 |         6 |
        |---:|---------:|---------:|---------:|---------:|-----------:|----------:|
        |  0 | 0.213412 | 0.300457 | 0.930696 | 0.938597 | -0.0964506 | -0.137783 |


        > [!WARNING] 
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final` completed  
  - **Koopmans DSCF Snapshot 2 of 5**
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
      - ✅ `03-dft_dummy` completed  
      - ✅ `04-dft_init` completed  
    - **Calculate Screening Via DSCF**
      - **Iteration 1**
        - ✅ `01-ki` completed  
        - **Power Spectrum Decomposition**
          - **Convert Orbital Files To XML**
            - ✅ `01-bin2xml_total_density` completed  
            - ✅ `02-bin2xml_occ_spin_0_orb_1_density` completed  
            - ✅ `03-bin2xml_occ_spin_0_orb_2_density` completed  
            - ✅ `04-bin2xml_occ_spin_0_orb_3_density` completed  
            - ✅ `05-bin2xml_occ_spin_0_orb_4_density` completed  
            - ✅ `06-bin2xml_emp_spin_0_orb_5_density` completed  
            - ✅ `07-bin2xml_emp_spin_0_orb_6_density` completed  
          - ✅ `02-extract_coefficients_from_xml` completed  
          - ✅ `03-compute_power_spectrum_orbital_1` completed  
          - ✅ `04-compute_power_spectrum_orbital_2` completed  
          - ✅ `05-compute_power_spectrum_orbital_3` completed  
          - ✅ `06-compute_power_spectrum_orbital_4` completed  
          - ✅ `07-compute_power_spectrum_orbital_5` completed  
          - ✅ `08-compute_power_spectrum_orbital_6` completed  
        - **Orbital 1**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 2**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 3**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 4**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 5**
          - ✅ `01-pz_print` completed  
          - ✅ `02-dft_n+1_dummy` completed  
          - ✅ `03-dft_n+1` completed  
        - **Orbital 6**
          - ✅ `01-pz_print` completed  
          - ✅ `02-dft_n+1_dummy` completed  
          - ✅ `03-dft_n+1` completed  

        **α**
        |    |        1 |        2 |        3 |       4 |        5 |        6 |
        |---:|---------:|---------:|---------:|--------:|---------:|---------:|
        |  0 | 0.6      | 0.6      | 0.6      | 0.6     | 0.6      | 0.6      |
        |  1 | 0.569337 | 0.582594 | 0.508371 | 0.50811 | 0.519749 | 0.570484 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |        3 |        4 |         5 |          6 |
        |---:|---------:|---------:|---------:|---------:|----------:|-----------:|
        |  0 | 0.278618 | 0.148526 | 0.850251 | 0.856718 | -0.165157 | -0.0940282 |


        > [!WARNING] 
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final` completed  
  - **Koopmans DSCF Snapshot 3 of 5**
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
      - ✅ `03-dft_dummy` completed  
      - ✅ `04-dft_init` completed  
    - **Calculate Screening Via DSCF**
      - **Iteration 1**
        - ✅ `01-ki` completed  
        - **Power Spectrum Decomposition**
          - **Convert Orbital Files To XML**
            - ✅ `01-bin2xml_total_density` completed  
            - ✅ `02-bin2xml_occ_spin_0_orb_1_density` completed  
            - ✅ `03-bin2xml_occ_spin_0_orb_2_density` completed  
            - ✅ `04-bin2xml_occ_spin_0_orb_3_density` completed  
            - ✅ `05-bin2xml_occ_spin_0_orb_4_density` completed  
            - ✅ `06-bin2xml_emp_spin_0_orb_5_density` completed  
            - ✅ `07-bin2xml_emp_spin_0_orb_6_density` completed  
          - ✅ `02-extract_coefficients_from_xml` completed  
          - ✅ `03-compute_power_spectrum_orbital_1` completed  
          - ✅ `04-compute_power_spectrum_orbital_2` completed  
          - ✅ `05-compute_power_spectrum_orbital_3` completed  
          - ✅ `06-compute_power_spectrum_orbital_4` completed  
          - ✅ `07-compute_power_spectrum_orbital_5` completed  
          - ✅ `08-compute_power_spectrum_orbital_6` completed  
        - **Orbital 1**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 2**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 3**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 4**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 5**
          - ✅ `01-pz_print` completed  
          - ✅ `02-dft_n+1_dummy` completed  
          - ✅ `03-dft_n+1` completed  
        - **Orbital 6**
          - ✅ `01-pz_print` completed  
          - ✅ `02-dft_n+1_dummy` completed  
          - ✅ `03-dft_n+1` completed  

        **α**
        |    |        1 |        2 |        3 |       4 |       5 |        6 |
        |---:|---------:|---------:|---------:|--------:|--------:|---------:|
        |  0 | 0.6      | 0.6      | 0.6      | 0.6     | 0.6     | 0.6      |
        |  1 | 0.578944 | 0.552147 | 0.492751 | 0.49322 | 0.58534 | 0.419386 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |        3 |        4 |          5 |         6 |
        |---:|---------:|---------:|---------:|---------:|-----------:|----------:|
        |  0 | 0.169357 | 0.434248 | 0.998794 | 0.991985 | -0.0596641 | -0.256215 |


        > [!WARNING] 
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final` completed  
  - **Koopmans DSCF Snapshot 4 of 5**
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
      - ✅ `03-dft_dummy` completed  
      - ✅ `04-dft_init` completed  
    - **Calculate Screening Via DSCF**
      - **Iteration 1**
        - ✅ `01-ki` completed  
        - **Power Spectrum Decomposition**
          - **Convert Orbital Files To XML**
            - ✅ `01-bin2xml_total_density` completed  
            - ✅ `02-bin2xml_occ_spin_0_orb_1_density` completed  
            - ✅ `03-bin2xml_occ_spin_0_orb_2_density` completed  
            - ✅ `04-bin2xml_occ_spin_0_orb_3_density` completed  
            - ✅ `05-bin2xml_occ_spin_0_orb_4_density` completed  
            - ✅ `06-bin2xml_emp_spin_0_orb_5_density` completed  
            - ✅ `07-bin2xml_emp_spin_0_orb_6_density` completed  
          - ✅ `02-extract_coefficients_from_xml` completed  
          - ✅ `03-compute_power_spectrum_orbital_1` completed  
          - ✅ `04-compute_power_spectrum_orbital_2` completed  
          - ✅ `05-compute_power_spectrum_orbital_3` completed  
          - ✅ `06-compute_power_spectrum_orbital_4` completed  
          - ✅ `07-compute_power_spectrum_orbital_5` completed  
          - ✅ `08-compute_power_spectrum_orbital_6` completed  
        - **Orbital 1**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 2**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 3**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 4**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 5**
          - ✅ `01-pz_print` completed  
          - ✅ `02-dft_n+1_dummy` completed  
          - ✅ `03-dft_n+1` completed  
        - **Orbital 6**
          - ✅ `01-pz_print` completed  
          - ✅ `02-dft_n+1_dummy` completed  
          - ✅ `03-dft_n+1` completed  

        **α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      |
        |  1 | 0.565378 | 0.589445 | 0.523162 | 0.522763 | 0.454254 | 0.565805 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |         2 |        3 |        4 |         5 |         6 |
        |---:|---------:|----------:|---------:|---------:|----------:|----------:|
        |  0 | 0.339028 | 0.0922657 | 0.708686 | 0.715767 | -0.233903 | -0.108281 |


        > [!WARNING] 
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final` completed  
  - **Koopmans DSCF Snapshot 5 of 5**
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
      - ✅ `03-dft_dummy` completed  
      - ✅ `04-dft_init` completed  
    - **Calculate Screening Via DSCF**
      - **Iteration 1**
        - ✅ `01-ki` completed  
        - **Power Spectrum Decomposition**
          - **Convert Orbital Files To XML**
            - ✅ `01-bin2xml_total_density` completed  
            - ✅ `02-bin2xml_occ_spin_0_orb_1_density` completed  
            - ✅ `03-bin2xml_occ_spin_0_orb_2_density` completed  
            - ✅ `04-bin2xml_occ_spin_0_orb_3_density` completed  
            - ✅ `05-bin2xml_occ_spin_0_orb_4_density` completed  
            - ✅ `06-bin2xml_emp_spin_0_orb_5_density` completed  
            - ✅ `07-bin2xml_emp_spin_0_orb_6_density` completed  
          - ✅ `02-extract_coefficients_from_xml` completed  
          - ✅ `03-compute_power_spectrum_orbital_1` completed  
          - ✅ `04-compute_power_spectrum_orbital_2` completed  
          - ✅ `05-compute_power_spectrum_orbital_3` completed  
          - ✅ `06-compute_power_spectrum_orbital_4` completed  
          - ✅ `07-compute_power_spectrum_orbital_5` completed  
          - ✅ `08-compute_power_spectrum_orbital_6` completed  
        - **Orbital 1**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 2**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 3**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 4**
          - ✅ `01-dft_n-1` completed  
        - **Orbital 5**
          - ✅ `01-pz_print` completed  
          - ✅ `02-dft_n+1_dummy` completed  
          - ✅ `03-dft_n+1` completed  
        - **Orbital 6**
          - ✅ `01-pz_print` completed  
          - ✅ `02-dft_n+1_dummy` completed  
          - ✅ `03-dft_n+1` completed  

        **α**
        |    |        1 |       2 |        3 |        4 |        5 |        6 |
        |---:|---------:|--------:|---------:|---------:|---------:|---------:|
        |  0 | 0.6      | 0.6     | 0.6      | 0.6      | 0.6      | 0.6      |
        |  1 | 0.583025 | 0.59729 | 0.559323 | 0.558505 | 0.382952 | 0.447094 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |         2 |       3 |       4 |         5 |         6 |
        |---:|---------:|----------:|--------:|--------:|----------:|----------:|
        |  0 | 0.181934 | 0.0270643 | 0.37397 | 0.38485 | -0.287828 | -0.308915 |


        > [!WARNING] 
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final` completed  

**Workflow complete** 🎉
