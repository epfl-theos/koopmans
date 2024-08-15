
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
  > Please cite the papers listed in `h2o_predict.bib` in work involving this calculation

  h2o_predict
  -----------

  > [!WARNING]
  > This system is not cubic and will therefore not have a uniform dielectric tensor. However, the image-correction
  > schemes that are currently implemented assume a uniform dielectric. Proceed with caution


  > [!WARNING]
  > Predicting screening parameters with machine-learning is an experimental feature; proceed with caution

  - **Koopmans DSCF Snapshot 1 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.491291 | 0.576154 | 0.560565 | 0.491663 | 0.447122 | 0.608699 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 2 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.566859 | 0.588328 | 0.522807 | 0.521998 | 0.467255 | 0.549119 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 3 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.572956 | 0.579031 | 0.511804 | 0.507079 | 0.522491 | 0.578475 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 4 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.579284 | 0.583308 | 0.518712 | 0.517387 | 0.500875 | 0.556827 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 5 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.586545 | 0.549919 | 0.511353 | 0.515493 | 0.623179 | 0.397015 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 6 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.577255 | 0.581919 | 0.521259 | 0.516971 | 0.530282 | 0.536405 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 7 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |      4 |        5 |        6 |
        |---:|---------:|---------:|---------:|-------:|---------:|---------:|
        |  0 | 0.561836 | 0.581848 | 0.516206 | 0.5205 | 0.451726 | 0.585901 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 8 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.572049 | 0.595366 | 0.543866 | 0.542408 | 0.405701 | 0.545905 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 9 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |       6 |
        |---:|---------:|---------:|---------:|---------:|---------:|--------:|
        |  0 | 0.548768 | 0.594412 | 0.508993 | 0.507261 | 0.451001 | 0.54979 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 10 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.557584 | 0.597334 | 0.520721 | 0.519211 | 0.429592 | 0.611618 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 11 of 15**
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

        **predicted α**
        |    |        1 |       2 |        3 |        4 |       5 |        6 |
        |---:|---------:|--------:|---------:|---------:|--------:|---------:|
        |  0 | 0.570493 | 0.59744 | 0.544795 | 0.542678 | 0.42195 | 0.519504 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 12 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.588503 | 0.589973 | 0.556208 | 0.554384 | 0.402689 | 0.450284 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 13 of 15**
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

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.595083 | 0.580196 | 0.540117 | 0.535208 | 0.548499 | 0.459095 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 14 of 15**
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

        **predicted α**
        |    |        1 |        2 |      3 |        4 |        5 |       6 |
        |---:|---------:|---------:|-------:|---------:|---------:|--------:|
        |  0 | 0.583543 | 0.579593 | 0.5423 | 0.541382 | 0.457817 | 0.47079 |

    - ✅ `03-ki_final_ml` completed  
  - **Koopmans DSCF Snapshot 15 of 15**
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

        **predicted α**
        |    |        1 |       2 |        3 |        4 |       5 |        6 |
        |---:|---------:|--------:|---------:|---------:|--------:|---------:|
        |  0 | 0.580167 | 0.57688 | 0.532816 | 0.532875 | 0.47138 | 0.499316 |

    - ✅ `03-ki_final_ml` completed  

**Workflow complete** 🎉
