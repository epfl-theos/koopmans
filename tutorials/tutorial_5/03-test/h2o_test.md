
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
  > Please cite the papers listed in `h2o_test.bib` in work involving this calculation

  h2o_test
  --------

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
        |  1 | 0.493417 | 0.576716 | 0.560552 | 0.492911 | 0.463213 | 0.578981 |

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.491291 | 0.576154 | 0.560565 | 0.491663 | 0.447122 | 0.608699 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |        3 |        4 |         5 |          6 |
        |---:|---------:|---------:|---------:|---------:|----------:|-----------:|
        |  0 | 0.991671 | 0.190919 | 0.346885 | 0.999876 | -0.223314 | -0.0793481 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |    |        1 |       2 |        3 |       4 |        5 |        6 |
        |---:|---------:|--------:|---------:|--------:|---------:|---------:|
        |  0 | 0.6      | 0.6     | 0.6      | 0.6     | 0.6      | 0.6      |
        |  1 | 0.567854 | 0.58555 | 0.521799 | 0.52209 | 0.464916 | 0.556765 |

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.566859 | 0.588328 | 0.522807 | 0.521998 | 0.467255 | 0.549119 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |       3 |        4 |         5 |         6 |
        |---:|---------:|---------:|--------:|---------:|----------:|----------:|
        |  0 | 0.308552 | 0.127731 | 0.72158 | 0.719919 | -0.221631 | -0.131173 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |    |        1 |        2 |        3 |        4 |        5 |       6 |
        |---:|---------:|---------:|---------:|---------:|---------:|--------:|
        |  0 | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6     |
        |  1 | 0.572237 | 0.578571 | 0.507954 | 0.508348 | 0.529859 | 0.56078 |

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.572956 | 0.579031 | 0.511804 | 0.507079 | 0.522491 | 0.578475 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |        3 |        4 |         5 |         6 |
        |---:|---------:|---------:|---------:|---------:|----------:|----------:|
        |  0 | 0.247779 | 0.186124 | 0.855863 | 0.851467 | -0.151585 | -0.116657 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |  1 | 0.574586 | 0.579708 | 0.518171 | 0.51738 | 0.497619 | 0.546833 |

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.579284 | 0.583308 | 0.518712 | 0.517387 | 0.500875 | 0.556827 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |        3 |        4 |         5 |         6 |
        |---:|---------:|---------:|---------:|---------:|----------:|----------:|
        |  0 | 0.232531 | 0.181479 | 0.755435 | 0.768346 | -0.188521 | -0.149484 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |  1 | 0.587597 | 0.536034 | 0.507366 | 0.507948 | 0.593472 | 0.348985 |

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.586545 | 0.549919 | 0.511353 | 0.515493 | 0.623179 | 0.397015 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |        3 |       4 |          5 |         6 |
        |---:|---------:|---------:|---------:|--------:|-----------:|----------:|
        |  0 | 0.100276 | 0.650265 | 0.856186 | 0.84596 | -0.0282292 | -0.303171 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |  1 | 0.576641 | 0.580238 | 0.517019 | 0.517468 | 0.517127 | 0.540264 |

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.577255 | 0.581919 | 0.521259 | 0.516971 | 0.530282 | 0.536405 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |       3 |        4 |        5 |         6 |
        |---:|---------:|---------:|--------:|---------:|---------:|----------:|
        |  0 | 0.212739 | 0.176823 | 0.76942 | 0.764724 | -0.17612 | -0.150683 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |  1 | 0.551894 | 0.582402 | 0.513221 | 0.513847 | 0.463958 | 0.573966 |

        **predicted α**
        |    |        1 |        2 |        3 |      4 |        5 |        6 |
        |---:|---------:|---------:|---------:|-------:|---------:|---------:|
        |  0 | 0.561836 | 0.581848 | 0.516206 | 0.5205 | 0.451726 | 0.585901 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |        3 |        4 |         5 |          6 |
        |---:|---------:|---------:|---------:|---------:|----------:|-----------:|
        |  0 | 0.466627 | 0.149535 | 0.800093 | 0.790631 | -0.222878 | -0.0919165 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |  1 | 0.566523 | 0.593375 | 0.541938 | 0.541537 | 0.403634 | 0.539545 |

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.572049 | 0.595366 | 0.543866 | 0.542408 | 0.405701 | 0.545905 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |         2 |        3 |        4 |         5 |         6 |
        |---:|---------:|----------:|---------:|---------:|----------:|----------:|
        |  0 | 0.351316 | 0.0612747 | 0.531791 | 0.537885 | -0.267762 | -0.170156 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |  1 | 0.549669 | 0.592365 | 0.508269 | 0.508768 | 0.331522 | 0.591193 |

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |       6 |
        |---:|---------:|---------:|---------:|---------:|---------:|--------:|
        |  0 | 0.548768 | 0.594412 | 0.508993 | 0.507261 | 0.451001 | 0.54979 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |         2 |        3 |        4 |        5 |          6 |
        |---:|---------:|----------:|---------:|---------:|---------:|-----------:|
        |  0 | 0.499557 | 0.0621687 | 0.851869 | 0.846048 | -0.31661 | -0.0362059 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |  1 | 0.552556 | 0.594246 | 0.520174 | 0.520856 | 0.444507 | 0.581726 |

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.557584 | 0.597334 | 0.520721 | 0.519211 | 0.429592 | 0.611618 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |         2 |        3 |        4 |         5 |          6 |
        |---:|---------:|----------:|---------:|---------:|----------:|-----------:|
        |  0 | 0.486761 | 0.0486094 | 0.738739 | 0.730438 | -0.240785 | -0.0669215 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |  1 | 0.578669 | 0.593759 | 0.542751 | 0.54324 | 0.417688 | 0.511266 |

        **predicted α**
        |    |        1 |       2 |        3 |        4 |       5 |        6 |
        |---:|---------:|--------:|---------:|---------:|--------:|---------:|
        |  0 | 0.570493 | 0.59744 | 0.544795 | 0.542678 | 0.42195 | 0.519504 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |         2 |        3 |        4 |         5 |         6 |
        |---:|---------:|----------:|---------:|---------:|----------:|----------:|
        |  0 | 0.216331 | 0.0588557 | 0.529442 | 0.524046 | -0.263908 | -0.216698 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |    |       1 |        2 |        3 |        4 |       5 |       6 |
        |---:|--------:|---------:|---------:|---------:|--------:|--------:|
        |  0 | 0.6     | 0.6      | 0.6      | 0.6      | 0.6     | 0.6     |
        |  1 | 0.58385 | 0.590943 | 0.554956 | 0.554785 | 0.39397 | 0.43811 |

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.588503 | 0.589973 | 0.556208 | 0.554384 | 0.402689 | 0.450284 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |         2 |        3 |        4 |         5 |         6 |
        |---:|---------:|----------:|---------:|---------:|----------:|----------:|
        |  0 | 0.167763 | 0.0908645 | 0.412826 | 0.416922 | -0.278471 | -0.320879 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |  1 | 0.591929 | 0.579028 | 0.532618 | 0.53386 | 0.531262 | 0.455302 |

        **predicted α**
        |    |        1 |        2 |        3 |        4 |        5 |        6 |
        |---:|---------:|---------:|---------:|---------:|---------:|---------:|
        |  0 | 0.595083 | 0.580196 | 0.540117 | 0.535208 | 0.548499 | 0.459095 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |         1 |        2 |        3 |        4 |         5 |         6 |
        |---:|----------:|---------:|---------:|---------:|----------:|----------:|
        |  0 | 0.0735648 | 0.204779 | 0.626319 | 0.611921 | -0.167729 | -0.246584 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |    |        1 |        2 |        3 |        4 |        5 |       6 |
        |---:|---------:|---------:|---------:|---------:|---------:|--------:|
        |  0 | 0.6      | 0.6      | 0.6      | 0.6      | 0.6      | 0.6     |
        |  1 | 0.588736 | 0.585594 | 0.542145 | 0.543084 | 0.460209 | 0.45727 |

        **predicted α**
        |    |        1 |        2 |      3 |        4 |        5 |       6 |
        |---:|---------:|---------:|-------:|---------:|---------:|--------:|
        |  0 | 0.583543 | 0.579593 | 0.5423 | 0.541382 | 0.457817 | 0.47079 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |        3 |        4 |         5 |         6 |
        |---:|---------:|---------:|---------:|---------:|----------:|----------:|
        |  0 | 0.108654 | 0.141421 | 0.535913 | 0.525147 | -0.243062 | -0.268694 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  
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
        |    |        1 |       2 |        3 |       4 |        5 |       6 |
        |---:|---------:|--------:|---------:|--------:|---------:|--------:|
        |  0 | 0.6      | 0.6     | 0.6      | 0.6     | 0.6      | 0.6     |
        |  1 | 0.577381 | 0.57514 | 0.529469 | 0.52982 | 0.457536 | 0.49218 |

        **predicted α**
        |    |        1 |       2 |        3 |        4 |       5 |        6 |
        |---:|---------:|--------:|---------:|---------:|--------:|---------:|
        |  0 | 0.580167 | 0.57688 | 0.532816 | 0.532875 | 0.47138 | 0.499316 |

        **ΔE<sub>i</sub> - λ<sub>ii</sub> (eV)**
        |    |        1 |        2 |        3 |       4 |         5 |         6 |
        |---:|---------:|---------:|---------:|--------:|----------:|----------:|
        |  0 | 0.211444 | 0.235669 | 0.646775 | 0.64367 | -0.228191 | -0.249771 |


        > [!WARNING] UserWarning
        > The screening parameters have been calculated but are not necessarily self-consistent. You may want to
        > increase `alpha_numsteps` to obtain a more accurate result.

    - ✅ `03-ki_final_ml` completed  
    - ✅ `04-ki_final` completed  

**Workflow complete** 🎉
