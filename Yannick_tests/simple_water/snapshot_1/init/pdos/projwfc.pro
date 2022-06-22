
     Program PROJWFC v.7.0 starts on 22Jun2022 at 13: 2:13 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
         "P. Giannozzi et al., J. Phys.:Condens. Matter 29 465901 (2017);
         "P. Giannozzi et al., J. Chem. Phys. 152 154105 (2020);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Parallel version (MPI), running on     8 processors

     MPI processes distributed on     1 nodes
     R & G space division:  proc/nbgrp/npool/nimage =       8
     25925 MiB available memory on the printing compute node when the environment starts
 

     Reading xml data from directory:

     /home/yshubert/code/koopmans/Yannick_tests/simple_water/snapshot_1/init/wannier/TMP/kc.save/

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation= PBE
                           (   1   4   3   4   0   0   0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want

 
     Parallelization info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Min         582     582    144                27708    27708    3460
     Max         583     583    145                27712    27712    3465
     Sum        4657    4657   1159               221683   221683   27705
 
     Using Slab Decomposition
 

     Gaussian broadening (read from input): ngauss,degauss=   0    0.010000


     Calling projwave .... 
     Subspace diagonalization in iterative solution of the eigenvalue problem:
     a serial algorithm will be used

 
  Problem Sizes 
  natomwfc =            6
  nbnd     =           12
  nkstot   =            1
  npwx     =         3465
  nkb      =           12
 

     Atomic states used for projection
     (read from pseudopotential files):

     state #   1: atom   1 (O  ), wfc  1 (l=0 m= 1)
     state #   2: atom   1 (O  ), wfc  2 (l=1 m= 1)
     state #   3: atom   1 (O  ), wfc  2 (l=1 m= 2)
     state #   4: atom   1 (O  ), wfc  2 (l=1 m= 3)
     state #   5: atom   2 (H  ), wfc  1 (l=0 m= 1)
     state #   6: atom   3 (H  ), wfc  1 (l=0 m= 1)

 k =   0.0000000000  0.0000000000  0.0000000000
==== e(   1) =   -24.82469 eV ==== 
     psi = 0.702*[#   1]+0.123*[#   5]+0.123*[#   6]+0.021*[#   3]+0.012*[#   4]
          +0.006*[#   2]
    |psi|^2 = 0.987
==== e(   2) =   -12.64222 eV ==== 
     psi = 0.380*[#   4]+0.280*[#   2]+0.143*[#   5]+0.142*[#   6]+0.032*[#   3]
          
    |psi|^2 = 0.977
==== e(   3) =    -8.86890 eV ==== 
     psi = 0.447*[#   3]+0.245*[#   4]+0.123*[#   1]+0.123*[#   2]+0.023*[#   5]
          +0.023*[#   6]
    |psi|^2 = 0.984
==== e(   4) =    -6.77102 eV ==== 
     psi = 0.441*[#   2]+0.401*[#   3]+0.148*[#   4]
    |psi|^2 = 0.990
==== e(   5) =    -1.07777 eV ==== 
     psi = 0.199*[#   5]+0.195*[#   6]+0.076*[#   1]+0.025*[#   3]+0.014*[#   4]
          +0.007*[#   2]
    |psi|^2 = 0.515
==== e(   6) =     1.80927 eV ==== 
     psi = 0.169*[#   6]+0.166*[#   5]+0.052*[#   4]+0.033*[#   2]+0.005*[#   3]
          
    |psi|^2 = 0.424
==== e(   7) =     2.13658 eV ==== 
     psi = 0.074*[#   6]+0.065*[#   5]+0.028*[#   1]+0.018*[#   3]+0.007*[#   4]
          +0.006*[#   2]
    |psi|^2 = 0.196
==== e(   8) =     2.63335 eV ==== 
     psi = 0.019*[#   6]+0.005*[#   1]+0.001*[#   3]
    |psi|^2 = 0.028
==== e(   9) =     2.75928 eV ==== 
     psi = 0.001*[#   2]
    |psi|^2 = 0.003
==== e(  10) =     3.08550 eV ==== 
     psi = 0.019*[#   5]+0.008*[#   6]+0.005*[#   2]+0.002*[#   4]
    |psi|^2 = 0.035
==== e(  11) =     3.19162 eV ==== 
     psi = 0.031*[#   5]+0.012*[#   1]+0.009*[#   6]+0.005*[#   4]+0.002*[#   3]
          
    |psi|^2 = 0.059
==== e(  12) =     4.63343 eV ==== 
     psi = 0.063*[#   6]+0.036*[#   5]+0.025*[#   4]+0.025*[#   2]+0.005*[#   3]
          
    |psi|^2 = 0.155

Lowdin Charges: 

     Atom #   1: total charge =   6.7201, s =  1.6512, 
     Atom #   1: total charge =   6.7201, p =  5.0689, pz=  1.6985, px=  1.8012, py=  1.5692, 
     Atom #   2: total charge =   0.5785, s =  0.5785, 
     Atom #   2: total charge =   0.5785, p =  0.0000, pz=  0.0000, px=  0.0000, py=  0.0000, 
     Atom #   3: total charge =   0.5775, s =  0.5775, 
     Atom #   3: total charge =   0.5775, p =  0.0000, pz=  0.0000, px=  0.0000, py=  0.0000, 
     Spilling Parameter:   0.0155
 
     PROJWFC      :      0.32s CPU      0.43s WALL

 
   This run was terminated on:  13: 2:14  22Jun2022            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
