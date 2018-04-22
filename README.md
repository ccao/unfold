# unfold
This code unfolds the electronic state of a folded BZ (due to supercell) to the primitive cell BZ. Requires a generated real-space Hamiltonian (using either Wannier90 or other code, in wannier90_hr.dat form). The code follows the idea from Phys. Rev. Lett. 109, 147003.

To compile:
  1) modify make.sys
  2) go into modules directory and make
  3) go into src directory and make
  
To perform calculations:
  1) construct supercell/primitive cell structure file
  2) construct orbital mapping definitions manually or using preprocess.x (from utils).
  3) construct unfold.inp, perform unfolding
  4) postprocess using unfold2plot.x (from utils).
  5) plot
