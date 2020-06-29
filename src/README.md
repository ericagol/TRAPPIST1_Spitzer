
6/18/2020

In this directory lies the source code which was
used in running the analyses presented in Agol et
al. (2020).  The analyses were run with Julia 1.0.1.
Here is a list of files and code used to carry out the
analyses:

1.  N-body model:
   * NbodyGradient/  This directory contains julia
     source code used for the following dynamical analyses
     (Agol & Hernandez, 2020, in prep.)

2.  Photometric model:
   * Limbdark/ Contains source code for computing exposure-time integrated
     transit model (based on Agol et al. 2020).

3.  N-body dynamical modeling:
   * v09_eps0.1_2k/  This directory contains julia scripts
     used to run the HMC analysis for the paper.

4.  Photodynamic modeling:
   * v06/ This directory contains julia scripts for running the photodynamic
     Markov chains for constraining stellar density, radius ratio, limb-darkening,
     and impact parameters assuming the best-fit N-body model for the transit
     times and sky velocities.

5.  Other code:
   * In this directory are other julia scripts used in making the figures
     (in ../tex/figures/julia/) and in running the source code in the sub-directories
     v06/ and v09_eps0.1_2k/
