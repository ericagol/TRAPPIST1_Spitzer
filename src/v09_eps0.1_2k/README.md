
2/11/2020

Okay, setting up a 10-day, 2000-step x 112 CPU
run.  Let's see how this goes!

This seems to be the run I'm going with.

6/18/2020

Here is some description of how to run an HMC chain, how to set up the slurm 
  runs on Hyak Mox, and what the output files contain:

0. First, you will need to install Julia v1.0.1 (although other Julia v1.*
  ought to work as well).  Then install Julia packages: ForwardDiff,
  DelimitedFiles, Printf, SpecialFunctions, LinearAlgebra, IterativeSolvers
  & JLD2.  The code was run with the following versions:
    "ForwardDiff"         => v"0.10.9"
    "DiffResults"         => v???
    "IterativeSolvers"    => v"0.8.1"
    "JLD2"                => v"0.1.11"
    "SpecialFunctions"    => v"0.8.0"

1. To run a single 10-step chain, from the prompt, type:

  v09_eps0.1_2k$ julia trappist1_run_hmc.jl 000 &> trappist1_hmc_000.txt &
  
  (Right now the 2000 step chains are commented out in trappist1_run_hmc.jl
  and an example of a 10-step chain will be run, which takes about one hour
  on my Mackbook.  To run a 2000 step chain, the commented lines should be 
  uncommented, and the 10-step lines should be commented; this took about
  9 days for each chain on Hyak.)
   
  To run the trappist1_run_hmc.jl script requires two input files:
  
  - ../T1_timings_20191203.txt        # Transit times and uncertainties
  - ../elements_noprior_students.txt  # Initial orbital elements
  
  Also required is source code:
  
    -trappist1_run_hmc.jl:  This is a script called by the slurm
      scripts which initializes one of the Markov chains.
    -run_hmc_background_student_rescale.jl: This is the workhorse
      script which calls NbodyGradient (via ttv.jl).
    -log_students_prob.jl: Student likelihood function (called
      in run_hmc_background_student_rescale.jl).
    
    -../compute_grad_num.jl:  
       Computes gradients numerically with finite differences to check that the code 
       is working properly (commented out - this takes a long time to run since it 
       computes finite differences at BigFloat precision).
    -../CGS.jl: 
       Defines some fundamental constants.
    -../extract_planet.jl: 
       Extracts transit times for individual planets.
    -../nlog_prior.jl:  
       Prior function which places bounds on parameters and corrects for eccentricity 
       vector bias with 1/e.
    -../loglinspace.jl: 
       Function to define a vector of linear or logarithmically spaced values.
    -../regress.jl:  
       Function which carries out linear regression.
    -../../NbodyGradient/src/ttv.jl:  
       This is the transit-time N-body integrator which computes gradients of each 
       transit time with respect to the initial orbital elements and mass ratios.
       It has several dependencies in that source directory (this has been copied
       from the NbodyGradient repository).

2. The 2000-step chains were run with slurm on Hyak Mox with the command:

  sbatch -p astro -A astro trappist1_run_hmc_05.slurm
  
  etc. for the four slurm scripts run on 4 Hyak nodes
  of 28 threads each, giving a total of 112 chains.  
  
  Note that the slurm scripts contain PATH definitions which need
  to point to the local Julia installation.
  
  Each chain was run for 2000 steps using epsilon0 = 0.1 (hence the name of this directory)
  and nleap0 = 20.
  
3. Each job run has an output file trappist1_hmc_***.txt
  which have been placed into a zip file trappist1_hmc_output.zip
  
  The results are output to JLD2 compressed HDF files
  which have been moved to the subdirectory ../../data/output_files/

The contensts of each jld2 file are:
```Hyak$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.0.1 (2018-09-29)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> using JLD2

julia> @load "T1_run_hmc_student_ecc_lgndof_V1exp2nuinv_nstep2000_eps0.1_nleap20_318.jld2"
26-element Array{Symbol,1}:
 :n            # Number of bodies = Number of planets + 1 (for the star).
 :cov_save     # Covariance matrix used in HMC.
 :hess_save    # Hessian matrix.
 :fname        # File name for the input data (T1_timings_20191203.txt  in this case).
 :felements    # File name for the initial orbital elements (elements_noprior_students.txt).
 :nparam       # Number of Markov chain parameters (5*(n-1)+2 in this case - 5 parameters
                 for each plane-parallel planet, and 2 Student's t distribution parameters).
 :t0           # Initial time of integration.
 :h            # Integration time step.
 :tmax         # Duration of integration (this needs increasing as more transit times are added).
 :ntrans       # Total number of observed transits (447 in this case).
 :iplanet      # Vector of bodies associated with each transit (1=star; 2=b; 3=c; ...; 8=h).
 :indx         # Index of each observed transit time associated in the array of computed transits.
 :tobs_tot     # Vector holding transit times.
 :sobs_tot     # Vector holding timing uncertainties.
 :data         # Array of input data (planet, transit ephemeris, transit time, uncertainty).
 :count1       # Vector length n of the total number of simulated transits for each planet (and 0 for the star).
 :state        # Results of the HMC chain.  Size is (10*(n-1)+6,nstep) (twice number of parameters for momenta
                 plus the log hamiltonian and log of the posterior function (i.e. without kinetic energy term).
 :hessian      # Hessian matrix.
 :cholh        # Cholesky decomposition of Hessian - used for choosing a momentum vector.
 :nleap0       # Maximum length of each leapfrog  HMC integration (chosen between 0.8*nleap0 and nleap0).
 :nacc         # Number of accepted steps.
 :nstep        # Number of steps in the Markov chain.
 :stats        # Five statistics which are saved for each leapfrog step: [epsilon,nleap,alpha,uni,accepted],
                 where epsilon is leapfrog timestep (chosen from absolute value of a Gaussian of width epsilon0).
                 nleap is the number of integration steps for leapfrog integration, alpha is the ratio of
                 the log posterios of the prior and trial steps, uni is uniform deviate for metropolis
                 rejection, and accepted is 1 if accepted, zero if not.
 :elements     # Initial orbital element array (n x 7) contains
                  mass, period, t0, e*cos(omega), e*sin(omega), inclination, Omega for each planet in Jacobi coordinates.
                  First line is for the star, which only has a mass value specified (usually 1.0).
 :chi_opt      # Optimized initial log likelihood (without log prior).
 :x_opt        # Values of initial optimized likelihood.
```
(Forgot to save the value of epsilon0=0.1)

The results given in the paper were taken from the state variable for the
112 chains run for 2000 leapfrog steps.  The state variables were combined
into a large array "state_total" which is saved in the path ../../tex/figures/julia/
for use in creating figures and tables for the paper.
