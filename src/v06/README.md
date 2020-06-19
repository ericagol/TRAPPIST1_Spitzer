
2/25/2020

Photodynamics run with 8 different assumptions, and three
chains each.

6/11/2020

Realized that there is a mistake in run_pd_hyak.jl.  The polynomial
fit was divided by the model with contamination equal to zero.
Need to rerun simulations which allow epsilon to vary.

Creating v20 to fix this (I've also moved the buggy one to
run_pd_hyak_bug.jl, and the corrected version is run_pd_hyak.jl).

6/18/2020

Here is some description of the directory contents, which contain
scripts for running the "photodynamical" model of the Spitzer data
on TRAPPIST-1 from Agol et al. (2020), but with the dynamics fixed.

* As described in Agol et al., there are three different assumptions
in the photodynamic analysis:

1. With or without limb-darkening prior (actually, the prior assumption
   was dropped from the paper as this did not affect the results
   materially).

2. With or without a contamination parameter, \varepsilon.

3. With or without a Gaussian inclination prior.

These are toggled with three binary parameters: limbprior, contamination,
  and incprior, whose default values are zero (for no prior/contamination).

We ran eight sets of Markov chains with these three parameters toggled
on or off, which are called from the scripts run_pd_background_01.jl 
to run_pd_background_08.jl.

The final results in the paper were reported without any priors or contamination, 
although the inclination and contamination results are discussed in the paper.

* Requirements to run the photodynamic code:

  1. Install Julia v1.0.1 (although should work with any v1.*)
  2. Install Julia packages: JLD2, DelimitedFiles, Statistics, PyPlot, AffineInvariantMCMC.
   The code was run with the following versions: (using Pkg;  Pkg.installed())
    "AffineInvariantMCMC"     => v"0.6.0"
    "PyPlot"                  => v"2.8.1"
    "JLD2"                    => v"0.1.11"

  3. Input files:
     - ../elements_noprior_students.txt:  Orbital elements used to compute the dynamics.
     - T1_Spitzer_data.jld2: Contains the Spitzer data used in the photodynamic analysis.

  4. Required scripts:
     - ../Limbdark/src/integrate_lightcurve.jl (This is borrowed from Limbdark.jl repository)

     - ../CGS.jl:
         Defines some fundamental constants.
     - ../regress.jl:
         Function which carries out linear regression.
     - ../../NbodyGradient/src/ttv.jl:
         This is the transit-time N-body integrator which computes gradients of each
         transit time with respect to the initial orbital elements and mass ratios.
         It has several dependencies in that source directory (this has been copied
         from the NbodyGradient repository).

* Example running in background (this only takes 50 steps):
    julia run_pd_background.jl 000 &> run_pd_000.out &
   This takes about 5 minutes on Macbook, and leads to an output file  
     T1_pd_MCMC_nw50_ns50_limbprior_001.jld2 which is defined in run_pd_background.jl


* The Hyak Mox Markov chain runs were carried out with the following command:

   sbatch -p vsm -A vsm  trappist1_pd_01.slurm

   which set up 24 sets of Markov chains to be run on a single VPL (vsm) node on
   Hyak Mox.  This set of runs took about 2 days of wall time for 50,000 steps each.
   These were checked with squeue -u agol to see which nodes were being used.

*  In the case of the .slurm runs, the results are saved to a file, for example, 
   T1_pd_MCMC_run_001_000.jld2, where "001" means that the chains were run with no 
   limb-prior (0), no contamination (0), and with the inclination prior (1).  
   The 000 is a label passed through during the running of run_pd_hyak.jl by
   trappist1_pd_01.slurm.  In this case the output were stored in slurm-1820442.out


* The *.jld2 file contains the following variables:

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

julia> @load "T1_pd_MCMC_run_000_001.jld2"
26-element Array{Symbol,1}:
  :chain               # Markov chain with size 19-21  x 50 x 4000 for
                         50 walkers and 4000 steps (after burn-in).
                         There are 19 variables for n=8 with: n-1 radius-ratios, n-1
                         impact parameters, 1 stellar density (solar units),
                         4 quadratic limb-darkening params for 2 Spitzer bands.
                         One extra parameter if contamination is included (\varepsilon); another
                         if inclination prior is included (\sigma_\theta).
  :llhoodvals          # Log likelihood values
  :flatchain           # Flattenend chains
  :flatllhoodvals      # Flattened log likelihood values
  :tobs                # Observed transit times
  :phot                # Vector of Spitzer photometric data (already detrended and with flares removed)
  :nfile               # Number of Spitzer photometric files = number of transit windows.
  :texp                # "Exposure time" which is actually the binned Spitzer cadence of 2.15 minutes
  :ninc                # Vector with number of planets to include in each transit window (length nfile).
  :ip_inc              # Vector of indices of planets to include in each transit window (1=b; 2=c;...; 7=h).
  :ind_inc             # Vector of indices for each planet in the computed transit-timing model for each transit window.
  :tt1                 # The array of transit times, sky-velocities, and impact parameters (not used)
                        for each body (1=star; 2=b;...;8=h) for the NbodyGradient model (computed with ttv.jl).
  :x0                  # Initial set of parameters for the Markov chain.
  :n                   # Number of bodies = number of planets + 1
  :sig_file            # Uncertainty of data points in each transit window computed from scatter of residuals
                        with respect to an initial fit which forces reduced chi-square to be of order unity.
  :i1                  # Vector of initial index associated with the phot array for each transit window (size nfile).
  :i2                  # Vector of final index associated with the phot array for each transit window (size nfile).
  :thinning            # Thinning value of the chains (set to 10 - every 10th step is stored).
  :numwalkers          # Number of walkers in the affined invariant markov chain.
  :numdims             # number of dimensions (19-21, depending on whether priors/contamination included).
  :numsamples_perwalker # Number of steps for the Markov chain.
  :burnin              # Length of burn-in (in practice this wasn't long enough).
  :astep               # Affine-invariant stepsize factor (default of 2.0).
  :limbprior           # Bool value of whether limb darkening prior was applied.
  :incprior            # Bool value of whether inclination prior was applied.
  :contamination       # Bool value of whether contamination was applied.
```

* The final results reported in the paper were combined from
     @load "T1_pd_MCMC_run_000_001.jld2"
     @load "T1_pd_MCMC_run_000_002.jld2"
     @load "T1_pd_MCMC_run_000_003.jld2"
   which have no priors and no contamination.  The inclination prior and contamination results
   are also reported on briefly.
