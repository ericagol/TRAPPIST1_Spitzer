

6/21/2020

This directory contains HDF files (in the *.jld2 format) which
are used in making the paper figures (in tex/figures/julia).

This files are too large to host on github (the .gitignore file
has *jld2), so please contact Eric Agol if you would like these.

References for some of these files:

- planets_2020.02.26_11.12.22.csv: NeXSci Exoplanets Database
- apjab3a3bt2_ascii.txt: Dai et al. (2019), Dressing et al. (2014) for Kepler-93b
- MR_trappist_Solar.ddat: Fully differentiated models computed by Caroline Dorn
   based upon the methods of Dorn et al. (2016), with some updates.
- MR_trappist_corefree_Solar_v02.ddat: Core-free models computed by Caroline Dorn
   based upon the methods of Dorn et al. (2016).
- T1_masses_libration_min_max.txt: Masses from TTV analysis posterior with
  the minimum and maximum eccentricities, and a random value.
- T1_timings_20191203.txt:  The 447 transit times and uncertainties used
  in the dynamical analysis for this paper.
- elements_noprior_students.txt:  The maximum-likelihood transit-timing model
  parameters consisting of mass, period, t0, e*cos(omega),e*sin(omega),inclination,
  and longitude of ascending node.
- planets_2020.02.26_11.12.22.csv: Data from NeXScI
- residuals_nomixture_rescale.txt: Residuals of optimum transit-timing fit normalized
  to the transit-timing uncertainties.
- sources.txt: Telescopes used for the 447 transit time observations.
- times_obs_and_posterior.txt:  Observed times and model posterior in ascii format.
- file_list.txt: List of HDF files (*jld2) and large txt files contained in this 
  directory (these are not in the github version since they are too large - .gitignore 
  masks these), and a brief description of the files and where they are read.
- t1_mass_radius.txt: List of masses and radii results.
- T1_photdyn_chain_noprior.txt: Posterior from photodynamical model.  Columns are
   radius-ratios (1-7), impact parameters (8-14), stellar density (15, solar units),
   limb-darkening parameters (16-19; q1-q4). 597,000 samples.
- T1_photdyn_chain_noprior_short.txt:  Same, but with only 10,000 samples.
- optimum_values.txt: Optimum values of the TTV model.  (mu,P,t0,e*cos(omega),e*sin(omega) x 7
   then log(nu), V_1*exp(1/(2*nu)).
- state_total: Contains the posterior states of the markov chain (excluding the likelihood
  parameters and the momentum terms).  Columns are the same as the rows of optimum_values.txt.
- T1_mass_radius_posterior_10k.zip:  10,000 posterior samples of the masses and radii of
   the star and planets.  Header in unzipped .txt file explains the columns.
- jld2_list_nstep2000_eps0.1.txt: List of HDF/JLD2 files from transit-timing model
   contained in the directory output_files/
- file_list.txt:  Listing of the other large HDF (*jld2) files and large .txt files.
- outlier.jl: Prints out a list of outlier data points.

The following lists the directories and their contents:

- v13/  This contains data from a search for an eighth planet, and a list of the
   HDF data files which contain the results, file_list.txt
- Grimm/ This contains data for the histograms and panels of the long-term GLF
   angle plot.
- POSTERIOR_CMF/ This contains the iron mass fractions for the differentiated
   models computed by Turbet based on Dorn et al. (2018).  Uses 10^4 draws
   from the posterior.
- POSTERIOR_NORM_DENSITY/ The "normalized density" (relative to a 20%iron/80% rock
   model computed by Turbet from Dorn et al. (2018) computed for 10^4 draws from the posterior.
- output_files/ This contains the results from the HMC transit-timing analysis
   based on code in v09_eps0.1_2k/
