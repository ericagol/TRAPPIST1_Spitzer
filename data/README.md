
6/21/2020

This directory contains HDF files (in the *.jld2 format) which
are used in making the paper figures (in tex/figures/julia).

This files are too large to host on github (the .gitignore file
has *jld2), so please contact Eric Agol if you would like these.

References for some of these files:

- massradiusEarthlikeRocky.txt:  Zeng et al. (2019)
- planets_2020.02.26_11.12.22.csv: NeXSci Exoplanets Database
- apjab3a3bt2_ascii.txt: Dai et al. (2019), Dressing et al. (2014) for Kepler-93b
- MR_trappist_corefree_Solar.ddat: Core-free model computed by Caroline Dorn
   based upon the models of Dorn et al. (2016).
-T1_masses_libration_min_max.txt: Masses from TTV analysis posterior with
  the minimum and maximum eccentricities, and a random value.
-T1_timings_20191203.txt:  The 447 transit times and uncertainties used
  in the dynamical analysis for this paper.
- elements_noprior_students.txt:  The maximum-likelihood transit-timing model
  parameters consisting of mass, period, t0, e*cos(omega),e*sin(omega),inclination,
  and longitude of ascending node.
- planets_2020.02.26_11.12.22.csv: Data from NeXScI
- residuals_nomixture_rescale.txt: Residuals of optimum transit-timing fit normalized
  to the transit-timing uncertainties.
- sources.txt: Telescopes used for the 447 transit time observations.
- times_obs_and_posterior.txt:  Observed times and model posterior in ascii format.

The following lists the directories and their contents:

- v13/  This contains data from a search for an eighth planet, and a list of the
   HDF data files which contain the results, file_list.txt
- Grimm/ This contains data for the histograms and panels of the long-term GLF
   angle plot.
- POSTERIOR_CMF/ This contains the iron mass fractions for the differentiated
   models computed by Turbet based on Zeng et al. (2019).  Uses 10^4 draws
   from the posterior.
- POSTERIOR_NORM_DENSITY/ The "normalized density" (relative to a 20%iron/80% rock
   model) computed for 10^4 draws from the posterior.
