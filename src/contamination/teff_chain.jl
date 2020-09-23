
# Run a chain with:
#  1). T_1, T_2, f with the constraint T_1^4(1-f)+f T_2^4 = T_eff^4 within the uncertainties.
#  2). Fit these to the transit transmission data (assuming an achromatic atmosphere).
#  3). See how much this affects the Spitzer IRAC band inferred depths.

# Load in requisite Julia packages:
using DelimitedFiles
using Statistics
using PyPlot

# Read in the Flux vs. Teff:
star_grid = readdlm("phot_teff_v02.csv",',',skipstart=1)
# Indices of the "surface" fluxes:
jsurf = collect(7:2:31)  # These bands are: Kepler,IRAC Ch1, IRAC Ch2, SPECULOOS
      # g, r, i, z, y, J, H, Ks, W1, W2, W3, W4.
# Create an effective temperature grid:
tgrid = collect(2000.0:20.0:2980.0)

# Now, read in the transit transmission spectrum for each
# planet from Ducrot et al. (2020):
data_b = readdlm("transmission_spectrum_b.csv",',',skipstart=1)
data_c = readdlm("transmission_spectrum_c.csv",',',skipstart=1)
data_d = readdlm("transmission_spectrum_d.csv",',',skipstart=1)
data_e = readdlm("transmission_spectrum_e.csv",',',skipstart=1)
data_f = readdlm("transmission_spectrum_f.csv",',',skipstart=1)
data_g = readdlm("transmission_spectrum_g.csv",',',skipstart=1)
data_h = readdlm("transmission_spectrum_h.csv",',',skipstart=1)
# Define the spectral indices from Adam's grid for each planet.
# These are the bands for each planet in the transit transmission
# spectra from Elsa's paper:
jb = [1,8,4,10,12,2,3] # K2, LT (z'), SPECULOOS, J, NB2090 (~Ks), Ch1, Ch2
jc = [1,8,4,10,12,2,3]
jd = [1,8,4,10,2,3]
je = [1,8,4,10,2,3]
jf = [1,4,10,2,3] 
jg = [1,4,10,2,3] 
jh = [1,8,4,2,3]

# Define a log likelihood function:
function log_likelihood(param)
  # Ten parameters: two temperatures, spot covering fraction, and seven depths.
  T1,T2,f,db,dc,dd,de,df,dg,dh = param
  # Okay, now interpolate effective temperatures. First check to see if in bounds:
  if (T1 < 2000.0) || (T1 > 2980.0) || (T2 < 2000.0) || (T2 > 2980.0) || (minimum(param) <= 0.0) || (f >= 1.0)
    return -Inf
  end
  # Interpolate the spectra for T1 and T2:
  it1 = ceil(Int64,(T1 - 2000)/980*49)
  it2 = ceil(Int64,(T2 - 2000)/980*49)
  # Now, find interpolated fluxes:
  f1 = star_grid[it1,jsurf]*(tgrid[it1+1]-T1)/20.0+star_grid[it1+1,jsurf]*(T1-tgrid[it1])/20.0
  f2 = star_grid[it2,jsurf]*(tgrid[it2+1]-T2)/20.0+star_grid[it2+1,jsurf]*(T2-tgrid[it2])/20.0
  # Compute the impact on transit depth (equation 1 in Rackham et al.).
  # We are assuming that *all* planets transit the same temperature,
  # so we require that f > 50%; otherwise the spots and surface flip roles:
  if f < 0.5
    tlam = 1 ./(1 .-  f .*(1 .- f2 ./ f1))
  else
    tlam = 1 ./(1 .- (1 - f) .*(1 .- f1 ./ f2))
  end
  # Compute the chi-squares:
  chisquare  = sum((data_b[:,4] .- tlam[jb] .* db).^2 ./ data_b[:,5].^2)
  chisquare += sum((data_c[:,4] .- tlam[jc] .* dc).^2 ./ data_c[:,5].^2)
  chisquare += sum((data_d[:,4] .- tlam[jd] .* dd).^2 ./ data_d[:,5].^2)
  chisquare += sum((data_e[:,4] .- tlam[je] .* de).^2 ./ data_e[:,5].^2)
  chisquare += sum((data_f[:,4] .- tlam[jf] .* df).^2 ./ data_f[:,5].^2)
  chisquare += sum((data_g[:,4] .- tlam[jg] .* dg).^2 ./ data_g[:,5].^2)
  chisquare += sum((data_h[:,4] .- tlam[jh] .* dh).^2 ./ data_h[:,5].^2)
  # That's it!
  # Add to it the constraint on Teff:
  chisquare += ((T1^4*(1-f)+T2^4*f)^0.25 - 2566)^2/25^2
  # Return the log likelihood:
  return -0.5*chisquare
end


# Compute the correction factors in the two IRAC bands:
function irac_correction(param)
  T1,T2,f,_ = param
# Computes the correction in the IRAC bands:
  if (T1 < 2000.0) || (T1 > 2980.0) || (T2 < 2000.0) || (T2 > 2980.0) || (minimum(param) <= 0.0) || (f >= 1.0)
    return -Inf
  end
  # Interpolate the spectra for T1 and T2:
  it1 = ceil(Int64,(T1 - 2000)/980*49)
  it2 = ceil(Int64,(T2 - 2000)/980*49)
  # Now, find interpolated fluxes:
  f1 = star_grid[it1,jsurf[2:3]]*(tgrid[it1+1]-T1)/20.0+star_grid[it1+1,jsurf[2:3]]*(T1-tgrid[it1])/20.0
  f2 = star_grid[it2,jsurf[2:3]]*(tgrid[it2+1]-T2)/20.0+star_grid[it2+1,jsurf[2:3]]*(T2-tgrid[it2])/20.0
  if f < 0.5
    return 1 ./(1 .- f .*(1 .- f2 ./ f1))
  else
    return 1 ./(1 .- (1-f) .*(1 .- f1 ./ f2))
  end
end

# Okay, now run a Markov chain:
using AffineInvariantMCMC

# Set up initial conditions for the chain:
nparam  = 10
nwalker = 30
astep   = 2.0
xstart = zeros(nparam,nwalker)
for i=1:nwalker
  param = zeros(nparam)
  param[1] = 2000.0 + 980.0 * rand()
  param[2] = 2000.0 + 980.0 * rand()
  param[3:10] = rand(8)
  xstart[:,i] = param
end

# Run the chain:
nsample = 100000
nthin   = 100
# Okay, now run the chain:
chain,llhoodvals = AffineInvariantMCMC.sample(log_likelihood,nwalker,xstart,nsample,nthin,astep)

# Now, compute the IRAC band corrections:
irac_grid = zeros(2,nwalker,div(nsample,nthin))
for i=1:nwalker
  for j=1:div(nsample,nthin)
    irac_grid[:,i,j] = irac_correction(chain[:,i,j])
  end
end

# Print the results:
println("IRAC Ch1: ",@sprintf("%6.4f",mean(sqrt.(irac_grid[1,:,150:end]))),"+-",@sprintf("%6.4f",std(sqrt.(irac_grid[1,:,150:end]))))
println("IRAC Ch2: ",@sprintf("%6.4f",mean(sqrt.(irac_grid[2,:,150:end]))),"+-",@sprintf("%6.4f",std(sqrt.(irac_grid[2,:,150:end]))))
