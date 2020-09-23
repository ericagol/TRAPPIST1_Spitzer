
# Run a chain with:
#  1). T_1, T_2, f with the constraint T_1^4(1-f)+f T_2^4 = T_eff^4 within the uncertainties.
#  2). Fit these to the transit transmission data (assuming an achromatic atmosphere).
#  3). See how much this affects the Spitzer IRAC band inferred depths.

# Load in requisite Julia packages:
using DelimitedFiles
using Statistics
using PyPlot
using Printf

# Read in the Flux vs. Teff:
star_grid = readdlm("phot_teff_v02.csv",',',skipstart=1)
# Indices of the "surface" fluxes:
jsurf = collect(7:2:37)  # These bands are: 1 Kepler,2 IRAC Ch1,3 IRAC Ch2,4 SPECULOOS
      # 5 g, 6 r, 7 i, 8 z, 9 y, 10 J, 11 H, 12 Ks, 13 W1, 14 W2, 15 W3, 16 W4.
lam_mod = [ 0.807, 3.6, 4.5, 0.9102, 0.4811, 0.6156, 0.7504, 0.8668, 0.9613, 
            1.235, 1.662, 2.159, 3.368, 4.618, 12.082, 22.194]
isort = sortperm(lam_mod)
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

function compute_tlam(param)
  # Computes the transit wavelength dependence:
  T1,T2,f,db,dc,dd,de,df,dg,dh = param
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
  return tlam
end

# Define a log likelihood function:
function log_likelihood(param)
  # Ten parameters: two temperatures, spot covering fraction, and seven depths.
  T1,T2,f,db,dc,dd,de,df,dg,dh = param
  # Okay, now interpolate effective temperatures. First check to see if in bounds:
  if (T1 < 2000.0) || (T1 > 2980.0) || (T2 < 2000.0) || (T2 > 2980.0) || (minimum(param) <= 0.0) || (f >= 1.0)
    return -Inf
  end
  # Compute the transit wavelength depedence:
  tlam = compute_tlam(param)
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
#nsample = 100000
#nthin   = 100
nsample = 10000
nthin   = 10
# First run burn-in for 2000 steps:
chain,llhoodvals = AffineInvariantMCMC.sample(log_likelihood,nwalker,xstart,2000,nthin,astep)
# Okay, now run the chain:
chain,llhoodvals = AffineInvariantMCMC.sample(log_likelihood,nwalker,chain[:,:,end],nsample,nthin,astep)
# Now flatten the chains:
flatchain, flatllhoodvals = AffineInvariantMCMC.flattenmcmcarray(chain, llhoodvals)
nflat = nwalker*div(nsample,nthin)
# Now, compute the IRAC band corrections:
irac_grid = zeros(2,nflat)
for i=1:nflat
  irac_grid[:,i] = irac_correction(flatchain[:,i])
end

# Print the results:
println("IRAC Ch1: ",@sprintf("%6.4f",mean(sqrt.(irac_grid[1,:]))),"+-",@sprintf("%6.4f",std(sqrt.(irac_grid[1,:]))))
println("IRAC Ch2: ",@sprintf("%6.4f",mean(sqrt.(irac_grid[2,:]))),"+-",@sprintf("%6.4f",std(sqrt.(irac_grid[2,:]))))

# Now, plot some results:
fig,axes = subplots(2,4)
ax = axes[1]
ax.errorbar(data_b[:,3],data_b[:,4],yerr=data_b[:,5],fmt="o",label="b")
ax.legend()
ax = axes[2]
ax.errorbar(data_c[:,3],data_c[:,4],yerr=data_c[:,5],fmt="o",label="c")
ax.legend()
ax = axes[3]
ax.errorbar(data_d[:,3],data_d[:,4],yerr=data_d[:,5],fmt="o",label="d")
ax.legend()
ax = axes[4]
ax.errorbar(data_e[:,3],data_e[:,4],yerr=data_e[:,5],fmt="o",label="e")
ax.legend()
ax = axes[5]
ax.errorbar(data_f[:,3],data_f[:,4],yerr=data_f[:,5],fmt="o",label="f")
ax.legend()
ax = axes[6]
ax.errorbar(data_g[:,3],data_g[:,4],yerr=data_g[:,5],fmt="o",label="g")
ax.legend()
ax = axes[7]
ax.errorbar(data_h[:,3],data_h[:,4],yerr=data_h[:,5],fmt="o",label="h")
ax.legend()

dsamp = zeros(nflat,7,16)
for i=1:nflat
  # Randomly choose posterior values:
  param = flatchain[:,i]
  tlam_samp = compute_tlam(param)
  for ip=1:7
     dsamp[i,ip,:] =  tlam_samp .* param[3+ip]
  end
end
for ip=1:7
  ax = axes[ip]
  depth_mean = zeros(16)
  depth_sig = zeros(16)
  for j=1:16
    depth_mean[j] = mean(dsamp[:,ip,j])
    depth_sig[j] = std(dsamp[:,ip,j])
  end
  ax.semilogx(lam_mod[isort],depth_mean[isort],color="g")
  ax.semilogx(lam_mod[isort],depth_mean[isort]-depth_sig[isort],color="g",linestyle=":")
  ax.semilogx(lam_mod[isort],depth_mean[isort]+depth_sig[isort],color="g",linestyle=":")
  ax.semilogx(lam_mod[isort],zeros(16) .+ mean(flatchain[3+ip,:]),color="r")
  ax.semilogx(lam_mod[isort],zeros(16) .+ (mean(flatchain[3+ip,:]) - std(flatchain[3+ip,:])),color="r",linestyle=":")
  ax.semilogx(lam_mod[isort],zeros(16) .+ (mean(flatchain[3+ip,:]) + std(flatchain[3+ip,:])),color="r",linestyle=":")
  ax.set_xlabel(L"$\lambda$ [micron]")
  ax.set_ylabel("Depth [%]")
end
