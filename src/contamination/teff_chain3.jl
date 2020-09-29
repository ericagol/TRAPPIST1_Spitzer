
# Run a chain with:
#  1). T_1, T_2, T_3, f_2, f_3 with the constraint T_1^4(1-f_2-f_3)+f_2 T_2^4 + f_3 T_3^4 = T_eff^4 within the uncertainties.
#  2). Fit these to the transit transmission data (assuming an achromatic atmosphere).
#  3). See how much this affects the Spitzer IRAC band inferred depths.

# Load in requisite Julia packages:
using DelimitedFiles
using Statistics
using PyPlot
using Printf
using Optim

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
  T1,T2,T3,f2,f3,db,dc,dd,de,df,dg,dh = param
  # Interpolate the spectra for T1, T2 & T3:
  it1 = ceil(Int64,(T1 - 2000)/980*49)
  it2 = ceil(Int64,(T2 - 2000)/980*49)
  it3 = ceil(Int64,(T3 - 2000)/980*49)
  # Now, find interpolated fluxes:
  flux1 = star_grid[it1,jsurf]*(tgrid[it1+1]-T1)/20.0+star_grid[it1+1,jsurf]*(T1-tgrid[it1])/20.0
  flux2 = star_grid[it2,jsurf]*(tgrid[it2+1]-T2)/20.0+star_grid[it2+1,jsurf]*(T2-tgrid[it2])/20.0
  flux3 = star_grid[it3,jsurf]*(tgrid[it3+1]-T3)/20.0+star_grid[it3+1,jsurf]*(T3-tgrid[it3])/20.0
  # Compute the impact on transit depth (equation 1 in Rackham et al.).
  # We are assuming that *all* planets transit the same temperature,
  # so we require that f > 50%; otherwise the spots and surface flip roles:
  if (f2+f3) < 2//3
    tlam = 1 ./((1 - f2 - f3) .+ f2 .*(flux2 ./ flux1) .+ f3 .*(flux3 ./ flux1))
  elseif (1-f2) < 2//3
    tlam = 1 ./(f2 .+ (1 - f2 - f3) .*(flux1 ./ flux2) .+ f3 .*(flux3 ./ flux2))
  else
    tlam = 1 ./(f3 .+ (1 - f2 - f3) .*(flux1 ./ flux3) .+ f2 .*(flux2 ./ flux3))
  end
  return tlam
end

# Define a log likelihood function:
function log_likelihood(param)
  # Ten parameters: two temperatures, spot covering fraction, and seven depths.
  T1,T2,T3,f2,f3,db,dc,dd,de,df,dg,dh = param
  # Okay, now interpolate effective temperatures. First check to see if in bounds:
  if (T1 < 2000.0) || (T1 > 2980.0) || (T2 < 2000.0) || (T2 > 2980.0) || (T3 < 2000.0) || (T3 > 2980.0) || (minimum(param) <= 0.0) || ((f2+f3) >= 1.0)
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
  chisquare += ((T1^4*(1-f2-f3)+T2^4*f2+T3^4+f3)^0.25 - 2566)^2/25^2
  # Return the log likelihood:
  return -0.5*chisquare
end


# Compute the correction factors in the two IRAC bands:
function irac_correction(param)
  T1,T2,T3,f2,f3,_ = param
# Computes the correction in the IRAC bands:
  if (T1 < 2000.0) || (T1 > 2980.0) || (T2 < 2000.0) || (T2 > 2980.0) || (T3 < 2000.0) || (T3 > 2980.0) || (minimum(param) <= 0.0) || ((f2+f3) >= 1.0)
    return -Inf
  end
  # Interpolate the spectra for T1, T2 & T3:
  it1 = ceil(Int64,(T1 - 2000)/980*49)
  it2 = ceil(Int64,(T2 - 2000)/980*49)
  it3 = ceil(Int64,(T3 - 2000)/980*49)
  # Now, find interpolated fluxes:
  flux1 = star_grid[it1,jsurf[2:3]]*(tgrid[it1+1]-T1)/20.0+star_grid[it1+1,jsurf[2:3]]*(T1-tgrid[it1])/20.0
  flux2 = star_grid[it2,jsurf[2:3]]*(tgrid[it2+1]-T2)/20.0+star_grid[it2+1,jsurf[2:3]]*(T2-tgrid[it2])/20.0
  flux3 = star_grid[it3,jsurf[2:3]]*(tgrid[it3+1]-T3)/20.0+star_grid[it3+1,jsurf[2:3]]*(T3-tgrid[it3])/20.0
  if (f2+f3) < 2//3
    return 1 ./((1 - f2 - f3) .+ f2 .*(flux2 ./ flux1) .+ f3 .*(flux3 ./ flux1))
  elseif (1-f2) < 2//3
    return 1 ./(f2 .+ (1 - f2 - f3) .*(flux1 ./ flux2) .+ f3 .*(flux3 ./ flux2))
  else
    return 1 ./(f3 .+ (1 - f2 - f3) .*(flux1 ./ flux3) .+ f2 .*(flux2 ./ flux3))
  end
end

# Okay, now run a Markov chain:
using AffineInvariantMCMC

# Set up initial conditions for the chain:
nparam  = 12
nwalker = 40
astep   = 2.0
xstart = zeros(nparam,nwalker)
for i=1:nwalker
  param = zeros(nparam)
  param[1] = 2000.0 + 980.0 * rand()
  param[2] = 2000.0 + 980.0 * rand()
  param[3] = 2000.0 + 980.0 * rand()
  param[4] = rand()
  param[5] = rand()*(1-param[4])
  param[6:12] = rand(7)
  xstart[:,i] = param
end

# Optimize the model:
nllmin = Inf
pmin = zeros(12)
for i=1:nwalker
  xtmp = zeros(12)
  xprior = xstart[:,i]
  while maximum(abs.(xtmp-xprior)) > 1e-4
    xtmp = xprior
    result = optimize((x)-> -log_likelihood(x),xtmp)
    xprior = result.minimizer
  end
  nlltmp = -log_likelihood(xtmp)
  if nlltmp < nllmin
    println(nlltmp," ",xtmp)
    pmin = xtmp
    global nllmin = nlltmp
  end
end

# Run the chain:
nsample = 100000
nthin   = 100
#nsample = 10000
#nthin   = 10
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
     dsamp[i,ip,:] =  tlam_samp .* param[5+ip]
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
  ax.semilogx(lam_mod[isort],zeros(16) .+ mean(flatchain[5+ip,:]),color="r")
  ax.semilogx(lam_mod[isort],zeros(16) .+ (mean(flatchain[5+ip,:]) - std(flatchain[5+ip,:])),color="r",linestyle=":")
  ax.semilogx(lam_mod[isort],zeros(16) .+ (mean(flatchain[5+ip,:]) + std(flatchain[5+ip,:])),color="r",linestyle=":")
  ax.set_xlabel(L"$\lambda$ [micron]")
  ax.set_ylabel("Depth [%]")
end

ibest = argmax(flatllhoodvals)
pbest = flatchain[:,ibest]
tlambest = compute_tlam(pbest)
for ip=1:7
  ax = axes[ip]
  ax.semilogx(lam_mod[isort],pbest[5+ip] .*tlambest[isort],color="b")
end

t1 = zeros(nflat)
t2 = zeros(nflat)
t3 = zeros(nflat)
f1 = zeros(nflat)
f2 = zeros(nflat)
f3 = zeros(nflat)
for i=1:nflat
  if flatchain[1,i] > flatchain[2,i] && flatchain[1,i] > flatchain[3,i]
    t1[i] = flatchain[1,i]; f1[i] = 1-flatchain[4,i]-flatchain[5,i]
    if flatchain[2,i] > flatchain[3,i]
      t2[i] = flatchain[2,i]; f2[i] = flatchain[4,i]; t3[i] = flatchain[3,i]; f3[i]=flatchain[5,i]
    else
      t2[i] = flatchain[3,i]; f2[i] = flatchain[5,i]; t3[i] = flatchain[2,i]; f3[i]=flatchain[4,i]
    end
  elseif flatchain[2,i] > flatchain[1,i] && flatchain[2,i] > flatchain[3,i]
    t1[i] = flatchain[2,i]; f1[i] = flatchain[4,i]
    if flatchain[1,i] > flatchain[3,i]
      t2[i] = flatchain[1,i]; f2[i] = 1-flatchain[4,i]-flatchain[5,i]; t3[i] = flatchain[3,i]; f3[i]=flatchain[5,i]
    else
      t2[i] = flatchain[3,i]; f2[i] = flatchain[5,i]; t3[i] = flatchain[1,i]; f3[i]=1-flatchain[4,i]-flatchain[5,i]
    end
  else
    t1[i] = flatchain[3,i]; f1[i] = flatchain[5,i]
    if flatchain[1,i] > flatchain[2,i]
      t2[i] = flatchain[1,i]; f2[i] = 1-flatchain[4,i]-flatchain[5,i]; t3[i] = flatchain[2,i]; f3[i]=flatchain[4,i]
    else
      t2[i] = flatchain[2,i]; f2[i] = flatchain[4,i]; t3[i] = flatchain[1,i]; f3[i]=1-flatchain[4,i]-flatchain[5,i]
    end
  end
end
