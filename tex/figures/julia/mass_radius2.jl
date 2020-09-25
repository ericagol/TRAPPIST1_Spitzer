

using PyPlot
if !@isdefined(CGS)
  include("../../../src/CGS.jl")
  using Main.CGS
end
include("../../../src/seager_mr.jl")

include("../../../src/loglinspace.jl")

#function mass_radius2(mass_ratio_total,mstar,smstar,flatchain,logplot)

nplanet = 7
# Give names to planets:
planet_name=["b","c","d","e","f","g","h"]
#cname = ["b","g","r","c","m","y","k"]
cname = ["C0","C1","C2","C3","C4","C5","C6","C7"]

r1 = 0.6
r2 = 1.2
nr = 200
if logplot
 m1 = 0.01
 m2 = 4.0
else
 m1 = 0.20
 m2 = 1.6
end
nm = 600
mrgrid = zeros(nplanet,nr,nm)

rgrid = zeros(nr)
mgrid = zeros(nm)
for i=1:nr
  rgrid[i] = (i-0.5)*(r2-r1)/nr+r1
end
for i=1:nm
  if logplot
    mgrid[i] = 10 .^((i-0.5)*(log10(m2)-log10(m1))/nm)*m1
  else
    mgrid[i] = (i-0.5)*(m2-m1)/nm+m1
  end
end

# Set up plotting window in Matplotlib (called by PyPlot):
#fig,axes = subplots()

# Number of posterior samples to select:
nsamp = 10000000
npoints = [100,100,100,100,100,100,1000]
i0 = 20001
# Factor to convert from solar masses to Earth masses:
fac_mass = MSUN/MEARTH
# Factor to convert from solar radii to Earth radii:
fac_rad = RSUN/REARTH

# Average mass of star in solar masses:
#mstar = 0.09 # Mass of star
# Stand deviation of mass of star in solar masses:
#smstar  = 0.01 # Uncertainty on stellar mass

# Loop over planets:
#conf = [0.997, 0.95, 0.683]
conf = [0.95, 0.683]
nconf = length(conf)
level_mr = zeros(nplanet,nconf)
mavg = zeros(nplanet)
merror = zeros(nplanet)
ravg = zeros(nplanet)
rerror = zeros(nplanet)
# Photodynamic size:
pdsize = size(flatchain)[2]
# N-body HMC size:
nhmc = size(mass_ratio_total)[2]
# Use latest markov chain radius ratios:
# Set up vectors for samples:
msamp = zeros(nplanet,nsamp)
rsamp = zeros(nplanet,nsamp)
mstar_samp = zeros(nsamp)
rho_samp = zeros(nsamp)
rstar_samp = zeros(nsamp)
# Loop over number of samples:
for isamp=1:nsamp
  # Draw from posterior distribution:
  mplanet = mass_ratio_total[:,ceil(Int64,rand()*nhmc)]
  # Select a stellar mass:
  mstar_samp[isamp] = mstar+randn()*smstar
  # Compute the mass of planet in Earth mass units:
  msamp[:,isamp] = mplanet .*(mstar_samp[isamp]*fac_mass)
  # Select a density (relative to Solar):
  ipd = ceil(Int64,rand()*pdsize)
  rho_samp[isamp] = flatchain[15,ipd]
  # Compute the radius of the star in solar radii:
  rstar_samp[isamp] = (mstar_samp[isamp]/rho_samp[isamp])^(1/3)
  # Compute the radius of planet in Earth radii units:
  # Randomly select a radius ratio:
  rsamp[:,isamp] = flatchain[1:nplanet,ipd] .*(rstar_samp[isamp]*fac_rad)
  # Now compute density distribution:
  for j=1:nplanet
    if msamp[j,isamp] > m1 && msamp[j,isamp] < m2 && rsamp[j,isamp] > r1 && rsamp[j,isamp] < r2
      ir = ceil(Int64,(rsamp[j,isamp]-r1)/(r2-r1)*nr)
      if logplot
        im = ceil(Int64,(log(msamp[j,isamp])-log(m1))/(log(m2)-log(m1))*nm)
      else
        im = ceil(Int64,(msamp[j,isamp]-m1)/(m2-m1)*nm)
      end
      mrgrid[j,ir,im] +=1.0 
    end
  end
end
for j=1:nplanet
  # Set levels for each planet:
  mrgrid[j,:,:] /= float(sum(mrgrid[j,:,:]))
  mrgrid_sort = sort(vec(mrgrid[j,:,:]))
  i0 = nm*nr
  mrgrid_cum = mrgrid_sort[i0]
  while mrgrid_cum < conf[1] && i0 > 1
    i0 -= 1
    for k=1:nconf
      if mrgrid_cum < conf[k] && (mrgrid_cum+mrgrid_sort[i0]) > conf[k]
        level_mr[j,k] = 0.5*(mrgrid_sort[i0]+mrgrid_sort[i0+1])
      end
    end
    mrgrid_cum += mrgrid_sort[i0]
  end
  println(j," Mass:   ",mean(msamp[j,:])," ",std(msamp[j,:])," fractional error: ",std(msamp[j,:])/mean(msamp[j,:])*100.0,"%")
  println(j," Radius: ",mean(rsamp[j,:])," ",std(rsamp[j,:])," fractional error: ",std(rsamp[j,:])/mean(rsamp[j,:])*100.0,"%")
  merror[j] = std(msamp[j,:])
  mavg[j] = mean(msamp[j,:])
  rerror[j] = std(rsamp[j,:])
  ravg[j] = mean(rsamp[j,:])
end
# Plot maximum likelihood values:
#mopt = readdlm("optimum_values.txt")
#ropt = [0.08581987985664834,0.08461786201718988,0.060328647716537474,0.07106532992445118, 0.08035101963995364, 0.08570949401350199, 0.057994985982196025, 53.833783736830554]
for j=1:nplanet
#  cs = axes.contour(mgrid,rgrid,mrgrid[j,:,:],levels = level_mr[j,:],colors=cname[j],linewidth=2.0,label=planet_name[j])
  cs = contour(mgrid,rgrid,mrgrid[j,:,:],levels = level_mr[j,:],colors=cname[j],linewidth=2.0)
#  plot(mopt[(j-1)*5+1]*fac_mass*mstar,ropt[j]*(mstar/ropt[8])^(1/3)*fac_rad,"o",label=planet_name[j],color=cname[j])
  plot(mavg[j],ravg[j],"o",label=planet_name[j],color=cname[j])
#  axes[:imshow](mrgrid[j,:,:],alpha=0.1)
end
# Make some plots:
#axes.set_xlabel(L"Mass $[M_\oplus]$")
xlabel(L"Mass $[M_\oplus]$")
#axes.set_ylabel(L"Radius $[R_\oplus]$")
ylabel(L"Radius $[R_\oplus]$")
#axes.axis([m1,m2,r1,r2])
axis([m1,m2,r1,r2])
text(1.0,1.00,"Earth")
text(0.85,0.93,"Venus")
if logplot
  #axes.set_xscale("log")
  set_xscale("log")
end
#mm = logarithmspace(-2.5,0.2,1000)
#plot(mm,seager_mr(mm,"H2O"),label="Water",linewidth=3,linestyle="--")
#plot(mm,seager_mr(mm,"Fe"),label="Iron",linewidth=3,linestyle="--")
#plot(mm,seager_mr(mm,"Earth"),label="Earth",linewidth=3,linestyle="--")
#plot(mm,seager_mr(mm,"MgSiO3"),label="Rock",linewidth=3,linestyle="--")
#axes[:legend](loc = "upper left",fontsize=10)
#axes.legend(loc = "upper left",fontsize=10)
#legend(loc = "upper left",fontsize=10)
#return merror
#end
