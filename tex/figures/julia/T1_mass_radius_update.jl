

if !@isdefined(CGS)
  include("../../../src/CGS.jl")
  using Main.CGS
end

using JLD2; using PyPlot
using DelimitedFiles
using Statistics

#rhostar = 51.14; sig_rho=0.67 # In units of solar density
#rhostar = 53.7; sig_rho=0.85 # In units of solar density
#rhostar = 53.308; sig_rho=0.0208 # In units of solar density
#rhostar = 51.14; sig_rho=0.001 # In units of solar density

# Artificially increase density of star to see what it
# would take to make an Earth-like density:
#rhostar *= 1.23

# Mass of star from Mann et al. (2019):  https://github.com/awmann/M_-M_K-
mstar = 0.0898; smstar = 0.0023  # Solar units
#mstar = 0.0898; smstar = 0.0001  # Solar units

# Latest masses of planets assuming a stellar mass of 0.09 M_sun.
# Take from plot_masses.jl in
# /Users/ericagol/Observing/Spitzer/Cycle13/TRAPPIST1/campaign03/April2019
# Includes delta-omega prior, includes mixture model for handling outliers.

#mp= [1.37024,1.31417,0.364033,0.652497,0.997529,1.28074,0.340946]
#sp= [0.0848279,0.0602858,0.022669,0.0230492,0.0202153,0.022538,0.0355592]
#mp = [1.3785, 1.3201, 0.3911, 0.7042, 1.0463, 1.3237, 0.3283]
#sp = [0.0530, 0.0397, 0.0076, 0.0116, 0.0128, 0.0146, 0.0190]

#mp = [1.37927,1.31014,0.388479,0.693422,1.04092,1.32341,0.325985]
#sp = [0.0596173,0.0465207,0.00748253,0.0127809,0.0154444,0.0168958,0.0185088]

# Load in the posterior from the transit-timing analysis:
#@load "/Users/ericagol/Observing/Spitzer/DDT2019/Campaign04/v09/Hyak/T1_hmc_total_02132020.jld2" state_total
@load "../../../data/T1_hmc_total_02212020.jld2" state_total

# Load in the posterior from the photodynamic analysis:
@load "../../../data/T1_pd_MCMC_run_000_001.jld2" flatchain
flatchain1 = flatchain[:,35001:200000]
#@load "/Users/ericagol/Observing/Spitzer/DDT2019/Campaign04/v06/Hyak/T1_pd_MCMC_run_000_002.jld2" flatchain
@load "../../../data/T1_pd_MCMC_run_000_002.jld2" flatchain
flatchain12 = hcat(flatchain1,flatchain[:,35001:200000])
#@load "/Users/ericagol/Observing/Spitzer/DDT2019/Campaign04/v06/Hyak/T1_pd_MCMC_run_000_003.jld2" flatchain
@load "../../../data/T1_pd_MCMC_run_000_003.jld2" flatchain
flatchain = hcat(flatchain12,flatchain[:,35001:200000])

# Make an array to hold posterior of mass-ratios:
nplanet = 7
mass_ratio_total = zeros(nplanet,size(state_total)[2])
for ip=1:nplanet
  mass_ratio_total[ip,:] .= state_total[(ip-1)*5+1,:]
end

periods = [1.51087637,2.42180746,4.049959,6.099043,9.205585,12.354473,18.767953]
t0 = [7322.51654,7282.80879,7670.14227,7660.37910,7671.39470,7665.35084,7662.55467]

logplot = false
include("mass_radius2.jl")

#merror = mass_radius2(mass_ratio_total,mstar,smstar,flatchain,false)

# Now, plot mass-radius relations from
# Zeng https://www.cfa.harvard.edu/~lzeng/tables/mrtable3.txt

dmr = readdlm("../../../data/zeng_mrtable3.txt",skipstart=2)

#plot(dmr[:,1],dmr[:,16],label="30% Fe",linestyle="--")  # 30% Fe
#plot(dmr[:,1],dmr[:,20],label="10% Fe",linestyle="--")  # 10% Fe
#plot(dmr[:,1],dmr[:,23],label=" 0% Fe",linestyle="--")  #  0% Fe
scatter([1.0,0.815],[1.0,0.95])

# Plot new Zeng rocky model:
# (from https://www.cfa.harvard.edu/~lzeng/planetmodels.html#mrtables)
dmr1 = readdlm("../../../data/zeng_rock.txt")
#plot(dmr1[1:18,1],dmr1[1:18,2],label="Rocky (Zeng et al. 2019)",linestyle="--")

# Plot new Zeng Earth-like mode:
dmr2 = readdlm("../../../data/massradiusEarthlikeRocky.txt")
plot(dmr2[1:23,1],dmr2[1:23,2],label="Earth-like (Zeng et al. 2019)",linestyle="--")

# Read in Carbon-rich planet model (for kicks):

#dmr_carbon = readdlm("../../../data/M_R_carbon-rich_Eric.ddat")
#plot(dmr_carbon[:,1],dmr_carbon[:,2],label="Carbon-rich (Miozzi et al. 2018)")

# Interpolate between the two Zeng models:
nm = 100
mm = logarithmspace(-2.5,0.2,nm)

function interpolate_models(iron,dmr1,dmr2,mm,nm)
  # Loop over each model and interpolate/extrapolate:
  r1 = zeros(nm)
  r2 = zeros(nm)
  r = zeros(nm)
  for i=1:nm
    if mm[i] < dmr1[1,1]
      i1 = 1; i2 = 2
    else
      for j=1:size(dmr1)[1]
        if mm[i] > dmr1[j,1] && mm[i] <= dmr1[j+1,1]
           i1 = j; i2 = j+1
        end
      end
    end
    r1[i] = dmr1[i1,2]*(dmr1[i2,2]/dmr1[i1,2])^((log(mm[i]/dmr1[i1,1]))/(log(dmr1[i2,1]/dmr1[i1,1])))
    if mm[i] < dmr2[1,1]
      i1 = 1; i2 = 2
    else
      for j=1:size(dmr2)[1]
        if mm[i] > dmr2[j,1] && mm[i] <= dmr2[j+1,1]
           i1 = j; i2 = j+1
        end
      end
    end
    r2[i] = dmr2[i1,2]*(dmr2[i2,2]/dmr2[i1,2])^(log(mm[i]/dmr2[i1,1])/log(dmr2[i2,1]/dmr2[i1,1]))
    # Now interpolate to the right iron content:
    r[i] = r1[i]+(r2[i]-r1[i])*iron/0.325
  end
  return r
end

#plot(mm,interpolate_models(0.10,dmr1,dmr2,mm,nm),label="10% iron (Zeng et al. 2019)",linestyle="--")
#plot(mm,interpolate_models(0.25,dmr1,dmr2,mm,nm),label="25% iron (Zeng et al. 2019)",linestyle="--")

legend(loc = "lower right",fontsize=10)

tight_layout()
#savefig("T1_mass_radius_Earth_carbon.pdf",bbox_inches="tight")
#read(stdin,Char)

## Semi-amplitudes:
fac = (2*pi*GRAV/(24*3600)/MSUN^2)^(1//3)*MEARTH
Kp = fac.*mavg./(mstar.^2 .*periods).^(1//3)
# Plot RV amplitude:
duration = 300
nrv = duration * 24
RV = zeros(nplanet+1,nrv)
time = zeros(nrv)
t0 .-= 7600.0
for i=1:nrv
  time[i] = i/24.0
  for j=1:nplanet
    RV[j,i] = Kp[j].*sin.((time[i].-t0[j])./periods[j])
  end
  RV[nplanet+1,i] = sum(RV[1:nplanet,i])
#    RV[j,i] = Kp[j].*sin.((time[i].-t0[j])./periods[j]))
end

clf()
plot(time,RV[nplanet+1,:]./100,label="Sum",color="k")
#plot(time,RV[nplanet+1,:],label="Sum",color="k")
xlabel("Time [days]")
#ylabel("RV [cm/sec]")
ylabel("RV [m/sec]")
#axis([0,duration,-1200,1200])
axis([0,duration,-12.00,12.00])
for i=1:7
#  plot(time,RV[i,:],label=planet_name[i],color=cname[i])
  plot(time,RV[i,:]./100,label=planet_name[i],color=cname[i])
#  errorbar([120+10*i],[-1000],Kp[i]*merror[i]/mavg[i],color=cname[i])
  errorbar([90+10*i],[-10],Kp[i]*merror[i]/mavg[i]/100,color=cname[i])
end
#text(40,-1020,"Equivalent RV precision:",fontsize=10)
text(10,-10.20,"Equivalent RV precision:",fontsize=10)
text(180,-10.2,"Measured RV precision: ")
errorbar([280],[-9.0],2.5,label="Hirano et al. (2020)",fmt="o")
legend(ncol=3)
tight_layout()
savefig("../Equivalent_RV_precision.pdf",bbox_inches="tight")

println("Semi-amplitudes: ",Kp," cm/s")
println("Precision: ",merror./mavg.*Kp," cm/s")
Kpsun = fac/(365.25)^(1//3)
println("Earth-Sun amplitude: ",Kpsun)
println("Equivalent Solar precision: ",Kpsun.*merror)

#include("/Users/ericagol/JWST/Exomoon/plot_moons.jl")

#@save "T1_mass_radius_cmf_000.jld2" nsamp msamp rsamp cmf rho_samp mstar_samp rstar_samp nr nm mrgrid rgrid mgrid level_mr conf m1 m2 r1 r2 mavg ravg merror rerror
