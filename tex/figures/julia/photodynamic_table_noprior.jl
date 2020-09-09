

# Create plots and a LaTeX table for the parameters derived from the photodynamic
# model.

using Statistics
using Printf
using JLD2; using PyPlot
using MCMCDiagnostics

if ~@isdefined(CGS)
  include("../../../src/CGS.jl")
  using Main.CGS
end

# This is run from T1_photodynamics_MCMC_02.ipynb
# (200,000 steps; 50 chains took ~3 days to run):
@load "../../../data/T1_pd_MCMC_run_000_001.jld2"
chain001 = chain
flatchain001 = flatchain
@load "../../../data/T1_pd_MCMC_run_000_002.jld2"
chain002 = chain
flatchain002 = flatchain
@load "../../../data/T1_pd_MCMC_run_000_003.jld2"
chain003 = chain
flatchain003 = flatchain
# (This doesn't have variable contamination factor).

#istart = 15000; iend = 500000
istart = 1001; iend = 200000
flatchain = hcat(flatchain001[:,istart:iend],flatchain002[:,istart:iend],flatchain003[:,istart:iend])
nsamples = size(flatchain)[2]

iburn = 201; ilength = 4000
nwalker = 50
# number of effective samples:
for i=1:19
 neff_tot = 0
 for j=1:nwalker
   neff_tot += effective_sample_size(chain001[i,j,iburn:ilength])
   neff_tot += effective_sample_size(chain002[i,j,iburn:ilength])
   neff_tot += effective_sample_size(chain003[i,j,iburn:ilength])
 end
 println(i," ",neff_tot)
end 

# First, make histograms of some parameters:

function histogram(param,nbin)
  p1 = minimum(param)-1e-15; p2 = maximum(param)+1e-15
  pbin_square = [p1]
  hist_square = [0.0]
  pbin = zeros(nbin)
  hist = zeros(nbin)
  psort = sort(param)
  i1 = 1; np = size(param)[1]
  #println(p1," ",psort[i1]," ",p2," ",psort[np])
  for i=1:nbin
    while psort[i1] <= p1+(p2-p1)/nbin*i
      hist[i] += 1.0
      i1 += 1
      if i1 == np
        break
      end
    end
    push!(hist_square,hist[i]); push!(hist_square,hist[i])
    pbin[i] = p1+(p2-p1)/nbin*(i-0.5)
    push!(pbin_square,p1+(p2-p1)/nbin*(i-1))
    push!(pbin_square,p1+(p2-p1)/nbin*i)
    if i1 == np
      break
    end
  end
  push!(hist_square,0.0)
  push!(pbin_square,p2)
  return pbin,hist,pbin_square,hist_square
end

# plot histogram of each parameter:
pname = ["b","c","d","e","f","g","h"]

clf()
nplanet = n-1
rp_mean = zeros(nplanet); rp_sig = zeros(nplanet)
depth_mean = zeros(nplanet); depth_sig = zeros(nplanet)
dur_mean = zeros(nplanet); dur_sig = zeros(nplanet)
tau_mean = zeros(nplanet); tau_sig = zeros(nplanet)
med_b = zeros(nplanet); b_sig = zeros(nplanet)
# Normalized velocity from the dynamical model (edge-on) computed
# at mid-transit for trappist-1 in units of AU/day:
vnorm = [0.10676978639803744, 0.09139150615416661, 0.07710600470245187, 0.06732285130838131, 0.0586281758248469, 
         0.05318938396566236, 0.04623982742725911]
rpstring = ""
depthstring = ""
durstring = ""
bstring = ""
taustring = ""
for i=1:nplanet
  # Impact parameter:
  bcurr = sort(flatchain[i+nplanet,:])
  b_bin, b_hist, b_bin_square,b_hist_square = histogram(bcurr,50)
  plot(b_bin,b_hist./maximum(b_hist),label=pname[i],linewidth=2)
  rcurr= flatchain[i,:]
  rp_mean[i] = mean(rcurr)
  rp_sig[i] = std(rcurr)
  depth_mean[i] = mean(flatchain[i,:].^2)
  depth_sig[i] = std(flatchain[i,:].^2)
  med_b[i] = median(bcurr)
  sig_b_low = med_b[i] - bcurr[ceil(Int64,0.158655*nsamples)]
  sig_b_high = bcurr[ceil(Int64,0.841345*nsamples)]-med_b[i]
  # Compute the transit duration:
  fac = flatchain[15,:].^(1//3)*AU/RSUN
  dur = 2sqrt.((1 .+flatchain[i,:]).^2 .- flatchain[i+nplanet,:].^2)./(vnorm[i]*fac)*(24*60)
  dur_mean[i] = mean(dur)
  dur_sig[i] = std(dur)
  tau = (sqrt.((1 .+flatchain[i,:]).^2 .- flatchain[i+nplanet,:].^2) - 
          sqrt.((1 .-flatchain[i,:]).^2 .- flatchain[i+nplanet,:].^2))./(vnorm[i]*fac)*(24*60)
  tau_mean[i] = mean(tau)
  tau_sig[i] = std(tau)
  global rpstring = string(rpstring," \$",@sprintf("%6.4f",rp_mean[i])," \\pm ",@sprintf("%6.4f",rp_sig[i]),"\$")
  global depthstring = string(depthstring," \$",@sprintf("%6.4f",depth_mean[i]*100)," \\pm ",@sprintf("%6.4f",depth_sig[i]*100),"\$")
  global durstring = string(durstring," \$",@sprintf("%6.4f",dur_mean[i])," \\pm ",@sprintf("%6.4f",dur_sig[i]),"\$")
  global taustring = string(taustring," \$",@sprintf("%6.4f",tau_mean[i])," \\pm ",@sprintf("%6.4f",tau_sig[i]),"\$")
  if i < nplanet
    global rpstring = string(rpstring," &")
    global depthstring = string(depthstring," &")
    global durstring = string(durstring," &")
    global taustring = string(taustring," &")
  else
    global rpstring = string(rpstring," \\cr")
    global depthstring = string(depthstring," \\cr")
    global durstring = string(durstring," \\cr")
    global taustring = string(taustring," \\cr")
  end
#  global bstring = string(bstring," \$",@sprintf("%6.4f",b_mean[i])," \\pm ",@sprintf("%6.4f",b_sig[i]),"\$")
  global bstring = string(bstring," \$",@sprintf("%6.4f",med_b[i]),"_{- ",@sprintf("%6.4f",sig_b_low),"}^{+",@sprintf("%6.4f",sig_b_high),"}\$")
  if i < nplanet
    global bstring = string(bstring," &")
  else
    global bstring = string(bstring," \\cr")
  end
end
legend()
axis([0.0,0.5,0,1.1])
#xlabel("Radius ratio"); ylabel("Probability")
xlabel("Impact parameter",fontsize=15); ylabel("Probability",fontsize=15)
savefig("../impact_parameter_noprior.pdf",bbox_inches="tight")

# Now, create a table of values:


# Compute 68% confidence interval and median of density:
density = sort(flatchain[15,:])
med_dens = median(density); sig_low = med_dens - density[ceil(Int64,0.158655*nsamples)]
sig_high = density[ceil(Int64,0.841345*nsamples)]-med_dens
#println("\$\\rho_*/\\rho_\\odot\$ & \$",@sprintf("%6.4f",med_dens),"\_\{- ",@sprintf("%6.4f",sig_low)"\}\^\{+",@sprintf("%6.4f",sig_high),"\$ & & & & & &\\cr")

# Also derive a/R_* and inclinations of the orbits, and add to the table. [ ]

q11 = flatchain[16,:]
q21 = flatchain[17,:]
q12 = flatchain[18,:]
q22 = flatchain[19,:]
#println("Parameter & \$\\rho_*/\\rho_\\odot\$ & \$q_{1,Ch 1}\$ & \$q_{2,Ch 1}\$ & \$ q_{1,Ch 2}\$ & \$q_{2,Ch 2}\$ & \$ & &\\cr")
#println("Value & \$",@sprintf("%6.4f",med_dens),"_{- ",@sprintf("%6.4f",sig_low),"}^{+",@sprintf("%6.4f",sig_high),"}\$ & \$",@sprintf("%6.4f",mean(q11))," \\pm ",@sprintf("%6.4f",std(q11))," \$ & \$",@sprintf("%6.4f",mean(q21))," \\pm ",@sprintf("%6.4f",std(q21))," \$ & \$",@sprintf("%6.4f",mean(q12))," \\pm ",@sprintf("%6.4f",std(q12))," \$ & \$",@sprintf("%6.4f",mean(q22))," \\pm ",@sprintf("%6.4f",std(q22))," \$ & &\\cr")
#println("\\hline")
#println(" Planet & b & c & d & e & f & g \\cr")
#println("\$R_p/R_*\$ & ",rpstring)
#println("\$b/R_*\$ & ",bstring)


#read(stdin,Char)
clf()
@load "../../../data/T1_density_incprior.jld2"

rho_bin, rho_hist, rho_bin_square,rho_hist_square = histogram(density,100)
plot(rho_bin,rho_hist./maximum(rho_hist),linewidth=2,label="No prior")
plot(rho_bin_incprior[1:end-2],rho_hist_incprior[1:end-2]./maximum(rho_hist_incprior),linestyle="--",linewidth=2,label="Inclination prior")
xlabel(L"$\rho_*/\rho_\odot$",fontsize=15)
ylabel("Probability",fontsize=15)
legend(fontsize=15)
axis([47,56,0,1.2],fontsize=15)
savefig("../stellar_density_noprior.pdf",bbox_inches="tight")


#plot(log.(rho_bin_square./54.0),log.(rho_hist_square))

#loglog(rho_bin_square./54,(rho_bin_square/54.0).^61 .* exp(12)./(exp.(rho_bin_square/54) .-1))

#loglog(rho_bin_square./54,(rho_bin_square ./54).^3 ./(exp.(rho_bin_square ./54).^3 .-1))

#loglog(rho_bin_square./54,(rho_bin_square ./54).^60 .*exp(11))
#loglog(rho_bin_square./54,exp(-(rho_bin_square./54))*1e5)


#rho_bin_min = copy(rho_bin_square)
#rhoc = 53.7
#plot(rho_bin_square ./ rhoc,rho_hist_square)
#rho_bin_min[rho_bin_min .< rhoc] .= rhoc
#semilogy(rho_bin_square./54,exp.((rho_bin_square .-45) .*1.25))
#plot(rho_bin_square ./ 54,exp.(-((rho_bin_min ./ 54) .-1).^2 .* 0.5 .*1e4)*6e4)
# This doesn't do too bad of a job in reproducing the distribution:
#plot(rho_bin_square./rhoc,exp.((rho_bin_square .-45) .*1.25).*exp.(-((rho_bin_min ./ rhoc) .-1).^2 .* 0.5 .*1e4))

# Plot histograms of a/R_* for all planets:
#read(stdin,Char)
#clf()

P = [1.510826,2.421938,4.049218,6.101013,9.207541,12.352445,18.772863]  # Orbit periods

aonr0 = (GRAV*MSUN*(24*3600)^2/(4pi^2*RSUN^3))^(1//3) 

med_dens = median(density); dens_low = density[ceil(Int64,0.158655*nsamples)]; dens_high = density[ceil(Int64,0.841345*nsamples)]
aonr_string = "\$a/R_*\$ & \$"
for i=1:7
  aonr_range = [dens_low,med_dens,dens_high].^(1//3) .* aonr0 .* P[i]^(2//3)
  println("aon_range: ",aonr_range)
  global aonr_string = string(aonr_string,@sprintf("%5.3f",aonr_range[2]),"_{- ",@sprintf("%5.3f",aonr_range[2]-aonr_range[1]),"}^{+",@sprintf("%5.3f",aonr_range[3]-aonr_range[2]),"}\$")
  if i < 7
    aonr_string = string(aonr_string," & \$")
  end
end 
aonr_string = string(aonr_string," \\cr")

#  aonr = density.^(1//3) .* (aonr0 * P[i]^(2//3))
#  aonr_bin,aonr_hist,aonr_bin_sq,aonr_hist_sq = histogram(aonr,50)
#  plot(aonr_bin_sq,aonr_hist_sq,label=pname[i])
#end
#xlabel(L"$a/R_*$")
#ylabel(L"Probability")
#legend()

# Derive the inclinations of the planets:
#read(stdin,Char)
clf()
nplanet = n-1
inc_mean = zeros(nplanet); inc_sig = zeros(nplanet)
incstring = " \$ I \$(\$^\\circ\$) & "
cname = ["C0","C1","C2","C3","C4","C5","C6","C7"]
nb = [30,30,30,40,40,40,40,40]
for i=1:nplanet
  density = sort(flatchain[15,:])
  aonr = density.^(1//3) .* (aonr0 * P[i]^(2//3))
  bcurr = flatchain[i+nplanet,:]
  inc_planet = acos.(bcurr ./ aonr) .* (180/pi)
  inc_bin, inc_hist, inc_bin_square,inc_hist_square = histogram(inc_planet,nb[i])
#  plot(inc_bin_square,inc_hist_square./maximum(inc_hist_square),label=pname[i],color=cname[i])
  plot([inc_bin;reverse(180 .-inc_bin)],[inc_hist;reverse(inc_hist)]./maximum(inc_hist),label=pname[i],color=cname[i],linewidth=2)
#  plot(180 .-inc_bin_square,inc_hist_square./maximum(inc_hist_square),color=cname[i])
  inc_mean[i] = mean(inc_planet)
  inc_sig[i] = std(inc_planet)
  global incstring = string(incstring," \$",@sprintf("%6.3f",inc_mean[i])," \\pm ",@sprintf("%6.3f",inc_sig[i]),"\$")
  if i < nplanet
    global incstring = string(incstring," & ")
  else
    global incstring = string(incstring," \\cr")
  end
end
legend(fontsize=15)
axis([89.25,90.75,0,1.1],fontsize=15)
#xlabel("Radius ratio"); ylabel("Probability")
xlabel("Inclination [deg]",fontsize=15); ylabel("Probability",fontsize=15)
tight_layout()
savefig("../inclination_noprior.pdf",bbox_inches="tight")
#read(stdin,Char)

# Limb-darkening:
clf()
u21 = sqrt.(q11).*(1 .-2q21)
u11 = 2*sqrt.(q11).*q21
u22 = sqrt.(q12).*(1 .-2q22)
u12 = 2*sqrt.(q12).*q22
#plot(u11,u21,color="r",".",alpha=0.005)
#plot(u12,u22,color="g",".",alpha=0.005)
errorbar([0.1633],[0.2549],xerr=[0.0364],yerr=[0.057],color="k",linewidth=3)
plot([0.1633],[0.2549],color="r","o",label="Ch 1")
errorbar([0.1442],[0.2173],xerr=[0.0324],yerr=[0.0482],color="k",linewidth=3)
plot([0.1442],[0.2173],color="g","o",label="Ch 2")
xlabel(L"$u_1$",fontsize=15); ylabel(L"$u_2$",fontsize=15)
legend()
include("../../../src/contour_limb.jl")
axis([u1i,u1f,u2i,u2f])
tight_layout()
savefig("../limb_darkening_nouprior.pdf",bbox_inches="tight")
rho_sun = MSUN/(4pi/3*RSUN^3)


println("Parameter & \$\\rho_*/\\rho_\\odot\$ & \$q_{1,Ch 1}\$ & \$q_{2,Ch 1}\$ & \$ q_{1,Ch 2}\$ & \$q_{2,Ch 2}\$ &  &\\cr")
println("Value & \$",@sprintf("%6.4f",med_dens),"_{- ",@sprintf("%6.4f",sig_low),"}^{+",@sprintf("%6.4f",sig_high),"}\$ & \$",@sprintf("%6.4f",mean(q11))," \\pm ",@sprintf("%6.4f",std(q11))," \$ & \$",@sprintf("%6.4f",mean(q21))," \\pm ",@sprintf("%6.4f",std(q21))," \$ & \$",@sprintf("%6.4f",mean(q12))," \\pm ",@sprintf("%6.4f",std(q12))," \$ & \$",@sprintf("%6.4f",mean(q22))," \\pm ",@sprintf("%6.4f",std(q22))," \$ & &\\cr")
#println("Parameter & \$ \\rho_*\$ \[g/cm\$\^2\$\] & \$u_\\mathrm{1,Ch 1}\$ & \$u_\\mathrm{2,Ch 1}\$ & \$u_\\mathrm{1,Ch 2}\$ & \$u_\\mathrm{2,Ch 2}\$ &  &\\cr")
println("Parameter & \$ \\rho_*\$ [g / cm\$^3\$] & \$u_\\mathrm{1,Ch 1}\$ & \$u_\\mathrm{2,Ch 1}\$ & \$u_\\mathrm{1,Ch 2}\$ & \$u_\\mathrm{2,Ch 2}\$ &  & \\cr")

println("Value & \$",@sprintf("%6.4f",med_dens*rho_sun),"_{- ",@sprintf("%6.4f",sig_low*rho_sun),"}^{+",@sprintf("%6.4f",sig_high*rho_sun),"}\$ & \$",@sprintf("%6.4f",mean(u11))," \\pm ",@sprintf("%6.4f",std(u11))," \$ & \$",@sprintf("%6.4f",mean(u21))," \\pm ",@sprintf("%6.4f",std(u21))," \$ & \$",@sprintf("%6.4f",mean(u12))," \\pm ",@sprintf("%6.4f",std(u12))," \$ & \$",@sprintf("%6.4f",mean(u22))," \\pm ",@sprintf("%6.4f",std(u22))," \$ & & \\cr")
println("\\hline")
println(" Planet & b & c & d & e & f & g \\cr")
println("\$R_p/R_*\$ & ",rpstring)
println("Depth [\\%] & ",depthstring)
println("T [min] & ",durstring)
println("\$\\tau\$ [min] & ",taustring)
println("\$b/R_*\$ & ",bstring)
println(aonr_string)
println(incstring)
