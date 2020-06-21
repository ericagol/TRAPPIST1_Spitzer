

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

#@load "T1_pd_AIMC_001.jld2"
# This is run from T1_photodynamics_MCMC_02.ipynb
# (50,000 steps; 100 chains took ~3 days to run):
# This needs to be updated with the corrected chains:
#@load "T1_pd_MCMC_004.jld2"
@load "../../../data/T1_pd_MCMC_run_001_001.jld2"
chain001 = chain
flatchain001 = flatchain
@load "../../../data/T1_pd_MCMC_run_001_002.jld2"
chain002 = chain
flatchain002 = flatchain
@load "../../../data/T1_pd_MCMC_run_001_003.jld2"
chain003 = chain
flatchain003 = flatchain

istart = 35001; iend = 200000
flatchain = hcat(flatchain001[:,istart:iend],flatchain002[:,istart:iend],flatchain003[:,istart:iend])

# (This doesn't have variable contamination factor).

#istart = 15000; iend = 500000
iburn = 701; ilength = 4000
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
med_b = zeros(nplanet); b_sig = zeros(nplanet)
rpstring = ""
bstring = ""
for i=1:nplanet
#  rp_bin, rp_hist, rp_bin_square,rp_hist_square = histogram(flatchain[i,:],40)
#  plot(rp_bin_square,rp_hist_square./maximum(rp_hist_square),label=pname[i])
  bcurr = sort(flatchain[i+nplanet,:])
  b_bin, b_hist, b_bin_square,b_hist_square = histogram(bcurr,50)
  plot(b_bin,b_hist./maximum(b_hist),label=pname[i])
  rcurr= flatchain[i,:]
  #rp_mean[i] = mean(flatchain[i,istart:iend])
  rp_mean[i] = mean(rcurr)
  #rp_sig[i] = std(flatchain[i,istart:iend])
  rp_sig[i] = std(rcurr)
#  b_mean[i] = mean(flatchain[i+nplanet,istart:iend])
  med_b[i] = median(bcurr)
  sig_b_low = med_b[i] - bcurr[ceil(Int64,0.158655*3*(iend-istart+1))]
  sig_b_high = bcurr[ceil(Int64,0.841345*3*(iend-istart+1))]-med_b[i]
#  b_sig[i] = std(flatchain[i+nplanet,istart:iend])
#  b_sig[i] = std(bcurr)
  global rpstring = string(rpstring," \$",@sprintf("%6.4f",rp_mean[i])," \\pm ",@sprintf("%6.4f",rp_sig[i]),"\$")
  if i < nplanet
    global rpstring = string(rpstring," &")
  else
    global rpstring = string(rpstring," \\cr")
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
axis([0.05,0.5,0,1.1])
#xlabel("Radius ratio"); ylabel("Probability")
xlabel("Impact parameter"); ylabel("Probability")

# Now, create a table of values:


#println("\$\\rho_*/\\rho_\\odot\$ & \$",@sprintf("%6.4f",mean(flatchain[15,istart:iend]))," \\pm ",@sprintf("%6.4f",std(flatchain[15,istart:iend])),"\$ & & & & & &\\cr")
# Compute 68% confidence interval and median of density:
density = sort(flatchain[15,:])
med_dens = median(density); sig_low = med_dens - density[ceil(Int64,0.158655*3*(iend-istart+1))]
sig_high = density[ceil(Int64,0.841345*3*(iend-istart+1))]-med_dens
#println("\$\\rho_*/\\rho_\\odot\$ & \$",@sprintf("%6.4f",med_dens),"\_\{- ",@sprintf("%6.4f",sig_low)"\}\^\{+",@sprintf("%6.4f",sig_high),"\$ & & & & & &\\cr")

# Also derive a/R_* and inclinations of the orbits, and add to the table. [ ]

println("\$\\rho_*/\\rho_\\odot\$ & \$",@sprintf("%6.4f",med_dens),"_{- ",@sprintf("%6.4f",sig_low),"}^{+",@sprintf("%6.4f",sig_high),"}\$ & & & & & &\\cr")
q11 = flatchain[16,:]
q21 = flatchain[17,:]
q12 = flatchain[18,:]
q22 = flatchain[19,:]
println("\$q_{1,Ch 1}\$ & \$",@sprintf("%6.4f",mean(q11))," \\pm ",@sprintf("%6.4f",std(q11))," \$ & & & & & &\\cr")
println("\$q_{2,Ch 1}\$ & \$",@sprintf("%6.4f",mean(q21))," \\pm ",@sprintf("%6.4f",std(q21))," \$ & & & & & &\\cr")
println("\$q_{1,Ch 2}\$ & \$",@sprintf("%6.4f",mean(q12))," \\pm ",@sprintf("%6.4f",std(q12))," \$ & & & & & &\\cr")
println("\$q_{2,Ch 2}\$ & \$",@sprintf("%6.4f",mean(q22))," \\pm ",@sprintf("%6.4f",std(q22))," \$ & & & & & &\\cr")
println("\\hline")
println(" Planet & b & c & d & e & f & g \\cr")
println("\$R_p/R_*\$ & ",rpstring)
println("\$b/R_*\$ & ",bstring)


#read(stdin,Char)
clf()
#rho_bin, rho_hist, rho_bin_square,rho_hist_square = histogram(flatchain[15,istart:iend],30)
rho_bin, rho_hist, rho_bin_square,rho_hist_square = histogram(density,100)
plot(rho_bin[1:end-2],rho_hist[1:end-2]./maximum(rho_hist))
xlabel(L"$\rho_*/\rho_\odot$")
ylabel("Probability")

rho_bin_incprior =  rho_bin
rho_hist_incprior = rho_hist
#@save "T1_density_incprior.jld2" rho_bin_incprior rho_hist_incprior

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

med_dens = median(density); dens_low = density[ceil(Int64,0.158655*(iend-istart+1))]
dens_high = density[ceil(Int64,0.841345*(iend-istart+1))]
aonr_string = "\$a/R_*\$ & \$"
for i=1:7
  aonr_range = [dens_low,med_dens,dens_high].^(1//3) .* aonr0 .* P[i]^(2//3)
  global aonr_string = string(aonr_string,@sprintf("%5.3f",aonr_range[2]),"_{- ",@sprintf("%5.3f",aonr_range[2]-aonr_range[1]),"}^{+",@sprintf("%5.3f",aonr_range[3]-aonr_range[2]),"}\$")
  if i < 7
    aonr_string = string(aonr_string," & \$")
  end
end 
aonr_string = string(aonr_string," \\cr")
println(aonr_string)

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
incstring = ""
cname = ["C0","C1","C2","C3","C4","C5","C6","C7"]
nb = [70,70,70,70,70,70,70,70]
for i=1:nplanet
  density = flatchain[15,:];
  aonr = density.^(1//3) .* (aonr0 * P[i]^(2//3))
  bcurr = flatchain[i+nplanet,:]
  inc_planet = acos.(bcurr ./ aonr) .* (180/pi)
  inc_bin, inc_hist, inc_bin_square,inc_hist_square = histogram(inc_planet,nb[i])
#  plot(inc_bin_square,inc_hist_square./maximum(inc_hist_square),label=pname[i],color=cname[i])
#  plot([inc_bin;reverse(180 .-inc_bin)],[inc_hist;reverse(inc_hist)]./maximum(inc_hist),label=pname[i],color=cname[i],linewidth=2)
  plot(inc_bin,inc_hist./maximum(inc_hist),label=pname[i],color=cname[i],linewidth=2)
#  plot(180 .-inc_bin_square,inc_hist_square./maximum(inc_hist_square),color=cname[i])
  inc_mean[i] = mean(inc_planet)
  inc_sig[i] = std(inc_planet)
  global incstring = string(incstring," \$",@sprintf("%6.4f",inc_mean[i])," \\pm ",@sprintf("%6.4f",inc_sig[i]),"\$")
  if i < nplanet
    global rpstring = string(incstring," &")
  else
    global rpstring = string(incstring," \\cr")
  end
end
legend(fontsize=15)
axis([89.5,90.1,0,1.1],fontsize=15)
#xlabel("Radius ratio"); ylabel("Probability")
xlabel("Inclination [deg] ",fontsize=15); ylabel("Probability",fontsize=15)
tight_layout()
savefig("../inclination_nouprior_incprior.pdf",bbox_inches="tight")
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


println("\$u_{1,Ch 1}\$ & \$",@sprintf("%6.4f",mean(u11))," \\pm ",@sprintf("%6.4f",std(u11))," \$ & & & & & &\\cr")
println("\$u_{2,Ch 1}\$ & \$",@sprintf("%6.4f",mean(u21))," \\pm ",@sprintf("%6.4f",std(u21))," \$ & & & & & &\\cr")
println("\$u_{1,Ch 2}\$ & \$",@sprintf("%6.4f",mean(u12))," \\pm ",@sprintf("%6.4f",std(u12))," \$ & & & & & &\\cr")
println("\$u_{2,Ch 2}\$ & \$",@sprintf("%6.4f",mean(u22))," \\pm ",@sprintf("%6.4f",std(u22))," \$ & & & & & &\\cr")

# Plot histogram of the inclination angle prior, \sigma_\theta:

st_bin, st_hist, st_bin_square,st_hist_square = histogram(log10.(flatchain[20,:]),49)
semilogx(10.0.^st_bin,st_hist./maximum(st_hist),linewidth=2)
xlabel(L"$\sigma_\theta$",fontsize=15)
ylabel(L"Probability")

sigtheta = sort(flatchain[20,:])
nsamp = size(sigtheta)[1]
i1 = ceil(Int64,0.158655*nsamp)
i2 = floor(Int64,0.841345*nsamp)
println("\$ \\sigma_\\theta = ",@sprintf("%5.3f",median(sigtheta)),"_{-",@sprintf("%5.3f",median(sigtheta)-sigtheta[i1]),"}^{+",@sprintf("%5.3f",sigtheta[i2]-median(sigtheta)),"}^\\circ \$")
