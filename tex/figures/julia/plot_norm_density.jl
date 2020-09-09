

# Recreate plot with normalized density vs. orbital period:

include("../../../src/regress.jl")
using Statistics
using PyPlot
using DelimitedFiles

period = [1.510826, 2.421937,4.049219, 6.101013, 9.207540,12.352446,18.772866]

dir = "../../../data/POSTERIOR_NORM_DENSITY/"
fname = ["out_T1b_norm","out_T1c_norm","out_T1d_norm","out_T1e_norm","out_T1f_norm","out_T1g_norm","out_T1h_norm"]

np = 7; nsamp = 10000


dens_norm = zeros(np,nsamp)
dens_mean = zeros(np)
dens_sig = zeros(np)
for ip=1:np
  data = readdlm(string(dir,fname[ip]))
  dens_norm[ip,:] .= data[:,2]
  dens_mean[ip] = mean(dens_norm[ip,:])
  dens_sig[ip] = std(dens_norm[ip,:])
end

# Now, carry out regressions for each of these:
coeff_samp = zeros(2,nsamp)
fn = zeros(2,np)
fn[1,:] .= 1.0
fn[2,:] .= period
sig = ones(np)
dens_avg = 0.0
for isamp=1:nsamp
  global dens_avg += mean(dens_norm[:,isamp])
  coeff,cov = regress(fn,dens_norm[:,isamp],sig)
  coeff_samp[:,isamp] = coeff
end
dens_avg /= nsamp
# Now, compute standard deviation:
dens_var = 0.0
for isamp=1:nsamp
  global dens_var += (mean(dens_norm[:,isamp]) - dens_avg)^2
end
dens_sig = sqrt(dens_var/nsamp)


# Now, draw 1-sigma and 2-sigma confidence intervals:
period_trial = collect(0.0:0.1:20.0)
conf = zeros(5,length(period_trial))
for (i,p) in enumerate(period_trial)
  # Predict samples at this period
  csamp = coeff_samp[1,:] .+ coeff_samp[2,:] .* p
  # Now, sort these and pick out 68.3% & 95% confidence limits:
  csamp = sort(csamp)
  conf[:,i] = [csamp[250],csamp[1415],csamp[5000],csamp[8415],csamp[9750]]
end
# Now, plot median, 1-sigma and 2-sigma confidence intervals:

clf()
#plot(period_trial,conf[1,:],color="C1",alpha=0.5,linestyle="--")
errorbar(period,dens_mean,yerr=dens_sig,fmt="o")
#plot(period_trial,conf[5,:],color="C1",alpha=0.5,linestyle="--")
plot(period_trial,conf[3,:],color="C1",linewidth=2)
plot(period_trial,conf[2,:],color="C1",linewidth=2,linestyle=":")
plot(period_trial,conf[4,:],color="C1",linewidth=2,linestyle=":")
#plot(period_trial,dens_avg  .+ 0.0 .* period_trial,color="C0",linewidth=2)
fill_between([0,20],[dens_avg-dens_sig,dens_avg-dens_sig],[dens_avg+dens_sig,dens_avg+dens_sig],alpha=0.2,color="C0")
axis([0,20,0.8,1.2])
xlabel("Orbital Period [d]")
ylabel("Normalized density")
tight_layout()
savefig("../Norm_dens_vs_period.pdf",bbox_inches = "tight")
