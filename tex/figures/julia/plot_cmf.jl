


# Recreate plot with CMF vs. orbital period:

include("../../../src/regress.jl")
using Statistics
using PyPlot
using DelimitedFiles
using Printf

period = [1.510826, 2.421937,4.049219, 6.101013, 9.207540,12.352446,18.772866]

dir = "../../../data/POSTERIOR_CMF/"
fname = ["T1b_IRON_FRAC.out","T1c_IRON_FRAC.out","T1d_IRON_FRAC.out","T1e_IRON_FRAC.out","T1f_IRON_FRAC.out","T1g_IRON_FRAC.out","T1h_IRON_FRAC.out"]

np = 7; nsamp = 10000

feonmg(cmf) = 1.778/(100/cmf-1)

cmf_norm = zeros(np,nsamp)
cmf_mean = zeros(np)
cmf_median = zeros(np)
cmf_sig = zeros(np)
cmf_sig1 = zeros(np)
cmf_sig2 = zeros(np)
feonmg_median = zeros(np)
feonmg_sig1 = zeros(np)
feonmg_sig2 = zeros(np)

for ip=1:np
  data = readdlm(string(dir,fname[ip]))
  cmf_norm[ip,:] .= data[:,2]
  cmf_mean[ip] = mean(cmf_norm[ip,:])
  cmf_median[ip] = median(cmf_norm[ip,:])
  cmf_sig[ip] = std(cmf_norm[ip,:])
  cmf_sig1[ip] = cmf_median[ip]-sort(cmf_norm[ip,:])[1587]
  cmf_sig2[ip] = sort(cmf_norm[ip,:])[8413]-cmf_median[ip]
  feonmg_median[ip] = median(feonmg.(cmf_norm[ip,:]))
  feonmg_sig1[ip] = feonmg_median[ip]-sort(feonmg.(cmf_norm[ip,:]))[1587]
  feonmg_sig2[ip] = sort(feonmg.(cmf_norm[ip,:]))[8413]-feonmg_median[ip]
end

# Now, carry out regressions for each of these:
coeff_samp = zeros(2,nsamp)
fn = zeros(2,np)
fn[1,:] .= 1.0
fn[2,:] .= period ./6
sig = ones(np)
cmf_avg = 0.0
for isamp=1:nsamp
  global cmf_avg += mean(cmf_norm[:,isamp])
  coeff,cov = regress(fn,cmf_norm[:,isamp],sig)
  coeff_samp[:,isamp] = coeff
end

cmf_avg /= nsamp

# Now, compute standard deviation:
cmf_var = 0.0
for isamp=1:nsamp
  global cmf_var += (mean(cmf_norm[:,isamp]) - cmf_avg)^2
end
cmf_rms = sqrt(cmf_var/nsamp)

# Now, draw 1-sigma and 2-sigma confidence intervals:
period_trial = collect(0.0:0.1:20.0)
conf = zeros(5,length(period_trial))
for (i,p) in enumerate(period_trial)
  # Predict samples at this period
  csamp = coeff_samp[1,:] .+ coeff_samp[2,:] .* p /6
  # Now, sort these and pick out 68.3% & 95% confidence limits:
  csamp = sort(csamp)
  conf[:,i] = [csamp[250],csamp[1415],csamp[5000],csamp[8415],csamp[9750]]
end
# Now, plot median, 1-sigma and 2-sigma confidence intervals:

clf()
#plot(period_trial,conf[1,:],color="C1",alpha=0.5,linestyle="--")
errorbar(period,cmf_mean,yerr=cmf_sig,fmt="o")
#errorbar(period,cmf_median,yerr=[cmf_sig1';cmf_sig2'],fmt="o")
#plot(period_trial,conf[5,:],color="C1",alpha=0.5,linestyle="--")
plot(period_trial,conf[3,:],color="C1",linewidth=2)
plot(period_trial,conf[2,:],color="C1",linewidth=2,linestyle=":")
plot(period_trial,conf[4,:],color="C1",linewidth=2,linestyle=":")
#plot(period_trial,cmf_avg  .+ 0.0 .* period_trial,color="C0",linewidth=2)
fill_between([0,20],[cmf_avg-cmf_rms,cmf_avg-cmf_rms],[cmf_avg+cmf_rms,cmf_avg+cmf_rms],alpha=0.2,color="C0")
xlabel("Orbital Period [d]")
ylabel("CMF/Iron mass fraction (%)")
# Overplot solar system bodies:
text([16],[28],"Mars")
plot([0,20],[30.0,30.0],"--",color="k")
text([1],[8],"Moon")
plot([0,20],[10.0,10.0],"--",color="k")
text([18],[33],"Earth")
plot([0,20],[32.0,32.0],"--",color="k")
#text([17],[31.5],"Venus")
#plot([0,20],[31.,31.],"--",color="k")
axis([0,20,0,40])
text([15],[20],"Chondrites (L/CI)")
#plot([0,20],[22.0,22.0],"--",color="k")
plot([0,20],[18.7,18.7],"--",color="k")
text([1],[20.6],"b")
text([2.7],[22.6],"c")
text([3.5],[13],"d")
text([6.3],[20],"e")
text([8.7],[15],"f")
text([11.8],[12],"g")
text([18.1],[11.8],"h")


tight_layout()
savefig("../Norm_cmf_vs_period.pdf",bbox_inches = "tight")

println("Coefficients of fit: ",mean(coeff_samp[1,:]),"+-",std(coeff_samp[1,:])," ",mean(coeff_samp[2,:]),"+-",std(coeff_samp[2,:]))

# Now, print out first line of Table 9:
tab09_01 = "CMF [wt\\%] & "
for i=1:7
  global tab09_01 = string(tab09_01," \$",@sprintf("%4.1f",cmf_median[i]),"_{-",@sprintf("%4.1f",cmf_sig1[i]),"}^{+",@sprintf("%4.1f",cmf_sig2[i]),"}\$ &")
end

tab09_01 = string(tab09_01," \$",@sprintf("%4.1f",cmf_avg),"\\pm",@sprintf("%4.1f",cmf_rms),"\$ \\cr")
println(tab09_01)

# Now, print out second line of Table 9:
tab09_02 = "Fe/Mg molar ratio & "
for i=1:7
  global tab09_02 = string(tab09_02," \$",@sprintf("%4.2f",feonmg_median[i]),"_{-",@sprintf("%4.2f",feonmg_sig1[i]),"}^{+",@sprintf("%4.2f",feonmg_sig2[i]),"}\$ &")
end

tab09_02 = string(tab09_02," \$",@sprintf("%4.2f",feonmg.(cmf_avg)),"\\pm",@sprintf("%4.2f",feonmg.(cmf_rms)),"\$ \\cr")
println(tab09_02)
