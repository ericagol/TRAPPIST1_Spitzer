

include("../../../src/histogram_code.jl")

if !@isdefined(CGS)
  include("../../../src/CGS.jl")
  using Main.CGS
end
using PyPlot
using Printf
using JLD2

# Load in the likelihood profile dataset:
@load "../../../data/T1_likelihood_profile_student_all_3.0sig.jld2"

# Load in the Markov chain dataset:
@load "../../../data/T1_hmc_total_02212020.jld2"

using Optim

function normal_data(m0,sm0,mgrid,dchi)
  function fit_normal(x)
    return sum((exp.(-0.5*dchi) .-exp.(-0.5*((mgrid .-x[1])/x[2]).^2)).^2)
  end
  result = Optim.optimize(fit_normal,[m0,sm0])
  return Optim.minimizer(result)
end

;clf()
fig,axes = subplots(3,1,figsize=(8,8))
fac = 0.09*MSUN/MEARTH; dchi = chi_grid_all .-minimum(chi_grid_all)
# Overax.plot Simon's values:
#mgrimm = [1.310,1.307,0.367,0.642,0.950,1.245,0.373]
mgrimm = [1.017,1.156,0.297,0.772,0.934,1.148,0.331]
#sgrimm = [0.082,0.073,0.022,0.027,0.026,0.027,0.034]
sgrimm1 = [0.143,0.131,0.035,0.075,0.078,0.095,0.049]
sgrimm2 = [0.154,0.142,0.039,0.079,0.080,0.098,0.056]
#cplanet = [:blue,:green,cplanet[i],cplanet[i],cplanet[i],:yellow,:black]
cplanet = ["C0","C1","C2","C3","C4","C5","C6"]
planet =["b","c","d","e","f","g","h"]
# Fit Gaussians to each set of parameters:
nplanet = n-1
xfit = zeros(2,nplanet)
for i=1:nplanet
  xfit[:,i] = normal_data(elements_grid_all[5i-4,i+1,1,16]*fac,0.01,elements_grid_all[5i-4,i+1,1,:]*fac,dchi[5i-4,:])
  println("planet: ",i," mass: ",@sprintf("%6.4f",xfit[1,i]),"+-",@sprintf("%6.4f",xfit[2,i]))
end
#i = 5
#xfit[:,i] = normal_data(1.0,sgrimm1[i],elements_grid_all[5i-4,i+1,1,:]*fac,dchi[5i-4,:])

# First, ax.plot planets d, h:
ax = axes[1]
i=3
mgrid = elements_grid_all[5i-4,i+1,1,:].*fac
ax.plot(mgrid,exp.(-.5*dchi[5i-4,:]),label=planet[i],linewidth=1,color=cplanet[i])
#ax.plot(mgrid,exp.(-.5*(mgrid.-xfit[1,i]).^2/xfit[2,i]^2),linestyle=":",color=:black)
ax.plot(mgrid,exp.(-.5*(mgrid.-xfit[1,i]).^2/(cov_save[(i-1)*5+1,(i-1)*5+1]*fac^2)),linestyle=":",color=cplanet[i],linewidth=2,alpha=0.5)
ax.plot(mgrimm[i]*[1,1],[0.4,0.7],color=cplanet[i])
ax.plot((mgrimm[i]-sgrimm1[i])*[1,1],[0.4,0.7],color=cplanet[i],linewidth=2)
ax.plot((mgrimm[i]+sgrimm2[i])*[1,1],[0.4,0.7],color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i] .+[-sgrimm1[i],sgrimm2[i]],[0.55,0.55],color=cplanet[i],linewidth=2)
#ax.errorbar(mgrimm[i],0.5,xerr=[sgrimm1[i],sgrimm2[i]],color=cplanet[i])
# Now, make a histogram:
md_bin,md_hist,md_bin_square,md_hist_square = histogram(state_total[(i-1)*5+1,:].*fac,50)
#ax.plot(md_bin,md_hist ./maximum(md_hist),color=cplanet[i],linewidth=2,linestyle="--")
ax.plot(md_bin_square,md_hist_square ./maximum(md_hist_square),color=cplanet[i],linewidth=2)
i=7
mgrid = elements_grid_all[5i-4,i+1,1,:].*fac
ax.plot(mgrid,exp.(-.5*dchi[5i-4,:]),label=planet[i],linewidth=1,color=cplanet[i])
#ax.plot(mgrid,exp.(-.5*(mgrid .-xfit[1,i]).^2/xfit[2,i]^2),linestyle=":",color=:black)
ax.plot(mgrid,exp.(-.5*(mgrid .-xfit[1,i]).^2/(cov_save[(i-1)*5+1,(i-1)*5+1]*fac^2)),linestyle=":",color=cplanet[i],linewidth=2,alpha=0.5)
ax.plot(mgrimm[i]*[1,1],[0.35,0.65],color=cplanet[i])
ax.plot((mgrimm[i]-sgrimm1[i])*[1,1],[0.35,0.65],color=cplanet[i],linewidth=2)
ax.plot((mgrimm[i]+sgrimm2[i])*[1,1],[0.35,0.65],color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i] .+[-sgrimm1[i],sgrimm2[i]],[0.5,0.5],color=cplanet[i],linewidth=2)
# Now plot HMC posterior:
mh_bin,mh_hist,mh_bin_square,mh_hist_square = histogram(state_total[(i-1)*5+1,:].*fac,50)
#ax.plot(mh_bin,mh_hist./maximum(mh_hist),color=cplanet[i],linewidth=2,linestyle="--")
ax.plot(mh_bin_square,mh_hist_square./maximum(mh_hist_square),color=cplanet[i],linewidth=2)
ax.set_xlabel(L"$(M/M_\oplus)(M_*/0.09 M_\odot)^{-1}$",fontsize=10)
#ax.set_xlabel(L"$\frac{M}{M_\oplus}\frac{0.09 M_\odot}{M_*}$",fontsize=10)
#ax.set_ylabel(L"$\exp(-\frac{1}{2}\Delta \chi^2)$",fontsize=10)
ax.set_ylabel("Probability",fontsize=10)
ax.legend(loc = "upper left")
ax.axis([0.25,0.45,0,1.1])
#ax.plot([0.2,0.5],exp(-0.5)*[1,1],color="k",linestyle=":")
# Next, e, f:
ax = axes[2]
i=4
mgrid = elements_grid_all[5i-4,i+1,1,:].*fac
ax.plot(mgrid,exp.(-.5*dchi[5i-4,:]),label=planet[i],linewidth=1,color=cplanet[i])
#ax.plot(mgrid,exp.(-.5*(mgrid.-xfit[1,i]).^2/xfit[2,i]^2),linestyle=":",color=cplanet[i],linewidth=2)
ax.plot(mgrid,exp.(-.5*(mgrid.-xfit[1,i]).^2/(cov_save[(i-1)*5+1,(i-1)*5+1]*fac^2)),linestyle=":",color=cplanet[i],alpha=0.5,linewidth=2)
ax.plot(mgrimm[i]*[1,1],[0.45,0.75],color=cplanet[i])
ax.plot((mgrimm[i]-sgrimm1[i])*[1,1],[0.45,0.75],color=cplanet[i],linewidth=2)
ax.plot((mgrimm[i]+sgrimm2[i])*[1,1],[0.45,0.75],color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i] .+[-sgrimm1[i],sgrimm2[i]],[0.60,0.60],color=cplanet[i],linewidth=2)
# Now plot HMC posterior:
me_bin,me_hist,me_bin_square,me_hist_square = histogram(state_total[(i-1)*5+1,:].*fac,50)
#ax.plot(me_bin,me_hist./maximum(me_hist),color=cplanet[i],linewidth=2,linestyle="--")
ax.plot(me_bin_square,me_hist_square./maximum(me_hist_square),color=cplanet[i],linewidth=2)
i=5
mgrid = elements_grid_all[5i-4,i+1,1,:].*fac
ax.plot(mgrid,exp.(-.5*dchi[5i-4,:]),label=planet[i],linewidth=1,color=cplanet[i])
#ax.plot(mgrid,exp.(-.5*(mgrid .-xfit[1,i]).^2/xfit[2,i]^2),linestyle=":",color=:black)
ax.plot(mgrid,exp.(-.5*(mgrid .-xfit[1,i]).^2/(cov_save[(i-1)*5+1,(i-1)*5+1]*fac^2)),linestyle=":",color=cplanet[i],alpha=0.5,linewidth=2)
# Now plot HMC posterior:
mf_bin,mf_hist,mf_bin_square,mf_hist_square = histogram(state_total[(i-1)*5+1,:].*fac,50)
ax.plot(mf_bin_square,mf_hist_square./maximum(mf_hist_square),color=cplanet[i],linewidth=2)
#ax.plot(mf_bin,mf_hist./maximum(mf_hist),color=cplanet[i],linewidth=2,linestyle="--")
ax.plot(mgrimm[i]*[1,1],[0.4,0.7],color=cplanet[i])
ax.plot((mgrimm[i]-sgrimm1[i])*[1,1],[0.4,0.7],color=cplanet[i],linewidth=2)
ax.plot((mgrimm[i]+sgrimm2[i])*[1,1],[0.4,0.7],color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i] .+[-sgrimm1[i],sgrimm2[i]],[0.55,0.55],color=cplanet[i],linewidth=2)
ax.set_xlabel(L"$(M/M_\oplus)(M_*/0.09 M_\odot)^{-1}$",fontsize=10)
#ax.set_ylabel(L"$\exp(-\frac{1}{2}\Delta \chi^2)$",fontsize=10)
ax.set_ylabel("Probability",fontsize=10)
ax.legend(loc = "upper left")
ax.axis([0.6,1.1,0,1.1])
#ax.plot([0.5,1.1],exp(-0.5)*[1,1],color="k",linestyle=":")
# Extend errors from b/c/g to this panel:
i=1
ax.plot(mgrimm[i]*[1,1],[0.45,0.75],color=cplanet[i])
ax.plot((mgrimm[i]-sgrimm1[i])*[1,1],[0.45,0.75],color=cplanet[i],linewidth=2)
ax.plot((mgrimm[i]+sgrimm2[i])*[1,1],[0.45,0.75],color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i] .+[-sgrimm1[i],sgrimm2[i]],[0.6,0.6],color=cplanet[i],linewidth=2)
i=2
ax.plot(mgrimm[i]*[1,1],[0.4,0.7],color=cplanet[i])
ax.plot((mgrimm[i]-sgrimm1[i])*[1,1],[0.4,0.7],color=cplanet[i],linewidth=2)
ax.plot((mgrimm[i]+sgrimm2[i])*[1,1],[0.4,0.7],color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i] .+[-sgrimm1[i],sgrimm2[i]],[0.55,0.55],color=cplanet[i],linewidth=2)
i=6
ax.plot(mgrimm[i]*[1,1],[0.35,0.65],color=cplanet[i])
ax.plot((mgrimm[i]-sgrimm1[i])*[1,1],[0.35,0.65],color=cplanet[i],linewidth=2)
ax.plot((mgrimm[i]+sgrimm2[i])*[1,1],[0.35,0.65],color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i] .+[-sgrimm1[i],sgrimm2[i]],[0.5,0.5],color=cplanet[i],linewidth=2)
# Next, b, c, g:
ax = axes[3]
i=1
mgrid = elements_grid_all[5i-4,i+1,1,:].*fac
ax.plot(mgrid,exp.(-.5*dchi[5i-4,:]),label=planet[i],linewidth=1,color=cplanet[i])
#ax.plot(mgrid,exp.(-.5*(mgrid.-xfit[1,i]).^2/xfit[2,i]^2),linestyle=":",color=:black)
ax.plot(mgrid,exp.(-.5*(mgrid.-xfit[1,i]).^2/(cov_save[(i-1)*5+1,(i-1)*5+1]*fac^2)),linestyle=":",color=cplanet[i],alpha=0.5,linewidth=2)
ax.plot(mgrimm[i]*[1,1],[0.45,0.75],color=cplanet[i])
ax.plot((mgrimm[i]-sgrimm1[i])*[1,1],[0.45,0.75],color=cplanet[i],linewidth=2)
ax.plot((mgrimm[i]+sgrimm2[i])*[1,1],[0.45,0.75],color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i] .+[-sgrimm1[i],sgrimm2[i]],[0.6,0.6],color=cplanet[i],linewidth=2)
# Now plot HMC posterior:
mb_bin,mb_hist,mb_bin_square,mb_hist_square = histogram(state_total[(i-1)*5+1,:].*fac,50)
#ax.plot(mb_bin,mb_hist./maximum(mb_hist),color=cplanet[i],linewidth=2,linestyle="--")
ax.plot(mb_bin_square,mb_hist_square./maximum(mb_hist_square),color=cplanet[i],linewidth=2)
i=2
mgrid = elements_grid_all[5i-4,i+1,1,:].*fac
ax.plot(mgrid,exp.(-.5*dchi[5i-4,:]),label=planet[i],linewidth=1,color=cplanet[i])
#ax.plot(mgrid,exp.(-.5*(mgrid.-xfit[1,i]).^2/xfit[2,i]^2),linestyle=":",color=:black)
ax.plot(mgrid,exp.(-.5*(mgrid.-xfit[1,i]).^2/(cov_save[(i-1)*5+1,(i-1)*5+1]*fac^2)),linestyle=":",color=cplanet[i],alpha=0.5,linewidth=2)
# Now plot HMC posterior:
mc_bin,mc_hist,mc_bin_square,mc_hist_square = histogram(state_total[(i-1)*5+1,:].*fac,50)
#ax.plot(mc_bin,mc_hist./maximum(mc_hist),color=cplanet[i],linewidth=2,linestyle="--")
ax.plot(mc_bin_square,mc_hist_square./maximum(mc_hist_square),color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i]*[1,1],[0.4,0.7],color=cplanet[i])
ax.plot((mgrimm[i]-sgrimm1[i])*[1,1],[0.4,0.7],color=cplanet[i],linewidth=2)
ax.plot((mgrimm[i]+sgrimm2[i])*[1,1],[0.4,0.7],color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i] .+[-sgrimm1[i],sgrimm2[i]],[0.55,0.55],color=cplanet[i],linewidth=2)
i=6
mgrid = elements_grid_all[5i-4,i+1,1,:].*fac
ax.plot(mgrid,exp.(-.5*dchi[5i-4,:]),label=planet[i],linewidth=1,color=cplanet[i])
#ax.plot(mgrid,exp.(-.5*(mgrid.-xfit[1,i]).^2/xfit[2,i]^2),linestyle=":",color=:black)
ax.plot(mgrid,exp.(-.5*(mgrid.-xfit[1,i]).^2/(cov_save[(i-1)*5+1,(i-1)*5+1]*fac^2)),linestyle=":",color=cplanet[i],alpha=0.5,linewidth=2)
# Now plot HMC posterior:
mg_bin,mg_hist,mg_bin_square,mg_hist_square = histogram(state_total[(i-1)*5+1,:].*fac,50)
#ax.plot(mg_bin,mg_hist./maximum(mg_hist),color=cplanet[i],linewidth=2,linestyle="--")
ax.plot(mg_bin_square,mg_hist_square./maximum(mg_hist_square),color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i]*[1,1],[0.35,0.65],color=cplanet[i])
ax.plot((mgrimm[i]-sgrimm1[i])*[1,1],[0.35,0.65],color=cplanet[i],linewidth=2)
ax.plot((mgrimm[i]+sgrimm2[i])*[1,1],[0.35,0.65],color=cplanet[i],linewidth=2)
ax.plot(mgrimm[i] .+[-sgrimm1[i],sgrimm2[i]],[0.5,0.5],color=cplanet[i],linewidth=2)
ax.set_xlabel(L"$(M/M_\oplus)(M_*/0.09 M_\odot)^{-1}$",fontsize=10)
#ax.set_ylabel(L"$\exp(-\frac{1}{2}\Delta \chi^2)$",fontsize=10)
ax.set_ylabel("Probability",fontsize=10)
ax.legend()
ax.axis([1.1,1.6,0,1.1])
#ax.plot([1.1,1.6],exp(-0.5)*[1,1],color="k",linestyle=":")

tight_layout()
savefig("../T1_masses_03312020.pdf", bbox_inches="tight")
