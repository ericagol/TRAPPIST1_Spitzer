

using JLD2; using PyPlot

include("../../../src/CGS.jl")
using Main.CGS

@load "../../../data/T1_JWST_all_transits_mass_error.jld2" JWST_sims nplanet nsim mstar planet

clf()
for ip=1:nplanet
  mp_in = zeros(Float64,nsim); mp_out = zeros(Float64,nsim); sp_out = zeros(Float64,nsim);
  for i=1:nsim
    mp_in[i]  = JWST_sims[i].elements_in[ip+1,1]*mstar*MSUN/MEARTH
    mp_out[i] = JWST_sims[i].elements_best[ip+1,1]*mstar*MSUN/MEARTH
    sp_out[i] = sqrt(JWST_sims[i].cov[(ip-1)*5+1,(ip-1)*5+1])*mstar*MSUN/MEARTH
  end
  errorbar(mp_in,mp_out-mp_in,sp_out,fmt="o",label=string("T1",planet[ip]),linewidth=2)
end
plot([0,1.5],[0,0])
axis([0,1.5,-0.002,0.002])
legend(loc="upper left",fontsize=15,ncol=2)
xlabel(L"$M_{in}/M_\oplus$",fontsize=15)
ylabel(L"$(M_{out}-M_{in})/M_\oplus$",fontsize=15)
tight_layout()
savefig("../Recovered_masses_JWST_5yr_all_transits_NIRSPEC.pdf",bbox_inches="tight")

