
# Creates a plot of the simulated JWST transit times for TRAPPIST-1:

using JLD2
include("plot_ttv.jl")
include("../../../src/regress.jl")

@load "../../../data/T1_JWST_all_transits_mass_error.jld2" JWST_sims data_obs

nplanet = JWST_sims[1].nplanet
count1 = JWST_sims[1].count
data_obs = JWST_sims[1].data_obs
elements_sim= JWST_sims[1].elements_in
tt_sim = JWST_sims[1].tt_actual
plot_ttv(data_obs,elements_sim,tt_sim)
subplots_adjust(hspace=0)
savefig("../T1_JWST_all_transits_tight.pdf",bbox_inches="tight")
