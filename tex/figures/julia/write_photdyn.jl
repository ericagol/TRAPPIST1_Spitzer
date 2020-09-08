

# Here we take the data from the photodynamic chains
# and write it to a text file so it can be read in
# Python:

using JLD2

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

using DelimitedFiles
writedlm("../../../data/T1_photdyn_chain_noprior.txt",flatchain')
