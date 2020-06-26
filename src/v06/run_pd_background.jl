

include("run_pd_hyak.jl")
datafile = "../../data/T1_Spitzer_data.jld2"
foutput = "T1_pd_MCMC_nw50_ns50_limbprior_001.jld2"
numwalkers = 50
burnin = 1
thinning = 1
astep = 2.0
nsteps = 50
run_pd(datafile,foutput,numwalkers,burnin,thinning,astep,nsteps;nout=10000,limbprior=true)
