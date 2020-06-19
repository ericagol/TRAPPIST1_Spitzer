

function show_args(args)
  @show args
end

include("run_pd_hyak.jl")
datafile = "T1_Spitzer_data.jld2"
#foutput = "T1_pd_MCMC_run_06.jld2"
foutput = string("T1_pd_MCMC_run_101_",show_args(ARGS)[1],".jld2")
numwalkers = 50
burnin = 1
thinning = 10
astep = 2.0
nsteps = 40000
run_pd(datafile,foutput,numwalkers,burnin,thinning,astep,nsteps;nout=10000,limbprior=true,incprior=true)
