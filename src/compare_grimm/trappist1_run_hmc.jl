

function show_args(args)
  @show args
end

# Example of running this routine:

include("run_hmc_background_student_rescale.jl")



# File containing the transit times:
#fname = "../../data/T1_timings_20191203.txt"
fname = "timing_simon.txt"

felements = "../../data/elements_noprior_students.txt"

# Name an output file:
#foutput = string("T1_run_hmc_student_ecc_lgndof_V1exp2nuinv_nstep2000_eps0.1_nleap20_",show_args(ARGS)[1],".jld2")
foutput = string("T1_run_hmc_student_ecc_lgndof_V1exp2nuinv_nstep10_eps0.1_nleap20_",show_args(ARGS)[1],".jld2")

# Call the likelihood profile function:
#run_hmc(fname,felements,foutput,2000,0.1,20,false)
run_hmc(fname,felements,foutput,10,0.1,20,false)
