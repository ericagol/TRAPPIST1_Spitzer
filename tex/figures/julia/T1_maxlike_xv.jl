
# Compute the starting x & v positions:

using DelimitedFiles
using JLD2

path = "/Users/ericagol/Software/TRAPPIST1_Spitzer"
# Relevant code:
include(string(path,"/src/regress.jl"))
include(string(path,"/src/NbodyGradient/src/ttv.jl"))

# Load in parameters and results from the likelihood profile:
#grid_file=string(path,"/data/T1_likelihood_profile_student_planetc_3.0sig.jld2")
@load "/Users/ericagol/Software/TRAPPIST1_Spitzer//data/T1_likelihood_profile_student_planetc_3.0sig.jld2"

# Choose the maximum likelihood elements:
elements = readdlm(string(path,"/data/elements_noprior_students.txt"),',')
x,v = init_nbody(elements,t0,n)

println("x: ",x)
println("v: ",v)

writedlm("T1_maxlike_xv.txt",[x;v])

open("T1_maxlike_tt.txt","w") do io
for i=2:8
  for j=1:count1[i]
    println(io,string(i," ",j," ",tt_grid[1,i,j,16]))
  end
end
end
