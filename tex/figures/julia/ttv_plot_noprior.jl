
# Optimized model and plotted TTVs:

using DelimitedFiles
using JLD2

include("../../../src/regress.jl")

@load "../../../data/T1_likelihood_profile_student_planetc_3.0sig.jld2"

data = readdlm("../../../data/T1_timings_20191203.txt",',')

#include("plot_ttv_6panel.jl")
#include("plot_ttv_4panel.jl")
include("plot_ttv_4panel_stacked.jl")
#include("plot_ttv_vertical.jl")
#include("plot_ttv_vertical_highres.jl")
plot_ttv(data,elements_grid[1,:,:,16],tt_grid[1,:,:,16],count1,7)

include("plot_ttv_diff.jl")
plot_ttv_diff(data,elements_grid[1,:,:,16],tt_grid[1,:,:,16],count1,7,20)
