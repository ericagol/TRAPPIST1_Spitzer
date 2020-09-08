
# Creates a tex file with 35 columns of state_total
# for use in plotCorner_ttv.py:

using JLD2
using DelimitedFiles

@load "../../../data/T1_hmc_total_02212020.jld2" state_total

writedlm("../../../data/state_total.txt",state_total[1:35,:]')
