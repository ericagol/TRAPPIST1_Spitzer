

# Creates a table of forecast times based on 1000 posterior draws
# run for 2600 days created in ttv_posterior.jl

using JLD2; using DelimitedFiles; using Printf; using Statistics

@load "../../data/T1_timing_posterior.jld2" n nplanet t0 h rstar tt1 count1 ngrid elements_grid tt_grid


# Now, try to improve the fit to the data:
# First, find the indices that match for each planet:
data_out = zeros(sum(count1[2:end]),4)
# Create a timings table:
io = open("timings_forecast.tex","w")
pname = ["b","c","d","e","f","g","h"]
itrans = 1
for ip=1:nplanet
  ittv = 1
  for j=1:count1[ip+1]
    # First, match the observation with a transit time:
    # Now compute the mean and standard deviation of the posterior:
    data_out[itrans,1] = ip
    data_out[itrans,2] = j
    data_out[itrans,3] = mean(tt_grid[ip+1,j,:])
    data_out[itrans,4] = std(tt_grid[ip+1,j,:])
    dstring = string(pname[ip]," & ",@sprintf("%i",j)," & ",@sprintf("%12.6f",data_out[itrans,3])," & ",@sprintf("%9.6f",data_out[itrans,4]),"\\cr")
    println(io,dstring)
    # Now, create output for timings table:
    global itrans += 1
  end
end
close(io)
writedlm("times_forecast.txt",data_out,',')
