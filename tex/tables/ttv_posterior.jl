
using JLD2; using DelimitedFiles; using Printf; using PyPlot
using Statistics

include("../../src/regress.jl")
include("../../src/NbodyGradient/src/ttv.jl")


# Load in the posterior from the transit-timing analysis:
@load "../../data/T1_hmc_total_02212020.jld2"

# Make a larger tmax to forecast to first ~2 years of JWST
# (given the current launch date of October 31, 2021):
tmax = 3000.0

# Boolean that determines if the recompute the posterior transit
# times:
recompute = true

# Set up parameters needed for time-integration:
n = 8; nplanet=7; t0 = 7257.93115525; h = 0.06; rstar =1e12

pname = ["b","c","d","e","f","g","h"]
# Read in orbital elements
elements = readdlm("../../data/elements_noprior_students.txt",',')

# Make an array, tti1,  to hold transit times:
ntt = zeros(Int64,n);
for j=2:n
  ntt[j] = ceil(Int64,tmax/elements[j,2])+1
end

# Run 10^3 realizations of the posterior:
ngrid = 1000
#ngrid = 10
tt1 = zeros(n,maximum(ntt))
tt_grid = zeros(n,maximum(ntt),ngrid)
count1 = zeros(Int64,n)
elements_grid = zeros(n,7,ngrid)


nsamples = size(state_total)[2]
if recompute
# The following code computes the posterior samples.
# However, this takes some time, so these have been saved to a file,
# and can be restored instead.
  for igrid=1:ngrid
    # Select a state from the transit-timing analysis:
    x = state_total[1:35,ceil(Int64,rand()*nsamples)]
    # Insert into the orbital element array:
    for i=2:n, j=1:5
      elements[i,j] = x[(i-2)*5+j]
    end
    @time dq0 = ttv_elements!(n,t0,h,tmax,elements,tt1,count1,0.0,0,0,rstar)
    elements_grid[:,:,igrid] = elements
    tt_grid[:,:,igrid] = tt1
  end
  @save "../../data/T1_timing_posterior.jld2" n nplanet t0 h rstar tt1 count1 ngrid elements_grid tt_grid
else
  # Here is where the timing results are restored:
  @load "../../data/T1_timing_posterior.jld2"
end

# Now, create two tables to save these, one with the measured times
# and posterior, and one with the forecast times:

include("../../src/extract_planet.jl")
data = readdlm("../../data/T1_timings_20191203.txt",',')

# Read in the sources of each of the transit times:
sources = readdlm("../../data/sources.txt",',')

# Now, try to improve the fit to the data:
# First, find the indices that match for each planet:
data_out = zeros(size(data)[1],6)
# Create a timings table:
io = open("timings_table_with_sources.tex","w")
itrans = 1
for ip=1:nplanet
  eobs,tobs,sobs,nobs = extract_planet(data,ip)
  ittv = 1
  for j=1:nobs
    ittv = 1; dt = Inf
    data_out[itrans,1:4] = data[itrans,:]
    # First, match the observation with a transit time:
    for jttv=1:count1[ip+1]
      if abs(tt1[ip+1,jttv]-tobs[j]) < dt
        ittv = jttv
        dt = abs(tt1[ip+1,jttv]-tobs[j])
      end
    end
    # Now compute the mean and standard deviation of the posterior:
    data_out[itrans,5] = mean(tt_grid[ip+1,ittv,:])
    data_out[itrans,6] = std(tt_grid[ip+1,ittv,:])
    dstring = string(pname[ip]," & ",@sprintf("%i",data[itrans,2])," & ",@sprintf("%12.6f",data[itrans,3]), " & ",@sprintf("%9.6f",data[itrans,4])," & ",@sprintf("%12.6f",data_out[itrans,5])," & ",@sprintf("%9.6f",data_out[itrans,6])," &  ",sources[itrans],"  \\cr")
    println(io,dstring)
    println(dstring)
    # Now, create output for timings table:
    global itrans += 1
  end
end
close(io)
writedlm("../../data/times_obs_and_posterior.txt",data_out,',')
