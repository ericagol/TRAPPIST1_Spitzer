

# From data (which is an array with planet number,
# epoch, time of observation, error), return
# the epochs, times, errors, and number of observations
# for the ith planet:
function extract_planet(data,i)
eobs = zeros(0)
tobs = zeros(0)
sobs = zeros(0)
nobs = 0
for j=1:size(data,1)
  if data[j,1] == i
# Extract data:
    if data[j,4] != 0.0
      push!(eobs,data[j,2])
      push!(tobs,data[j,3])
      push!(sobs,data[j,4])
      nobs += 1
    end
  end
end
return eobs,tobs,sobs,nobs
end
