using PyPlot

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

function plot_ttv(data,elements,tt1)

# Now make some plots:
fig,axes = subplots(3,2,sharex="col",figsize=(12,8))
plabel = ["T1b","T1c","T1d","T1e","T1f","T1g","T1h"]
range = [1.0,1.0,1.0,2.0,10.0,8.0,3.0]
# Now make some plots:
#data = readdlm(fname,',')
for ip=1:nplanet
  if ip == 1
    ax = axes[ip]
  else
    ax = axes[ip-1]
  end
  coeff = zeros(Float64,2)
  # Plot each "Mode"
  eobs,tobs,sobs,nobs = extract_planet(data,ip)
  if nobs > 1
    fn = zeros(2,nobs)
    fn[1,:] .= 1.0
    fn[2,:] .= eobs
    coeff,cov = regress(fn,tobs,sobs)
    println(ip," ",coeff)
  end
  # Now plot TTVFast:
  fn = zeros(Float64,2,count1[ip+1])
  sig = ones(count1[ip+1])
  tti1 = tt1[ip+1,1:count1[ip+1]]
  tt_ref1 = zeros(count1[ip+1])
  for j=1:count1[ip+1]
    fn[1,j] = 1.0
    fn[2,j] = round(Int64,(tti1[j]-elements[ip+1,3])/elements[ip+1,2])
    epoch = round((tti1[j]-coeff[1])/coeff[2])
    tt_ref1[j] = coeff[1]+coeff[2]*epoch
  end
  coeff,cov = regress(fn,tti1,sig)
  tt_ref1 = coeff[1] .+coeff[2]*fn[2,:]
  ttv1 = (tti1 .-tt_ref1)*24*60.
  ax.plot(tti1,ttv1,label=plabel[ip],zorder=-32)
  if nobs > 0
    ttv_obs = (tobs .- coeff[1] .- coeff[2]*eobs)*24*60.
#  ax[:plot](tobs,ttv_obs)
    ax.errorbar(tobs, ttv_obs, sobs*24*60.,fmt=".")
  end
  if ip == 4 || ip == 7
    ax.set_xlabel(L"BJD$_\mathrm{TDB}$-2,450,000 [d]")
  end
  ax.set_ylabel("TTV [min]")
#  ax.axis([7260,8860,-range[ip],range[ip]])
#  if ip != 2 && ip != 3 && ip != 6
#    ax.legend(loc="upper right")
#  else
#    ax.legend(loc="lower right")
    ax.legend()
#  end
end
tight_layout()
return
end
