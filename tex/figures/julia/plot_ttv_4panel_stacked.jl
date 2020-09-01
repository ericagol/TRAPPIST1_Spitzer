

using PyPlot
include("../../../src/extract_planet.jl")

function plot_ttv(data,elements,tt1,count1,nplanet)

# Now make some plots:
fig,axes = subplots(4,1,sharex="col",figsize=(6,8))
plabel = ["T1b","T1c","T1d","T1e","T1f","T1g","T1h"]
range = [1.0,1.0,1.0,2.0,12.0,8.0,3.0]

# Now make some plots:
eoffset = [43,10,75,9,7,3,22]
ord = [5,5,18,24,18,18,18]
i1 = 1
i2 = 0
for ip=1:nplanet
  if ip == 1 || ip == 2
    ax = axes[1]
  elseif ip == 3
    ax = axes[2]
  elseif ip == 4 || ip == 5
    ax = axes[3]
  else
    ax = axes[4]
  end
  coeff = zeros(Float64,2)
  # Plot each "Mode"
  eobs,tobs,sobs,nobs = extract_planet(data,ip)
  println("ip: ",ip," nobs: ",nobs)
  fn = zeros(2,nobs)
  fn[1,:] .= 1.0
  fn[2,:] .= eobs
  coeff,cov = regress(fn,tobs,sobs)
  println(ip," ",coeff)
  # Now plot TTVFast:
  fn = zeros(Float64,2,count1[ip+1])
  sig = ones(count1[ip+1])
  tti1 = tt1[ip+1,1:count1[ip+1]]
  tt_ref1 = zeros(count1[ip+1])
  epoch = zeros(count1[ip+1])
  for j=1:count1[ip+1]
    fn[1,j] = 1.0
    fn[2,j] = round(Int64,(tti1[j]-elements[ip+1,3])/elements[ip+1,2])
    epoch[j] = round((tti1[j]-coeff[1])/coeff[2])
    tt_ref1[j] = coeff[1]+coeff[2]*epoch[j]
  end
  coeff,cov = regress(fn,tti1,sig)
  println(ip," ",coeff)
  tt_ref1 = coeff[1] .+coeff[2] .*fn[2,:]
  ttv1 = (tti1 .-tt_ref1) .*(24*60)
  ax.plot(tti1,ttv1,label=plabel[ip],linewidth=1.5)
  ttv_obs = (tobs .- coeff[1] .- coeff[2] .*(eobs .+eoffset[ip]))*(24*60)
  ax.errorbar(tobs, ttv_obs, sobs*(24*60),fmt=".",linewidth=1.5,markersize=4.0)
  i2 = i1 + count1[ip+1] - 1
  resid = (tobs .- tt1[ip+1,indx[iplanet .== (ip+1)]])./sobs
  i1 = i2 + 1
  iout = abs.(resid) .> 3.0
  ax.errorbar(tobs[iout], ttv_obs[iout], sobs[iout]*(24*60),fmt=".",linewidth=1.5,markersize=4.0,color="C5")
  if ip == 7
    ax.set_xlabel(L"BJD$_\mathrm{TDB}$-2,450,000")
  end
  ax.grid(linestyle=":")
  ax.set_ylabel("TTV [min]")
  ax.legend()
  # Now, fit a high-order polynomial to TTVs, and overplot it:
  fn = zeros(ord[ip]+1,count1[ip+1])
  fn[1,:] .= 1.0
  for j=1:ord[ip]
    fn[j+1,:] .= epoch.^j
  end
  coeff,cov = regress(fn,tti1,sig)
  mod_poly = ones(count1[ip+1])*coeff[1]
  for j=1:ord[ip]
     mod_poly .+= fn[j+1,:] .*coeff[j+1]
  end
  ax.plot(tti1,(mod_poly .-tt_ref1) .*(24*60),linewidth=1.5,alpha=0.5)
end
#ax = axes[8]
#ax.axis("off")
tight_layout()
subplots_adjust(hspace=0)
savefig("../T1_ttvs_4panel_stacked.pdf",bbox_inches="tight")
return
end

