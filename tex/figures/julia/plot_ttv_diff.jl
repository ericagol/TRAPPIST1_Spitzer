

using PyPlot
using Statistics
include("../../../src/extract_planet.jl")
include("../../../src/laplace_wisdom.jl")

function chop_coeff_inner(alpha,j)
# Computes f_1^(j) from equation (10) in Deck & Agol (2015).
beta = j*(1-alpha^1.5)
# Equation (11) in Deck & Agol (2015):
f1 = 2*beta*laplace_wisdom(1//2,j,1,alpha)+
    j*(3+beta^2)*laplace_wisdom(1//2,j,0,alpha)
if j == 1
  f1 -= alpha*(beta^2+2*beta+3)
end
f1 *= alpha/(beta^2*(1-beta^2))
return f1
end

function chop_coeff_outer(alpha,j)
# Computes f_2^(j) from equation (14) in Deck & Agol (2015).
kappa = j*(1/alpha^1.5-1)
# Equation (15) in Deck & Agol (2015):
f2 = (j*(kappa^2+3)+2*kappa)*laplace_wisdom(1//2,j,0,alpha)+
    2*kappa*laplace_wisdom(1//2,j,1,alpha)
if j == 1
  f2 -= (kappa^2-2*kappa+3)/alpha^2
end
f2 /= (kappa^2*(kappa^2-1))
return f2
end


function plot_ttv_diff(data,elements,tt1,count1,nplanet,jmax)

# Now make some plots:
#fig,axes = subplots(4,2,sharex="col",constrained_layout="True")
fig,axes = subplots(7,1,sharex="col",constrained_layout="True",figsize=(8,8))
#fig,axes = subplots(7,1,sharex="col",figsize=(24,24))
#fig,axes = subplots(7,1,sharex="col")
plabel = ["T1b","T1c","T1d","T1e","T1f","T1g","T1h"]
range = [1.5,1.5,1.5,3.0,12.0,8.0,4.0]

# First, go through the planets and compute their
# ephemerides:
coeff_planet = zeros(2,nplanet)
for ip=1:nplanet
  eobs,tobs,sobs,nobs = extract_planet(data,ip)
  fn = zeros(2,nobs); fn[1,:] .= 1.0; fn[2,:] .= eobs
  coeff,cov = regress(fn,tobs,sobs)
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
  coeff_planet[:,ip] .= coeff
end

# Now make some plots:
eoffset = [43,10,75,9,7,3,22]
ord = [5,5,18,30,18,18,18]
for ip=1:nplanet
  ax = axes[ip]
  coeff = zeros(Float64,2)
  # Plot each "Mode"
  eobs,tobs,sobs,nobs = extract_planet(data,ip)
#  println("ip: ",ip," nobs: ",nobs)
  fn = zeros(2,nobs)
  fn[1,:] .= 1.0
  fn[2,:] .= eobs
  coeff,cov = regress(fn,tobs,sobs)
#  println(ip," ",coeff)
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
#  println(ip," ",coeff)
  tt_ref1 = coeff[1] .+coeff[2] .*fn[2,:]
  ttv1 = (tti1 .-tt_ref1) .*(24*60)
#  ax.plot(tti1,ttv1,label=plabel[ip])
  # Now, fit a high-order polynomial to TTVs, and overplot it:
  fn = zeros(ord[ip]+1,count1[ip+1])
  fn[1,:] .= 1.0
  for j=1:ord[ip]
    fn[j+1,:] .= epoch.^j
  end
  coeff_poly,cov = regress(fn,tti1,sig)
  mod_poly = ones(count1[ip+1])*coeff_poly[1]
  for j=1:ord[ip]
     mod_poly .+= fn[j+1,:] .*coeff_poly[j+1]
  end
#  ax.plot(tti1,(mod_poly .-tt_ref1) .*(24*60))
  ttv_mod = (tti1 .- mod_poly) .*(24*60)
  ax.plot(tti1,ttv_mod,label=plabel[ip])
#  ttv_obs = tobs .- coeff_poly[1] .- coeff_poly[2] .*(eobs .+eoffset[ip])
  ttv_obs = tobs .- coeff_poly[1] .- coeff_poly[2] .*eobs
  for j=2:ord[ip]
#    ttv_obs .-= coeff[j+1] .*(eobs .+eoffset[ip]).^j
    ttv_obs .-= coeff_poly[j+1] .*eobs.^j
  end
  ttv_obs .*= 24*60
  sobs_min = sobs .*(24*60)
  if ip < 5
    igood = sobs_min .< 0.5*(maximum(ttv_mod)-minimum(ttv_mod))
  else
    igood = ones(Bool,size(ttv_obs)[1])
  end
  println("Fraction of points: ",sum(igood)/size(sobs)[1])
  ax.errorbar(tobs[igood], ttv_obs[igood], sobs_min[igood],fmt=".")
  if ip == 7
    ax.set_xlabel(L"BJD$_\mathrm{TDB}$-2,450,000")
  end
  if ip == 4
    ax.set_ylabel("TTV - Polynomial [min]")
  end
#  ax.legend(loc="upper right",fontsize=10)
  ax.legend(fontsize=8)
  ax.axis([7200,8900,-range[ip],range[ip]])
#  ax.axis([7960,8860,-range[ip],range[ip]])
  # Now, overplot chopping due to all adjacent planets:
  dt_chop = zeros(length(tti1))
  if ip > 1
    for jp=ip-1:-1:1
      # Compute the observed mean longitudes at the time of transit:
      lam1 = 2*pi*(tti1 .-coeff_planet[1,jp])./coeff_planet[2,jp]
      lam2 = 2*pi*(tti1 .-coeff_planet[1,ip])./coeff_planet[2,ip]
      dlambda = lam1-lam2
      alpha = (coeff_planet[2,jp]/coeff_planet[2,ip])^(2//3)
      for j=1:jmax
        dt_chop .+= elements[1+jp,1]*chop_coeff_outer(alpha,j)*sin.(j*dlambda)
      end
    end
  end
  if ip < nplanet
    for jp = ip+1:nplanet
      # Compute the observed mean longitudes at the time of transit:
      lam1 = 2*pi*(tti1 .-coeff_planet[1,ip])./coeff_planet[2,ip]
      lam2 = 2*pi*(tti1 .-coeff_planet[1,jp])./coeff_planet[2,jp]
      dlambda = lam1-lam2
      alpha = (coeff_planet[2,ip]/coeff_planet[2,jp])^(2//3)
      for j=1:jmax
        dt_chop .+= elements[1+jp,1]*chop_coeff_inner(alpha,j)*sin.(j*dlambda)
      end
    end
  end
  dt_chop .*= coeff_planet[2,ip]/(2pi)
  println(ip," range of chopping: ",minimum(dt_chop*(24*60))," ",maximum(dt_chop*(24*60)))
  # Remove a polynomial:
  coeff_poly,cov = regress(fn,dt_chop,sig)
  fill!(mod_poly,coeff_poly[1])
  for j=1:ord[ip]
     mod_poly .+= fn[j+1,:] .*coeff_poly[j+1]
  end
  ax.plot(tti1,(dt_chop-mod_poly)*(24*60),".",markersize=2.0)
#  ax.plot(tti1,(dt_chop-mod_poly)*(24*60) .- ttv_mod)
end
#ax = axes[8]
#ax.axis("off")
tight_layout()
subplots_adjust(hspace=0)
#println("Ephemerides: ",coeff_planet)
savefig("../T1_chopping.pdf",bbox_inches="tight")
return
end
