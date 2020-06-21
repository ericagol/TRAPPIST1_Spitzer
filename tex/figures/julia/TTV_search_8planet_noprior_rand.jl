# Examining results from the search for the eighth planet
# without the Earth-like core-mass-fraction constraint, with
# random initial conditions.

using JLD2; using PyPlot; using DelimitedFiles

elements_7planet = readdlm("../../../data/elements_noprior_students.txt",',')

fnames = readdlm("../../../data/v13/file_list.txt");

size(fnames)

# Minimum chi-square with 7 planets:
chi7 = 1396.205; chimin = 1396.205
chi_all = Float64[]; p8_all = Float64[]; m8_all = Float64[]; t08_all = Float64[]; ecc8_all = Float64[]
m8_nonzero = Float64[]; chi_nonzero = Float64[]
mb_all = Float64[]; mc_all = Float64[]; md_all = Float64[]; me_all = Float64[]; mf_all = Float64[]; mg_all = Float64[]; mh_all = Float64[]
count = 0
elements_all = zeros(9,7,1651)
iall = 1
for fname in fnames
    @load string("../../../data/v13/",fname) chi_grid elements_grid
    for j=1:100
        if chi_grid[j] != 0.0 && chi_grid[j] < chi7 && elements_grid[9,1,j] > 0
            push!(m8_nonzero,elements_grid[9,1,j])
            push!(chi_nonzero,chi_grid[j])
            #println(fname," ",j," ",elements_grid[2:9,1:5,j]," ",chi_grid[j])
        end
        if chi_grid[j] != 0.0
            global count +=1
        end
        ecc8 = sqrt(elements_grid[9,4,j]^2+elements_grid[9,5,j]^2)
        if chi_grid[j] != 0.0 && chi_grid[j] < chi7 && elements_grid[9,1,j] > 0.0 # && ecc8 < 0.01
          #println(fname," ",j," ",chi_grid[j]," ",chi7-chi_grid[j]," ",elements_grid[2:9,1:5,j])
          if chi_grid[j] < chimin
                global chimin = chi_grid[j]
                println(fname," j: ",j," chi_min: ",chimin," ",elements_grid[2:9,1:5,j])
          end 
          push!(chi_all,chi_grid[j])
          push!(p8_all,elements_grid[9,2,j])          
          push!(t08_all,elements_grid[9,3,j])
          #push!(ecc8_all,sqrt(elements_grid[9,4,j]^2+elements_grid[9,5,j]^2))
          push!(ecc8_all,ecc8)
          push!(m8_all,elements_grid[9,1,j]) 
          push!(mb_all,elements_grid[2,1,j]);  push!(mc_all,elements_grid[3,1,j]);  
          push!(md_all,elements_grid[4,1,j]);  push!(me_all,elements_grid[5,1,j])
          push!(mf_all,elements_grid[6,1,j]);  push!(mg_all,elements_grid[7,1,j])
          push!(mh_all,elements_grid[8,1,j])
          elements_all[:,:,iall] .= elements_grid[:,:,j]
          global iall += 1
        end
    end
end
nall = iall - 1

# Okay, so the results show that the minimum chi-square
# is 1362.7, which is about ~33 less than a fit with
# seven planets with CMF-earth constraint.  So, not a
# huge improvement, and this is close to a 2:1 resonance,
# but not near a Laplace resonance.
# which breaks the pattern seen with the inner 7.  Also, the 
# eccentricity is pretty large (~0.08), and the mass is fairly
# small (1.24x10^{-5} ~ 0.37 M_earth.  BIC expected is ~5 log(447) 
# ~ 30.5, while AIC is ~10, so either way, not a great improvement.
# It's interesting that the best-fit values are similar to the
# CMF-earth-like constraint case.  But, it would be more
# convincing if the delta-chi-square were significant.  Either
# way we can use this for a photometric search for transits,
# as well as an estimate of the systematic errors.

pygui(false)

if ~@isdefined(CGS)
    include("/Users/ericagol/Computer/Julia/CGS.jl")
    using Main.CGS
end

fac_mass = 0.09*MSUN/MEARTH

fig,axes = subplots(1,3,sharey="row",constrained_layout="True",figsize=(6,4))
ax = axes[1]
ax.plot(p8_all,chi7 .-chi_all,".",alpha=0.5)
ax.plot(p8_all[ecc8_all .< 0.01],chi7 .-chi_all[ecc8_all .< 0.01],"o")

ax.set_xlabel(L"$P_i$ [days]")
ax.plot([20.0,45.0],[30.5,30.5])
ax.axis([20,45,0,40])
ax.text([30],[29],L"$\Delta$ BIC = 0")
ax.set_ylabel(L"$\Delta 2 \log\mathcal{L}$")

ax = axes[2]
ax.plot(ecc8_all,chi7 .-chi_all,".",alpha=0.5)
ax.plot(ecc8_all[ecc8_all .< 0.01],chi7 .-chi_all[ecc8_all .< 0.01],"o")

ax.set_xlabel(L"$e_i$")
ax.plot([0.0,0.35],[30.5,30.5])
ax.axis([0,0.35,0,40])
ax.text([0.1],[29],L"$\Delta$ BIC = 0")
#ax.set_ylabel(L"$\Delta 2\log\mathcal{L}$")


ax = axes[3]
ax.plot(m8_all*fac_mass,chi7 .-chi_all,".",alpha=0.5)
ax.plot(m8_all[ecc8_all .< 0.01]*fac_mass,chi7 .-chi_all[ecc8_all .< 0.01],"o")

ax.plot([0.0,1.0],[30.5,30.5])
ax.text([0.5],[28],L"$\Delta$ BIC = 0")
ax.axis([0.0,1.0,0,40])
ax.set_xlabel(L"$M_i/M_\oplus$")
#ax.set_ylabel(L"$\Delta 2\log\mathcal{L}$")
subplots_adjust(wspace=0)
#tight_layout()
savefig("../Planet_i_properties.pdf",bbox_inches="tight")





#savefig("Planet_i_mass.pdf",bbox_inches="tight")

using Statistics
(median(chi_all)-minimum(chi_all))/std(chi_all)

using Statistics; using Printf
m7 =[1.3774,1.3109,0.3885, 0.6933, 1.0410, 1.3237, 0.3264]
m7_opt = elements_7planet[2:8,1]*fac_mass
s7 =[0.0593,0.0450, 0.0074, 0.0128, 0.0155, 0.0171, 0.0186]


println("b ",@sprintf("%6.3f",median(mb_all)*fac_mass-m7[1]),"+-",@sprintf("%5.3f",s7[1]),"+-",@sprintf("%5.3f",std(mb_all)*fac_mass))
println("c ",@sprintf("%6.3f",median(mc_all)*fac_mass-m7[2]),"+-",@sprintf("%5.3f",s7[2]),"+-",@sprintf("%5.3f",std(mc_all)*fac_mass))
println("d ",@sprintf("%6.3f",median(md_all)*fac_mass-m7[3]),"+-",@sprintf("%5.3f",s7[3]),"+-",@sprintf("%5.3f",std(md_all)*fac_mass))
println("e ",@sprintf("%6.3f",median(me_all)*fac_mass-m7[4]),"+-",@sprintf("%5.3f",s7[4]),"+-",@sprintf("%5.3f",std(me_all)*fac_mass))
println("f ",@sprintf("%6.3f",median(mf_all)*fac_mass-m7[5]),"+-",@sprintf("%5.3f",s7[5]),"+-",@sprintf("%5.3f",std(mf_all)*fac_mass))
println("g ",@sprintf("%6.3f",median(mg_all)*fac_mass-m7[6]),"+-",@sprintf("%5.3f",s7[6]),"+-",@sprintf("%5.3f",std(mg_all)*fac_mass))
println("h ",@sprintf("%6.3f",median(mh_all)*fac_mass-m7[7]),"+-",@sprintf("%5.3f",s7[7]),"+-",@sprintf("%5.3f",std(mh_all)*fac_mass))

fig,axes = subplots(2,4,sharey="row",figsize=(10,10))
ax = axes[1]
ax.plot(mb_all.*fac_mass,chi7 .-chi_all,".",label="b")
ax.plot(mb_all[ecc8_all .< 0.01].*fac_mass,chi7 .-chi_all[ecc8_all .< 0.01],".",label="low ecc")
ax.fill_between([m7[1]-s7[1],m7[1]+s7[1]],[0,0],[40,40],alpha=0.2,color="C1")
#ax.fill_between([m7_opt[1]-s7[1],m7_opt[1]+s7[1]],[0,0],[40,40],alpha=0.2,color="C2")
ax.set_xlabel(L"$M_b/M_\oplus$")
ax.set_ylabel(L"$\Delta 2\log \mathcal{L}$")
ax.legend()

ax = axes[2]
ax.plot(mc_all.*fac_mass,chi7 .-chi_all,".",label="c")
ax.plot(mc_all[ecc8_all .< 0.01].*fac_mass,chi7 .-chi_all[ecc8_all .< 0.01],".",label="low ecc")
ax.fill_between([m7[2]-s7[2],m7[2]+s7[2]],[0,0],[40,40],alpha=0.2,color="C1")
#ax.fill_between([m7_opt[2]-s7[2],m7_opt[2]+s7[2]],[0,0],[40,40],alpha=0.2,color="C2")

ax.legend()
ax.set_xlabel(L"$M_c/M_\oplus$")
ax.set_ylabel(L"$\Delta 2\log \mathcal{L}$")

ax = axes[3]
ax.plot(md_all.*fac_mass,chi7 .-chi_all,".",label="d")
ax.plot(md_all[ecc8_all .< 0.01].*fac_mass,chi7 .-chi_all[ecc8_all .< 0.01],".",label="low ecc")
ax.fill_between([m7[3]-s7[3],m7[3]+s7[3]],[0,0],[40,40],alpha=0.2,color="C1")
#ax.fill_between([m7_opt[3]-s7[3],m7_opt[3]+s7[3]],[0,0],[40,40],alpha=0.2,color="C2")

ax.legend()
ax.set_xlabel(L"$M_d/M_\oplus$")
#ax.set_ylabel(L"$\Delta 2\log \mathcal{L}$")

ax = axes[4]
ax.plot(me_all.*fac_mass,chi7 .-chi_all,".",label="e")
ax.plot(me_all[ecc8_all .< 0.01].*fac_mass,chi7 .-chi_all[ecc8_all .< 0.01],".",label="low ecc")
ax.fill_between([m7[4]-s7[4],m7[4]+s7[4]],[0,0],[40,40],alpha=0.2,color="C1")
#ax.fill_between([m7_opt[4]-s7[4],m7_opt[4]+s7[4]],[0,0],[40,40],alpha=0.2,color="C2")

ax.legend()
ax.set_xlabel(L"$M_e/M_\oplus$")
#ax.set_ylabel(L"$\Delta 2\log \mathcal{L}$")

ax = axes[5]
ax.plot(mf_all.*fac_mass,chi7 .-chi_all,".",label="f")
ax.plot(mf_all[ecc8_all .< 0.01].*fac_mass,chi7 .-chi_all[ecc8_all .< 0.01],".",label="low ecc")
ax.fill_between([m7[5]-s7[5],m7[5]+s7[5]],[0,0],[40,40],alpha=0.2,color="C1")
#ax.fill_between([m7_opt[5]-s7[5],m7_opt[5]+s7[5]],[0,0],[40,40],alpha=0.2,color="C2")

ax.legend()
ax.set_xlabel(L"$M_f/M_\oplus$")
#ax.set_ylabel(L"$\Delta 2\log \mathcal{L}$")

ax = axes[6]
ax.plot(mg_all.*fac_mass,chi7 .-chi_all,".",label="g")
ax.plot(mg_all[ecc8_all .< 0.01].*fac_mass,chi7 .-chi_all[ecc8_all .< 0.01],".",label="low ecc")
ax.fill_between([m7[6]-s7[6],m7[6]+s7[6]],[0,0],[40,40],alpha=0.2,color="C1")
#ax.fill_between([m7_opt[6]-s7[6],m7_opt[6]+s7[6]],[0,0],[40,40],alpha=0.2,color="C2")

ax.legend()
ax.set_xlabel(L"$M_g/M_\oplus$")
#ax.set_ylabel(L"$\Delta 2\log \mathcal{L}$")

ax = axes[7]
ax.plot(mh_all.*fac_mass,chi7 .-chi_all,".",label="h")
ax.plot(mh_all[ecc8_all .< 0.01].*fac_mass,chi7 .-chi_all[ecc8_all .< 0.01],".",label="low ecc")
ax.fill_between([m7[7]-s7[7],m7[7]+s7[7]],[0,0],[40,40],alpha=0.2,color="C1")
#ax.fill_between([m7_opt[7]-s7[7],m7_opt[7]+s7[7]],[0,0],[40,40],alpha=0.2,color="C2")

ax.legend()
ax.set_xlabel(L"$M_h/M_\oplus$")
#ax.set_ylabel(L"$\Delta 2\log \mathcal{L}$")
ax = axes[8]; ax.axis("off"); tight_layout()
subplots_adjust(wspace=0)
savefig("T1_8planet_mass_comarison.png",bbox_inches="tight")

m7 .- m7_opt

plot(mc_all.*fac_mass,chi_all,"."); elements_7planet[3,1]*fac_mass

plot(mg_all.*fac_mass,chi_all,"."); elements_7planet[7,1]*fac_mass

# Wow, this case has a *very* small mass for planet i,
# yet improves on the chi-square from our 7-planet fit.
# That's not good...  but it's a very high eccentricity.
@load "../../../data/T1_CMFearth_run094.jld2" chi_grid elements_grid

# Note to future self:  this cell below (between the multi-line commens)
# took about 5 hours to run.
# The results are saved in T1_8planet_noprior_efree.jld2, but I should 
# have saved specific variables. To get back the interesting quantities 
# in the future without having to rerun these simulations (which took about 5 hours),
# just carry out:

@load "../../../data/T1_8planet_noprior_efree.jld2" ecosfree esinfree

# The rest of the necessary quantities were computed in the cells above,
# which are quick to run.

# Now, let's compute the free eccentricities of these planets
# and see if we can find a relation to the initial eccentricities:

#=
include("/Users/ericagol/Observing/Spitzer/DDT2019/Campaign04/v09/compute_elements.jl")
include("/Users/ericagol/Software/NbodyGradient/src/ttv.jl")

using Statistics

function normalize(phi)
while phi < 0
  phi += 2pi
end
while phi > 2pi
  phi -= 2pi
end
return phi
end

function offset(phi,nt)
for i=1:nt-1
  while phi[i+1] > (phi[i]+3.0)
    phi[i+1] -= pi
  end
  while phi[i+1] < (phi[i]-3.0)
    phi[i+1] += pi
  end
end
return phi
end

@load "../../../data/T1_CMFearth_run001.jld2"
const NDIM  = 3
const third = 1.0/3.0
const alpha0 = 0.0


ecosfree = zeros(n-1,nall)
esinfree = zeros(n-1,nall)
tt1 = zeros(n,1160)
rstar = 1e12
tmax = 1000.0
i1 = 1; i2 = 16668
elements  = zeros(n,7)

# I'm going to compute the free eccentricities in the
# same manner as I did for the v09/Hyak/stability_10k/ case
# /Users/ericagol/Observing/Spitzer/DDT2019/Campaign04/v09/Hyak/stability_10k/compute_free_eccentricity_10k.jl
# so that I can see if I can use the same transformation in
# the eight planet case as I did in the seven planet case:

for iall=1:nall
  # Insert into the orbital element array:
  elements .= elements_all[:,:,iall]   
  fileout = string("T1_xv_8planet_",lpad(iall,5,"0"),".txt")
  rm(fileout,force=true)
  @time dq0 = ttv_elements!(n,t0,h,tmax,elements,tt1,count1,0.0,0,0,rstar;fout=fileout,iout=1)
  # Now, compute the time-dependent elements:
  t,eloft = compute_elements(fileout,elements[:,1],n)
  # Now, measure free eccentricities of the planets:
  for ip=1:n-1
    ecosfree[ip,iall] = mean(eloft[i1:i2,3,ip].*cos.(eloft[i1:i2,4,ip]))
    esinfree[ip,iall] = mean(eloft[i1:i2,3,ip].*sin.(eloft[i1:i2,4,ip]))
  end
  # Need to remove the file since these are large:
  rm(fileout,force=true)
end

# Now save the results:

#@save "T1_8planet_noprior_efree.jld2"

=#



ecosfree[:,1]

esinfree[:,1]

#@load "/Users/ericagol/Observing/Spitzer/DDT2019/Campaign04/v09/Hyak/stability_10k/efree_osc_coeff.jld2"
@load "../../../data/efree_osc_coeff.jld2"

# Convert free eccentricities to elements, and compare:
ecos_osc = zeros(7,nall); esin_osc = zeros(7,nall)
for iall = 1:nall, ip=1:7
  ecos_osc[ip,iall] = c_ecos_f2o[ip,1] + c_ecos_f2o[ip,2]*ecosfree[ip,iall] + c_ecos_f2o[ip,3]*esinfree[ip,iall]
  esin_osc[ip,iall] = c_esin_f2o[ip,1] + c_esin_f2o[ip,2]*ecosfree[ip,iall] + c_esin_f2o[ip,3]*esinfree[ip,iall]
  #println(iall," ",ip," ",elements_all[ip+1,4:5,iall]," ",ecos_osc[ip,iall]," ",esin_osc[ip,iall]," ",
  #      elements_all[ip+1,4,iall]-ecos_osc[ip,iall]," ",elements_all[ip+1,5,iall]-esin_osc[ip,iall])
end

using PyPlot

# Make plots of the residuals from the prediction for the osculating elements from
# the free elements (based on the fit to seven planets from 
# v09/Hyak/stability_10k/plot_eccentricity_correlations.jl
# efree_osc_coeff.jld2 in that same directory.  Looks pretty good
# - predicted values agree to ~10^{-4} rms.  Worst agreement is for h,
# not surprisingly.

pname = ["b","c","d","e","f","g","h"]
fig,axes = subplots(1,2)
ax = axes[1]
for ip=1:7
  ax.plot(elements_all[ip+1,4,:],ecos_osc[ip,:] .-elements_all[ip+1,4,:],".",label=pname[ip])
  println(ip," ecos: ",std(ecos_osc[ip,:] .-elements_all[ip+1,4,:]))
end
ax.legend()
ax.axis([-0.015,0.015,-0.015,0.015])
ax = axes[2]
for ip=1:7
  ax.plot(elements_all[ip+1,5,:],esin_osc[ip,:] .-elements_all[ip+1,5,:],".",label=pname[ip])
  println(ip," esin: ",std(esin_osc[ip,:] .-elements_all[ip+1,5,:]))
end
ax.legend()
ax.axis([-0.015,0.015,-0.015,0.015])

# Plot chi-square of 8-planet fit vs. residuals of forecast
# osculating eccentricity vectors (based on free eccentricity
# vector).  Some of the lower values of chi-square have the
# significant outliers.

for ip=7:7
  plot(ecos_osc[ip,:] .-elements_all[ip+1,4,:],chi_all,".")
  plot(esin_osc[ip,:] .-elements_all[ip+1,5,:],chi_all,".",label=pname[ip])
end
legend()

# Plot residuals of the prediction vs. orbital period of outer planet.
# Looks most significant near ~39 days (~2:1 period ratio, or p=q=1):
for ip=7:7
  plot(ecos_osc[ip,:] .-elements_all[ip+1,4,:],elements_all[9,2,:]./elements_all[8,2,:],".")
  plot(esin_osc[ip,:] .-elements_all[ip+1,5,:],elements_all[9,2,:]./elements_all[8,2,:],".",label=pname[ip])
end
legend()

# Plot residuals of the prediction vs. eccentricity of outer planet.
# Looks most significant near small eccentricity:
for ip=7:7
  plot(ecos_osc[ip,:] .-elements_all[ip+1,4,:],elements_all[9,4,:],".")
  plot(esin_osc[ip,:] .-elements_all[ip+1,5,:],elements_all[9,5,:],".",label=pname[ip])
end
legend()


