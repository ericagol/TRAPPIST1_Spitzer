
# Prints out a list of outliers for best-fit model:
using JLD2
using PyPlot
using Printf

include("../src/NbodyGradient/src/ttv.jl")

#@load "output_files/T1_run_hmc_student_ecc_lgndof_V1exp2nuinv_nstep2000_eps0.1_nleap20_201.jld2"
@load "T1_timing_posterior.jld2"
@load "T1_likelihood_profile_student_planetc_3.0sig.jld2"


# Print deviations of each transit:
chi_planet = zeros(Float64,n-1)
nobs_planet = zeros(Int64,n-1)
dev = zeros(Float64,ntrans)
elements = elements_grid[1,:,:,jgrid[1]]
dq0 = ttv_elements!(n,t0,h,tmax,elements,tt1,count1,0.0,0,0,rstar)
println("| Transit number | Planet | Model time | Observed time | Obs - mod | (Obs-mod)/sig | (Obs-mod)^2/sig^2 |")
println("| ---| --- | --- | --- | ---- | --- | --- |")
for j=1:ntrans
  dev[j] = (tobs_tot[j]-tt1[iplanet[j],indx[j]])/sobs_tot[j]
  chi_planet[iplanet[j]-1] += dev[j]^2
  nobs_planet[iplanet[j]-1] += 1
  if abs(dev[j]) > 3
    println("|",@sprintf("%4i",j),"| ",@sprintf("%2i",iplanet[j]),"| ",@sprintf("%4i",indx[j]),"| ",
        @sprintf("%12.5f",tt1[iplanet[j],indx[j]]),"| ",@sprintf("%12.5f",tobs_tot[j]),"| ",
        @sprintf("%8.5f",tobs_tot[j]-tt1[iplanet[j],indx[j]]),"| ",
        @sprintf("%8.5f",sobs_tot[j]),"| ",@sprintf("%5.1f",dev[j]),"| ",@sprintf("%5.1f",dev[j]^2),"|")
  end
end
