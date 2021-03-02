
#  phifh_1 = normalize.(2eloft[i1:i2,2,7] .- eloft[i1:i2,2,5] .- eloft[i1:i2,4,5])
#  phifh_1[phifh_1 .> pi] .-= 2pi
#  phifh_2 = normalize.(2eloft[i1:i2,2,7] .- eloft[i1:i2,2,5] .- eloft[i1:i2,4,7])
#  phifh_2[phifh_2 .> pi] .-= 2pi
#  phifh_2[phifh_2 .< -pi/2] .+= 2pi
#  phifh_3 = normalize.(eloft[i1:i2,4,5] .- eloft[i1:i2,4,7])
#  ax = axes[11]
#  ax.plot(t[i1:i2],phifh_1[i1:i2].*rtd,label=L"$\phi_{fh,1}$",color=cp[1],".");
#  ax.plot(t[i1:i2],phifh_2[i1:i2].*rtd,label=L"$\phi_{fh,2}$",color=cp[2],".");
#  ax.plot(t[i1:i2],phifh_3[i1:i2].*rtd,label=L"$\phi_{fh,3}$",color=cp[3],".");
#  ax.legend(ncol=3,fontsize=4);
#  ax.set_xlabel("Time [days]");
#  ax.set_ylabel(L"$\phi$ [$^\circ$]");
#  ax = axes[12]
#  ax.plot(eloft[i1:i2,3,6] .*cos.(eloft[i1:i2,4,6]),eloft[i1:i2,3,6] .*sin.(eloft[i1:i2,4,6]),label="g")
#  ax.plot(eloft[i1:i2,3,7] .*cos.(eloft[i1:i2,4,7]),eloft[i1:i2,3,7] .*sin.(eloft[i1:i2,4,7]),label="h")
#  ax.set_xlabel(L"$e_\mathrm{free} \cos{\omega_\mathrm{free}}$")
#  ax.set_ylabel(L"$e_\mathrm{free} \sin{\omega_\mathrm{free}}$")
#  ax.legend(fontsize=4)
#  ax.plot([0],[0],"o")



using DelimitedFiles
using PyPlot
using Statistics
#include("../../compute_elements.jl")
include("../../../src/NbodyGradient/src/ttv.jl")
include("compute_elements.jl")

# Re-run the maximum likelihood model for 1600 days, and print out the
# state to a file:
elements = readdlm("../../../data/elements_noprior_students.txt",',')
fileout =   "T1_maxlike_xv_1600.txt"

rm(fileout,force=true)
n = 8
nplanet = n-1
t0 = 7257.93115525
h = 0.06
tmax = 1600.0

# Set up array for holding times of transit:
ntt = zeros(Int64,n);
for i=2:n
  ntt[i] = ceil(Int64,tmax/elements[i,2])+1
end
tt1 = zeros(n,maximum(ntt));
# Save a counter for the actual number of transit times of each planet:
count1 = zeros(Int64,n);
# Set the "size" of the star (only checks that close to transit):
rstar = 1e11
@time dq0 = ttv_elements!(n,t0,h,tmax,elements,tt1,count1,0.0,0,0,rstar;fout=fileout,iout=1)

#mass_grid = readdlm("T1_masses_grimm_10k.txt")



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

nskip = 100
rtd = 180/pi
cp = ["C0","C1","C2","C3","C4","C5","C6","C7"]
i1 = 1
i2 = 16668
# T1_xv_grimm_10k_03991.txt              T1_grimm_100_efree_phi_all.jld2         T1_xv_grimm_10k_06564.txt 
  fname = fileout

  t,eloft = compute_elements(fname,elements[:,1],n)
  nt = length(t)

  fig,axes = subplots(4,4)

#  phibc_1 = normalize.(3eloft[i1:i2,2,2] .- 2eloft[i1:i2,2,1] .- eloft[i1:i2,4,1])
  phibc_1 = normalize.(8eloft[i1:i2,2,2] .- 5eloft[i1:i2,2,1] .- eloft[i1:i2,4,1])
  phibc_1[phibc_1 .> pi] .-= 2pi
#  phibc_2 = normalize.(3eloft[i1:i2,2,2] .- 2eloft[i1:i2,2,1] .- eloft[i1:i2,4,2])
  phibc_2 = normalize.(8eloft[i1:i2,2,2] .- 5eloft[i1:i2,2,1] .- eloft[i1:i2,4,2])
  phibc_2[phibc_2 .> pi] .-= 2pi
  phibc_2[phibc_2 .< -pi/2] .+= 2pi
  phibc_3 = normalize.(eloft[i1:i2,4,1] .- eloft[i1:i2,4,2])
  ax = axes[1]
  ax.plot(t[i1:i2],phibc_1[i1:i2].*rtd,label=L"$\phi_{bc,1}$",color=cp[1],".");
  ax.plot(t[i1:i2],phibc_2[i1:i2].*rtd,label=L"$\phi_{bc,2}$",color=cp[2],".");
  ax.plot(t[i1:i2],phibc_3[i1:i2].*rtd,label=L"$\phi_{bc,3}$",color=cp[3],".");
  ax.legend(ncol=3,fontsize=4);
  ax.set_xlabel("Time [days]");
  ax.set_ylabel(L"$\phi$ [$^\circ$]");
  ax = axes[2]
  ax.plot(eloft[i1:i2,3,1] .*cos.(eloft[i1:i2,4,1]),eloft[i1:i2,3,1] .*sin.(eloft[i1:i2,4,1]),label="b")
  ax.plot(eloft[i1:i2,3,2] .*cos.(eloft[i1:i2,4,2]),eloft[i1:i2,3,2] .*sin.(eloft[i1:i2,4,2]),label="c")
  ax.plot([0],[0],"o")
  ax.set_xlabel(L"$e_\mathrm{free} \cos{\omega_\mathrm{free}}$")
  ax.set_ylabel(L"$e_\mathrm{free} \sin{\omega_\mathrm{free}}$")
  ax.legend(fontsize=4)
  println(" Range, 1: ",maximum(phibc_1)-minimum(phibc_1)," Range 2: ",maximum(phibc_1)-minimum(phibc_1))

#  phicd_1 = normalize.(3eloft[i1:i2,2,3] .- 2eloft[i1:i2,2,2] .- eloft[i1:i2,4,2])
  phicd_1 = normalize.(5eloft[i1:i2,2,3] .- 3eloft[i1:i2,2,2] .- 2offset(eloft[i1:i2,4,2],i2-i1+1))
  phicd_1[phicd_1 .> pi] .-= 2pi
#  phicd_2 = normalize.(3eloft[i1:i2,2,3] .- 2eloft[i1:i2,2,2] .- offset(eloft[i1:i2,4,3],i2-i1+1))
  phicd_2 = normalize.(5eloft[i1:i2,2,3] .- 3eloft[i1:i2,2,2] .- 2offset(eloft[i1:i2,4,3],i2-i1+1))
  phicd_2[phicd_2 .> pi] .-= 2pi
  phicd_2[phicd_2 .< -pi/2] .+= 2pi
  phicd_3 = normalize.(eloft[i1:i2,4,2] .- eloft[i1:i2,4,3])
  ax = axes[3]
  ax.plot(t[i1:i2],phicd_1[i1:i2].*rtd,label=L"$\phi_{cd,1}$",color=cp[1],".");
  ax.plot(t[i1:i2],phicd_2[i1:i2].*rtd,label=L"$\phi_{cd,2}$",color=cp[2],".");
  ax.plot(t[i1:i2],phicd_3[i1:i2].*rtd,label=L"$\phi_{cd,3}$",color=cp[3],".");
  ax.legend(ncol=3,fontsize=4);
  ax.set_xlabel("Time [days]");
  ax.set_ylabel(L"$\phi$ [$^\circ$]");
  ax = axes[4]
  ax.plot(eloft[i1:i2,3,2] .*cos.(eloft[i1:i2,4,2]),eloft[i1:i2,3,2] .*sin.(eloft[i1:i2,4,2]),label="c")
  ax.plot(eloft[i1:i2,3,3] .*cos.(eloft[i1:i2,4,3]),eloft[i1:i2,3,3] .*sin.(eloft[i1:i2,4,3]),label="d")
  ax.set_xlabel(L"$e_\mathrm{free} \cos{\omega_\mathrm{free}}$")
  ax.set_ylabel(L"$e_\mathrm{free} \sin{\omega_\mathrm{free}}$")
  ax.legend(fontsize=4)
  ax.plot([0],[0],"o")

#i2 = 83000
i2 = 16668
  phide_1 = normalize.(3eloft[i1:i2,2,4] .- 2eloft[i1:i2,2,3] .- eloft[i1:i2,4,3])
  phide_1[phide_1 .> pi] .-= 2pi
  phide_2 = normalize.(3eloft[i1:i2,2,4] .- 2eloft[i1:i2,2,3] .- eloft[i1:i2,4,4])
  phide_2[phide_2 .> pi] .-= 2pi
  phide_2[phide_2 .< -pi/2] .+= 2pi
  phide_3 = normalize.(eloft[i1:i2,4,3] .- eloft[i1:i2,4,4])
  ax = axes[5]
  ax.plot(t[i1:i2],phide_1[i1:i2].*rtd,label=L"$\phi_{de,1}$",color=cp[1],".");
  ax.plot(t[i1:i2],phide_2[i1:i2].*rtd,label=L"$\phi_{de,2}$",color=cp[2],".");
  ax.plot(t[i1:i2],phide_3[i1:i2].*rtd,label=L"$\phi_{de,3}$",color=cp[3],".");
  ax.legend(ncol=3,fontsize=4);
  ax.set_xlabel("Time [days]");
  ax.set_ylabel(L"$\phi$ [$^\circ$]");
  ax = axes[6]
  ax.plot(eloft[i1:i2,3,3] .*cos.(eloft[i1:i2,4,3]),eloft[i1:i2,3,3] .*sin.(eloft[i1:i2,4,3]),label="d")
  ax.plot(eloft[i1:i2,3,4] .*cos.(eloft[i1:i2,4,4]),eloft[i1:i2,3,4] .*sin.(eloft[i1:i2,4,4]),label="e")
  ax.set_xlabel(L"$e_\mathrm{free} \cos{\omega_\mathrm{free}}$")
  ax.set_ylabel(L"$e_\mathrm{free} \sin{\omega_\mathrm{free}}$")
  ax.legend(fontsize=4)
  ax.plot([0],[0],"o")

  phief_1 = normalize.(3eloft[i1:i2,2,5] .- 2eloft[i1:i2,2,4] .- eloft[i1:i2,4,4])
  phief_1[phief_1 .> pi] .-= 2pi
  phief_2 = normalize.(3eloft[i1:i2,2,5] .- 2eloft[i1:i2,2,4] .- eloft[i1:i2,4,5])
  phief_2[phief_2 .> pi] .-= 2pi
  phief_2[phief_2 .< -pi/2] .+= 2pi
  phief_3 = normalize.(eloft[i1:i2,4,4] .- eloft[i1:i2,4,5])
  ax = axes[7]
  ax.plot(t[i1:i2],phief_1[i1:i2].*rtd,label=L"$\phi_{ef,1}$",color=cp[1],".");
  ax.plot(t[i1:i2],phief_2[i1:i2].*rtd,label=L"$\phi_{ef,2}$",color=cp[2],".");
  ax.plot(t[i1:i2],phief_3[i1:i2].*rtd,label=L"$\phi_{ef,3}$",color=cp[3],".");
  ax.legend(ncol=3,fontsize=4);
  ax.set_xlabel("Time [days]");
  ax.set_ylabel(L"$\phi$ [$^\circ$]");
  ax = axes[8]
  ax.plot(eloft[i1:i2,3,4] .*cos.(eloft[i1:i2,4,4]),eloft[i1:i2,3,4] .*sin.(eloft[i1:i2,4,4]),label="e")
  ax.plot(eloft[i1:i2,3,5] .*cos.(eloft[i1:i2,4,5]),eloft[i1:i2,3,5] .*sin.(eloft[i1:i2,4,5]),label="f")
  ax.set_xlabel(L"$e_\mathrm{free} \cos{\omega_\mathrm{free}}$")
  ax.set_ylabel(L"$e_\mathrm{free} \sin{\omega_\mathrm{free}}$")
  ax.legend(fontsize=4)
  ax.plot([0],[0],"o")

  phifg_1 = normalize.(4eloft[i1:i2,2,6] .- 3eloft[i1:i2,2,5] .- eloft[i1:i2,4,5])
  phifg_1[phifg_1 .> pi] .-= 2pi
  phifg_2 = normalize.(4eloft[i1:i2,2,6] .- 3eloft[i1:i2,2,5] .- eloft[i1:i2,4,6])
  phifg_2[phifg_2 .> pi] .-= 2pi
  phifg_2[phifg_2 .< -pi/2] .+= 2pi
  phifg_3 = normalize.(eloft[i1:i2,4,5] .- eloft[i1:i2,4,6])
  ax = axes[9]
  ax.plot(t[i1:i2],phifg_1[i1:i2].*rtd,label=L"$\phi_{fg,1}$",color=cp[1],".");
  ax.plot(t[i1:i2],phifg_2[i1:i2].*rtd,label=L"$\phi_{fg,2}$",color=cp[2],".");
  ax.plot(t[i1:i2],phifg_3[i1:i2].*rtd,label=L"$\phi_{fg,3}$",color=cp[3],".");
  ax.legend(ncol=3,fontsize=4);
  ax.set_xlabel("Time [days]");
  ax.set_ylabel(L"$\phi$ [$^\circ$]");
  ax = axes[10]
  ax.plot(eloft[i1:i2,3,5] .*cos.(eloft[i1:i2,4,5]),eloft[i1:i2,3,5] .*sin.(eloft[i1:i2,4,5]),label="f")
  ax.plot(eloft[i1:i2,3,6] .*cos.(eloft[i1:i2,4,6]),eloft[i1:i2,3,6] .*sin.(eloft[i1:i2,4,6]),label="g")
  ax.set_xlabel(L"$e_\mathrm{free} \cos{\omega_\mathrm{free}}$")
  ax.set_ylabel(L"$e_\mathrm{free} \sin{\omega_\mathrm{free}}$")
  ax.legend(fontsize=4)
  ax.plot([0],[0],"o")

  phigh_1 = normalize.(3eloft[i1:i2,2,7] .- 2eloft[i1:i2,2,6] .- eloft[i1:i2,4,6])
  phigh_1[phigh_1 .> pi] .-= 2pi
  phigh_2 = normalize.(3eloft[i1:i2,2,7] .- 2eloft[i1:i2,2,6] .- eloft[i1:i2,4,7])
  phigh_2[phigh_2 .> pi] .-= 2pi
  phigh_2[phigh_2 .< -pi/2] .+= 2pi
  phigh_3 = normalize.(eloft[i1:i2,4,6] .- eloft[i1:i2,4,7])
  ax = axes[11]
  ax.plot(t[i1:i2],phigh_1[i1:i2].*rtd,label=L"$\phi_{gh,1}$",color=cp[1],".");
  ax.plot(t[i1:i2],phigh_2[i1:i2].*rtd,label=L"$\phi_{gh,2}$",color=cp[2],".");
  ax.plot(t[i1:i2],phigh_3[i1:i2].*rtd,label=L"$\phi_{gh,3}$",color=cp[3],".");
  ax.legend(ncol=3,fontsize=4);
  ax.set_xlabel("Time [days]");
  ax.set_ylabel(L"$\phi$ [$^\circ$]");
  ax = axes[12]
  ax.plot(eloft[i1:i2,3,6] .*cos.(eloft[i1:i2,4,6]),eloft[i1:i2,3,6] .*sin.(eloft[i1:i2,4,6]),label="g")
  ax.plot(eloft[i1:i2,3,7] .*cos.(eloft[i1:i2,4,7]),eloft[i1:i2,3,7] .*sin.(eloft[i1:i2,4,7]),label="h")
  ax.set_xlabel(L"$e_\mathrm{free} \cos{\omega_\mathrm{free}}$")
  ax.set_ylabel(L"$e_\mathrm{free} \sin{\omega_\mathrm{free}}$")
  ax.legend(fontsize=4)
  ax.plot([0],[0],"o")

  phieg_1 = normalize.(2eloft[i1:i2,2,6] .- eloft[i1:i2,2,4] .- eloft[i1:i2,4,4])
  phieg_1[phieg_1 .> pi] .-= 2pi
  phieg_2 = normalize.(2eloft[i1:i2,2,6] .- eloft[i1:i2,2,4] .- eloft[i1:i2,4,6])
  phieg_2[phieg_2 .< 0.0] .+= 2pi
  phieg_3 = normalize.(eloft[i1:i2,4,4] .- eloft[i1:i2,4,6])
  ax = axes[13]
  ax.plot(t[i1:i2],phieg_1[i1:i2].*rtd,label=L"$\phi_{eg,1}$",color=cp[1],".");
  ax.plot(t[i1:i2],phieg_2[i1:i2].*rtd,label=L"$\phi_{eg,2}$",color=cp[2],".");
  ax.plot(t[i1:i2],phieg_3[i1:i2].*rtd,label=L"$\phi_{eg,3}$",color=cp[3],".");
  ax.legend(ncol=3,fontsize=4);
  ax.set_xlabel("Time [days]");
  ax.set_ylabel(L"$\phi$ [$^\circ$]");
  ax = axes[14]
  ax.plot(eloft[i1:i2,3,4] .*cos.(eloft[i1:i2,4,4]),eloft[i1:i2,3,4] .*sin.(eloft[i1:i2,4,4]),label="e")
  ax.plot(eloft[i1:i2,3,6] .*cos.(eloft[i1:i2,4,6]),eloft[i1:i2,3,6] .*sin.(eloft[i1:i2,4,6]),label="g")
  ax.set_xlabel(L"$e_\mathrm{free} \cos{\omega_\mathrm{free}}$")
  ax.set_ylabel(L"$e_\mathrm{free} \sin{\omega_\mathrm{free}}$")
  ax.legend(fontsize=4)
  ax.plot([0],[0],"o")

  phifh_1 = normalize.(2eloft[i1:i2,2,7] .- eloft[i1:i2,2,5] .- eloft[i1:i2,4,5])
  phifh_1[phifh_1 .> pi] .-= 2pi
  phifh_2 = normalize.(2eloft[i1:i2,2,7] .- eloft[i1:i2,2,5] .- eloft[i1:i2,4,7])
  phifh_2[phifh_2 .> pi] .-= 2pi
  phifh_2[phifh_2 .< -pi/2] .+= 2pi
  phifh_3 = normalize.(eloft[i1:i2,4,5] .- eloft[i1:i2,4,7])
  phifh_3[phifh_3 .> pi] .-= 2pi
  ax = axes[15]
  ax.plot(t[i1:i2],phifh_1[i1:i2].*rtd,label=L"$\phi_{fh,1}$",color=cp[1],".");
  ax.plot(t[i1:i2],phifh_2[i1:i2].*rtd,label=L"$\phi_{fh,2}$",color=cp[2],".");
  ax.plot(t[i1:i2],phifh_3[i1:i2].*rtd,label=L"$\phi_{fh,3}$",color=cp[3],".");
  ax.legend(ncol=3,fontsize=4);
  ax.set_xlabel("Time [days]");
  ax.set_ylabel(L"$\phi$ [$^\circ$]");
  ax = axes[16]
  ax.plot(eloft[i1:i2,3,5] .*cos.(eloft[i1:i2,4,5]),eloft[i1:i2,3,5] .*sin.(eloft[i1:i2,4,5]),label="f")
  ax.plot(eloft[i1:i2,3,7] .*cos.(eloft[i1:i2,4,7]),eloft[i1:i2,3,7] .*sin.(eloft[i1:i2,4,7]),label="h")
  ax.set_xlabel(L"$e_\mathrm{free} \cos{\omega_\mathrm{free}}$")
  ax.set_ylabel(L"$e_\mathrm{free} \sin{\omega_\mathrm{free}}$")
  ax.legend(fontsize=4)
  ax.plot([0],[0],"o")
