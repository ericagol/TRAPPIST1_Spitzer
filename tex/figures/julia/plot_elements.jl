

using DelimitedFiles
using PyPlot
using Statistics
include("compute_elements.jl")

n = 8
mass_grid = readdlm("../../../data/T1_masses_libration_min_max.txt")
fnames = ["T1_xv_libration_min_1.txt","T1_xv_libration_max_2.txt","T1_positions_velocities0001.txt"]

ifile = 1
fname = string("../../../data/",fnames[ifile])
t,eloft1 = compute_elements(fname,[1.0;mass_grid[ifile,:]],n)
ifile = 2
fname = string("../../../data/",fnames[ifile])
t,eloft2 = compute_elements(fname,[1.0;mass_grid[ifile,:]],n)
ifile = 3
fname = string("../../../data/",fnames[ifile])
t,eloft3 = compute_elements(fname,[1.0;mass_grid[ifile,:]],n)
#mass_grid = [ 4.5651630946375855e-5, 4.509775588066495e-5, 1.2793574298379903e-5, 2.2636788066718058e-5, 3.4259988884100445e-5, 4.390338504975531e-5, 9.898177265368541e-6]
#t,eloft = compute_elements("Hyak/Fabrycky/T1_positions_velocities0001.txt",[1.0;mass_grid],n)

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
  if phi[i+1] > (phi[i]+3.0)
    phi[i+1] -= pi
  end
end
return phi
end

nskip = 100
rtd = 180/pi
cp = ["C0","C1","C2","C3","C4","C5","C6","C7"]
nt = length(t)
phi11 = normalize.(2eloft1[:,2,1] .- 5eloft1[:,2,2] + 3eloft1[:,2,3])
phi11 = offset(phi11,nt)
phi12 = normalize.(2eloft2[:,2,1] .- 5eloft2[:,2,2] + 3eloft2[:,2,3])
phi12 = offset(phi12,nt)
phi13 = normalize.(2eloft3[:,2,1] .- 5eloft3[:,2,2] + 3eloft3[:,2,3])
phi13 = offset(phi13,nt)

clf();
figure(figsize=(8,8))
#plot(t[1:nskip:end],360 .- phi1[1:nskip:end].*rtd,label=L"360-$\phi_{bcd}$",color=cp[1],alpha=0.5);
plot(t[1:nskip:end], phi11[1:nskip:end].*rtd,label=L"$\phi_{bcd}$",color=cp[1],alpha=0.5);
plot(t[1:nskip:end], phi12[1:nskip:end].*rtd,color=cp[1],alpha=0.5);
plot(t[1:nskip:end], phi13[1:nskip:end].*rtd,color=cp[1],alpha=0.5);
plot([t[1],t[end]],[0.0,0.0] .+ (360 - 196.5),color=cp[1],linestyle="--");
phi21 = normalize.(eloft1[:,2,2] .- 3eloft1[:,2,3] + 2eloft1[:,2,4] .+8pi);
phi21 = offset(phi21,nt);
phi22 = normalize.(eloft2[:,2,2] .- 3eloft2[:,2,3] + 2eloft2[:,2,4] .+8pi);
phi22 = offset(phi22,nt);
phi23 = normalize.(eloft3[:,2,2] .- 3eloft3[:,2,3] + 2eloft3[:,2,4] .+8pi);
phi23 = offset(phi23,nt);
#plot(t[1:nskip:end],360 .- phi2[1:nskip:end].*rtd,label=L"360-$\phi_{cde}$",color=cp[2],alpha=0.5);
plot(t[1:nskip:end],phi21[1:nskip:end].*rtd,label=L"$\phi_{cde}$",color=cp[2],alpha=0.5);
plot(t[1:nskip:end],phi22[1:nskip:end].*rtd,color=cp[2],alpha=0.5);
plot(t[1:nskip:end],phi23[1:nskip:end].*rtd,color=cp[2],alpha=0.5);
plot([t[1],t[end]],[0.0,0.0] .+ (360-308.8),color=cp[2],linestyle="--");
phi31 = normalize.(2eloft1[:,2,3] .- 5eloft1[:,2,4] + 3eloft1[:,2,5] .-4pi);
phi31 = offset(phi31,nt);
phi32 = normalize.(2eloft2[:,2,3] .- 5eloft2[:,2,4] + 3eloft2[:,2,5] .-4pi);
phi32 = offset(phi32,nt);
phi33 = normalize.(2eloft3[:,2,3] .- 5eloft3[:,2,4] + 3eloft3[:,2,5] .-4pi);
phi33 = offset(phi33,nt);
plot(t[1:nskip:end],phi31[1:nskip:end].*rtd,label=L"$\phi_{def}$",color=cp[3],alpha=0.5);
plot(t[1:nskip:end],phi32[1:nskip:end].*rtd,color=cp[3],alpha=0.5);
plot(t[1:nskip:end],phi33[1:nskip:end].*rtd,color=cp[3],alpha=0.5);
plot([t[1],t[end]],[0.0,0.0] .+ 219.7,color=cp[3],linestyle="--");
phi41 = normalize.(eloft1[:,2,4] .- 3eloft1[:,2,5] + 2eloft1[:,2,6] .+2pi);
phi41 = offset(phi41,nt);
phi42 = normalize.(eloft2[:,2,4] .- 3eloft2[:,2,5] + 2eloft2[:,2,6] .+2pi);
phi42 = offset(phi42,nt);
phi43 = normalize.(eloft3[:,2,4] .- 3eloft3[:,2,5] + 2eloft3[:,2,6] .+2pi);
phi43 = offset(phi43,nt);
plot(t[1:nskip:end],phi41[1:nskip:end].*rtd,label=L"$\phi_{efg}$",color=cp[4],alpha=0.5);
plot(t[1:nskip:end],phi42[1:nskip:end].*rtd,color=cp[4],alpha=0.5);
plot(t[1:nskip:end],phi43[1:nskip:end].*rtd,color=cp[4],alpha=0.5);
plot([t[1],t[end]],[0.0,0.0].+ 282.3,color=cp[4],linestyle="--");
phi51 = normalize.(eloft1[:,2,5] .- 2eloft1[:,2,6] +  eloft1[:,2,7] .+2pi);
phi51 = offset(phi51,nt);
phi52 = normalize.(eloft2[:,2,5] .- 2eloft2[:,2,6] +  eloft2[:,2,7] .+2pi);
phi52 = offset(phi52,nt);
phi53 = normalize.(eloft3[:,2,5] .- 2eloft3[:,2,6] +  eloft3[:,2,7] .+2pi);
phi53 = offset(phi53,nt);
plot(t[1:nskip:end],phi51[1:nskip:end].*rtd,label=L"$\phi_{fgh}$",color=cp[5],alpha=0.5);
plot(t[1:nskip:end],phi52[1:nskip:end].*rtd,color=cp[5],alpha=0.5);
plot(t[1:nskip:end],phi53[1:nskip:end].*rtd,color=cp[5],alpha=0.5);
plot([t[1],t[end]],[0.0,0.0] .+180.4,color=cp[5],linestyle="--");
legend(ncol=3,fontsize=20);
xlabel("Time [days]",fontsize=15);
ylabel(L"$\phi$ [$^\circ$]",fontsize=15);
;
println("bcd: ",mean(phi11)*rtd," ",mean(phi12)*rtd," ",mean(phi13)*rtd);
println("cde: ",mean(phi21)*rtd," ",mean(phi22)*rtd," ",mean(phi23)*rtd);
println("def: ",mean(phi31)*rtd," ",mean(phi32)*rtd," ",mean(phi33)*rtd);
println("efg: ",mean(phi41)*rtd," ",mean(phi42)*rtd," ",mean(phi43)*rtd);
println("fgh: ",mean(phi51)*rtd," ",mean(phi52)*rtd," ",mean(phi53)*rtd);
axis([t[1],t[end],0,360]);
tight_layout()
savefig("../Laplace_angle.pdf",bbox_inches="tight")
