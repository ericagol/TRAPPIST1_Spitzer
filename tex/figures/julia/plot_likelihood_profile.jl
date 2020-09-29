

using JLD2
using PyPlot
using Statistics
using Printf
nfile = 8
include("../../../src/loglinspace.jl")
include("../../../src/histogram_code.jl")

# Load in likelihood profile:
#@load "T1_likelihood_profile_student_all_3.0sig.jld2"
@load "../../../data/T1_likelihood_profile_student_all_3.0sig_v02.jld2"

# Load in the Markov chain dataset:
@load "../../../data/T1_hmc_total_02212020.jld2"
#@load "/Users/ericagol/Observing/Spitzer/DDT2019/Campaign04/v15/Hyak/T1_hmc_total_02212020.jld2"

## Now, make plots of each parameter versus the others:
#planet =["b","c","d","e","f","g","h"]
#var = ["m","P","t0","ecos","esin","log(nu)","V1 e^{1/(2nu)}"]
#
#for i=1:36
#  if i < 36
#    ip1 = ceil(Int64,i/5); ip2 = mod(i-1,5)+1
#    x = vec(elements_grid_all[:,ip1+1,ip2,:])
#    xname = string(planet[ip1]," ",var[ip2])
#  else
#    x = log.(vec(ndof_grid_all))
#    xname = var[6]
#  end
#  for j=i+1:37
##  for j=36:37
#    if j < 36
#      jp1 = ceil(Int64,j/5); jp2 = mod(j-1,5)+1
#      y = vec(elements_grid_all[:,jp1+1,jp2,:])
#      yname = string(planet[jp1]," ",var[jp2])
#    elseif j < 37
#      y = log.(vec(ndof_grid_all))
#      yname = var[6]
#    else
#      y = exp.(vec(lnV1_grid_all) .+ 1 ./(2vec(ndof_grid_all)))
#      yname = var[7]
#    end
#    if(abs(cor(x,y)) > 0.5)
#      clf()
#      cmap = exp.(-0.5 .*(vec(chi_grid_all) .- minimum(chi_grid_all)))
#      scatter(x[cmap .> 0.6],y[cmap .> 0.6],c=cmap[cmap .> 0.6])
#      println(i," ",xname," ",j," ",yname," ",cor(x,y))
#      read(stdin,Char)
#    end
#  end
#end

planet =["b","c","d","e","f","g","h"]
cp = ["C0","C1","C2","C3","C4","C5","C6","C7","C8","C9"]
x2=collect(linearspace(-0.015,0.015,1000))
fig,axes = subplots(4,2,figsize=(10,8),sharex="col",sharey="row")
for i=1:7
  ax = axes[i]
  x= elements_grid_all[(i-1)*5+4,i+1,4,:]; ecos0 = elements_grid_all[(i-1)*5+4,i+1,4,16]
  ecc = sqrt.(x.^2 .+ elements_grid_all[(i-1)*5+5,i+1,5,:].^2)
  prob = exp.(-0.5*(chi_grid_all[(i-1)*5+4,:] .-chi_grid_all[(i-1)*5+4,16]))
  ax.plot(x,prob./maximum(prob),color=cp[i],linewidth=1)
  prob2 = exp.(-0.5 .*(x2 .-ecos0).^2 ./cov_save[(i-1)*5+4,(i-1)*5+4])
  ax.plot(x2,prob2,color=cp[i],alpha=0.3,linewidth=1)
  # Plot histogram:
  ecos_bin,ecos_hist,ecos_bin_square,ecos_hist_square = histogram(state_total[(i-1)*5+4,:],50)
#  ax.plot(ecos_bin_square,ecos_hist_square ./maximum(ecos_hist_square),color=cp[i],linewidth=3,label=L"$e\cos{\omega}$)
  ax.plot(ecos_bin_square,ecos_hist_square ./maximum(ecos_hist_square),color=cp[i],linewidth=3,label=L"$e\cos \omega$")
  x = elements_grid_all[(i-1)*5+5,i+1,5,:]; esin0 = elements_grid_all[i*5,i+1,5,16]
  ecc = sqrt.(x.^2 + elements_grid_all[5i-4,i+1,4,:].^2)
  prob = exp.(-0.5*(chi_grid_all[(i-1)*5+5,:] .-chi_grid_all[(i-1)*5+5,16]))
  ax.plot(x,prob./maximum(prob),linestyle=":",color=cp[i],linewidth=1)
  prob2 = exp.(-0.5 .*(x2 .-esin0).^2 ./cov_save[i*5,i*5])
  ax.plot(x2,prob2,color=cp[i],linestyle=":",alpha=0.3,linewidth=1)
  # Plot histogram:
  esin_bin,esin_hist,esin_bin_square,esin_hist_square = histogram(state_total[(i-1)*5+5,:],50)
  ax.plot(esin_bin_square,esin_hist_square ./maximum(esin_hist_square),color=cp[i],linestyle=":",linewidth=3,label=L"$e\sin{\omega}$")
#  ax.plot(esin_bin_square,esin_hist_square ./maximum(esin_hist_square),color=cp[i],linestyle=":",linewidth=3)
  ax.plot([0,0],[0,1.05],linestyle="--",color=cp[i],linewidth=2)
  ax.legend(); 
  ax.axis([-0.0175,0.0175,0,1.05]); ax.annotate(string("(",planet[i],")"),xy=[-0.014;0.8])
  println(planet[i]," ",@sprintf("%6.4f",ecos0),"+-",@sprintf("%6.4f",sqrt(cov_save[5i-1,5i-1])),
     " ",@sprintf("%6.4f",esin0),"+-",@sprintf("%6.4f",sqrt(cov_save[5i,5i])))
  ax.grid(linestyle=":")
end
ax = axes[8]
#ax.axis("off")
# Plot the prior:
include("compute_ecc_prior.jl")
ax.grid(linestyle=":")
tight_layout()
subplots_adjust(hspace = 0,wspace=0)
savefig("../T1_eccentricity_vectors_likelihood_profile_hmc.pdf",bbox_inches="tight")
#read(stdin,Char)

# Make a plot of eccentricity histogram:

clf()
for i=1:7
  ecc_bin,ecc_hist,ecc_bin_square,ecc_hist_square = histogram(sqrt.(state_total[(i-1)*5+4,:].^2+state_total[(i-1)*5+5,:].^2),50)
  plot(ecc_bin_square,ecc_hist_square./maximum(ecc_hist_square),color=cp[i],linewidth=3,label=planet[i])
end
xlabel("Eccentricity",fontsize=15)
ylabel("Probability",fontsize=15)
legend(fontsize=15)
axis([0,0.015,0,1.05])
xticks(fontsize=15)
yticks(fontsize=15)
#read(stdin,Char)

savefig("../eccentricity_posterior.pdf",bbox_inches="tight")

#read(stdin,Char)

fig,axes = subplots(4,2)
for i=1:7
  ax = axes[i]
  # Plot period:
  P0 = elements_grid_all[(i-1)*5+2,i+1,2,16]; x= elements_grid_all[(i-1)*5+2,i+1,2,:] .- P0
  prob = exp.(-0.5*(chi_grid_all[(i-1)*5+4,:] .-chi_grid_all[(i-1)*5+4,16]))
  sigP = sqrt(cov_save[(i-1)*5+2,(i-1)*5+2])
  x2=collect(linearspace(-4sigP,4sigP,1000))
  ax.plot(x,prob,label="Period",color=cp[i])
  prob2 = exp.(-0.5 .*x2.^2 ./cov_save[(i-1)*5+2,(i-1)*5+2])
  ax.plot(x2,prob2,color=cp[i],alpha=0.3)
  t0 = elements_grid_all[(i-1)*5+3,i+1,3,16]
  prob = exp.(-0.5*(chi_grid_all[(i-1)*5+3,:] .-chi_grid_all[(i-1)*5+3,16]))
  x = elements_grid_all[(i-1)*5+3,i+1,3,:] .- t0
  ax.plot(x,prob,label=L"$t_0$",linestyle="--",color=cp[i])
  sigt0 = sqrt(cov_save[(i-1)*5+3,(i-1)*5+3])
  x2=collect(linearspace(-4sigt0,4sigt0,1000))
  prob2 = exp.(-0.5 .*x2.^2 ./cov_save[(i-1)*5+3,(i-1)*5+3])
  ax.plot(x2,prob2,color=cp[i],linestyle="--",alpha=0.3)
  ax.plot([0,0],[0,1],linestyle=":",color=cp[i])
  ax.legend(); 
  ax.axis([-0.015,0.015,0,1]); ax.annotate(string("(",planet[i],")"),xy=[-0.014;0.8])
  println(planet[i]," ",@sprintf("%6.4f",P0),"+-",@sprintf("%6.4f",sigP),
     " ",@sprintf("%6.4f",t0),"+-",@sprintf("%6.4f",sigt0))
end

# Finally plot log(ndof) & V1e^{1/(2nu)}:

fig,axes = subplots(1,2,sharey="row")

ax = axes[1]

#ndof0 = ndof_grid_all[nparam-1,16]; x= ndof_grid_all[nparam-1,:]
#prob = exp.(-0.5*(chi_grid_all[nparam-1,:] .-chi_grid_all[nparam-1,16]))
#signdof = sqrt(cov_save[nparam-1,nparam-1])
#x2=collect(linearspace(-4signdof,4signdof,1000)) 
#ax.plot(x,prob,label=L"$\nu$",color=cp[n])
#prob2 = exp.(-0.5 .*x2.^2 ./cov_save[nparam-1,nparam-1])
#x2 .+= ndof0
#ax.plot(x2,prob2,color=cp[n],alpha=0.3)
#lgndof0 = log(ndof_grid_all[nparam-1,16]); x= log.(ndof_grid_all[nparam-1,:])
lgndof0 = lgndof_grid_all[nparam-1,16]; x= lgndof_grid_all[nparam-1,:]
prob = exp.(-0.5*(chi_grid_all[nparam-1,:] .-chi_grid_all[nparam-1,16]))
signdof = sqrt(cov_save[nparam-1,nparam-1])/exp(lgndof0)
x2=collect(linearspace(-4signdof,4signdof,1000)) 
ax.plot(x,prob,label=L"$\nu$",color=cp[n],linewidth=3)
prob2 = exp.(-0.5 .*x2.^2 ./(cov_save[nparam-1,nparam-1]/exp(2lgndof0)))
x2 .+= lgndof0
ax.plot(x2,prob2,color=cp[n],alpha=0.3,linewidth=3)
# Plot histogram of log(ndof) parameter:
lndof_bin,lndof_hist,lndof_bin_square,lndof_hist_square = histogram(state_total[36,:],50)
ax.plot(lndof_bin_square,lndof_hist_square./maximum(lndof_hist_square))
ax.set_xlabel(L"$\log{\nu}$",fontsize=15)
ax.set_ylabel("Probability",fontsize=15)

ax = axes[2]
V1exp2nuinv0 = V1exp2nuinv_grid_all[nparam,16]; x= V1exp2nuinv_grid_all[nparam,:] 
prob = exp.(-0.5*(chi_grid_all[nparam,:] .-chi_grid_all[nparam,16]))
#sigV1exp2nuinv = sqrt(cov_save[nparam,nparam])
sigV1exp2nuinv =  0.09
x2=collect(linearspace(-4sigV1exp2nuinv,4sigV1exp2nuinv,1000))
ax.plot(x,prob,label=L"$\ln(V_1)$",color=cp[n+1],linewidth=3)
#prob2 = exp.(-0.5 .*x2.^2 ./cov_save[nparam,nparam])
prob2 = exp.(-0.5 .*x2.^2 ./sigV1exp2nuinv^2)
x2 .+= V1exp2nuinv0
ax.plot(x2,prob2,color=cp[n+1],alpha=0.3,linewidth=3)
# Plot histogram of this parameter:
V1expinv2nu_bin,V1expinv2nu_hist,V1expinv2nu_bin_square,V1expinv2nu_hist_square = histogram(state_total[37,:],50)
ax.plot(V1expinv2nu_bin_square,V1expinv2nu_hist_square./maximum(V1expinv2nu_hist_square))
ax.set_xlabel(L"$V_1 e^{1/(2\nu)}$",fontsize=15)
subplots_adjust(wspace=0)
savefig("../T1_students_params_transformed.pdf",bbox_inches="tight")
