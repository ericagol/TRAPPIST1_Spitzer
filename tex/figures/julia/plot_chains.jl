


# Load in, combine, and plot chains:
using DelimitedFiles
using JLD2
using Statistics
using Printf
using MCMCDiagnostics
using PyPlot
include("../nlog_prior.jl")
if !@isdefined(CGS)
  include("/Users/ericagol/Computer/Julia/CGS.jl")
  using Main.CGS
end

n = 8
nplanet = n-1
nparam = (n-1)*5+2
nstep = 2000

#fnames = readdlm("jld2_list.txt")
#fnames = readdlm("jld2_list_eps0.1.txt")
fnames = readdlm("jld2_list_nstep2000_eps0.1.txt")

nfile = size(fnames)[1]
#ncurr = 250
#ncurr = 350
ncurr = 2000
state_total = zeros(nparam*2+2,ncurr*nfile)

nacc_all = zeros(Int64,nfile)
neff_tot = zeros(nfile,37)
for i=1:nfile
  i1 = (i-1)*ncurr+1
  i2 = i*ncurr
  @load fnames[i] state nacc
  nacc_all[i] = nacc
  state_total[:,i1:i2] = state[:,1:ncurr]
  for j=1:37
    neff_tot[i,j] += effective_sample_size(state[j,1:ncurr])
  end
end
for j=1:37
  println(j," ",sum(neff_tot[:,j]))
end

clf()
plot(nacc_all)
read(stdin,Char)

@load "T1_run_hmc_student_ecc_lgndof_V1exp2nuinv_nstep2000_eps0.1_nleap20_201.jld2"

x_opt = [4.598791344049799e-5, 1.5108213441190745, 7257.550484415569, -0.00495015358026971, 0.004671261296998595, 
         4.40588192304326e-5, 2.4219505838480635, 7258.587230337616, -0.001188890950465861, 0.0014051713506332184, 
         1.305145694403463e-5, 4.049193371144702, 7257.06806273418, -0.0064708229088599216, 0.003529480254615964, 
         2.3551037452505444e-5, 6.101020579519802, 7257.827698633612, 0.003217590595460454, -0.003935834596847923,
         3.4925310671129774e-5, 9.207561405004297, 7257.073658382135, -0.009581969253521241, -0.00016501113996023978,
         4.417497904782645e-5, 12.352405680378052, 7257.715102776943, 0.00288626733667452, 0.0016434273040851228, 
         1.0931051156428014e-5, 18.772977044007867, 7249.605686691679, -0.004439910317349172, 0.0002702876058394804, 
         1.2119892906238703, 0.8222963482043892]

planet = ["b","c","d","e","f","g","h"]
param = ["mu","P","t0","ecos","esin"]
prior = zeros(nfile*ncurr)
for j=1:nfile*ncurr
  prior[j] = nlog_prior(state_total[1:37,j])
end

# compute f0:
#f0total = zeros(nfile*ncurr)
#for j=1:nfile*ncurr
#  f0total[j] = f0(state_total[1:37,j])
#end

for iparam=1:nparam
  clf()
  cparam = state_total[iparam,:]
#  plot(cparam[cparam .== circshift(cparam,100)],".")
  plot(cparam,".")
  plot([0,ncurr*nfile],x_opt[iparam] .+ [0,0])
  plot([0,ncurr*nfile],x_opt[iparam] + sqrt(cov_save[iparam,iparam]) .+ [0,0],linestyle="--")
  plot([0,ncurr*nfile],x_opt[iparam] - sqrt(cov_save[iparam,iparam]) .+ [0,0],linestyle="--")
#  axis([0,700,minimum(cparam),maximum(cparam)])
  if iparam < nparam-1
    pname = string(planet[div(iparam-1,5)+1]," ",param[mod(iparam-1,5)+1])
  elseif iparam == nparam-1
    pname = "log(nu)"
  else
    pname = "V1 e^{1/(2nu)}"
  end
  println(iparam," ",pname," chains: ",@sprintf("%6.3g",mean(state_total[iparam,:])),"+-",@sprintf("%6.3g",std(state_total[iparam,:]))," Fisher: ",@sprintf("%6.3g",x_opt[iparam]),"+-",@sprintf("%6.3g",sqrt(cov_save[iparam,iparam]))," ratio: ",@sprintf("%6.3g",std(state_total[iparam,:])/sqrt(cov_save[iparam,iparam])))
  read(stdin,Char)
end

clf()
# Probability distribution of e*sin(omega) vs. e*cos(omega)
plot(state_total[4,:],state_total[5,:],".",alpha=0.02)
plot(state_total[9,:],state_total[10,:],".",alpha=0.02)
plot([0,0],[0,0],".")
xlabel(L"$e_b \cos{\omega_b}$")
ylabel(L"$e_b \sin{\omega_b}$")
read(stdin,Char)

include("histogram_code.jl")

# plot histogram of each

esinb_bin, esinb_hist, esinb_bin_square,esinb_hist_square = histogram(state_total[4,:],50)
ecosb_bin, ecosb_hist, ecosb_bin_square,ecosb_hist_square = histogram(state_total[5,:],50)
esinc_bin, esinc_hist, esinc_bin_square,esinc_hist_square = histogram(state_total[9,:],50)
ecosc_bin, ecosc_hist, ecosc_bin_square,ecosc_hist_square = histogram(state_total[10,:],50)

clf(); 
plot(esinb_bin_square,esinb_hist_square ./ maximum(esinb_hist_square),label=L"$e_b \sin{\omega_b}$")
plot(ecosb_bin_square,ecosb_hist_square ./ maximum(ecosb_hist_square),label=L"$e_b \cos{\omega_b}$")
plot(esinc_bin_square,esinc_hist_square ./ maximum(esinc_hist_square),label=L"$e_c \sin{\omega_c}$")
plot(ecosc_bin_square,ecosc_hist_square ./ maximum(ecosc_hist_square),label=L"$e_c \cos{\omega_c}$")
axis([-0.01,0.01,0,1.05])
legend()
read(stdin,Char)

# Now plot histogram of longitudes of inner two planets:
 
#delta_omega = atan.(state_total[10,:],state_total[9,:]) .- atan.(state_total[5,:],state_total[4,:]) .-pi
ecosc = state_total[9,:]
esinc = state_total[10,:]
ecosb = state_total[4,:]
esinb = state_total[5,:]
delta_omega = atan.(ecosc .*esinb .-ecosb .*esinc, -ecosb .*ecosc .-esinb .*esinc)

delta_omega[delta_omega .< -pi] .+= 2\pi
delta_omega[delta_omega .< -pi] .+= 2\pi
delta_omega[delta_omega .>  pi] .-= 2\pi
delta_omega[delta_omega .>  pi] .-= 2\pi

domega_bin, domega_hist, domega_bin_square, domega_hist_square = histogram(delta_omega,50)

#domega_hist_square[domega_hist_square .== 0.0] .+= 200.0
clf(); plot(domega_bin_square*180/pi,domega_hist_square./maximum(domega_hist_square),linewidth=2)

xlabel(L"$\omega_c - \omega_b - 180^\circ$",fontsize=20)
ylabel("Probability",fontsize=20)
axis([-180,180,0,1.05])

# Now print mass ratios:
xfit = zeros(2,nplanet)
fac = 0.09*MSUN/MEARTH
for i=1:nplanet
  xfit[:,i] = [mean(state_total[(i-1)*5+1,:])*fac,std(state_total[(i-1)*5+1,:])*fac]
  println("planet: ",i," mass: ",@sprintf("%6.4f",xfit[1,i]),"+-",@sprintf("%6.4f",xfit[2,i]))
end

@save "T1_hmc_total_02212020.jld2" state_total
