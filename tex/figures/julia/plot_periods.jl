
# Make a plot of the periods of the planets versus integer values:
using PyPlot; using DelimitedFiles
using Statistics

include("../../../src/regress.jl")

data = readdlm("../../tables/times_forecast.txt",',')

# Now compute the periods and uncertainties:
period = Float64[]
sigp   = Float64[]
nplanet = 7
for i=1:nplanet
  indx = data[:,1] .== i
  nt = sum(indx)
  fn = zeros(2,nt)
  fn[1,:] .= 1.0
  fn[2,:] = data[indx,2]
  coeff,cov = regress(fn,data[indx,3],data[indx,4])
  push!(period,coeff[2])
  push!(sigp,sqrt(cov[2,2]))
end

#period = [ 1.510826, 2.421937, 4.049219, 6.101013, 9.207540, 12.352446, 18.772866]
#sigp = [0.000006, 0.000018, 0.000026, 0.000035, 0.000032, 0.000054, 0.000214]
i0 = [24,15,9,6,4,3,2]
n_i = 1 ./ period
sign = sigp .* n_i ./period
# Now, regress:
#fn = zeros(3,nplanet)
fn = zeros(2,nplanet)
fn[1,:] .= 1.0
fn[2,:] = i0
#fn[3,:] = i0.^2
coeff,cov = regress(fn,n_i,sign)
#errorbar(i0,period,sigp,fmt=".")
#errorbar(i0,n_i,sigp,fmt=".")
#errorbar(i0,n_i,sigp,fmt=".")
rc("xtick",labelsize=15)
rc("ytick",labelsize=15)
fig,axes = subplots(2,1,sharex="col")
ax = axes[1]
#ax.errorbar(i0,n_i .- coeff[2] .* i0,sigp,fmt=".")
ax.errorbar(i0,n_i,sign,fmt=".",label=L"$P^{-1}$")
#errorbar(i0,n_i .- coeff[2] .* i0 .- coeff[3] .* i0.^2,sigp,fmt=".")
#plot(i0,1 ./(coeff[1] .+ coeff[2] .* i0))
ax.plot([0;i0;25],[coeff[1];(coeff[1] .+ coeff[2] .* i0);coeff[1]+coeff[2]*25],label=L"$nP_0^{-1} - P_{TTV}^{-1}$")
#ax.plot([0,25],[coeff[1],coeff[1]])
#plot([0;i0;25],[0;coeff[2] .*i0; coeff[2]*25],linestyle=":")
#plot([0;i0;25],[0;coeff[2] .*i0; coeff[2]*25],linestyle=":")
pname = ["b","c","d","e","f","g","h"]
for i=1:nplanet
#  text(i0[i]+1,period[i]+0.1,string(pname[i],": ",i0[i]))
  ax.text(i0[i],n_i[i]-0.07,string(i0[i]))
  ax.text(i0[i],n_i[i]+0.03,string(pname[i]))
end
#ax.set_xlabel(L"n: Integer multiple of frequency $P_0^{-1}$")
ax.set_ylabel("1/P [1/days]",fontsize=15)
ax.legend()
#axis([0,25,0,20])
ax.axis([0,25,-0.1,0.8])
ax.plot([0,25],[0,0],linestyle=":")
#ax.text(2,0.6,L"$P^{-1} = n P_0^{-1} - P_{TTV}^{-1}$",fontsize=15)
#ax.text(5,0.65,L"$P_0 = 36.15007\pm 0.00013\ \mathrm{d} $",fontsize=10)
#ax.text(5,0.60,L"$P_{TTV} = 492.3\pm 0.1\ \mathrm{d} $",fontsize=10)
resid = (period .- 1 ./(coeff[1] .+ coeff[2] .* i0)) ./period .* 100
println("Residuals: ",resid," % ",std(resid))

ax = axes[2]
ax.errorbar(i0,n_i .- coeff[2] .* i0,sign,fmt=".",label=L"$P^{-1} - nP_0^{-1}$")
#ax.errorbar(i0,n_i .- coeff[2] .* i0 .- coeff[3] .* i0.^2,sign,fmt=".",label=L"$P^{-1} - nP_0^{-1}$")
ax.plot([0,25],[coeff[1],coeff[1]],label=L"$-P_{TTV}^{-1}$")
#ax.plot([0,25],[coeff[1]-sqrt(cov[1,1]),coeff[1]-sqrt(cov[1,1])],linestyle=":")
#ax.plot([0,25],[coeff[1]+sqrt(cov[1,1]),coeff[1]+sqrt(cov[1,1])],linestyle=":")
ax.set_xlabel(L"n: Integer multiple of frequency $P_0^{-1}$",fontsize=15)
ax.set_ylabel(L"1/P [1/days] - $nP_0^{-1}$",fontsize=15)
ax.legend()
ax.plot([0,25],[0,0],linestyle=":")
ax.axis([0,25,-0.003,0.001])
tight_layout()
subplots_adjust(hspace = 0)
savefig("T1_orbital_frequency_relation.png",bbox_inches="tight")
