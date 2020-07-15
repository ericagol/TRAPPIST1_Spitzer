 


# Load in, combine, and plot chains:
using DelimitedFiles
using JLD2
using Statistics
using Printf
using PyPlot
include("nlog_prior.jl")
if !@isdefined(CGS)
  include("../../src/CGS.jl")
  using Main.CGS
end

n = 8
nplanet = n-1
nparam = (n-1)*5+2
nstep = 2000


#@load "T1_hmc_total_02132020.jld2"
@load "../../data/T1_hmc_total_02212020.jld2"

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


# Now print mass ratios:
xfit = zeros(2,nplanet)
fac = 0.09*MSUN/MEARTH
#println("\$M_p/M_* (M_\\oplus/0.09 M_\\odot)\$ & \$ M_p/M_*\$ & \$\\sigma(M_p/M_*)/(M_p/M_*)\$ & \$ P\$ [day] & \$ t0\$ [JD - 2,450,000] & \$ e\\cos{\\omega}\$ & \$ e\\sin{\\omega}\$  \cr")
for i=1:nplanet
  planetstring = string(planet[i]," &")
  xfit[:,i] = [mean(state_total[(i-1)*5+1,:])*fac,std(state_total[(i-1)*5+1,:])*fac]
  global planetstring = string(planetstring," \$",@sprintf("%6.4f",xfit[1,i])," \\pm ",@sprintf("%6.4f",xfit[2,i]),"\$ &")
  xfit[:,i] = [mean(state_total[(i-1)*5+1,:])*1e5,std(state_total[(i-1)*5+1,:])*1e5]
  global planetstring = string(planetstring," \$",@sprintf("%6.3f",xfit[1,i])," \\pm ",@sprintf("%6.3f",xfit[2,i]),"\$ &")
  planetstring = string(planetstring," \$",@sprintf("%3.2f",xfit[2,i]/xfit[1,i]*100),"\$ &")
  xfit[:,i] = [mean(state_total[(i-1)*5+2,:]),std(state_total[(i-1)*5+2,:])]
  planetstring= string(planetstring," \$",@sprintf("%8.6f",xfit[1,i])," \\pm ",@sprintf("%8.6f",xfit[2,i]),"\$ &")
  xfit[:,i] = [mean(state_total[(i-1)*5+3,:]),std(state_total[(i-1)*5+3,:])]
  planetstring= string(planetstring," \$",@sprintf("%7.5f",xfit[1,i])," \\pm ",@sprintf("%7.5f",xfit[2,i]),"\$ &")
  xfit[:,i] = [mean(state_total[(i-1)*5+4,:]),std(state_total[(i-1)*5+4,:])]
  planetstring= string(planetstring," \$",@sprintf("%6.5f",xfit[1,i])," \\pm ",@sprintf("%6.5f",xfit[2,i]),"\$ &")
  xfit[:,i] = [mean(state_total[(i-1)*5+5,:]),std(state_total[(i-1)*5+5,:])]
  planetstring= string(planetstring," \$",@sprintf("%6.5f",xfit[1,i])," \\pm ",@sprintf("%6.5f",xfit[2,i]),"\$ \\cr ")
  println(planetstring)
end

println("\$ \\log{\\nu} = ",@sprintf("%6.4f",mean(state_total[36,:]))," \\pm ",@sprintf("%6.4f",std(state_total[36,:])),"\$")
println("\$ V_1 e^{1/(2\\nu)} = ",@sprintf("%6.4f",mean(state_total[37,:]))," \\pm ",@sprintf("%6.4f",std(state_total[37,:])),"\$")
