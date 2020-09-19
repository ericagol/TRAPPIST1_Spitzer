2


# Read in NexSci planets with non-zero lower limits
# and both masses & radii measured.

if ~@isdefined(CGS)
  include("../../../src/CGS.jl")
  using Main.CGS
end

using DelimitedFiles; using PyPlot

clf()
figure(figsize=(6,6))
datadai = readdlm("../../../data/apjab3a3bt2_ascii.txt",',',comments=true,comment_char='#')
datadai[:,10] .*= -1
datadai[:,13] .*= -1
#datadaiflt = convert(Array{Float64,2},datadai[:,9:14])
#errorbar(datadai[:,9],datadai[:,12],xerr=permutedims(datadai[:,10:11]),yerr=permutedims(datadai[:,13:14]),fmt="o")
offsetx = zeros(13) .+ 0.35
#offsetx[13] *= 0.6  # K-93b
#offsetx[7] *= 1.1  # K2-141b
xerrs = reverse(permutedims(datadai[:,10:11]),dims=1)
yerrs = reverse(permutedims(datadai[:,13:14]),dims=1)
for i=1:13
  errorbar(datadai[i,9],datadai[i,12],xerr=reshape(xerrs[:,i],2,1),yerr=reshape(yerrs[:,i],2,1),fmt="o",alpha=0.18/minimum(xerrs[:,i]),color="C3",linewidth=2.0)
#  text([datadai[i,9]] .+ offsetx[i],[datadai[i,12]-0.01],datadai[i,1],fontsize=8)
#  text([1.0],[datadai[i,12]-0.01+offsety[i]],datadai[i,1],fontsize=4)
#  text([datadai[i,12]^3.7*offsetx[i]],[datadai[i,12]-0.01],datadai[i,1],fontsize=8)
#  text([datadai[i,12]^3.7*offsetx[i]],[datadai[i,12]-0.01],datadai[i,1],fontsize=5,rotation=90)
  println(i," ",offsetx[i]," ",datadai[i,1])
end
#dataeu = readdlm("exoplanet.eu_catalog.csv",',',skipstart=1,comments=true,comment_char='#')

#errorbar(dataeu[:,3].*(MJUPITER/MEARTH),dataeu[:,9].*(RJUPITER/REARTH),xerr=permutedims(dataeu[:,4:5].*(MJUPITER/MEARTH)),
#  yerr=permutedims(data[:,10:11].*(RJUPITER/REARTH)),fmt="o")

data = readdlm("../../../data/planets_2020.02.26_11.12.22.csv",',',skipstart=34,comments=true,comment_char='#')
data[:,14] .*= -1
data[:,18] .*= -1
#igood = (data[:,12] .> (5data[:,14])) .& (data[:,16] .> 12.5data[:,17]) .& (data[:,16] .> 12.5data[:,18])
igood = data[:,12] .> 5data[:,14]
#xerrs = reverse(permutedims(data[igood,13:14]),dims=1)
xerrs = reverse(permutedims(data[:,13:14]),dims=1)
#yerrs = reverse(permutedims(data[igood,17:18]),dims=1)
yerrs = reverse(permutedims(data[:,17:18]),dims=1)
#errorbar(data[igood,12],data[igood,16]./data[igood,12].^(1/3.7),xerr=xerrs,yerr=yerrs,fmt=".")
offset = zeros(31) .+0.35
#offset[15] *= 0.7 # Kepler-105c
#offset[13] *= 1.55 # GJ 357b
#offset[23] *= 0.75  # Kepler-36b
#offset[22] *= 0.8  # Kepler-10b
#offset[24] *= 1.2  # HD 219134 c
#offset[13] = -0.8; offset[18] = 0.6; offset[23] = -2.3; offset[25] = -2.4;  offset[24] = 0.5; offset[30] = -2.6; offset[27] = -2.4; offset[15] = -2.8; offset[31] = 0.7
for i=1:31
  if igood[i] 
    if i != 23
      errorbar(data[i,12],data[i,16],xerr=reshape(xerrs[:,i],2,1),yerr=reshape(yerrs[:,i],2,1),fmt="o",alpha=0.18/minimum(xerrs[:,i]),color="C3",linewidth=2.0)
    else
      errorbar(data[i,12],data[i,16],xerr=reshape(xerrs[:,i],2,1),yerr=reshape(yerrs[:,i],2,1),fmt="o",label="NexSci/Dai/Dressing",alpha=0.18/minimum(xerrs[:,i]),color="C3",linewidth=2.0)
    end
#    text([data[i,12]] .+ offset[i],[data[i,16]-0.01],string(data[i,1]," ",data[i,2]),fontsize=8)
#    text([1.0],[data[i,16]-0.01],string(data[i,1]," ",data[i,2]),fontsize=4)
#    text([data[i,16]^3.7*offset[i]],[data[i,16]-0.01],string(data[i,1]," ",data[i,2]),fontsize=5,rotation=90)
    println(i," ",offset[i]," ",data[i,1]," ",data[i,2])
  end
end
m = collect(0.05:0.01:10.0)
#semilogx(m,m.^(1/3.7)./m.^(1/3.7))
#semilogx(m,m.^(1/3.7),label="CMF=0.33")
#mrearth = readdlm("../../../data/massradiusEarthlikeRocky.txt")
mrearth = readdlm("../../../data/MR_trappist_Solar.ddat")
semilogx(mrearth[:,1],mrearth[:,2],label="Solar")
#mrearth_case5 = readdlm("zeng2013_tab02_case2.txt")
#plot(mrearth_case5[10:29,1],mrearth_case5[10:29,2],label="Case 2 CMF=0.1787 (Zeng & Sasselov 2013)")
# Read in the ternary diagram (Table 3) from Zeng & Sasselov (2013), but 
# without any of the water planet cases:
#mrearth_cmf = readdlm("zeng2013_tab03_dry.txt",skipstart=1)
#mrearth_cmf = reshape(mrearth_cmf,5,41,8)
## Find a 17% CMF case:
#mrearth_17pct = zeros(2,41)
#mrearth_17pct[1,:] = mrearth_cmf[1,:,1]
#pctfe = 0.17
#for i=1:41
#  for j=1:4
#    if (mrearth_cmf[j,i,6] > pctfe) && (mrearth_cmf[j+1,i,6] < pctfe)
#      x = (mrearth_cmf[j,i,6] - pctfe)/(mrearth_cmf[j,i,6] - mrearth_cmf[j+1,i,6])
#      mrearth_17pct[2,i] = mrearth_cmf[j,i,2] + (mrearth_cmf[j+1,i,2] - mrearth_cmf[j,i,2])*x
#    end
#  end
#end
#plot(mrearth_17pct[1,1:30],mrearth_17pct[2,1:30],label="Zeng & Sasselov (2013); 17% CMF")
#semilogx(m,m.^(1/3.7)*(1.07-0.21*0.25)./m.^(1/3.7))
#semilogx(m,m.^(1/3.7)*(1.07-0.21*0.25),label="CMF=0.25")

# Now plot trappist:

mrt1 = readdlm("../../../data/t1_mass_radius.txt",',')
#errorbar(mrt1[:,2],mrt1[:,4]./mrt1[:,2].^(1/3.7),xerr=mrt1[:,3],yerr=mrt1[:,5],fmt=".")
errorbar(mrt1[:,2],mrt1[:,4],xerr=mrt1[:,3],yerr=mrt1[:,5],fmt=".",label="Trappist-1",linewidth=2.0)

# Plot Earth, Venus, Mars, Mercury(?):
mss = [5.9736e27,0.3302e27,4.8685e27,0.64185e27]./5.9736e27
rss = [6371.00e5,2439.7e5,6051.8e5,3389.5e5]./6371.0e5

#plot(mss,rss./mss.^(1/3.7),"o")
plot(mss,rss,"o",label="Solar System")

#dmr_carbon = readdlm("../Figures/Dorn/M_R_carbon-rich_Eric.ddat")
dmr_nocore_solar = readdlm("../../../data/MR_trappist_corefree_Solar_v02.ddat",skipstart=1)
#dmr_nocore_unter = readdlm("../Figures/Dorn/MR_trappist_corefree_Unterborn.ddat",skipstart=1)

#plot(dmr_carbon[:,1],dmr_carbon[:,2],label="Carbon-rich (Miozzi et al. 2018)")
plot(dmr_nocore_solar[:,1],dmr_nocore_solar[:,2],label="No-core, Solar",color="C4")
#plot(dmr_nocore_unter[:,1],dmr_nocore_unter[:,2],label="No-core, U17")

axis([0.04,10.0,0.35,2.0])

xlabel(L"$M/M_\oplus$")
#ylabel(L"$(R/R_\oplus)(M_\oplus/M)^{0.27}$")
ylabel(L"$R/R_\oplus$")
legend()
tight_layout()
text([0.28],[0.77],"h")
text([0.35],[0.81],"d") 
text([0.64],[0.94],"e")
text([0.96],[1.05],"f")
text([1.22],[1.15],"g")
text([1.15],[1.08],"c")
text([1.40],[1.125],"b")

text([1.07],[0.98],"Earth")
text([0.88],[0.93],"Venus")
text([0.12],[0.52],"Mars")
text([0.064],[0.38],"Mercury")
savefig("../mass_radius_relation_comparison.pdf",bbox_inches="tight")
