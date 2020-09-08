

# Make corner plots for the paper:
using JLD2; using PyPlot; using Statistics; using DelimitedFiles

# Load in the posteriors:
@load "../../../data/T1_hmc_total_02212020.jld2" state_total

clf()

nplanet = 7 
nx = 30
ny = 30
#nsamp = 224000
nsamp = 1000

conf = [0.95, 0.683]
nconf = length(conf)
level_conf = zeros(nconf)
planet_name=["b","c","d","e","f","g","h"]

cname = ["C0","C1","C2","C3","C4","C5","C6","C7"]
mopt = readdlm("../../../data/optimum_values.txt")

fig,axes = subplots(35,35)
# Loop over the 35 parameters:
for i=1:34
  for j=i+1:35
    histgrid = zeros(nx,ny)
    xgrid = zeros(nx)
    x1 = minimum(state_total[i,:])
    x2 = maximum(state_total[i,:])
    y1 = minimum(state_total[j,:])
    y2 = maximum(state_total[j,:])
    ygrid = zeros(ny)
    for i=1:nx
      xgrid[i] = (i-0.5)*(x2-x1)/nx+x1
    end
    for i=1:ny
      ygrid[i] = (i-0.5)*(y2-y1)/ny+y1
    end
    # Set up vectors for samples:
    xsamp= zeros(nsamp)
    ysamp = zeros(nsamp)
    # Loop over number of samples:
    for isamp=1:nsamp
      # Draw from posterior distribution:
      # Compute the mass of planet in Earth mass units:
      xsamp[isamp] = state_total[i,isamp]
      ysamp[isamp] = state_total[j,isamp]
      # Now compute density distribution:
      if xsamp[isamp] > x1 && xsamp[isamp] < x2 && ysamp[isamp] > y1 && ysamp[isamp] < y2
        ir = ceil(Int64,(xsamp[isamp]-x1)/(x2-x1)*nx)
        im = ceil(Int64,(ysamp[isamp]-y1)/(y2-y1)*ny)
        histgrid[ir,im] +=1.0
      end
    end
    histgrid ./= float(sum(histgrid))
    histgrid_sort = sort(vec(histgrid))
    # Set levels:
    i0 = nx*ny
    histgrid_cum = histgrid_sort[i0]
    while histgrid_cum < conf[1] && i0 > 1
      i0 -= 1
      for k=1:nconf
        if histgrid_cum < conf[k] && (histgrid_cum+histgrid_sort[i0]) > conf[k]
          level_conf[k] = 0.5*(histgrid_sort[i0]+histgrid_sort[i0+1])
        end
      end
      histgrid_cum += histgrid_sort[i0]
    end
    ax = axes[i,j]
    cs = ax.contour(xgrid,ygrid,histgrid,levels = level_conf)
    ax.plot(mopt[i],mopt[j],"o")
    ax.set_xticks([])
    ax.set_yticks([])
    println(j," ",mean(xsamp)," ",std(xsamp)," fractional error: ",std(xsamp)/mean(xsamp)*100.0,"%")
  end
end
#figure(figsize=(6,6))
# Plot maximum likelihood values:
#for j=1:nplanet
#  text(mopt[(j-1)*5+4]+0.0005,mopt[(j-1)*5+5]-0.0005,planet_name[j],color=cname[j],fontsize=15)
#end
#legend(fontsize=15)
#xlabel(L"$e \cos{\omega}$",fontsize=15)
#ylabel(L"$e \sin{\omega}$",fontsize=15)
#axis([-0.0125,0.0125,-0.0125,0.0125])
#plot([-0.015,0.015],[0,0],linestyle=":",color="k")
#plot([0,0],[-0.015,0.015],linestyle=":",color="k")
#tight_layout()
#savefig("../esin_vs_ecos.pdf",bbox_inches="tight")
