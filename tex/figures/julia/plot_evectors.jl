

using JLD2; using PyPlot; using Statistics; using DelimitedFiles

# Load in the posteriors:
@load "../../../data/T1_hmc_total_02212020.jld2" state_total

clf()
# Make a confidence plot of eccentricity vectors, and
# show the longitudes of periastron.

nplanet = 7 

e2 = 0.015; e1 = -0.015
ncos = 30
nsin = 30
evecgrid = zeros(nplanet,ncos,nsin)
ecosgrid = zeros(ncos)
esingrid = zeros(nsin)
for i=1:ncos
  ecosgrid[i] = (i-0.5)*(e2-e1)/ncos+e1
  esingrid[i] = (i-0.5)*(e2-e1)/nsin+e1
end

nsamp = 224000

conf = [0.95, 0.683]
nconf = length(conf)
level_ecc = zeros(nplanet,nconf)
planet_name=["b","c","d","e","f","g","h"]

cname = ["C0","C1","C2","C3","C4","C5","C6","C7"]

for j=1:nplanet
  ncosp = 10000
  # Set up vectors for samples:
  ecossamp= zeros(nsamp)
  esinsamp = zeros(nsamp)
  # Loop over number of samples:
  for isamp=1:nsamp
    # Draw from posterior distribution:
    # Compute the mass of planet in Earth mass units:
    ecossamp[isamp] = state_total[(j-1)*5+4,isamp]
    esinsamp[isamp] = state_total[(j-1)*5+5,isamp]
    # Now compute density distribution:
    if ecossamp[isamp] > e1 && ecossamp[isamp] < e2 && esinsamp[isamp] > e1 && esinsamp[isamp] < e2
      ir = ceil(Int64,(esinsamp[isamp]-e1)/(e2-e1)*ncos)
      im = ceil(Int64,(ecossamp[isamp]-e1)/(e2-e1)*nsin)
      evecgrid[j,ir,im] +=1.0
    end
  end
  # Set levels:
  evecgrid[j,:,:] /= float(sum(evecgrid[j,:,:]))
  evecgrid_sort = sort(vec(evecgrid[j,:,:]))
  i0 = nsin*ncos
  evecgrid_cum = evecgrid_sort[i0]
  while evecgrid_cum < conf[1] && i0 > 1
    i0 -= 1
    for k=1:nconf
      if evecgrid_cum < conf[k] && (evecgrid_cum+evecgrid_sort[i0]) > conf[k]
        level_ecc[j,k] = 0.5*(evecgrid_sort[i0]+evecgrid_sort[i0+1])
      end
    end
    evecgrid_cum += evecgrid_sort[i0]
  end
  println(j," ",mean(ecossamp)," ",std(ecossamp)," fractional error: ",std(ecossamp)/mean(ecossamp)*100.0,"%")
end
mopt = readdlm("../../../data/optimum_values.txt")
figure(figsize=(6,6))
# Plot maximum likelihood values:
for j=1:nplanet
  cs = contour(ecosgrid,esingrid,evecgrid[j,:,:],levels = level_ecc[j,:],colors=cname[j],linewidth=2.0)
  plot(mopt[(j-1)*5+4],mopt[(j-1)*5+5],"o",label=planet_name[j],color=cname[j])
  text(mopt[(j-1)*5+4]+0.0005,mopt[(j-1)*5+5]-0.0005,planet_name[j],color=cname[j],fontsize=15)
end
#legend(fontsize=15)
xlabel(L"$e \cos{\omega}$",fontsize=15)
ylabel(L"$e \sin{\omega}$",fontsize=15)
axis([-0.0125,0.0125,-0.0125,0.0125])
plot([-0.015,0.015],[0,0],linestyle=":",color="k")
plot([0,0],[-0.015,0.015],linestyle=":",color="k")
tight_layout()
savefig("../esin_vs_ecos.pdf",bbox_inches="tight")
