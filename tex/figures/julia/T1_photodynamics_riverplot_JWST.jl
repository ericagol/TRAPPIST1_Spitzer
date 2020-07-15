

# Creating a photodynamical model of TRAPPIST-1 based on the Spitzer data:
# Create a function to read in the Spitzer datasets:

using DelimitedFiles
using Statistics
using PyPlot
using Optim
if ~@isdefined(CGS)
  include("../../../src/CGS.jl")
  using Main.CGS
end

# Now, call ttv to get a list of transit times and pair these with the data:

include("../../../src/NbodyGradient/src/ttv.jl")

# Read in initial orbital elements:
elements = readdlm("../../../data/elements_noprior_students.txt",',')

# Set integration parameters:
n=8
nplanet = n-1
t0 = 7257.93115525
h  = 0.06
tmax = 3000.0

# Set up array for holding times of transit:
ntt = zeros(Int64,n);
for i=2:n
  ntt[i] = ceil(Int64,tmax/elements[i,2])+1
end
tt1 = zeros(3,n,maximum(ntt));
# Save a counter for the actual number of transit times of each planet:
count1 = zeros(Int64,n);
# Set the "size" of the star (only checks that close to transit):
rstar = 1e11
# Call the ttv function:
@time dq=ttvbv_elements!(n,t0,h,tmax,elements,tt1,count1,0.0,0,0,rstar)


# Now, call Limbdark to create a model for each planet:

using Limbdark

# Take the radius ratios and limb-darkening parameters from
# Delrez et al. (2018):

depth_ch2  = [0.7277,0.6940,0.3566,0.4802,0.634,0.764,0.346] .* 0.01

rp = sqrt.(depth_ch2) # Ratio of planet radius to stellar radius

# Fractional radius errors from Delrez et al.:
srp = [0.028,0.028,0.020,0.025,0.026,0.029,0.025]./[1.127,1.100,0.788,0.915,1.052,1.154,0.777].*rp

# Set limb-darkening parameters:

u_ch1 = [0.168,0.244]
u_ch2 = [0.161,0.208]

# Set transit durations:
T = [36.19,42.31,49.33,55.92,63.14,68.53,76.92] ./(60*24)

# Set impact parameters:
b = [0.157,0.148,0.08,0.240,0.337,0.406,0.392]

# Set the stellar density:

include("../../../src/regress.jl")

function compute_model_all(tt1,count,rho_T1,texp,b,rp,u_ch2,neach)
      # Create a photometric model for all data.
      # Number of points per transit:
      # Create a matrix to hold the simulations:
      phot_mod = zeros(sum(count),neach)
      # Create a time array for each transit:
      t_mod = zeros(sum(count),neach)
      # Create array to tell me which transit belongs to which planet:
      ip_trans = zeros(Int64,sum(count))

      # Density of T1 star in units of solar density:
      fac = rho_T1^(1//3)*AU/RSUN

      # Exposure time in days:
      texp_day = texp/(60*24)

      # Two limb-darkening coefficients:
      nu = 2

      itrans = 0
      # Initialize transit structure (no derivatives for now):
      trans =  Limbdark.transit_init(0.1,0.5,u_ch2,false)
      # Loop over each planet:
        
      for ip=1:nplanet
      # Carry out a fit to the simulated transits:
        tt_ref1 = zeros(count1[ip+1])
        fn = zeros(Float64,2,count1[ip+1])
        sig = ones(count1[ip+1])
        tti1 = tt1[1,ip+1,1:count1[ip+1]]
        for j=1:count1[ip+1]
          fn[1,j] = 1.0
          fn[2,j] = round(Int64,(tti1[j]-elements[ip+1,3])/elements[ip+1,2])
          epoch = round((tti1[j]-elements[ip+1,3])/elements[ip+1,2])
        end
        coeff,cov = regress(fn,tti1,sig)
        println(ip," ",coeff)
        tt_ref1 = coeff[1] .+coeff[2] .*fn[2,:]
      # Now, Loop over each transit:
      for ifile=1:count[ip+1]
        # Create an time array for the current transit:
        tcurr = tt_ref1[ifile] .+ (collect(-div(neach,2):div(neach,2)-1) .+0.5) .*texp_day
        # First set the limb-darkening parameters
        # for the current light curve
        
        # Number of exposure times for this case:
        nt = neach
        # Set up array to hold the time-integrated flux for each planet:
        favg1 = zeros(5+nu,nt)
        # Set up arrays to hold statistics on integrations:
        neval1 = zeros(Int64,nt)
        depthmax1 = zeros(Int64,nt)
        # Loop over each transit:
        itrans += 1
        ip_trans[itrans] = ip
        ttbv = tt1[:,ip+1,ifile]
        # Now, create a model.  Convert impact parameter and
        # sky velocity to units of radius of the star:
        param = [ttbv[1],ttbv[2]*fac,ttbv[3]*fac]
        # Set the radius ratio and impact parameter to the current planet:
        param[3] = b[ip]
        trans.r = rp[ip]
        Limbdark.integrate_lightcurve!(trans,param,tcurr,texp_day,favg1,nt,1e-11,40,neval1,depthmax1)
        phot_mod[itrans,:] .= favg1[1,:]
        t_mod[itrans,:] .= tcurr
      end
    end
    return t_mod,ip_trans,phot_mod
end

rho_T1 = 52.5363620668614

neach = 400

t_mod_all,ip_trans, phot_mod_all = compute_model_all(tt1,count1,rho_T1,0.5,b,rp,u_ch2,neach)

#pygui(false)
# Now, emplace into an image:
pname = ["b","c","d","e","f","g","h"]
maxdepth = minimum(phot_mod_all)
fig,axes = subplots(1,7,figsize=[12,5])
# Range of the transit x-axis:
texp_plot = 0.5
t_xaxis = (collect(-div(neach,2):div(neach,2)-1) .+0.5) .*texp_plot
tmin = (t_xaxis[1]-0.5*texp_plot)/60
tmax = (t_xaxis[neach]+0.5*texp_plot)/60
stretch = 1
for i=1:nplanet
  ax = axes[i]
  if i > 2
    img = zeros(count1[i+1]*stretch,neach)
    i0 = sum(count1[1:i])
    for j=1:count1[i+1]
      for k=1:stretch
        img[(j-1)*stretch+k,:] = phot_mod_all[i0+j,:]
      end
    end
  else 
    img = zeros(count1[i+1]*stretch,div(neach,2))
    i0 = sum(count1[1:i])
    for j=1:count1[i+1]
      for k=1:stretch
        img[(j-1)*stretch+k,:] = phot_mod_all[i0+j,div(neach,4)+1:neach-div(neach,4)]
      end
    end   
  end
  
  # Set one pixel equal to maximum depth:
  img[1,175]=maxdepth
  print(i," ",size(img))
  if i > 2
    #ax.imshow(img,aspect="auto",cmap="winter",interpolation="nearest",origin="lower",extent = [tmin,tmax,0,3000])
    ax.imshow(img,aspect="auto",cmap="winter",interpolation="nearest",origin="lower",extent = [tmin,tmax,0,3000])
  else
    #ax.imshow(img,aspect="auto",cmap="winter",interpolation="nearest",origin="lower",extent = [tmin/2,tmax/2,0,3000])
    ax.imshow(img,aspect="auto",cmap="winter",interpolation = "nearest",origin="lower",extent = [tmin/2,tmax/2,0,3000])
  end
  ax.title.set_text(pname[i])
  #ax.set_xticks(fontsize=8)
  #ax.title(pname[i])
  #ax.tick_params(axis="x",labelsize=6)
  if i == 4 
    ax.set_xlabel(L"Time [hr] - $\mathdefault{t_n}$")
  end
  if i == 1
    ax.set_ylabel(L"Time [days] - $\mathdefault{t_0}$")
  else
        ax.set_yticks([])
  end
end
subplots_adjust(hspace=0,wspace=0.05)
savefig("../T1_riverplot_JWST.pdf",bbox_inches="tight")

i0 = 0
pc = ["C0","C1","C2","C3","C4","C5","C6","C7"]
for i=1:7
  for j=1:count1[i+1]
  plot(t_mod_all[i0+j,:],phot_mod_all[i0+j,:],color=pc[i])
  end
  global i0 += count1[i+1]
end
