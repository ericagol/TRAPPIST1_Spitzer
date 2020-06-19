

# Runs a "photodynamic" model (in background).
# (Dynamics is currently held fixed).

using JLD2
#using Limbdark
include("/gscratch/astro/agol/julia-1.0.5/bin/dev/Limbdark/src/integrate_lightcurve.jl")
using DelimitedFiles
using Statistics
using PyPlot
using AffineInvariantMCMC

if ~@isdefined(CGS)
#  include("/Users/ericagol/Computer/Julia/CGS.jl")
  include("/gscratch/astro/agol/src/CGS.jl")
  using Main.CGS
end

#include("/Users/ericagol/Computer/Julia/regress.jl")
 include("/gscratch/astro/agol/src/regress.jl")
#include("/Users/ericagol/Software/NbodyGradient/src/ttv.jl")
 include("/gscratch/astro/agol/NbodyGradient/src/ttv.jl")

function run_pd(datafile,foutput,numwalkers,burnin,thinning,astep,numsamples_perwalker;nout=1000,
    limbprior=false,incprior=false,contamination=false)
  # Load in the Spitzer data:
  @load datafile tobs phot error i1 i2 nfile channel nexp nspitzer fnames texp

  # Now, call ttv to get a list of transit times and pair these with the data:
  # Change this path to point to ttv.jl:  
  
  # Read in optimum orbital elements:
  elements = readdlm("../../data/elements_noprior_students.txt",',')
  
  # Set integration parameters:
  n=8
  nplanet = n-1
  t0 = 7257.93115525
  h  = 0.06
  tmax = 1600.0
  
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

  # Match the data to the model:
  # Take the radius ratios and limb-darkening parameters from
  # Delrez et al. (2018):
  
  depth_ch2  = [0.7277,0.6940,0.3566,0.4802,0.634,0.764,0.346] .* 0.01
  
  rp = sqrt.(depth_ch2) # Ratio of planet radius to stellar radius
  
  # Fractional radius errors from Delrez et al.:
  srp = [0.028,0.028,0.020,0.025,0.026,0.029,0.025]./[1.127,1.100,0.788,0.915,1.052,1.154,0.777].*rp
  
  # Set limb-darkening parameters:
  
#  u_ch1 = [0.168,0.244]
#  u_ch2 = [0.161,0.208]
  
  # Set transit durations:
  T = [36.19,42.31,49.33,55.92,63.14,68.53,76.92] ./(60*24)
  
  # Set impact parameters:
  b = [0.157,0.148,0.08,0.240,0.337,0.406,0.392]

  # Set orbital periods:
  P = [1.510826,2.421938,4.049218,6.101013,9.207541,12.352445,18.772863]  # Orbit periods in days

  # Set variable for computing semi-major axis to star ratio:
  aonr0 = (GRAV*MSUN*(24*3600)^2/(4pi^2*RSUN^3))^(1//3)
 
  # Now, loop over the datasets, figure out which planet(s) are
  # transiting in each, and compute the light curves:
  ip_inc = Int64[]   # Which planet for each transit
  ind_inc = Int64[]  # Transit indices for each planet
  ninc = zeros(Int64,nfile)  # Number of planet transits to include for each file
  for i=1:nfile
    # See which planets to include:
    for ip=1:nplanet
      tobs_start = tobs[i1[i]] - T[ip]/2
      tobs_end = tobs[i2[i]] + T[ip]/2
      for j=1:count1[ip+1]
       # See if the transit mid-point is less than half of a transit duration
        if (tt1[1,ip+1,j] > tobs_start) && (tt1[1,ip+1,j] < tobs_end)
          push!(ip_inc,ip)
          push!(ind_inc,j)
          ninc[i] +=1
        end
      end
    end
  end

  # Create a function which fits a polynomial and subtracts it:
  function fit_poly(t,f,ord)
      nt = size(t)[1]
      # create the basis functions:
      x = zeros(ord+1,nt)
      x[1,:] .= 1.0
      t0 = mean(t)
      sig = ones(nt)
      for i=1:ord
        x[i+1,:] .= (t .- t0).^i
      end
      coeff,cov = regress(x,f,sig)
      mod = zeros(nt) .+ coeff[1]
      for i=1:ord
        mod .+= coeff[i+1] .* x[i+1,:]
      end
      return mod
  end


#  function compute_model(time,phot,nfile,x0[15],texp,x0[8:14],x0[1:7],x0[16:17],x0[18:19],ninc,ip_inc,ind_inc,tt1)
  function compute_model(time,phot,nfile,rho_T1,texp,b,rp,u_ch1,u_ch2,ninc,ip_inc,ind_inc,tt1)
        # Create a photometric model with detrending:
        phot_mod = copy(phot)
        fill!(phot_mod,0.0)
        # Create a vector for a model of stellar/instrumental variabilty:
        phot_vary = copy(phot)
        fill!(phot_vary,0.0)
  
        # Density of T1 star in units of solar density:
        fac = rho_T1^(1//3)*AU/RSUN
  
        # Exposure time in days:
        texp_day = texp/(60*24)
  
        # Two limb-darkening coefficients:
        nu = 2
  
        itrans = 0
        # Initialize transit structure (no derivatives for now):
        #trans =  Limbdark.transit_init(0.1,0.5,[0.2,0.2],false)
        #trans =  transit_init(0.1,0.5,[0.2,0.2],false)
        # Loop over each light curve:
          
        for ifile=1:nfile
          tcurr = tobs[i1[ifile]:i2[ifile]]
          # First set the limb-darkening parameters
          # for the current light curve
          if channel[ifile] == "c1"
            #trans.u_n = u_ch1
            trans =  transit_init(0.1,0.5,u_ch1,false)
          else
            #trans.u_n = u_ch2
            trans =  transit_init(0.1,0.5,u_ch2,false)
          end
          # Number of exposure times for this case:
          nt = i2[ifile]-i1[ifile]+1
          # Set up array to hold the time-integrated flux for this light curve:
          #favg = zeros(5+nu,nt)
          favg = ones(5+nu,nt)
          # Set up array to hold the time-integrated flux for each planet:
          favg1 = zeros(5+nu,nt)
          # Set up arrays to hold statistics on integrations:
          neval1 = zeros(Int64,nt)
          depthmax1 = zeros(Int64,nt)
          # Loop over each transit:
          for it = 1:ninc[ifile]
            fill!(favg1,0.0)
           itrans += 1
            ttbv = tt1[:,ip_inc[itrans]+1,ind_inc[itrans]]
            # Now, create a model.  Convert impact parameter and
            # sky velocity to units of radius of the star:
            param = [ttbv[1],ttbv[2]*fac,ttbv[3]*fac]
            # Set the radius ratio and impact parameter to the current planet:
            param[3] = b[ip_inc[itrans]]
            trans.r = rp[ip_inc[itrans]]
            #Limbdark.integrate_lightcurve!(trans,param,tcurr,texp_day,favg1,nt,1e-11,40,neval1,depthmax1)
            integrate_lightcurve!(trans,param,tcurr,texp_day,favg1,nt,1e-11,40,neval1,depthmax1)
            favg .+= favg1
            # Subtract off 1 for additional planets
            #if it > 1
            #  favg[1,:] .-= 1
            #end
          end
          # Now, remove a polynomial:
          phot_mod[i1[ifile]:i2[ifile]] .= favg[1,:]
          phot_vary[i1[ifile]:i2[ifile]] .= fit_poly(tcurr,phot[i1[ifile]:i2[ifile]]./favg[1,:],3)
        end
      return phot_mod,phot_vary
  end

  # Compute the scatter:
  # Previously optimized parameters:
#  x0 = [ 0.08510408205855542, 0.08351749544695007, 0.060057015907927426 , 0.07003499168340532, 0.07938534728488905, 0.08571785728160501, 0.05767832213757894, 0.06777531110636684, 0.07268527380271149, 0.0003227708869417413, 0.17817318017323253, 0.2926207489330042, 0.36533945405158097, 0.3735584459205608, 52.46814828992533, 0.0961562252101869, 0.5059574689387107, 0.30961490912899814, 0.29644507311377927]
  x0 = [0.08511305133271845, 0.08352547037122816, 0.06005344100378423, 0.0700421002759512, 0.07962187102240074, 0.08611264006593941, 0.05767381066632817, 0.07876438173347886, 0.08351738136283536, 0.0008769148946956045, 0.18225051020178037, 0.30070386193030363, 0.36817942906674744, 0.3756424617589329, 52.34377029283971, 0.11235270319764341, 0.42037661035916857, 0.352424321959808, 0.2864053200404355]

  # Convert limb-darkening:
  q = x0[16:19]
  u_ch1_start = [2*sqrt(q[1])*q[2],sqrt(q[1])*(1 -2q[2])]
  u_ch2_start = [2*sqrt(q[3])*q[4],sqrt(q[3])*(1 -2q[4])]

  # This function computes an initial model in order to derive the scatter parameters (with no contamination or priors):
  phot_mod,phot_vary = compute_model(tobs,phot,nfile,x0[15],texp,x0[8:14],x0[1:7],u_ch1_start,u_ch2_start,ninc,ip_inc,ind_inc,tt1)
  #  clf(); plot(phot ./ phot_vary,"."); plot(phot_mod); plot(phot ./ (phot_vary .* phot_mod), "o")
 
  # Now compute the standard deviation of residuals for each time series: 
  sig_file = zeros(nfile)
  resid = phot .- (phot_vary .* phot_mod)
  for ifile=1:nfile
      i10 = i1[ifile]; i20 = i2[ifile]
      sig_file[ifile] = std(resid[i10:i20])
  end

  # Now compute a light curve, and run a Markov chain:
  
  # using  AdvancedHMC
  
  function compute_chain(tobs,phot,nfile,texp,ninc,ip_inc,ind_inc,tt1,x0,srp,
          n,sig_file,numdims,numwalkers,thinning,numsamples_perwalker,burnin,astep)
      tstart = time()
      # Keep track of how many times the likelihood function is called:
      nlike = 0
      nreject = 0
      # Creates a closure for sampling from the photometric model:
      np = n-1

      # Set up function to compute the log prior:
      function log_prior(x::Array{T,1}) where {T <: Real}
        # Returns log(prior) and its derivative:
        rp = x[1:np]; b = x[np+1:2*np]; rho_T1=x[2np+1]; q = x[2np+2:2np+5]
        # Check that impact parameters are in bounds:
        if maximum(abs.(b)) >= 1 || (!incprior && minimum(b) <= 0)
           return -Inf::T
        end
        # Check that limb-darkening parameters are in bounds:
        if minimum(q) <= 0 || maximum(q) >= 1 
          return -Inf::T
        end 
        # Check that radius ratios are in bounds:
        if minimum(rp) <= 0 || maximum(rp) > 0.2
          return -Inf::T
        end
        # Check that density is in bounds:
        if rho_T1 <= 0 || rho_T1 >= 100
          return -Inf::T
        end
        # Check that contamination parameter is in bounds:
        if contamination
          if abs(x[2np+6]) >= 1
            return -Inf::T
          end
        end
        # Add in prior on the limb-darkening parameters:
        logp = 0.0
        if limbprior
          u_ch1 = [2*sqrt(q[1])*q[2],sqrt(q[1])*(1 -2q[2])]
          u_ch2 = [2*sqrt(q[3])*q[4],sqrt(q[3])*(1 -2q[4])]
          # Priors from Elsa's analysis:
          logp -= 0.5*((u_ch1[1]-0.1633)/0.0364)^2
          logp -= 0.5*((u_ch1[2]-0.2549)/0.0570)^2
          logp -= 0.5*((u_ch2[1]-0.1442)/0.0324)^2
          logp -= 0.5*((u_ch2[2]-0.2173)/0.0482)^2
          #return logp::T
        end
        # Add in prior on the inclinations:
        if incprior
          # Now compute prior - sigtheta is the last parameter:
          sigtheta = x[end]
          # Check if sigtheta is in bounds:
          if sigtheta <= 0.0 || sigtheta >= 5.0  # Safe bet than inclination scatter is < 5 degrees
            return -Inf::T
          end
          # Compute a/R_* from density, period (ignore small uncertainty on orbital period):
          aonr = rho_T1^(1//3) .* (aonr0 * P.^(2//3))
          # Compute inclinations of each planet, ignore small eccentricity:
          inc_planet = acos.(b ./ aonr) .* (180/pi)  # in degrees
          # Compute mean inclination:
          inc_mean = sum(inc_planet)/np
          # Now, compute prior:
          for ip=1:np
            logp -= 0.5*((inc_planet[ip]-inc_mean)/sigtheta)^2
          end
          # Add in integral normalization
          logp -= np*log(sigtheta)
        end
        return logp::T 
      end
          
      # Convert to Kipping's parameterization:
      #q_ch1 = [(u_ch1[1]+u_ch1[2])^2,0.5*u_ch1[1]/(u_ch1[1]+u_ch1[2])]
      #q_ch2 = [(u_ch2[1]+u_ch2[2])^2,0.5*u_ch2[1]/(u_ch2[1]+u_ch2[2])]
      
      # For now we are holding the dynamical model fixed, and simply varying the model parameters.
      #x0 = permutedims([rp;b;[rho_T1];q_ch1;q_ch2])  # Radius-ratios, impact paramters, density, limb-darkening
    
      function generate_walkers(x0::Array{T,1},numwalkers::Int64) where {T <: Real}
          # Generates an array of walkers:
          xstart = zeros(size(x0)[1],numwalkers)
          println("size of xstart: ",size(xstart))
          llike_start = zeros(numwalkers)
          for i=1:numwalkers
              llike_start[i] = -Inf
              while llike_start[i] == -Inf
                # rp drawn from Delrez et al. uncertainties:
                xstart[1:np,i] .= abs.(x0[1:np] .+ randn(np) .* srp)
                # b drawn from Gaussian of width 0.1:
                xstart[np+1:2np,i] .= abs.(x0[np+1:2np] .+ 0.1 .*rand(np))
                # rho drawn from a Gaussian:
                xstart[2np+1,i] = x0[2np+1] + 1.2 *randn()
                # q's drawn with some scatter about optimum:
                xstart[2np+2:2np+5,i] .= abs.(x0[16:19] .+ 0.25 .* randn(4))
                # Draw contamination from N(0,0.05):
                if contamination
                  xstart[2np+6,i] = 0.05*randn()
                end
                # Draw sigtheta from 0.01 to 1.0 degrees:
                if incprior
                  xstart[end,i] = 10^(rand()*2-2)
                end
                # Compute likelihood:
                llike_start[i] = compute_loglike(xstart[:,i])
              end
          end
          # That's it!
      return xstart, llike_start
      end
  
      # Exposure time in days:
      texp_day = texp/(60*24)
      # Create a vector to hold the photometric model...
      phot_mod = copy(phot)
      # Small-planet approximation:
      #phot_small = copy(phot)
      # Create a vector for a model of stellar/instrumental variabilty:
      phot_vary = copy(phot)
  
   #   function compute_model(time,phot,nfile,rho_T1,texp,b,rp,u_ch1,u_ch2,ninc,ip_inc,ind_inc,tt1)
      function compute_loglike(x::Array{T,1}) where {T <: Real}
        # Create a photometric model with detrending.
        # Check that we are not out of bounds:
        lprior  = log_prior(x)
        if lprior == (-Inf)
           nreject += 1
           return lprior
        end
        # Unpack the parameters:
        rp = x[1:np]; b = x[np+1:2*np]; rho_T1,q1_ch1,q2_ch1,q1_ch2,q2_ch2 = x[2np+1:2np+5]
        if contamination
          epsilon = x[2np+6]
        end
        
        # Density of T1 star in units of solar density for conversion
        # of N-body code units to Trappist-1 units:
        fac = rho_T1^(1//3)*AU/RSUN
          
        # Fill photometric model with zeros:
        fill!(phot_mod,0.0)
        # ... and detrending model:
        fill!(phot_vary,0.0)
  
        # Two limb-darkening coefficients:
        nu = 2
        # Convert from Kipping parameterization back to u_1, u_2:
        u_ch1_curr = [2*sqrt(q1_ch1)*q2_ch1,sqrt(q1_ch1)*(1 -2q2_ch1)]
        u_ch2_curr = [2*sqrt(q1_ch2)*q2_ch2,sqrt(q1_ch2)*(1 -2q2_ch2)]
      
        itrans = 0
        # Initialize transit structure (no derivatives for now):
        #trans =  Limbdark.transit_init(0.1,0.5,[0.2,0.2],false)
        # Loop over each light curve:
        for ifile=1:nfile
          # An array of the current Spitzer observation times being modeled:
          tcurr = tobs[i1[ifile]:i2[ifile]]
          # First set the limb-darkening parameters
          # for the current light curve    
          if channel[ifile] == "c1"
            #trans =  Limbdark.transit_init(0.1,0.5,u_ch1_curr,false)
            trans =  transit_init(0.1,0.5,u_ch1_curr,false)
            #trans.u_n .= u_ch1_curr
            #trans.u_n[1] = u_ch1_curr[1]
            #trans.u_n[2] = u_ch1_curr[2]
          else
            #trans =  Limbdark.transit_init(0.1,0.5,u_ch2_curr,false)
            trans =  transit_init(0.1,0.5,u_ch2_curr,false)
            #trans.u_n .= u_ch2_curr
            #trans.u_n[1] = u_ch2_curr[1]
            #trans.u_n[2] = u_ch2_curr[2]
          end
          # Number of exposure times for this case:
          nt = i2[ifile]-i1[ifile]+1
          # Set up array to hold the time-integrated flux for this light curve:
          #favg = zeros(5+nu,nt)
          favg = ones(5+nu,nt)
          # Small-planet approximation:
          #fsmall = ones(nt)
          # Set up array to hold the time-integrated flux for each planet:
          favg1 = zeros(5+nu,nt)
          # Set up arrays to hold statistics on time integrations:
          neval1 = zeros(Int64,nt)
          depthmax1 = zeros(Int64,nt)
          # Loop over each transit:
          for it = 1:ninc[ifile]
            fill!(favg1,0.0)
            itrans += 1
            # ttbv holds the transit time, sky-velocity, and impact parameter from N-body model:
            ttbv = tt1[:,ip_inc[itrans]+1,ind_inc[itrans]]
            # Now, create a photometric model.  Convert impact parameter and
            # sky velocity to units of radius of the star (although we're not
            # using impact-parameter variations at the moment):
            # param = [ttbv[1],ttbv[2]*fac,ttbv[3]*fac]
            # Set the radius ratio and impact parameter to the current planet:
            param = [ttbv[1],ttbv[2]*fac,b[ip_inc[itrans]]]
            # Set radius in transit model structure to radius of current planet:
            trans.r = rp[ip_inc[itrans]]
            # Evaluate the adaptively time-integrated quadratic limb-darkened transit model:
            #Limbdark.integrate_lightcurve!(trans,param,tcurr,texp_day,favg1,nt,1e-11,40,neval1,depthmax1)
            integrate_lightcurve!(trans,param,tcurr,texp_day,favg1,nt,1e-11,40,neval1,depthmax1)
            # Small planet approximation:
            #boft2 = ((tcurr .- param[1]) .*param[2]).^2 .+ param[3]^2; mu = sqrt.(1 .- boft2[boft2 .<= 1])
            #fsmall[boft2 .<= 1] .-= trans.r^2 .* (1 .- trans.u_n[1] .* (1 .- mu) - trans.u_n[2] .* (1 .- mu).^2) ./ 
            #        (1 - trans.u_n[1]/3 - trans.u_n[2]/6)
            # Add to the current model:
            favg .+= favg1
            # Subtract off 1 for additional planets
            #if it > 1
            #  favg[1,:] .-= 1
            #end
          end
          # Now, remove a polynomial:
          if contamination
            phot_mod[i1[ifile]:i2[ifile]] .= favg[1,:] .* (1-epsilon) .+ epsilon
            #phot_small[i1[ifile]:i2[ifile]] .= fsmall .* (1-epsilon) .+ epsilon
          else
            phot_mod[i1[ifile]:i2[ifile]] .= favg[1,:]
            #phot_small[i1[ifile]:i2[ifile]] .= fsmall
          end
#          if channel[ifile] == "c1"
#            println("Limb-darkening: ",trans.u_n," ",[2sqrt(x[2np+2])*x[2np+3],sqrt(x[2np+2])*(1-2x[2np+3])])
#            clf(); plot(tcurr,phot_mod[i1[ifile]:i2[ifile]]); plot(tcurr,phot_small[i1[ifile]:i2[ifile]])
#            read(stdin,Char)
#          else
#            println("Limb-darkening: ",trans.u_n," ",[2sqrt(x[2np+4])*x[2np+5],sqrt(x[2np+4])*(1-2x[2np+5])])
#            clf(); plot(tcurr,phot_mod[i1[ifile]:i2[ifile]]); plot(tcurr,phot_small[i1[ifile]:i2[ifile]])
#            read(stdin,Char)
#          end
          phot_vary[i1[ifile]:i2[ifile]] .= fit_poly(tcurr,phot[i1[ifile]:i2[ifile]]./favg[1,:],3)
        end
        # Compute the chi_square:
        chi_square = 0.0
        resid = phot .- (phot_vary .* phot_mod)
        for ifile=1:nfile
          i10 = i1[ifile]; i20 = i2[ifile]
          chi_square += sum(resid[i10:i20].^2) / sig_file[ifile]^2
        end
        # Return log posterior = log(prior) + log(likelihood):
        llike = -0.5*chi_square + lprior
        nlike += 1
        if mod(nlike,nout) == 0
          println(nreject," ",nlike," log like: ",llike," time [hr]: ",(time()-tstart)/3600)
        end
        return llike
      end
  
      # Now sample.  Following example on AffineInvariantMCMC page:
      xstart, llstart = generate_walkers(x0,numwalkers)
      println("llstart: ",llstart)
      chain, llhoodvals = AffineInvariantMCMC.sample(compute_loglike,numwalkers,xstart,burnin,1,astep)
      chain, llhoodvals = AffineInvariantMCMC.sample(compute_loglike,numwalkers,chain[:, :, end], numsamples_perwalker, thinning,astep)
      flatchain, flatllhoodvals = AffineInvariantMCMC.flattenmcmcarray(chain, llhoodvals)
      return chain,llhoodvals,flatchain,flatllhoodvals
  end 

  # Start with optimized values:
#  x0 = [ 0.08510408205855542, 0.08351749544695007, 0.060057015907927426 , 0.07003499168340532, 0.07938534728488905, 0.08571785728160501, 0.05767832213757894, 0.06777531110636684, 0.07268527380271149, 0.0003227708869417413, 0.17817318017323253, 0.2926207489330042, 0.36533945405158097, 0.3735584459205608, 52.46814828992533, 0.0961562252101869, 0.5059574689387107, 0.30961490912899814, 0.29644507311377927]
  x0 =[0.08511305133271845, 0.08352547037122816, 0.06005344100378423, 0.0700421002759512, 0.07962187102240074, 0.08611264006593941, 0.05767381066632817, 0.07876438173347886, 0.08351738136283536, 0.0008769148946956045, 0.18225051020178037, 0.30070386193030363, 0.36817942906674744, 0.3756424617589329, 52.34377029283971, 0.11235270319764341, 0.42037661035916857, 0.352424321959808, 0.2864053200404355]
  numdims = 19
  if contamination
    # Add in contamination parameter:
    numdims += 1
    push!(x0,0.0)
  end
  if incprior
    numdims += 1
    # Add in sigtheta parameter:
    push!(x0,0.5)
  end
  # numwalkers = 50
  # thinning = 10
  # This is going to run for ~2 days:
  # numsamples_perwalker = 500
  # burnin = 10
  # Take smaller steps to increase the acceptance rate:
  # astep = 1.2
  chain,llhoodvals,flatchain,flatllhoodvals = 
  compute_chain(tobs,phot,nfile,texp,ninc,ip_inc,ind_inc,tt1,x0,srp,n,
      sig_file,numdims,numwalkers,thinning,numsamples_perwalker,burnin,astep)
  @save foutput chain llhoodvals flatchain flatllhoodvals tobs phot nfile texp ninc ip_inc ind_inc tt1 x0 n sig_file i1 i2 thinning numwalkers numdims numsamples_perwalker burnin astep limbprior incprior contamination

  # Now save to output file
return
end

#= Example:
 datafile = "T1_Spitzer_data.jld2"; foutput = "T1_pd_MCMC_v01.jld2"; numwalkers = 50; burnin = 10; thinning = 1; astep = 2.0; nsteps = 50
 run_pd(datafile,foutput,numwalkers,burnin,thinning,astep,nsteps;nout=100)
=#
