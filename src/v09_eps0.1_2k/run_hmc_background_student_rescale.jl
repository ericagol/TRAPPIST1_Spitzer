

# Carries out a fit to TRAPPIST-1 data from December 2019
# Optimizes, runs HMC.  No prior on longitude difference.
# Eccentricity vectors used without the 1/e prior for now.
# Student's-t distribution for likelihood.
# Uses full data, a larger timestep.

use_student = true

# Carry out a numerical check of the gradient once (NOTE: this is slow!)
test_gradient = false
if test_gradient
  include("../compute_grad_num.jl")
end
include("../CGS.jl")
include("log_students_prob.jl")
include("../extract_planet.jl")
include("../nlog_prior.jl")
include("../loglinspace.jl")
include("../regress.jl")
# Use NbodyGradient code:
include("../../src/NbodyGradient/ttv.jl")

using Main.CGS
#using PyPlot
using ForwardDiff
using DelimitedFiles
using Printf
using SpecialFunctions
using LinearAlgebra
using JLD2

# Uses e\cos{\omega} & e\sin{\omega} as parameters for Markov chain.

# Read in elements_grid from prior optimization that showed a jump
# in the minimum:

function run_hmc(fname::String,felements::String,foutput::String,nstep::Int64,epsilon0::Float64,nleap0::Int64,showplots::Bool)
# Number of bodies:
n = 8
nplanet = n-1

# Prior on masses:
sigm2 = 1e-8^2

# Number of grid points for each mass - this will give ~0.12-sigma spacing:
ngrid = 1
# Mass grid in standard deviations:
sig_grid = zeros(ngrid);
if ngrid > 1
  sig_grid .= linearspace(-3.0,3.0,ngrid)
end
mass_grid = copy(sig_grid)
# Use Simon's values as a guide; convert these to mass-ratio:
sgrimm = [0.082,0.073,0.022,0.027,0.026,0.027,0.034]*MEARTH/(0.09*MSUN)



t0 = 7257.93115525
h  = 0.06
tmax = 1600.0

# Read in initial conditions from the first optimization:
#elements_minimum = readdlm("elements_dompi.txt",',')
#elements_minimum = readdlm("elements_noprior_mixture_rescale.txt",',')
#elements_minimum = readdlm("elements_noprior_students.txt",',')
elements_minimum = readdlm(felements,',')
elements_minimum = elements_minimum[1:n,:]
elements = copy(elements_minimum)
elements_trial = copy(elements)

# Make an array, tt,  to hold transit times:
# First, though, make sure it is large enough:
ntt = zeros(Int64,n);
for i=2:n
  ntt[i] = ceil(Int64,tmax/elements[i,2])+100
end
#println("ntt: ",ntt)
tt1 = zeros(n,maximum(ntt));
tt2 = zeros(n,maximum(ntt));
dtdq0 = zeros(n,maximum(ntt),7,n);
dtdelements = zeros(n,maximum(ntt),7,n);
# Save a counter for the actual number of transit times of each planet:
count1 = zeros(Int64,n);
# Set the "size" of the star (only checks that close to transit):
rstar = 1e11
# Call the ttv function:
@time dq=ttv_elements!(n,t0,h,tmax,elements,tt1,count1,0.0,0,0,rstar)
@time ttv_elements!(n,t0,h,tmax,elements,tt1,count1,dtdq0,rstar)

# Now call with half the timestep:
count2 = zeros(Int64,n)

# Set the initial outlier parameters: ndof =  number of degrees of freedom (\nu); 
#    lnV1 = log ratio by which variance increases

#ndof = 3.27; lnV1 =  -0.36

ndof = 3.36016234786184600
lnV1 = -0.34445676197382336

# Convert to new parameterization:
lgndof = log(ndof)
V1exp2nuinv = exp(lnV1+0.5/ndof)


# Now, try to improve the fit to the data:
# First, find the indices that match for each planet:
iplanet = zeros(Int64,0)
indx = zeros(Int64,0)
tobs_tot = zeros(Float64,0)
sobs_tot = zeros(Float64,0)
ntrans = 0
chi_ttv = 0.0
chi_ttv_error = 0.0
chi_ttv += ntrans*log_norm_students_t([lgndof,V1exp2nuinv])
data = readdlm(fname,',')
for ip=1:nplanet
  eobs,tobs,sobs,nobs = extract_planet(data,ip)
  ittv = 1
  for j=1:nobs
    ittv = 1; dt = Inf
    for jttv=1:count1[ip+1]
      if abs(tt1[ip+1,jttv]-tobs[j]) < dt
        ittv = jttv
        dt = abs(tt1[ip+1,jttv]-tobs[j])
      end
    end
    push!(iplanet,ip+1)
    push!(indx,ittv)
    push!(tobs_tot,tobs[j])
    push!(sobs_tot,sobs[j])
    chi_ttv += lnprob([tobs[j],tt1[ip+1,ittv],sobs[j],lgndof,V1exp2nuinv])
    ntrans +=1
  end
end
# Now add normalization:
println("Total number of transits: ",ntrans," Chi-square: ",chi_ttv," error: ",chi_ttv_error)

# Now, optimize the model:

# First, vary just the period & the t0 values for all planets:
param_vary = [1,2,3,4,5]
nvary = length(param_vary)
dx = [1e-10,1e-5,1e-5,1e-5,1e-5]
#nparam = length(param_vary) * (n-1) # Exclude the star
if use_student
  nparam = length(param_vary) * (n-1) + 2 # Exclude the star, include outlier parameters
else
  nparam = length(param_vary) * (n-1)  # Exclude the star, exclude outlier/mixture parameters
end

# Now, try to implement a  Levenburg-Marquardt algorithm:
dchi = Inf
# Tolerance with which to converge to:
tol  = 3e-5
iter = 0
itmax = 100
hessian = zeros(Float64,nparam,nparam)
#lambda = [1e3,316.,100.,31.6,10.0,3.16,1.0,3.16e-1,1e-1,3.16e-2,1e-2,3.16e-3,1e-3,3.16e-4,1e-4,3.16e-5,1e-5,3.16e-6,1e-6]
#lambda = [1e3,316.,100.,31.6,10.0,3.16,1.0,3.16e-1,1e-1,1e-2,3.16e-3,1e-3,3.16e-4,1e-4]
#lambda = [1e4,1e3,100.,10.0,1.0,0.1,.01,.001,1e-4,1e-5,1e-6]
#lambda = [1e3,1e2,10.,1.0,0.1,.01,.001,1e-4,1e-5,1e-6,1e-7,0.]
#nlambda = 12
lambda = [1e3,1e2,10.,1.0,0.1,.01,.001,1e-4,0.]
nlambda = 9
#lambda = [1.0]
#nlambda = 1
#lambda = 5e-4
# Gradient contains derivatives with respect to mass, period, t0, e*cos(omega), e*sin(omega) for each planet:
gradf = zeros(Float64,nparam)
# Numerical gradient check:
gradf_num = zeros(Float64,nparam)


# Save gradient to compare with prior step when new chi-square is rejected:
# Set up a grids for looping over masses of each planet (there are nplanet=n-1 of these):
chi_grid = zeros(Float64,nplanet,ngrid)
prior_grid = zeros(Float64,nplanet,ngrid)
elements_grid = zeros(Float64,nplanet,n,7,ngrid)
lgndof_grid = zeros(Float64,nplanet,ngrid)
V1exp2nuinv_grid = zeros(Float64,nplanet,ngrid)
tt_grid = zeros(nplanet,n,maximum(ntt),ngrid)
#jgrid=[11,10,9,8,7,6,5,4,3,2,1,12,13,14,15,16,17,18,19,20,21]
#jgrid=[16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
jgrid=[1]
# Loop over the mass-ratios of each planet:
#for imass = 1:nplanet
# We're interested in the optimum for now, so just "iterate" over one planet:
for imass = 1:1
# Loop over the grid points for each planet:
for igrid in jgrid
  start = time()
  if igrid == jgrid[1]
    elements = copy(elements_minimum)
    # First, optimize fit, so no constraint upon mass:
    sigm2 = Inf
    if imass == 1
      compute_grad_num_q = true
    else
      compute_grad_num_q = false
    end
  elseif igrid == jgrid[1]+1
    elements = copy(elements_grid[imass,:,:,jgrid[1]])
    compute_grad_num_q = false
  end
  # Compute chi-square
  dq0 = ttv_elements!(n,t0,h,tmax,elements,tt2,count2,0.0,0,0,rstar)
  # Add in prior:
  prior_ttv = (elements[1+imass,1]-mass_grid[igrid])^2/sigm2
  lgndof_best = lgndof; V1exp2nuinv_best = V1exp2nuinv
  chi_ttv = prior_ttv
  for j=1:ntrans
    chi_ttv += lnprob([tobs_tot[j],tt2[iplanet[j],indx[j]],sobs_tot[j],lgndof_best,V1exp2nuinv_best])
  end
  dchi = Inf; iter = 0
  while abs(dchi) > tol && iter < itmax
    tstart = time()
    dtdelements = ttv_elements!(n,t0,h,tmax,elements,tt1,count1,dtdq0,rstar)
    # Compute the "chi-square", gradient & Hessian:
    fill!(hessian,0.0)
    fill!(gradf,0.0)
    i = 0; chi_ttv = 0.0
  # Iterate over planets & varying parameters:
    for ip=1:nplanet, ivary=1:nvary
      i +=1
      # If varying mass, this is at end of dtdq0 parameter list:
      if ivary == 1; iv = 7; else; iv= ivary-1; end
      # Gradient of "chi-square" with respect to each varying parameter:
      for it=1:ntrans
        # Set up parameter vector for computing likelihood of this data point:
        p = [tobs_tot[it],tt1[iplanet[it],indx[it]],sobs_tot[it],lgndof,V1exp2nuinv]
        if i == 1
          # Compute the chi-square (but only once):
          chi_ttv += lnprob(p)
        end
        # Compute gradient of -2*ln(p_i) with respect to these parameters:
        grad_tmp = grad_point(p)
        # Take derivative of -2*log likelihood with respect to model at this point, and then compute derivative of model with
        # respect to model parameters:
        gradf[i] += grad_tmp[2]*dtdelements[iplanet[it],indx[it],iv,ip+1]
        if use_student
          if i == 1
          # Derivative wrt ndof:
            gradf[nparam-1] += grad_tmp[4]
          # Derivative wrt V1exp2nuinv:
            gradf[nparam  ] += grad_tmp[5]
          end
        end
  #      gradf[i] -= dtdelements[iplanet[it],indx[it],iv,ip+1]*(tobs_tot[it]-tt1[iplanet[it],indx[it]])/sobs_tot[it]^2
      end
      # Now compute the Hessian:
      j = 0
      for jp=1:nplanet, jvary=1:nvary
        j +=1
        if jvary == 1; jv = 7; else; jv= jvary-1; end
        for it=1:ntrans
          # Set up parameter vector for computing likelihood of this data point:
          p = [tobs_tot[it],tt1[iplanet[it],indx[it]],sobs_tot[it],lgndof,V1exp2nuinv]
          # Compute hessian of -2ln(p_i) with respect to these parameters:
          hess_tmp = hess_point(p)
          hessian[i,j] += hess_tmp[2,2]*dtdelements[iplanet[it],indx[it],iv,ip+1]*dtdelements[iplanet[it],indx[it],jv,jp+1]
          if use_student
            if j == 1
              hessian[i,nparam-1] += hess_tmp[2,4]*dtdelements[iplanet[it],indx[it],iv,ip+1]
              hessian[i,nparam  ] += hess_tmp[2,5]*dtdelements[iplanet[it],indx[it],iv,ip+1]
            end
            if i == 1
              hessian[nparam-1,j] += hess_tmp[4,2]*dtdelements[iplanet[it],indx[it],jv,jp+1]
              hessian[nparam  ,j] += hess_tmp[5,2]*dtdelements[iplanet[it],indx[it],jv,jp+1]
            end
            if i == 1 && j == 1
              hessian[nparam-1,nparam-1] += hess_tmp[4,4]
              hessian[nparam-1,nparam  ] += hess_tmp[4,5]
              hessian[nparam  ,nparam-1] += hess_tmp[5,4]
              hessian[nparam  ,nparam  ] += hess_tmp[5,5]
            end
          end
  #        hessian[i,j] += dtdelements[iplanet[it],indx[it],iv,ip+1]*dtdelements[iplanet[it],indx[it],jv,jp+1]/sobs_tot[it]^2
        end
      end
    end
    # Add in grid constraint upon mass:
    chi_ttv += (elements[imass+1,1]-mass_grid[igrid])^2/sigm2
    # Add in normalization:
    chi_ttv += ntrans*log_norm_students_t([lgndof,V1exp2nuinv])
    # Add in prior gradient:
    gradf[(imass-1)*5+1] += 2.0*(elements[imass+1,1]-mass_grid[igrid])/sigm2
    # Add in gradient with respect to (lgndof,V1exp2nuinv):
    gradf[nparam-1:nparam] .+= ntrans*grad_norm([lgndof,V1exp2nuinv])
    # Compare with numerical gradient:
    if compute_grad_num_q && iter == 0 && test_gradient == true
       gradf_num = convert(Array{Float64,1},
        compute_grad_num(n,h,t0,tmax,elements,tt2,count2,rstar,tobs_tot,sobs_tot,lgndof,V1exp2nuinv,indx,nparam,nplanet,nvary))
      println("grad: ",gradf," gradf_num: ",gradf_num," diff: ",gradf-gradf_num)
    end
    # Add in Hessian with respect to prior:
    hessian[(imass-1)*5+1,(imass-1)*5+1] += 2.0/sigm2
    # Add in hessian of normalization:
    hessian[nparam-1:nparam,nparam-1:nparam] .+= ntrans*hess_norm([lgndof,V1exp2nuinv])
    # Now, carry out Newton step (can no longer use Levenberg-Marquardt!):
    chi_best = chi_ttv; elements_best = copy(elements); lgndof_best = copy(lgndof); V1exp2nuinv_best=copy(V1exp2nuinv);  prior_best = copy(prior_ttv)
    elements_trial = copy(elements); lgndof_trial = copy(lgndof); V1exp2nuinv_trial = copy(V1exp2nuinv); prior_trial = copy(prior_ttv)
    for ilambda = 1:nlambda
    # Update the parameters:
      mat = copy(hessian)
      for j=1:nparam
        mat[j,j] += lambda[ilambda]*hessian[j,j]
      end
      # Newton step:
      deltax = -\(mat,gradf)
      # deltax may take some parameters out of bounds - in which case
      # we'll shrink it until inbounds:
      out_of_bounds = true; factor = 1.0;  ecc2_trial = zeros(n-1)
      elements_trial = copy(elements); lgndof_trial = copy(lgndof); V1exp2nuinv_trial = copy(V1exp2nuinv)
      while out_of_bounds && factor > 1e-3
        elements_trial = copy(elements); lgndof_trial = copy(lgndof); V1exp2nuinv_trial = copy(V1exp2nuinv)
      # Loop over planets & varied parameters:
        i = 0
        for ip=2:n
          for ivary=1:nvary
            i += 1
            elements_trial[ip,param_vary[ivary]] += factor*deltax[i]
          end
          ecc2_trial[ip-1] = elements_trial[ip,4]^2+elements_trial[ip,5]^2
        end
        if use_student
          lgndof_trial+= factor*deltax[nparam-1]
          V1exp2nuinv_trial += factor*deltax[nparam]
        end
        if maximum(ecc2_trial) < 0.2 && log(V1exp2nuinv_trial) > -10.0 && exp(lgndof_trial) > 0.5 && exp(lgndof_trial) < 100.0
          out_of_bounds = false
        else
          # Since we're out of bounds, shrink deltax:
          factor *= 0.9
        end
      end
    # Compute chi-square
      if maximum(ecc2_trial) < 0.2 && log(V1exp2nuinv_trial) > -10.0 && exp(lgndof_trial) > 1.0 && exp(lgndof_trial) < 100.0
        chi_trial = Inf
        while chi_trial > chi_ttv && factor > 1e-2
          # Now, do a line search:
          dq0 = ttv_elements!(n,t0,h,tmax,elements_trial,tt2,count2,0.0,0,0,rstar)
          chi_trial = 0.0
          for j=1:ntrans
            chi_trial += lnprob([tobs_tot[j],tt2[iplanet[j],indx[j]],sobs_tot[j],lgndof_trial,V1exp2nuinv_trial])
          end
          # Add in normalization:
          chi_trial += ntrans*log_norm_students_t([lgndof_trial,V1exp2nuinv_trial])
          # Add in prior:
          prior_trial = (elements_trial[imass+1,1]-mass_grid[igrid])^2/sigm2
          chi_trial += prior_trial
          if chi_trial > chi_ttv
            factor *= 0.5
            elements_trial = copy(elements); lgndof_trial = copy(lgndof); V1exp2nuinv_trial = copy(V1exp2nuinv)
      # Loop over planets & varied parameters:
            i = 0
            for ip=2:n, ivary=1:nvary
              i += 1
              elements_trial[ip,param_vary[ivary]] += factor*deltax[i]
            end
            if use_student
              lgndof_trial += factor*deltax[nparam-1]
              V1exp2nuinv_trial += factor*deltax[nparam]
            end
          end
        end
      else
        chi_trial = Inf
      end
      if chi_trial < chi_ttv && chi_trial < chi_best
        println("New chi2: ",chi_trial," ",chi_trial," prior_trial: ",prior_trial,
           " target_mass: ",mass_grid[imass]," elapsed time: ",(time()-start)/60," min")
      # A better chi-square is found - update elements, transit times & chi-square:
        dchi = chi_ttv - chi_trial
        elements_best .= elements_trial 
        if use_student
          lgndof_best = copy(lgndof_trial)
          V1exp2nuinv_best = copy(V1exp2nuinv_trial)
        end
        tt1 .= tt2
        chi_best = copy(chi_trial)
        prior_best = copy(prior_trial)
      else
        println("Chi-square worse. ",chi_trial," best: ",chi_best," chi_ttv: ",chi_ttv)
      end
    end
    # Don't include prior:
    if chi_best < chi_ttv
      telapse = time()-tstart
      println("New chi2: ",chi_best," prior: ",prior_best," target m: ",mass_grid[igrid]," elapsed time: ",telapse/60," min")
      chi_grid[imass,igrid] = chi_best - prior_best
      if use_student
        lgndof_grid[imass,igrid] = lgndof_best
        V1exp2nuinv_grid[imass,igrid] = V1exp2nuinv_best
      end
      elements_grid[imass,:,:,igrid] .= elements_best
      tt_grid[imass,:,:,igrid] .= tt1
      prior_grid[imass,igrid] = prior_best
      chi_ttv = copy(chi_best); prior_ttv=copy(prior_best); elements .= elements_best
      if use_student
        lgndof = copy(lgndof_best); V1exp2nuinv = copy(V1exp2nuinv_best)
      end
      iter += 1
    else
#      println("Can't improve chi-square. elapsed time: ",toq()/60.," min")
      # Save the initial parameters:
      chi_grid[imass,igrid] = chi_best - prior_best
      if use_student
        lgndof_grid[imass,igrid] = lgndof_best
        V1exp2nuinv_grid[imass,igrid] = V1exp2nuinv_best
      end
      elements_grid[imass,:,:,igrid] .= elements_best
      tt_grid[imass,:,:,igrid] .= tt1
      prior_grid[imass,igrid] = prior_best
      chi_ttv = copy(chi_best); prior_ttv=copy(prior_best); elements .= elements_best
      if use_student
        lgndof = copy(lgndof_best); V1exp2nuinv = copy(V1exp2nuinv_best)
      end
      iter = itmax
    end
  end
  println("Planet ",imass," Finished grid point: ",igrid," ",elements," ",lgndof," ",V1exp2nuinv,
          " ",chi_ttv-prior_ttv," ",prior_ttv," Elapsed time: ",time()-start)
  # For first grid point, we optimized chi-square.  Now, we want to add grid constraints.
  if igrid == jgrid[1]
    # Now, create mass grid, spanning +-3-sigma of Simon Grimm's uncertainty:
    mass_grid = sig_grid.*sgrimm[imass]
    # Offset grid about optimum:
    mass_grid .+= elements_grid[imass,imass+1,1,jgrid[1]]
    # Now, cause grid to be reinstated:
    sigm2 = 1e-8^2
    itmax = 10
  end
end
end

cov_save = inv(0.5*hessian)

hess_save = 0.5*hessian

planet = ["b","c","d","e","f","g","h"]
sigm = zeros(nplanet)
mass = zeros(nplanet)
ecc = zeros(nplanet)
pomega = zeros(nplanet)
for i=1:nplanet
  mass[i] = elements[i+1,1]*.09*MSUN/MEARTH
  sigm[i] = sqrt(cov_save[5*(i-1)+1,5*(i-1)+1])*.09*MSUN/MEARTH
  ecc[i] = sqrt(elements[i+1,4]^2+elements[i+1,5]^2)
  pomega[i] = atan(elements[i+1,5],elements[i+1,4])*180/pi
  println(planet[i]," ",@sprintf("%6.4f",mass[i]),"+-",@sprintf("%6.4f",sigm[i])," ",
   @sprintf("%6.4f",elements[i+1,4]),"+-",@sprintf("%6.4f",sqrt(cov_save[5*(i-1)+4,5*(i-1)+4])),
   " ",@sprintf("%6.4f",elements[i+1,5]),"+-",@sprintf("%6.4f",sqrt(cov_save[5*(i-1)+5,5*(i-1)+5])),
   " ",@sprintf("%6.4f",ecc[i])," ",@sprintf("%7.2f",pomega[i]))
end

# Print deviations of each transit:
chi_planet = zeros(Float64,n-1)
nobs_planet = zeros(Int64,n-1)
dev = zeros(Float64,ntrans)
elements = elements_grid[1,:,:,jgrid[1]]
#        dq0 = ttv_elements!(n,t0,h,tmax,elements,tt2,count2,0.0,0,0,rstar)
dtdelements = ttv_elements!(n,t0,h,tmax,elements,tt1,count1,dtdq0,rstar)
for j=1:ntrans
  dev[j] = (tobs_tot[j]-tt1[iplanet[j],indx[j]])/sobs_tot[j]
  chi_planet[iplanet[j]-1] += dev[j]^2
  nobs_planet[iplanet[j]-1] += 1
  if abs(dev[j]) > 3
    println(j," ",iplanet[j]," ",indx[j]," ",tt1[iplanet[j],indx[j]]," ",tobs_tot[j]," ",sobs_tot[j]," ",dev[j]," ",dev[j]^2)
  end
end


hinv = 0.5*(inv(hess_save)+inv(hess_save)')
cholhinv = cholesky(hinv)

x = zeros(Float64,nparam)
# Insert elements into parameter vector:
for i=2:n, j=1:5
  x[(i-2)*5+j]= elements[i,j]
end
# Select values of mixture model at optimum parameters:
lgndof = lgndof_grid[1,jgrid[1]]
V1exp2nuinv = V1exp2nuinv_grid[1,jgrid[1]]
# Add in mixture model parameters to parameter vector:
if use_student
  # Convert to log(ndof) and V1*exp(1/(2ndof)):
  x[nparam-1] = lgndof
  x[nparam  ] = V1exp2nuinv
end

# Compute the negative log likelihood:
function f0(x)
# computes negative log likelihood:
for i=2:n, j=1:5
  elements_trial[i,j] = x[(i-2)*5+j]
end
dq0 = ttv_elements!(n,t0,h,tmax,elements_trial,tt2,count2,0.0,0,0,rstar)
# Ignoring priors:
chi_trial = 0.0
lgndof = x[nparam-1]; V1exp2nuinv = x[nparam]
for j=1:ntrans
  chi_trial += lnprob([tobs_tot[j],tt2[iplanet[j],indx[j]],sobs_tot[j],lgndof,V1exp2nuinv])
end
# Add in normalization:
chi_trial += ntrans*log_norm_students_t([lgndof,V1exp2nuinv])
# Convert from -2log(likelihood) to -log(likelihood)
chi_trial *= 0.5
# Add in prior:
chi_trial += nlog_prior(x)
return chi_trial
end

# Check that this routine works:
chi_opt = chi_grid[1,jgrid[1]]
println("f0(x): ",f0(x)," 0.5 chi_opt: ",0.5*chi_opt)

x_opt = x
writedlm("optimum_values.txt",x_opt,chi_opt)

#read(stdin,Char)

nrand = 10
dxrand = zeros(nparam,nrand)
nll_exact = zeros(nrand)
nll_approx = zeros(nrand)
for i=1:nrand
#for i=1:nparam
# Draw random parameters:
  dev = randn(nparam)/sqrt(nparam)
#  dev = zeros(nparam); dev[i]=1.0
  dx = cholhinv.L*dev
  x_rand = x + dx
  dxrand[:,i] = dx
  nll_exact[i]  = f0(x_rand);
  nll_approx[i] =  0.5*chi_opt + 0.5*dx'*hess_save*dx; #  Figure out the correct quadratic approximation
  # Compare:
  println(i," nll_exact: ",nll_exact[i]-0.5*chi_opt," nll_approx: ",nll_approx[i]-0.5*chi_opt," Difference: ",nll_approx[i]-nll_exact[i])
#  read(stdin,Char)
end



#--------------------------------------------------------------
# Now, run a Hamiltonian Markov chain:

# Create a parameter vector to hold orbital elements & mass ratios of planets:
npar = (n-1)*5
if use_student
  # Add two variables for Student's t model:
  npar += 2
end

# Recompute chi-square:
chi_ttv = 0.0
dq0 = ttv_elements!(n,t0,h,tmax,elements,tt2,count2,0.0,0,0,rstar)
chi_ttv = 0.0
for j=1:ntrans
  chi_ttv += lnprob([tobs_tot[j],tt2[iplanet[j],indx[j]],sobs_tot[j],lgndof,V1exp2nuinv])
end
# Add in normalization:
chi_ttv += ntrans*log_norm_students_t([lgndof,V1exp2nuinv])


# Check that prior is zero:
println("-Log prior: ",nlog_prior(x))

# The optimum chi-square:
chi_opt = chi_grid[1,jgrid[1]]

# Now invert to compute momentum "inverse mass" matrix:
#invmass = inv(hess_save)
invmass = 0.5*(inv(hess_save)+inv(hess_save)')

# "Time" step is 0.05:
#epsilon0 = 0.20

# Start out with a modest number of steps:
#nstep = 1
#nstep = 500
state = zeros(2*npar+2,nstep)

# Cholesky decompose hessian for computation of momentum:
#hessian = .5*(hessian + hessian')
hessian = .5*(hess_save + hess_save')
#cholh = transpose(chol(hessian))
#cholh = permutedims(chol(hessian))
cholh = permutedims(cholesky(hessian).U)

# Leapfrog (Neal 2012):
#nleap0 = 1000
#nleap0 = 10
if showplots
  substep = zeros(npar,nleap0)
end
nacc = 0
grad = zeros(Float64,npar)
hamnew = 0.0
#tic(); elapsed = 0.0
tstart = time(); elapsed = 0.0

function compute_gradient!(grad,x)
# Check that variables are defined:
#println("n: ",n," ntrans: ",ntrans," t0: ",t0," h: ",h," tmax: ",tmax)
# Begin gradient computation.
# Put the x parameters into the elements matrix:
for i=2:n
  for j=1:5
    elements_trial[i,j] = x[(i-2)*5+j]
  end
end
dtdelements = ttv_elements!(n,t0,h,tmax,elements_trial,tt1,count1,dtdq0,rstar)
# Compute the gradient & Hessian:
i = 0; fill!(grad,0.0)
lgndof = x[nparam-1]; V1exp2nuinv = x[nparam]
# Iterate over planets & varying parameters:
for ip=1:n-1, ivary=1:nvary
  i +=1
  # If varying mass, this is at end of dtdq0 parameter list:
  if ivary == 1; iv = 7; else; iv= ivary-1; end
  # Gradient of chi-square with respect to each varying parameter:
  for it=1:ntrans
    p = [tobs_tot[it],tt1[iplanet[it],indx[it]],sobs_tot[it],lgndof,V1exp2nuinv]
    # Take derivative of -2*log likelihood with respect to model at this point, and then compute derivative of model with
    # respect to model parameters:
    grad_tmp = grad_point(p)
    dtdxi = dtdelements[iplanet[it],indx[it],iv,ip+1]
    grad[i] += 0.5*grad_tmp[2]*dtdxi
    if use_student
      if i == 1
        # Derivative wrt lgndof:
        grad[nparam-1] += 0.5*grad_tmp[4]
        # Derivative wrt V1exp2nuinv:
        grad[nparam  ] += 0.5*grad_tmp[5]
      end
    end
  end
end
# Add in normalization:
grad[nparam-1:nparam] .+= 0.5*ntrans*grad_norm([lgndof,V1exp2nuinv])
# Add in prior:
grad .+= grad_nlog_prior(x)
# End gradient computation
return
end

function check_bounds(x)
# Check that we're within bounds:
check = true
for i=2:n
  # Check if eccentricity >= 1.0
  if x[(i-2)*5+4]^2+x[(i-1)*5]^2 >= 1.0
    check = false
  end
  # Check if masses are non-positive:
  if x[(i-2)*5+1] <= 0.0
    check = false
  end
end
# Check if ndof is in bounds:
if exp(x[nparam-1]) < 0.5 || exp(x[nparam-1]) > 100.0
  check = false
end
return check
end

# Create array to save data on Markov chain:
stats = zeros(Float64,5,nstep)  # save epsilon, nleap, alpha, accept/reject
fofx = f0(x)
fofxnew = 0.0

# Loop over steps in Hamiltonian Markov Chain:
for istep = 1:nstep
# Draw momentum from the kinetic energy distribution:
  p = cholh*randn(npar)
  ham = 0.5*p'*invmass*p + fofx
  pnew = p
  xnew = x
  nleap = round(Int64,nleap0*(0.8+0.2*rand()))
  if showplots
    # Save substeps to make plots of progress:
    substep = zeros(npar,nleap+1)
    substep[:,1]=xnew
  end
  # Compute gradient with respect to parameters:
  compute_gradient!(grad,xnew)
  # Choose a random value of epsilon:
#  epsilon = abs(randn())*epsilon0
# Choose between 0.03-0.04:
  epsilon = (0.75+0.25*rand())*epsilon0
  # Initial step of leapfrog integration:
  pnew -= 0.5*epsilon*grad
  # Loop over leapfrog steps:
  ileap = 1
  while ileap <= nleap
    # parameter update:
    xnew += epsilon*invmass*pnew
#    println("xnew: ",xnew)
    # check that we're within bounds for eccentricity:
    if check_bounds(xnew)
      if showplots
        substep[:,ileap+1]=xnew
      end
    # momentum update:
      compute_gradient!(grad,xnew)
      pnew -= epsilon*grad
      ileap += 1
    else
      println("Parameters out of bounds: ",xnew)
      ileap = nleap+1
    end
  end
  # Final momentum update:
  pnew += 0.5*epsilon*grad
  # Now, compute hamiltonian again:
  if check_bounds(xnew)
    fofxnew = f0(xnew)
    hamnew = 0.5*pnew'*invmass*pnew+ fofxnew 
  # Compute whether to accept or reject:
    alpha = exp(ham-hamnew)
  else
    alpha = 0.0
    hamnew = 0.0
  end
  uni = rand()
  if uni < alpha
    x .= xnew
    p .= pnew
    fofx = copy(fofxnew)
    println(istep," accepted, 1: ",ham," 2: ",hamnew," eps: ",epsilon," L: ",nleap," ",fofx)
    ham = copy(hamnew)
    nacc += 1
    # plot masses of inner two planets:
    if showplots
      plot(substep[1,:],substep[6,:])
    end
    stats[:,istep] = [epsilon;float(nleap);alpha;uni;1.0]
  else
    println(istep," rejected, 1: ",ham," 2: ",hamnew," eps: ",epsilon," L: ",nleap," ",fofx)
    stats[:,istep] = [epsilon;float(nleap);alpha;uni;0.0]
  end
  state[:,istep] = [x;-p;ham;fofx]
  if mod(istep,10) == 0
    elapsed += time()-tstart 
    println(istep," elapsed: ",elapsed/3600.," facc: ",nacc/istep," ",(nstep-istep)/istep*elapsed/3600.," hr remaining")
    tstart = time()
    # That's it!  Now save every variable to an HDF5 JLD2 file (checkpoint):
    @save foutput n cov_save hess_save fname felements nparam t0 h tmax ntrans iplanet indx tobs_tot sobs_tot data count1 state hessian cholh nleap0 nacc nstep stats elements chi_opt x_opt
  end
end
println("Fraction accepted: ",nacc/nstep)
#--------------------------------------------------------------

# End of run_hmc
return
end

# Example of running this routine:

# File containing the transit times:
#fname = "T1_timings_20191203.txt"

#felements = "elements_noprior_students.txt"

# Name an output file:
#foutput = "T1_run_hmc_student_ecc_lgndof_V1exp2nuinv_001.jld2"

# Call the likelihood profile function:
#run_hmc(fname,felements,foutput,nstep,epsilon0,nleap0,showplots)
#run_hmc(fname,felements,foutput,100,0.2,10,true)
