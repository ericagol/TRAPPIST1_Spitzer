
include("bounds.jl")

# Compute the negative log prior:

function nlog_prior_noecc(x::Array{T,1}) where {T <: Real}
# We will place a joint prior on eccentricity vector
# such that each planet has an eccentricity which lies between
# 0 and emax2, with a gradual decrease from emax1 to emax2:
# Eccentricity priors:
emax1 = 0.2; emax2 = 0.3
# Mass priors:
lgm1 = -8.0; lgm2 = -7.0
lgm3 = -3.0; lgm4 = -2.0
# Period priors:
P1 = [1.49,2.40,4.03,6.08,9.18,12.33,18.75]; P2 = [1.50,2.41,4.04,6.09,9.19,12.34,18.76]
P3 = [1.52,2.43,4.06,6.11,9.22,12.36,18.78]; P4 = [1.53,2.44,4.07,6.12,9.23,12.37,18.79]
# t0 priors:
t01 = [7257.53,7258.57,7257.05,7257.81,7257.05,7257.70,7249.58]
t02 = [7257.54,7258.58,7257.06,7257.82,7257.06,7257.71,7249.59]
t03 = [7257.56,7258.60,7257.08,7257.84,7257.08,7257.73,7249.61]
t04 = [7257.57,7258.61,7257.09,7257.85,7257.09,7257.74,7249.62]
# Define -log(prior):
nlprior = zero(T)
nparam = size(x)[1]
n = div(nparam-2,5)+1
# Loop over planets:
for i=2:n
  # Place prior on eccentricities:
  ecc = sqrt(x[(i-2)*5+4]^2+x[(i-2)*5+5]^2)
  nlprior_tmp,dpdx = log_bounds_upper(ecc,emax1,emax2)
  nlprior -= nlprior_tmp

  # We will place broad lower and upper limits on log_10(mass):
  log10mass = log10(x[(i-2)*5+1])
  nlprior_tmp,dpdx = log_bounds_lower(log10mass,lgm1,lgm2)
  nlprior -= nlprior_tmp
  nlprior_tmp,dpdx = log_bounds_upper(log10mass,lgm3,lgm4)
  nlprior -= nlprior_tmp
  # Broad priors on period:
  P = x[(i-2)*5+2]
  nlprior_tmp,dpdx = log_bounds_lower(P,P1[i-1],P2[i-1])
  nlprior -= nlprior_tmp
  nlprior_tmp,dpdx = log_bounds_upper(P,P3[i-1],P4[i-1])
  nlprior -= nlprior_tmp
  # Broad priors on t0:
  t0 = x[(i-2)*5+3]
  nlprior_tmp,dpdx = log_bounds_lower(t0,t01[i-1],t02[i-1])
  nlprior -= nlprior_tmp
  nlprior_tmp,dpdx = log_bounds_upper(t0,t03[i-1],t04[i-1])
  nlprior -= nlprior_tmp
end
# Add in prior on ndof & lnV1:
ndof1 = 0.5; ndof2 = 1.0
ndof3 = 50.0; ndof4 = 100.0
ndof = exp(x[nparam-1])
nlprior_tmp,dpdx = log_bounds_lower(ndof,ndof1,ndof2)
nlprior -= nlprior_tmp
nlprior_tmp,dpdx = log_bounds_upper(ndof,ndof3,ndof4)
nlprior -= nlprior_tmp
lnV11 = -2.0; lnV12 = -1.0
lnV13 = 5.0; lnV14 = 10.0
lnV1 = log(x[nparam]/exp(0.5/ndof))
nlprior_tmp,dpdx = log_bounds_lower(lnV1,lnV11,lnV12)
nlprior -= nlprior_tmp
nlprior_tmp,dpdx = log_bounds_upper(lnV1,lnV13,lnV14)
nlprior -= nlprior_tmp
# That's it!
return nlprior::T
end

# Take gradient and hessian:
grad_nlog_prior_noecc= x -> ForwardDiff.gradient(nlog_prior_noecc,x)
hess_nlog_prior_noecc= x -> ForwardDiff.hessian(nlog_prior_noecc,x)
