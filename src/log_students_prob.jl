# Define functions for the Student's t model:

function log_students_t_pdf(x::T,sig::T,nu::T) where {T <: Real}
# -2 log PDF of the Student's t distribution, without normalization:
# fac = gamma(alp)/sqrt(pi*nu)/gamma(0.5*nu)/sig
 # pdf[i] = fac/(1+(x[i]/sig)^2/nu)^alp
return (nu+1)*log(1+(x/sig)^2/nu)
#return pdf/sig
end

# Define -2 log normalization of Student's t distribution:
function log_norm_students_t(x::Array{T,1}) where {T <: Real}
nu,lnV1 = x
# Normalization is: Gamma((\nu+1)/2)/\sqrt{\nu \pi}/Gamma(\nu/2)/\sigma
# (sigma is factor by which the error bars are inflated - it equals sqrt(exp(lnV1)))
return -2*(log(gamma(0.5*(nu+1)))-log(gamma(0.5*nu)))+log(pi*nu)+lnV1
end


# Full log probabilistic model for a single datapoint (missing prior - should add this in [ ]):
function lnprob(p::Array{T,1}) where {T <: Real}
  # Returns -2 log(likelihood) for Student's t distribution.
  # y = observed value; model = model prediction; yerr = uncertainty (sqrt variance)
  # ndof = number of degrees of freedom (\nu) for Student's t distribution;
  # lnV1 = log of factor inflating variance foreground = log(sig^2)
  y, model, yerr, ndof, lnV1 = p

  ll = log_students_t_pdf((y-model)/yerr,exp(0.5*lnV1),ndof)

# return negative log likelihood (so optimize can minimize this,
# or mazimize likelihood).  If the likelihood were exp(-\chi^2/2), which
# it is for ndof = infinity and lnV1 = 0, then this would return \chi^2:
  return ll #+log_norm_students_t([ndof,lnV1])
end

# Compute gradients of this function:
grad_point = p -> ForwardDiff.gradient(lnprob,p)
hess_point = p -> ForwardDiff.hessian(lnprob,p)

# Compute gradients of normalization:
grad_norm = x -> ForwardDiff.gradient(log_norm_students_t,x)
hess_norm = x -> ForwardDiff.hessian(log_norm_students_t,x)

