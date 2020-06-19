
# Define PDF & CDF for Student's t distribution:

using SpecialFunctions
using GSL
using Optim
using DelimitedFiles

include("../../../src/loglinspace.jl")

function students_t_pdf(x,nx,sig,nu)
# PDF of the Student's t distribution.
pdf = zeros(nx)
alp = 0.5*(nu+1)
fac = gamma(alp)/sqrt(pi*nu)/gamma(0.5*nu)
for i=1:nx
  pdf[i] = fac/(1+(x[i]/sig)^2/nu)^alp
end
return pdf/sig
end

function students_t_cdf(x,nx,sig,nu)
# Computes the cumulative distribution function for Student's t distribution:
cdf = zeros(nx)
alp = 0.5*(nu+1)
fac = gamma(alp)/sqrt(pi*nu)/gamma(0.5*nu)
for i=1:nx
  y = (x[i]/sig)^2/nu
  if y < 1.0
    cdf[i] = 0.5+x[i]/sig*fac*hypergeom([0.5,alp],1.5,-y)
  else
    cdf[i] = 0.5+x[i]/sig*fac*hypergeom([0.5,1.5-alp],1.5,y/(1+y))/sqrt(1+y)
  end
end
return cdf
end

function fit_students_t(resid,s1,ndof)
# Fits a discrete distribution:
  # Fits a students-t distribution to the 
  x0 = [s1,ndof]
  nresid = size(resid)[1]
  function student_model(x)
    s1,ndof = x
    println("s1: ",s1," ndof: ",ndof)
    if ndof >= 0.5 && s1 > 0.0
      pdf = students_t_pdf(resid,nresid,s1,ndof)
    else
      return 0.0
    end
    # Return a "chi-square":
    return -sum(log.(pdf))
  end
  xbest = Optim.optimize(student_model, x0)
return xbest.minimizer
end

function fit_students_t_hist(resid,binhist,s1,ndof)
# Fits a binned distribution:
  norm_hist = sum(binhist)
  # Fits a students-t distribution to the 
  x0 = [s1,ndof,norm_hist]
  sigbin = copy(binhist)
  sigbin[sigbin .== 0.0] .= 1.0
  function student_model(x)
    s1,ndof,norm_hist = x
    println("s1: ",s1," ndof: ",ndof," norm: ",norm_hist)
    if ndof >= 0.0
      pdf = students_t_pdf(resid,size(resid)[1],s1,ndof)
    else
      return 0.0
    end
    # Return a "chi-square":
#    return sum((binhist[binhist .!= 0.0] .- pdf[binhist .!= 0.0].*norm_hist).^2 ./binhist[binhist .!= 0.0])
    return sum((binhist .- pdf.*norm_hist).^2 ./sigbin)
  end
  xbest = Optim.optimize(student_model, x0)
return xbest.minimizer
end

# Load in Julia wrapper of matplotlib:
using PyPlot
# Load Special Functions package for computing erf:
using SpecialFunctions

fig,axes = subplots(1,2,figsize=(10,5))

#clf()
# Read in the residuals:
dev = readdlm("../../../data/residuals_nomixture_rescale.txt")
# Plot the cumulative distribution function of the residuals:
ax = axes[1]
nobs = 447
ax.plot(sort(dev[:,1]),linearspace(0,1,nobs),linewidth=3,label="CDF")
# Create an array for plotting normalized residuals:
nx = 10000
x = linearspace(-7.5,7.5,nx)
# Plot CDF of a Gaussian with \sigma = 0.9:
#p = 0.8
p = 1.0 - 0.3150298440418438
s1 = 1.0
d1 = 0.5*(1 .+erf.(x./sqrt(2)./s1))
ax.plot(x,d1,label="Single Gaussian",linewidth=1)
# Plot CDF of a Gaussian with \sigma = 2.0:
#s1 = 0.87
s1 = 0.7312227481615532
d1 = 0.5*(1 .+erf.(x./sqrt(2)./s1))
#s2 = 2.0
s2 = 1.9192905742029427
d2 = 0.5*(1 .+erf.(x./sqrt(2)./s2))
#plot(x,d2)
# Sum of CDFs:
#ax.plot(x,p*d1 .+ (1-p)*d2,label="Double Gaussian",linewidth=2,linestyle="--")
#ax.set_title("Cumulative distribution of normalized residuals of best TTV model")
# Student's t distribution - first optimize the fit to the histogram:
sig_est = 0.835
dof = 3.27
ax.plot(x,students_t_cdf(x,nx,sig_est,dof),label=L"Student's-t",linewidth=2,linestyle=":")
ax.set_title("Cumulative distribution")
# Label the plot:
ax.set_xlabel(L"$z = (t_{\mathrm{obs},ij}-t_{ij}(\mathbf{x}))/\sigma_{ij}$")
ax.set_ylabel(L"$P(>z)$")
ax.legend(loc = "lower right",fontsize=10)
ax.axis([-8,8,-0.05,1.05])

# Now plot a histogram:
nbin = 50; x10 = -7.5; x20 = 7.5; dx = (x20-x10)/nbin
dev_hist = zeros(nbin)
for i=1:length(dev)
  ibin = ceil(Int64,(dev[i]-x10)/dx)
  if ibin > 0 && ibin <= nbin
    dev_hist[ibin] += 1.0
  end
end
res = []
resc = []
hbin = []
for ibin=1:nbin
  push!(res,x10 + (ibin-1)*dx)
  push!(resc,x10 + (ibin-0.5)*dx)
  push!(hbin,dev_hist[ibin])
  push!(res,x10 + ibin*dx)
  push!(hbin,dev_hist[ibin])
end
ax = axes[2]
ax.semilogy(res,hbin,linestyle="None")
errs = sqrt.(dev_hist)
#errs[errs .== 0.0] .= 1.0
ax.errorbar(resc,dev_hist,yerr=errs,fmt=".",color="#1f77b4")

# Optimize Student's-t distribution:
ndof = 2.8
s1 = 0.79
#xbest = fit_students_t_hist(resc,dev_hist,s1,ndof)
xbest = fit_students_t(dev,s1,ndof)
sig_est = xbest[1]
dof = xbest[2]
norm_hist = size(dev)[1]

# Now overplot single and double Gaussian models:
nx = 1000
devs = linearspace(-8,8,nx)
s1 = 1.0
g1 = dx*nobs/sqrt(2*pi*s1^2).*exp.(-0.5.*(devs/s1).^2)
ax.plot(devs,g1,linewidth=2)
ax.axis([-8,8,1e-1,1e2])
#ax.set_xlabel("x = (data-model)/sigma")
ax.set_xlabel(L"$z = (t_{\mathrm{obs},ij}-t_{ij}(\mathbf{x}))/\sigma_{ij}$")
ax.set_ylabel(L"$dP/dz$")
#s1 = 0.87
s1 = 0.7312227481615532
g2 = p*dx*nobs/sqrt(2*pi*s1^2).*exp.(-0.5.*(devs/s1).^2).+
     (1-p)*dx*nobs/sqrt(2*pi*s2^2).*exp.(-0.5.*(devs/s2).^2)
#ax.plot(devs,g2,linestyle="--",linewidth=2)
ax.set_title("Histogram")

# Overplot Student's t distribution:
stpdf = students_t_pdf(devs,nx,sig_est,dof)
ax.plot(devs,stpdf*norm_hist*0.3,linestyle=":",linewidth=3)
#ax.legend(ncol=3,fontsize=8)

savefig("../Students_t_optimized.pdf",bbox_inches="tight")
