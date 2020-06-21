

# Converts from cartesian coordinates to generalized 
# Jacobi Keplerian coordinates, and then converts to orbital
# elements:

using DelimitedFiles
using LinearAlgebra

if ~@isdefined YEAR
  const YEAR  = 365.242
  const GNEWT = 39.4845/YEAR^2
end

include("../../../src/NbodyGradient/src/init_nbody.jl")

function amatrix(mass::Array{T,1},n_body) where {T <: Real}
# Set up "A" matrix (Hamers & Portegies-Zwart 2016) which transforms from
# cartesian coordinates to Keplerian orbits (we are using opposite sign
# convention of HPZ16, for example, r_1 = R_2-R_1).
# Mass vector is "mass".  Returns A matrix, "amat":
amat = zeros(T,n_body,n_body)
# Sum of masses for each Keplerian:
mkep = zeros(T,n_body-1)
# Get the indices:
indices = get_indices_planetary(n_body)
# Fill in the A matrix & compute the Keplerian elements:
for i=1:n_body-1
  # Sums of masses for two components of Keplerian:
  m1 = zero(T)
  m2 = zero(T)
  for j=1:n_body
    if indices[i,j] == 1
      m1 += mass[j]
    end
    if indices[i,j] == -1
      m2 += mass[j]
    end
  end
  mkep[i] = m1+m2
  # Compute Kepler problem: r is a vector of positions of "body" 2 with respect to "body" 1; rdot is velocity vector
  # For now set inclination to Inclination = pi/2 and longitude of nodes to Omega = pi:
  # Fill in the A matrix
  for j=1:n_body
    if indices[i,j] == 1
      amat[i,j] = -mass[j]/m1
    end
    if indices[i,j] == -1
      amat[i,j] =  mass[j]/m2
    end
  end
end
mtot = sum(mass)
# Final row is for center-of-mass of system:
for j=1:n_body
  amat[n_body,j] = mass[j]/mtot
end
return amat,mkep
end

function cartesian_to_elements(t::T,M::T,r::Array{T,1},rdot::Array{T,1}) where {T <: Real}
# Takes 3-d vectors r and rdot, and converts
# these to [a,lambda,e,omega] (assuming y=0):
# Equations (2.126-2.130) in M&D:
r2 = dot(r,r) 
rabs = sqrt(r2)
v2 = dot(rdot,rdot)
rrdot = dot(r,rdot)
h = cross(r,rdot)
h2 = dot(h,h)
habs = sqrt(h2)
dotr = sign(rrdot)*sqrt(abs(v2-h2/r2))
# Now, compute elements:
semi = 1/(2/rabs-v2/(GNEWT*M))
ecc2 = 1-h2/(GNEWT*M*semi)
ecc  = sqrt(ecc2)
sinompf = -r[3]/rabs
#cosompf = -r[1]/rabs  # Do I want to assum \Omega = pi? [ ]
cosompf = r[1]/rabs  # Do I want to assum \Omega = pi? [ ]
ompf = atan(sinompf,cosompf)
sinf = semi*(1-ecc2)*dotr/(habs*ecc)
cosf = (semi*(1-ecc2)/rabs-1)/ecc
omega = atan(sinompf*cosf-cosompf*sinf,cosompf*cosf+sinompf*sinf)
# Eccentric anomaly (2.41): (does this have ambiguity? [ ])
#ean  = acos((1-rabs/semi)/ecc)
# Eccentric anomaly (2.46):
ean = 2 * atan(sqrt(1-ecc)*sinf,sqrt(1+ecc)*(1+cosf))
# Mean anomaly:
man  = ean - ecc*sin(ean)
# Now, we can compute lambda:
lam = man + omega
return [semi,lam,ecc,omega]
end

# Function which takes file name, reads in data, and computes
# the orbital elements:
function compute_elements(fname::String,masses::Array{T,1},n::Int64) where {T <: Real}
data = readdlm(fname)
# Fill out position and velocity arrays:
t = data[:,1]
nt = length(t)
x = reshape(data[:,2:n*3+1],nt,3,n)
v = reshape(data[:,n*3+2:end],nt,3,n)
# First, convert from Cartesian coordinates to
# Keplerian coordinates:
rkepler    = zeros(T,nt,n,3)
rdotkepler = zeros(T,nt,n,3)
# Construct the A matrix:
amat,mkep = amatrix(masses,n)
# Loop over time, and convert each set of cartesian coordinates
# to Keplerian coordinates:
for it=1:nt
  for i=1:n,j=1:3,k=1:n
    rkepler[it,i,j] += amat[i,k]*x[it,j,k]
    rdotkepler[it,i,j] += amat[i,k]*v[it,j,k]
  end 
end  
# Okay, now for conversion of Keplerian coordinates to
# orbital elements.  We will assume plane-parallel,
# and will compute [a,\lambda,e,\omega] at each time step:
elements = zeros(T,nt,4,n)
for it=1:nt
  # Compute orbital elements:
  for i=1:n-1
    elements[it,:,i] = cartesian_to_elements(t[it],mkep[i],rkepler[it,i,:],rdotkepler[it,i,:])
    if it > 1
      while elements[it,2,i] < (-pi/2 + elements[it-1,2,i])
        elements[it,2,i] += 2pi
      end
      while elements[it,2,i] > ( pi/2 + elements[it-1,2,i])
        elements[it,2,i] -= 2pi
      end
    end
  end
end
# Return array of elements as a function of time for each Keplerian
# in the Jacobi hierarchy:
return t,elements
end
