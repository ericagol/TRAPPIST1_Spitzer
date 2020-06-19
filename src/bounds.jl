

# Function for placing bounds on a variable:

function log_bounds_upper(x::T,x1::Float64,x2::Float64) where {T <: Real}
@assert(x2 > x1)
# Places bounds on variable, gradually decreasing probability
# from x1 to x2:
if x <= x1
  # Return prior probability and derivative:
  return zero(T),zero(T)
elseif x < x2
  # Cubic interpolation:
  y = (x2-x)/(x2-x1)
  return log(3y^2-2y^3),6/(x2-x1)*(y^2-y)
end
return -Inf,zero(T)
end


function log_bounds_lower(x::T,x1::Float64,x2::Float64) where {T <: Real}
@assert(x2 > x1)
# Places bounds on variable, gradually decreasing probability
# from x2 to x1:
if x <= x1
  # Return prior probability of zero below x1: 
  return -Inf,zero(T)
elseif x < x2
  # Cubic interpolation between x1 and x2:
  y = (x-x1)/(x2-x1)
  return log(3y^2-2y^3),6/(x2-x1)*(y-y^2)
end
# Above x2, return probability of unity:
return zero(T),zero(T)
end

