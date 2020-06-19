function compute_grad_num(n::Int64,h::Float64,t0::Float64,tmax::Float64,elements::Array{Float64,2},tt2::Array{Float64,2},
   count2::Array{Int64,1},rstar::Float64,tobs_tot::Array{Float64,1},sobs_tot::Array{Float64,1},lgndof::Float64,V1exp2nuinv::Float64,
   indx::Array{Int64,1},nparam::Int64,nplanet::Int64,nvary::Int64)
# compute_grad_num(n,h,t0,tmax,elements,tt2,count2,rstar,tobs_tot,sobs_tot,lgndof,V1exp2nuinv,indx,nparam,nplanet,nvary)
# Computes gradient numerically.
# First, compute model in BigFloat precision:
hbig = big(h); t0big = big(t0); tmaxbig = big(tmax); elements_big = big.(elements); rstarbig = big(rstar)
tt2big = big.(tt2)
zbig = zero(BigFloat)
@time dq0big = ttv_elements!(n,t0big,hbig,tmaxbig,elements_big,tt2big,count2,zbig,0,0,rstarbig)
lgndof_big = big(lgndof); V1exp2nuinv_big = big(V1exp2nuinv); tobs_tot_big = big.(tobs_tot); sobs_tot_big = big.(sobs_tot)
chi_ttv = zero(BigFloat)
for j=1:ntrans
  chi_ttv += lnprob([tobs_tot_big[j],tt2big[iplanet[j],indx[j]],sobs_tot_big[j],lgndof_big,V1exp2nuinv_big])
end
# Add in Student's t normalization:
chi_ttv += ntrans*log_norm_students_t([lgndof_big,V1exp2nuinv_big])
gradf_num = zeros(BigFloat,nparam)
i = 0
# Now perturb each element and compute numerical derivative:
println("Computing numerical gradient")
dq = big(1e-15)
for ip=1:nplanet, ivary=1:nvary
  tstart = time()
  i +=1
  # If varying mass, this is at end of dtdq0 parameter list:
  if ivary == 1; iv = 7; else; iv= ivary-1; end
  # Perturb the parameter associated with ivary:
  lgndof_1 = big(0.0); lgndof_2 = big(0.0); V1exp2nuinv_1 = big(0.0); V1exp2nuinv_2 = big(0.0)
  if ivary < nparam-1
    # Compute finite difference derivative with respect to elements[ip+1,ivary]
    elements_big_1 = copy(elements_big); tt2big_1 = copy(tt2big)
    elements_big_1[ip+1,ivary] = elements_big[ip+1,ivary]*(1-dq)
    elements_big_2 = copy(elements_big); tt2big_2 = copy(tt2big)
    elements_big_2[ip+1,ivary] = elements_big[ip+1,ivary]*(1+dq)
    dq0big1 = ttv_elements!(n,t0big,hbig,tmaxbig,elements_big_1,tt2big_1,count2,zbig,0,0,rstarbig)
    dq0big2 = ttv_elements!(n,t0big,hbig,tmaxbig,elements_big_2,tt2big_2,count2,zbig,0,0,rstarbig)
  end 
  # Gradient of "chi-square" with respect to each varying parameter:
  for it=1:ntrans
    # Set up parameter vector for computing likelihood of this data point:
    if ivary < nparam-1
      p1 = [tobs_tot_big[it],tt2big_1[iplanet[it],indx[it]],sobs_tot_big[it],lgndof_big,V1exp2nuinv_big]
      p2 = [tobs_tot_big[it],tt2big_2[iplanet[it],indx[it]],sobs_tot_big[it],lgndof_big,V1exp2nuinv_big]
      gradf_num[i] += (lnprob(p2)-lnprob(p1))/(2dq*elements_big[ip+1,ivary])
    end
    # gradf_num[i] += grad_tmp[2]*dtdelements[iplanet[it],indx[it],iv,ip+1]
    # Compute with finite difference
    if use_student
      if i == 1
      # Derivative wrt ndof:
        # gradf_num[nparam-1] += grad_tmp[4]
        lgndof_1 = lgndof_big *(1-dq)
        lgndof_2 = lgndof_big *(1+dq)
        p1 = [tobs_tot_big[it],tt2big[iplanet[it],indx[it]],sobs_tot_big[it],lgndof_1,V1exp2nuinv_big]
        p2 = [tobs_tot_big[it],tt2big[iplanet[it],indx[it]],sobs_tot_big[it],lgndof_2,V1exp2nuinv_big]
        deriv_num = (lnprob(p2)-lnprob(p1))/(2dq*lgndof_big)
        if isnan(deriv_num)
          println("NaN derivative.  p1: ",p1," p2: ",p2," dq: ",dq," lgndof_big: ",lgndof_big," lnprob(p1): ",lnprob(p1)," lnprob(p2): ",lnprob(p2))
        else
          gradf_num[nparam-1] += deriv_num
        end
        # Derivative wrt V1exp2nuinv:
        # gradf_num[nparam  ] += grad_tmp[5]
        V1exp2nuinv_1 = V1exp2nuinv_big * (1-dq)
        V1exp2nuinv_2 = V1exp2nuinv_big * (1+dq)
        p1 = [tobs_tot_big[it],tt2big[iplanet[it],indx[it]],sobs_tot_big[it],lgndof_big,V1exp2nuinv_1]
        p2 = [tobs_tot_big[it],tt2big[iplanet[it],indx[it]],sobs_tot_big[it],lgndof_big,V1exp2nuinv_2]
#        gradf_num[nparam] +=
        deriv_num = (lnprob(p2)-lnprob(p1))/(2dq*V1exp2nuinv_big)
        if isnan(deriv_num)
          println("NaN derivative.  p1: ",p1," p2: ",p2," dq: ",dq," V1exp2nuinv_big: ",V1exp2nuinv_big," lnprob(p1): ",lnprob(p1)," lnprob(p2): ",lnprob(p2))
        else
          gradf_num[nparam] += deriv_num
        end
      end
    end
  end
  telapse = time()- tstart
  println("Done with component: ",ip,", ",ivary," grad: ",gradf_num[i]," time: ",telapse/60," min" )
end
# Now add in derivative of normalization:
lgndof_1 = lgndof_big *(1-dq)
lgndof_2 = lgndof_big *(1+dq)
norm_1 = log_norm_students_t([lgndof_1,V1exp2nuinv_big])
norm_2 = log_norm_students_t([lgndof_2,V1exp2nuinv_big])
deriv_num = ntrans*(norm_2-norm_1)/(2dq*lgndof_big)
#gradf_num[nparam-1] += ntrans*(norm_2-norm_1)/(2*lgndof_big*dq)
if isnan(deriv_num)
  println("NaN derivative.  lgndof_1: ",lgndof_1," lgndof_2: ",lgndof_2," dq: ",dq," V1exp2nuinv_big: ",V1exp2nuinv_big," norm_1: ",norm_1," norm_2: ",norm_2," ntrans: ",ntrans)
else
  gradf_num[nparam-1] += deriv_num
end
V1exp2nuinv_1 = V1exp2nuinv_big * (1-dq)
V1exp2nuinv_2 = V1exp2nuinv_big * (1+dq)
norm_1 = log_norm_students_t([lgndof_big,V1exp2nuinv_1])
norm_2 = log_norm_students_t([lgndof_big,V1exp2nuinv_2])
deriv_num = ntrans*(norm_2-norm_1)/(2dq*V1exp2nuinv_big)
#gradf_num[nparam  ] += ntrans*(norm_2-norm_1)/(2*V1exp2nuinv_big*dq)
if isnan(deriv_num)
  println("NaN derivative.  V1exp2nuinv_1: ",V1exp2nuinv_1," V1exp2nuinv_2: ",V1exp2nuinv_2," dq: ",dq," lgndof_big: ",lgndof_big," norm_1: ",norm_1," norm_2: ",norm_2," ntrans: ",ntrans)
else
  gradf_num[nparam] += deriv_num
end
return gradf_num
end 
