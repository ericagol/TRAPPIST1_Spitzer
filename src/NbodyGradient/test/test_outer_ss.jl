

# Basing off of NbodyGradient/test/outer_ss/NBG_performance.jl
# Trying to profile to see how much time is spent on allocations
# for moving memory back and forth.

include("../src/ttv.jl")

n=5
xout = zeros(3,n)
# Positions at time September 5, 1994 at 0h00 in days (from Hairer, Lubich & Wanner
# 2006, Geometric Numerical Integration, 2nd Edition, Springer, pp. 13-14):

xout .= transpose([-2.079997415328555E-04  7.127853194812450E-03 -1.352450694676177E-05;
        -3.502576700516146E+00 -4.111754741095586E+00  9.546978009906396E-02;
         9.075323061767737E+00 -3.443060862268533E+00 -3.008002403885198E-01;
         8.309900066449559E+00 -1.782348877489204E+01 -1.738826162402036E-01;
         1.147049510166812E+01 -2.790203169301273E+01  3.102324955757055E-01]) #;
        #-1.553841709421204E+01 -2.440295115792555E+01  7.105854443660053E+00])
vout = zeros(3,n)
vout .= transpose([-6.227982601533108E-06  2.641634501527718E-06  1.564697381040213E-07;
         5.647185656190083E-03 -4.540768041260330E-03 -1.077099720398784E-04;
         1.677252499111402E-03  5.205044577942047E-03 -1.577215030049337E-04;
         3.535508197097127E-03  1.479452678720917E-03 -4.019422185567764E-05;
         2.882592399188369E-03  1.211095412047072E-03 -9.118527716949448E-05]) #;
         #2.754640676017983E-03 -2.105690992946069E-03 -5.607958889969929E-04]);
# Units of velocity are AU/day

# Specify masses, including terrestrial planets in the Sun:
m = [1.00000597682,0.000954786104043,0.000285583733151,
         0.0000437273164546,0.0000517759138449] #,6.58086572e-9];

# Compute the center-of-mass:
vcm = zeros(3);
xcm = zeros(3);
for j=1:n
    vcm .+= m[j]*vout[:,j];
    xcm .+= m[j]*xout[:,j];
end
vcm ./= sum(m);
xcm ./= sum(m);
# Adjust so CoM is stationary
for j=1:n
    vout[:,j] .-= vcm[:];
    xout[:,j] .-= xcm[:];
end

h = 50.0
n = 5
pair = zeros(Bool, n, n)  # Pairwise Keplerians for all bodies
xerror = similar(xout)
verror = similar(vout)
fill!(xerror,0.0)
fill!(verror,0.0)
nstep = 1000000
jac_step = Matrix{Float64}(I,7*n,7*n)
jac_error = zeros(Float64,7*n,7*n)

# Now, time it:
@time for i=1:nstep; ah18!(xout,vout,xerror,verror,h,m,n,jac_step,jac_error,pair); end

