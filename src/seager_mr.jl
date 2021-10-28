function seager_mr(mp,comp)
# computes the radius of a planet
# given the mass & composition of the planet
# according to the formulae in Seager et al. (2007)
k1=-0.20945
k2=0.0804
k3=0.394
if comp == "Fe"
  m1=4.34
  r1=2.23
end
if comp == "Earth"
  m1=6.41
  r1=3.19
#  r1=2.84
end
if comp == "MgSiO3"
  m1=7.38
  r1=3.58
end
if comp == "H2O"
  m1=8.16
  r1=4.73
end
if comp == "FeMgH2O"
#  m1 = 6.88
#  r1 = 4.02
  m1 = 6.41
  r1 = 3.63
end
ms=mp/m1
lgrp=k1 .+ log10.(ms)/3.0 .-k2*ms.^k3
return 10.0.^lgrp*r1
end
