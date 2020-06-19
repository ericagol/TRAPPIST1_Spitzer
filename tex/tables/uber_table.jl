


# Create the "Uber" table with the combined results from
# transit-timing, photodynamics & mass estimate of star:
if ~@isdefined(CGS)
  include("../../src/CGS.jl")
  using Main.CGS
end

using JLD2; using Statistics; using Printf

@load "T1_mass_radius_cmf_000.jld2"

# "T1_mass_radius_cmf_000.jld2": nsamp msamp rsamp cmf rho_samp mstar_samp rstar_samp nr nm mrgrid rgrid mgrid level_mr conf m1 m2 r1 r2 mavg ravg merror rerror

# Things to compute:
# Mass, radius, density, surface gravity (for planets and star).
# Semi-major axis, angular momentum.
# Angular momentum deficit.
# Transit duration, transit depth.
# Irradiation.
# Based on Ducrot et al. (2020) luminosity of L/L_sun = 0.000553±0.000019
# we can compute the updated effective temperature and irradiance of each
# of the planets, and log(g) of the star.

# Compute stellar parameters:

Tstar = zeros(nsamp)
# Draw from Ducrot et al. (2020) luminosity:
Lbol = 0.000553 .+ randn(nsamp)*0.000019
# Now, compute the temperature:
Tstar = sort(TSUN * Lbol.^0.25  ./ sqrt.(rstar_samp))

i1 = ceil(Int64,0.158655*nsamp)
i2 = floor(Int64,0.841345*nsamp)
# Compute confidence limits on the temperature:
#tmed = median(Tstar); sigt_lower = tmed - Tstar[i1]; sigt_upper = Tstar[i2]-tmed
tmed = median(Tstar); sigt = std(Tstar)
# Compute confidence limits on the stellar radius:
rsort = sort(rstar_samp)
#rmed = median(rsort); sigr_lower = rmed - rsort[i1]; sigr_upper = rsort[i2]-rmed
rmed = median(rsort); sigr= std(rstar_samp)
# Compute confidence limits on log(g):
logg = log10.(GRAV*MSUN/RSUN^2) .+ sort( log10.(mstar_samp) .- 2log10.(rstar_samp))
loggmed = median(logg); siglg_lower = loggmed - logg[i1]; siglg_upper = logg[i2]-loggmed

# Should we repeat density?  Should we add gravitational redshift?

println("\$ M \$ [\$ M_\\odot \$] & \$ 0.0898 \\pm 0.0023 \$ \\cr")
println("\$ R \$ [\$ R_\\odot \$] & \$ ",@sprintf("%6.4f",rmed), " \\pm ",@sprintf("%6.4f",sigr), "\$ \\cr")
println("\$ L \$ [\$ L_\\odot \$] & \$ 0.000553 \\pm 0.000019 \$ \\cr")
println("\$ T_{eff} \$ [K] & \$",@sprintf("%6.0f",tmed)," \\pm ",@sprintf("%6.0f",sigt), "\$ \\cr")
println("\$ \\log_{10} (g \\mathrm{[cm / s^2]}) \$ & \$",@sprintf("%6.4f",loggmed),"_{-",@sprintf("%6.4f",siglg_lower), "}^{+",@sprintf("%6.4f",siglg_upper),"} \$ \\cr")

# Now compute the radii of the planets:
rstring = " \$ R\$ [\$R_\\oplus\$] & "
rcm_string = " \$ R\$ [\$10^8\$ cm] & "
for i=1:7
  radius = sort(rsamp[i,:])
  rmed = median(radius); sigr_lower = rmed - radius[i1]; sigr_upper = radius[i2]-rmed
  global rstring = string(rstring," \$ ",@sprintf("%6.3f",rmed),"_{-",@sprintf("%6.3f",sigr_lower),"}^{+",@sprintf("%6.3f",sigr_upper),"} \$ ")
  global rcm_string = string(rcm_string," \$ ",@sprintf("%6.3f",rmed*REARTH*1e-8),"_{-",@sprintf("%6.3f",sigr_lower*REARTH*1e-8),"}^{+",@sprintf("%6.3f",sigr_upper*REARTH*1e-8),"} \$ ")
  if i < 7
     rstring = string(rstring," & ")
     rcm_string = string(rcm_string," & ")
  else
     rstring = string(rstring," \\cr ")
     rcm_string = string(rcm_string," \\cr ")
  end
end
println(rstring)

# Now compute the masses of the planets:
mstring = " \$ M\$ [\$M_\\oplus\$] & "
mgram_string = " \$ M\$ [\$10^{27}\$ g] & "
for i=1:7
  mass = sort(msamp[i,:])
#  mmed = median(mass); sigm_lower = mmed - mass[i1]; sigm_upper = mass[i2]-mmed
  mmed = mean(mass); sigm = std(mass)
#  global mstring = string(mstring," \$ ",@sprintf("%6.3f",mmed),"_{-",@sprintf("%6.3f",sigm_lower),"}^{+",@sprintf("%6.3f",sigm_upper),"} \$ ")
  global mstring = string(mstring," \$ ",@sprintf("%6.3f",mmed)," \\pm ",@sprintf("%6.3f",sigm)," \$ ")
#  global mgram_string = string(mgram_string," \$ ",@sprintf("%6.3f",mmed*MEARTH*1e-27),"_{-",@sprintf("%6.3f",sigm_lower*MEARTH*1e-27),"}^{+",@sprintf("%6.3f",sigm_upper*MEARTH*1e-27),"} \$ ")
  global mgram_string = string(mgram_string," \$ ",@sprintf("%6.3f",mmed*MEARTH*1e-27)," \\pm ",@sprintf("%6.3f",sigm*MEARTH*1e-27)," \$ ")
  if i < 7
     mstring = string(mstring," & ")
     mgram_string = string(mgram_string," & ")
  else
     mstring = string(mstring," \\cr ")
     mgram_string = string(mgram_string," \\cr ")
  end
end
println(mstring)

# Now compute the densities of the planets:
rhostring = " \$ \\rho\$ [\$\\rho_\\oplus\$] & "
rhocgs_string = " \$ \\rho \$ [ g cm\$^{-3}\$] & "
rhoearth=3MEARTH/(4pi*REARTH^3)
for i=1:7
  density = sort(msamp[i,:]./rsamp[i,:].^3)
  rhomed = median(density); sigrho_lower = rhomed - density[i1]; sigrho_upper = density[i2]-rhomed
  global rhostring = string(rhostring," \$ ",@sprintf("%6.3f",rhomed),"_{-",@sprintf("%6.3f",sigrho_lower),"}^{+",@sprintf("%6.3f",sigrho_upper),"} \$ ")
  global rhocgs_string = string(rhocgs_string," \$ ",@sprintf("%6.3f",rhomed*rhoearth),"_{-",@sprintf("%6.3f",sigrho_lower*rhoearth),"}^{+",@sprintf("%6.3f",sigrho_upper*rhoearth),"} \$ ")
  if i < 7
     rhostring = string(rhostring," & ")
     rhocgs_string = string(rhocgs_string," & ")
  else
     rhostring = string(rhostring," \\cr ")
     rhocgs_string = string(rhocgs_string," \\cr ")
  end
end
println(rhostring)

# Now compute the surface gravities of the planets:
gstring = " \$ g \$ [\$g_\\oplus\$] & "
gcgs_string = " \$ g\$ [\$10^3\$ cm s\$^{-2}\$] & "
gearth = GRAV*MEARTH/REARTH^2
for i=1:7
  gravity = sort(msamp[i,:]./rsamp[i,:].^2)
#  gmed = median(gravity); sigg_lower = gmed - gravity[i1]; sigg_upper = gravity[i2]-gmed
  gmed = median(gravity); sigg = std(gravity)
#  global gcgs_string = string(gcgs_string," \$ ",@sprintf("%6.3f",gmed*gearth*1e-3),"_{-",@sprintf("%6.3f",sigg_lower*gearth*1e-3),"}^{+",@sprintf("%6.3f",sigg_upper*gearth*1e-3),"} \$ ")
  global gcgs_string = string(gcgs_string," \$ ",@sprintf("%6.3f",gmed*gearth*1e-3)," \\pm ",@sprintf("%6.3f",sigg*gearth*1e-3)," \$ ")
#  global gstring = string(gstring," \$ ",@sprintf("%6.3f",gmed),"_{-",@sprintf("%6.3f",sigg_lower),"}^{+",@sprintf("%6.3f",sigg_upper),"} \$ ")
  global gstring = string(gstring," \$ ",@sprintf("%6.3f",gmed)," \\pm ",@sprintf("%6.3f",sigg)," \$ ")
  if i < 7
     gstring = string(gstring," & ")
     gcgs_string = string(gcgs_string," & ")
  else
     gstring = string(gstring," \\cr ")
     gcgs_string = string(gcgs_string," \\cr ")
  end
end
println(gstring)

# Escape velocities:
vescstring = " \$ v_\\mathrm{esc} \$ [\$v_\\mathrm{esc,\\oplus}\$] & "
vesccgs_string = " \$ v_\\mathrm{esc}\$ [km s\$^{-1}\$] & "
vescearth = sqrt(2*GRAV*MEARTH/REARTH)/1e5  # Escape velocity of Earth in km/sec
for i=1:7
  vesc = sort(sqrt.(msamp[i,:]./rsamp[i,:]))
  vescmed = median(vesc); sigvesc = std(vesc)
  global vesccgs_string = string(vesccgs_string," \$ ",@sprintf("%6.3f",vescmed*vescearth)," \\pm ",@sprintf("%6.3f",sigvesc*vescearth)," \$ ")
  global vescstring = string(vescstring," \$ ",@sprintf("%6.3f",vescmed)," \\pm ",@sprintf("%6.3f",sigvesc)," \$ ")
  if i < 7
     vescstring = string(vescstring," & ")
     vesccgs_string = string(vesccgs_string," & ")
  else
     vescstring = string(vescstring," \\cr ")
     vesccgs_string = string(vesccgs_string," \\cr ")
  end
end
println(vescstring)




# Incident flux:
periods = [1.51087637,2.42180746,4.049959,6.099043,9.205585,12.354473,18.767953]

sstring = " \$ S \$ [\$ S_\\oplus \$] & "
scgs_string = " \$ S \$ [\$10^{6}\\frac{\\mathrm{erg}}{\\mathrm{cm}\\mathrm{s}}\$] & "
year = 365.242
searth = LSUN/(4pi*AU^2)
for i=1:7
  s  = sort(Lbol .* (periods[i]/year)^(-4//3) ./ mstar_samp.^(2//3))
  smed = median(s); sigs_lower = smed - s[i1]; sigs_upper = s[i2]-smed
  global sstring = string(sstring," \$ ",@sprintf("%6.3f",smed),"_{-",@sprintf("%6.3f",sigs_lower),"}^{+",@sprintf("%6.3f",sigs_upper),"} \$ ")
  global scgs_string = string(scgs_string," \$ ",@sprintf("%6.3f",smed*searth*1e-6),"_{-",@sprintf("%6.3f",sigs_lower*searth*1e-6),"}^{+",@sprintf("%6.3f",sigs_upper*searth*1e-6),"} \$ ")
  if i < 7
     sstring = string(sstring," & ")
     scgs_string = string(scgs_string," & ")
  else
     sstring = string(sstring," \\cr ")
     scgs_string = string(scgs_string," \\cr ")
  end
end
println(sstring)

# Semi-major axis:
astring = " \$ a \$ [\$ 10^{-2} \\mathrm{AU} \$] & "
acm_string = " \$ a \$ [\$10^{11}\$ cm] & "

for i=1:7
  a  = sort((periods[i]/year)^(2//3) .* mstar_samp.^(1//3) .*100)
#  amed = median(a); siga_lower = amed - a[i1]; siga_upper = a[i2]-amed
  amed = mean(a); siga= std(a)
  global astring = string(astring," \$ ",@sprintf("%6.3f",amed)," \\pm ",@sprintf("%6.3f",siga)," \$ ")
#  global acm_string = string(acm_string," \$ ",@sprintf("%6.3f",amed*1e-13*AU),"_{-",@sprintf("%6.3f",siga_lower*1e-13*AU),"}^{+",@sprintf("%6.3f",siga_upper*1e-13*AU),"} \$ ")
  global acm_string = string(acm_string," \$ ",@sprintf("%6.3f",amed*1e-13*AU)," \\pm ",@sprintf("%6.3f",siga*1e-13*AU)," \$ ")
  if i < 7
     astring = string(astring," & ")
     acm_string = string(acm_string," & ")
  else
     astring = string(astring," \\cr ")
     acm_string = string(acm_string," \\cr ")
  end
end
println(astring)

println(rcm_string)
println(mgram_string)
println(rhocgs_string)
println(gcgs_string)
println(vesccgs_string)
println(scgs_string)
println(acm_string)
