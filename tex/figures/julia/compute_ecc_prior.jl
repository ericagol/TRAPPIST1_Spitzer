

# Plots the eccentricity prior:

include("/Users/ericagol/Observing/Spitzer/DDT2019/Campaign04/v09/Hyak/histogram_code.jl")
necc = 100000000
ecc = 0.1 .*rand(necc)
omega = 2pi .* rand(necc)

ecos = ecc .* cos.(omega)
esin = ecc .* sin.(omega)

ecos_bin,ecos_hist,ecos_bin_sq,ecos_hist_sq = histogram(ecos,333)
esin_bin,esin_hist,esin_bin_sq,esin_hist_sq = histogram(esin,333)
eboth_hist = 0.5*(ecos_hist_sq .+ esin_hist_sq)
ax.plot(ecos_bin_sq,eboth_hist./maximum(eboth_hist),label="Uniform eccentricity prior")
ax.axis([-0.0175,0.0175,0,1.05])
ax.legend()
