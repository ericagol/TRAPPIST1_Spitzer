import matplotlib.pyplot as plt

figures = [
    "figure_water_period_TTV_trim.png",
    "T1_JWST_all_transits_tight.png",
    "tlM.png",
    "Equivalent_RV_precision.pdf",
    "Laplace_angle.pdf",
    "Norm_cmf_vs_period.pdf",
    "Norm_dens_vs_period.pdf",
    "Planet_i_properties.pdf",
    "Recovered_masses_JWST_5yr_all_transits_NIRSPEC.pdf",
    "Students_t_optimized.pdf",
    "T1_chopping.pdf",
    "T1_eccentricity_vectors_likelihood_profile_hmc.pdf",
    "T1_evector.pdf",
    "T1_evector_forced.pdf",
    "T1_masses_03312020.pdf",
    "T1_riverplot.pdf",
    "T1_students_params_transformed.pdf",
    "T1_ttvs_4panel_stacked.pdf",
    "delta_omega_bc.pdf",
    "eccentricity_posterior.pdf",
    "esin_vs_ecos.pdf",
    "impact_parameter_noprior.pdf",
    "inclination_noprior.pdf",
    "inclination_nouprior_incprior.pdf",
    "limb_darkening_nouprior.pdf",
    "mass_radius_relation_comparison.pdf",
    "plot_normalized_density_agol_et_al.pdf",
    "stellar_density_noprior.pdf",
]

for figure in figures:
    fig, ax = plt.subplots(1)
    ax.axis("off")
    ax.annotate(
        "Figure failed to compile.",
        xy=(0.5, 0.5),
        xycoords="data",
        va="center",
        ha="center",
        fontsize=30,
        color="red",
    )
    fig.savefig(figure, bbox_inches="tight")
    plt.close()
