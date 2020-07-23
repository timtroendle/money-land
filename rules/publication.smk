localrules: figures, prepare_publication


rule figures:
    message: "Collect and rename all figures."
    input:
        "build/output/{resolution}/supply-shares.{plot_suffix}",
        "build/output/{resolution}/land-use/observations.{plot_suffix}",
        "build/output/{resolution}/land-use/ternary.{plot_suffix}",
        "build/output/{resolution}/land-use/technology.{plot_suffix}",
        "build/output/{resolution}/land-use/boxenplot-absolute.{plot_suffix}",
        "build/output/{resolution}/land-use/wind.{plot_suffix}",
        "build/output/{resolution}/flexibility.{plot_suffix}"
    output:
        "build/output/{resolution}/figures/Fig1.{plot_suffix}",
        "build/output/{resolution}/figures/Fig2.{plot_suffix}",
        "build/output/{resolution}/figures/Fig3.{plot_suffix}",
        "build/output/{resolution}/figures/Fig4.{plot_suffix}",
        "build/output/{resolution}/figures/Fig5.{plot_suffix}",
        "build/output/{resolution}/figures/Fig6.{plot_suffix}",
        "build/output/{resolution}/figures/Fig7.{plot_suffix}"
    run:
        from shutil import copyfile
        for i in range(len(output)):
            copyfile(input[i], output[i])


rule prepare_publication:
    message: "Prepare results for publication."
    input:
        "build/model/build-metadata.yaml",
        "build/model/demand-techs.yaml",
        "build/model/interest-rate.yaml",
        "build/model/link-techs.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/storage-techs.yaml",
        "build/model/national/capacityfactors-hydro-reservoir-inflow.csv",
        "build/model/national/capacityfactors-hydro-ror.csv",
        "build/model/national/capacityfactors-open-field-pv.csv",
        "build/model/national/capacityfactors-rooftop-pv-e-w.csv",
        "build/model/national/capacityfactors-rooftop-pv-n.csv",
        "build/model/national/capacityfactors-rooftop-pv-s-flat.csv",
        "build/model/national/capacityfactors-wind-offshore.csv",
        "build/model/national/capacityfactors-wind-onshore.csv",
        "build/model/national/directional-rooftop.yaml",
        "build/model/national/electricity-demand.csv",
        "build/model/national/import-restrictions.yaml",
        "build/model/national/link-all-neighbours.yaml",
        "build/model/national/locations.yaml",
        "build/model/national/model.yaml",
        "build/output/national/aggregation.nc",
        "build/output/national/land-use/xy.nc"
    output:
        "build/publication/model/build-metadata.yaml",
        "build/publication/model/demand-techs.yaml",
        "build/publication/model/interest-rate.yaml",
        "build/publication/model/link-techs.yaml",
        "build/publication/model/renewable-techs.yaml",
        "build/publication/model/storage-techs.yaml",
        "build/publication/model/national/capacityfactors-hydro-reservoir-inflow.csv",
        "build/publication/model/national/capacityfactors-hydro-ror.csv",
        "build/publication/model/national/capacityfactors-open-field-pv.csv",
        "build/publication/model/national/capacityfactors-rooftop-pv-e-w.csv",
        "build/publication/model/national/capacityfactors-rooftop-pv-n.csv",
        "build/publication/model/national/capacityfactors-rooftop-pv-s-flat.csv",
        "build/publication/model/national/capacityfactors-wind-offshore.csv",
        "build/publication/model/national/capacityfactors-wind-onshore.csv",
        "build/publication/model/national/directional-rooftop.yaml",
        "build/publication/model/national/electricity-demand.csv",
        "build/publication/model/national/import-restrictions.yaml",
        "build/publication/model/national/link-all-neighbours.yaml",
        "build/publication/model/national/locations.yaml",
        "build/publication/model/national/model.yaml",
        "build/publication/system-designs.nc",
        "build/publication/samples.nc"
    run:
        from shutil import copyfile
        for i in range(len(output)):
            copyfile(input[i], output[i])
