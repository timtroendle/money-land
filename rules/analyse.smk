RESAMPLE_OVERRIDE_DICT = {
    "model.time.function": "resample",
    "model.time.function_options": {"resolution": config["resolution"]["time"]}
}
CALLIOPE_OVERRIDE_DICT = {
    **config["calliope-parameters"],
    **RESAMPLE_OVERRIDE_DICT
}
CAPACITY_SHARE_SAMPLES = list(range(0, 110, 10))
ALL_SCENARIOS = (
    [
        f"build/output/{{resolution}}/runs/roof-{roof}-util-{util}-wind-{wind}-offshore-{offshore}.nc"
        for roof in CAPACITY_SHARE_SAMPLES
        for util in CAPACITY_SHARE_SAMPLES
        for wind in CAPACITY_SHARE_SAMPLES
        for offshore in CAPACITY_SHARE_SAMPLES
        if roof + util + wind + offshore == 100
    ]
)
wildcard_constraints:
    case = "((roof)|(offshore))",
    land = "((land-use)|(footprint-only))"


rule run_continental:
    message: "Run the model on continental resolution with {wildcards.roof}/{wildcards.util}/{wildcards.wind}/{wildcards.offshore}."
    input:
        src = "src/analyse/run.py",
        model = "build/model/continental/model.yaml"
    params:
        override_dict = CALLIOPE_OVERRIDE_DICT,
        resolution = "continental",
        no_shore = config["units-without-shore"]["continental"]
    output: "build/output/continental/runs/roof-{roof}-util-{util}-wind-{wind}-offshore-{offshore}.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/run.py"


rule run_national:
    message: "Run the model on national resolution with {wildcards.roof}/{wildcards.util}/{wildcards.wind}/{wildcards.offshore}."
    input:
        src = "src/analyse/run.py",
        model = "build/model/national/model.yaml"
    params:
        override_dict = CALLIOPE_OVERRIDE_DICT,
        resolution = "national",
        no_shore = config["units-without-shore"]["national"]
    output: "build/output/national/runs/roof-{roof}-util-{util}-wind-{wind}-offshore-{offshore}.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/run.py"


rule aggregated_results:
    message: "Create NetCDF overview over all results."
    input:
        src = "src/analyse/aggregation.py",
        scenarios = ALL_SCENARIOS,
        units = eurocalliope("build/data/{resolution}/units.csv")
    params: scaling_factors = config["scaling-factors"]
    conda: "../envs/calliope.yaml"
    output: "build/output/{resolution}/aggregation.nc"
    script: "../src/analyse/aggregation.py"


rule visualise_supply_shares:
    message: "Visualise supply share space."
    input:
        src = "src/analyse/supply_shares.py",
        results = rules.aggregated_results.output[0]
    output: "build/output/{resolution}/supply-shares.{plot_suffix}"
    conda: "../envs/default.yaml"
    script: "../src/analyse/supply_shares.py"


rule sample_x:
    message: "Create samples from input uncertainty."
    input:
        src = "src/analyse/x.py",
    params:
        runs = config["uncertainty"]["number-runs"],
        uncertain_parameters = config["uncertainty"]["parameters"]
    conda: "../envs/default.yaml"
    output:
        x = "build/output/{resolution}/{land}/x.csv",
    script: "../src/analyse/x.py"


rule xy:
    message: "Create samples from output uncertainty."
    input:
        src  = "src/analyse/y.py",
        x = rules.sample_x.output[0],
        model_output = rules.aggregated_results.output[0]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/{land}/xy.nc"
    script: "../src/analyse/y.py"


rule plot_observations:
    message: "Create overview plot of all observations."
    input:
        src = "src/analyse/observations.py",
        xy = rules.xy.output[0]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/{land}/observations.{plot_suffix}"
    script: "../src/analyse/observations.py"


rule ternary_plots:
    message: "Create ternary plots of results."
    input:
        src = "src/analyse/ternary.py",
        results = rules.xy.output[0]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/{land}/ternary.{plot_suffix}"
    script: "../src/analyse/ternary.py"


rule scatter_plot:
    message: "Create scatter plots of results."
    input:
        src = "src/analyse/scatter.py",
        results = rules.aggregated_results.output[0]
    params: land_factors = lambda wildcards: config["parameters"][wildcards["land"]]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/{land}/scatter-{case}.{plot_suffix}"
    script: "../src/analyse/scatter.py"


rule flexibility_plot:
    message: "Create plot of flexibility needs."
    input:
        src = "src/analyse/flexibility.py",
        results = rules.aggregated_results.output[0]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/flexibility.{plot_suffix}"
    script: "../src/analyse/flexibility.py"


rule technology_plot:
    message: "Create plot for single technologies."
    input:
        src = "src/analyse/technology_plot.py",
        results = rules.xy.output[0]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/{land}/technology.{plot_suffix}"
    script: "../src/analyse/technology_plot.py"


rule relative_boxenplot_uncertainty:
    message: "Create a boxenplot of uncertainty of relative land reduction."
    input:
        src = "src/analyse/boxenplot_relative.py",
        xy = rules.xy.output[0]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/{land}/boxenplot-relative.{plot_suffix}"
    script: "../src/analyse/boxenplot_relative.py"


rule absolute_boxenplot_uncertainty:
    message: "Create a boxenplot of uncertainty of absolute land reduction."
    input:
        src = "src/analyse/boxenplot_absolute.py",
        xy = rules.xy.output[0]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/{land}/boxenplot-absolute.{plot_suffix}"
    script: "../src/analyse/boxenplot_absolute.py"


rule wind_share_plot:
    message: "Create plot of max wind share."
    input:
        src = "src/analyse/wind.py",
        xy = rules.xy.output[0]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/{land}/wind.{plot_suffix}"
    script: "../src/analyse/wind.py"


rule overview_cost_assumptions:
    message: "Create table of key cost assumptions."
    input:
        src = "src/analyse/cost_assumptions.py",
        model = "build/output/{resolution}/runs/roof-0-util-0-wind-0-offshore-100.nc"
    params: scaling_factors = config["scaling-factors"]
    output: "build/output/{resolution}/overview-cost-assumptions.csv"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/cost_assumptions.py"


rule overview_uncertainty_parameters:
    message: "Create table of uncertainty parameters."
    input: src = "src/analyse/uncertain_parameters.py"
    params: parameters = config["uncertainty"]["parameters"]
    output: "build/output/{resolution}/overview-uncertain-parameters.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/uncertain_parameters.py"


rule test:
    message: "Run tests"
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        "tests/test_supply_shares.py",
        results = ALL_SCENARIOS,
        biofuel_potentials = eurocalliope("build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv".format(
            scenario=config["parameters"]["jrc-biofuel"]["scenario"]
        )),
        units = eurocalliope("build/data/{resolution}/units.csv")
    params: scaling_factors = config["scaling-factors"]
    output: "build/logs/{resolution}/test-report.html"
    conda: "../envs/test.yaml"
    script: "../src/analyse/test_runner.py"
