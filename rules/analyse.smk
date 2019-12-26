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
        f"build/output/{{resolution}}/runs/roof-{roof}-util-{util}-wind-{wind}-offshore-0.nc"
        for roof in CAPACITY_SHARE_SAMPLES
        for util in CAPACITY_SHARE_SAMPLES
        for wind in CAPACITY_SHARE_SAMPLES
        if roof + util + wind == 100
    ] +
    [
        f"build/output/{{resolution}}/runs/roof-0-util-{util}-wind-{wind}-offshore-{offshore}.nc"
        for offshore in CAPACITY_SHARE_SAMPLES
        for util in CAPACITY_SHARE_SAMPLES
        for wind in CAPACITY_SHARE_SAMPLES
        if offshore + util + wind == 100
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
    conda: "../envs/default.yaml"
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
    conda: "../envs/default.yaml"
    script: "../src/analyse/run.py"


rule aggregated_results:
    message: "Create NetCDF overview over all results."
    input:
        src = "src/analyse/aggregation.py",
        scenarios = ALL_SCENARIOS,
        units = eurocalliope("build/data/{resolution}/units.csv")
    params: scaling_factors = config["scaling-factors"]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/aggregation.nc"
    script: "../src/analyse/aggregation.py"


rule ternary_plots:
    message: "Create ternary plots of results."
    input:
        src = "src/analyse/ternary.py",
        results = rules.aggregated_results.output[0]
    params: land_factors = lambda wildcards: config["parameters"][wildcards["land"]]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/{land}/ternary-{case}.{plot_suffix}"
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


rule test:
    message: "Run tests"
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        results = ALL_SCENARIOS,
        biofuel_potentials = eurocalliope("build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv".format(
            scenario=config["parameters"]["jrc-biofuel"]["scenario"]
        )),
        units = eurocalliope("build/data/{resolution}/units.csv")
    params: scaling_factors = config["scaling-factors"]
    output: "build/logs/{resolution}/test-report.html"
    conda: "../envs/test.yaml"
    script: "../src/analyse/test_runner.py"
