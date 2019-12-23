RESAMPLE_OVERRIDE_DICT = {
    "model.time.function": "resample",
    "model.time.function_options": {"resolution": config["resolution"]["time"]}
}
CALLIOPE_OVERRIDE_DICT = {
    **config["calliope-parameters"],
    **RESAMPLE_OVERRIDE_DICT
}
CAPACITY_SHARE_SAMPLES = list(range(0, 110, 10))


rule run_continental:
    message: "Run the model on continental resolution with {wildcards.roof}/{wildcards.util}/{wildcards.wind}."
    input:
        src = "src/analyse/run.py",
        model = "build/model/continental/model.yaml"
    params:
        override_dict = CALLIOPE_OVERRIDE_DICT,
        resolution = "continental"
    output: "build/output/continental/runs/roof-{roof}-util-{util}-wind-{wind}.nc"
    conda: "../envs/default.yaml"
    script: "../src/analyse/run.py"


rule run_national:
    message: "Run the model on national resolution with {wildcards.roof}/{wildcards.util}/{wildcards.wind}."
    input:
        src = "src/analyse/run.py",
        model = "build/model/national/model.yaml"
    params:
        override_dict = CALLIOPE_OVERRIDE_DICT,
        resolution = "national"
    output: "build/output/national/runs/roof-{roof}-util-{util}-wind-{wind}.nc"
    conda: "../envs/default.yaml"
    script: "../src/analyse/run.py"


rule time_aggregated_results:
    message: "Create NetCDF overview over all results."
    input:
        src = "src/analyse/aggregation.py",
        scenarios = [f"build/output/{{resolution}}/runs/roof-{roof}-util-{util}-wind-{wind}.nc"
                     for roof in CAPACITY_SHARE_SAMPLES
                     for util in CAPACITY_SHARE_SAMPLES
                     for wind in CAPACITY_SHARE_SAMPLES
                     if roof + util + wind == 100],
        units = eurocalliope("build/data/{resolution}/units.csv")
    params: scaling_factors = config["scaling-factors"]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/aggregation.nc"
    script: "../src/analyse/aggregation.py"


rule ternary_plots:
    message: "Create ternary plots of results."
    input:
        src = "src/analyse/ternary.py",
        results = rules.time_aggregated_results.output[0]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/ternary.{plot_suffix}"
    script: "../src/analyse/ternary.py"


rule test:
    message: "Run tests"
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        results = [f"build/output/{{resolution}}/runs/roof-{roof}-util-{util}-wind-{wind}.nc"
                     for roof in CAPACITY_SHARE_SAMPLES
                     for util in CAPACITY_SHARE_SAMPLES
                     for wind in CAPACITY_SHARE_SAMPLES
                     if roof + util + wind == 100],
        biofuel_potentials = eurocalliope("build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv".format(
            scenario=config["parameters"]["jrc-biofuel"]["scenario"]
        )),
        units = eurocalliope("build/data/{resolution}/units.csv")
    params: scaling_factors = config["scaling-factors"]
    output: "build/logs/{resolution}/test-report.html"
    conda: "../envs/test.yaml"
    script: "../src/analyse/test_runner.py"
