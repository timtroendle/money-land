RESAMPLE_OVERRIDE_DICT = {
    "model.time.function": "resample",
    "model.time.function_options": {"resolution": config["resolution"]["time"]}
}
CALLIOPE_OVERRIDE_DICT = {
    **config["calliope-parameters"],
    **RESAMPLE_OVERRIDE_DICT
}

rule run_continental:
    message: "Run the model on continental resolution."
    input:
        model = "build/model/continental/model.yaml"
    params: override_dict = CALLIOPE_OVERRIDE_DICT
    output: "build/output/continental/results.nc"
    conda: "../envs/default.yaml"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} \
            --override_dict="{params.override_dict}"
        """


rule time_aggregated_results:
    message: "Create NetCDF overview over all results."
    input:
        src = "src/analyse/aggregation.py",
        scenarios = ["build/output/{resolution}/results.nc"],
        units = eurocalliope("build/data/{resolution}/units.csv")
    params: scaling_factors = config["scaling-factors"]
    conda: "../envs/default.yaml"
    output: "build/output/{resolution}/aggregation.nc"
    script: "../src/analyse/aggregation.py"


rule test:
    message: "Run tests"
    input:
        "src/analyse/test_runner.py",
        "tests/test_feasibility.py",
        results = ["build/output/{resolution}/results.nc"],
        biofuel_potentials = eurocalliope("build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv".format(
            scenario=config["parameters"]["jrc-biofuel"]["scenario"]
        )),
        units = eurocalliope("build/data/{resolution}/units.csv")
    params: scaling_factors = config["scaling-factors"]
    output: "build/logs/{resolution}/test-report.html"
    conda: "../envs/test.yaml"
    script: "../src/analyse/test_runner.py"
