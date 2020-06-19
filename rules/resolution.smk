# Some rules to analyse the impact of the temporal resolution.

configfile: "./config/default.yaml"
localrules: all
onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'money-land succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'money-land crashed' {config[email]}")

rule all:
    input:
        "build/output/national/runs/restest/roof-0-util-70-wind-30-offshore-0-4h.nc",
        "build/output/national/runs/restest/roof-0-util-70-wind-30-offshore-0-3h.nc",
        "build/output/national/runs/restest/roof-0-util-70-wind-30-offshore-0-2h.nc",
        "build/output/national/runs/restest/roof-0-util-70-wind-30-offshore-0-1h.nc"


def calliope_override(wildcards):
    resample_override_dict = {
        "model.time.function": "resample",
        "model.time.function_options": {"resolution": str(wildcards["resolution"]) + "H"}
    }
    return {
        **config["calliope-parameters"],
        **resample_override_dict
    }


rule run_national_test_temporalresolution:
    message: "Run the model with {wildcards.resolution}h resolution with {wildcards.roof}/{wildcards.util}/{wildcards.wind}/{wildcards.offshore}."
    input:
        src = "src/analyse/run.py",
        model = "build/model/national/model.yaml"
    params:
        override_dict = calliope_override,
        no_shore = config["units-without-shore"]["national"],
        overrides = ["no-hydro-costs", "stylised-storage", "directional-rooftop-pv",
                     "national-autarky-100-percent", "freeze-hydro-capacities"]
    output: "build/output/national/runs/restest/roof-{roof}-util-{util}-wind-{wind}-offshore-{offshore}-{resolution}h.nc"
    conda: "../envs/calliope.yaml"
    script: "../src/analyse/run.py"
