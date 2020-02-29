subworkflow eurocalliope:
    workdir: "./euro-calliope"
    snakefile: "./euro-calliope/Snakefile"
    configfile: "./config/default.yaml"

localrules: copy_euro_calliope, copy_resolution_specific_euro_calliope, model
ruleorder: model > import_restrictions > copy_euro_calliope > copy_resolution_specific_euro_calliope
wildcard_constraints:
    definition_file = "[^\/]*" # must not travers into directories


rule copy_euro_calliope:
    message: "Copy file {input[0]} from euro-calliope."
    input: eurocalliope("build/model/{definition_file}.{suffix}"),
    output: "build/model/{definition_file}.{suffix}"
    shell: "ln {input} {output}"


rule copy_resolution_specific_euro_calliope:
    message: "Copy file {input[0]} from euro-calliope."
    input: eurocalliope("build/model/{resolution}/{definition_file}.{suffix}"),
    output: "build/model/{resolution}/{definition_file}.{suffix}"
    shell: "ln {input} {output}"


rule import_restrictions:
    message: "Create import restriction overrides for {wildcards.resolution} resolution."
    input:
        src = "src/construct/import.py",
        units = eurocalliope("build/data/{resolution}/units.csv")
    params:
        restrictions = [0],
        connected_regions = config["connected-regions"]
    conda: "../envs/default.yaml"
    output: "build/model/{resolution}/import-restrictions.yaml"
    script: "../src/construct/import.py"


rule model:
    message: "Build entire model on resolution {wildcards.resolution}."
    input:
        "build/model/interest-rate.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/storage-techs.yaml",
        "build/model/link-techs.yaml",
        "build/model/{resolution}/locations.yaml",
        "build/model/{resolution}/link-all-neighbours.yaml",
        "build/model/{resolution}/electricity-demand.csv",
        "build/model/{resolution}/directional-rooftop.yaml",
        rules.import_restrictions.output,
        expand(
            "build/model/{{resolution}}/capacityfactors-{technology}.csv",
            technology=["open-field-pv", "rooftop-pv-n", "rooftop-pv-e-w", "rooftop-pv-s-flat",
                        "wind-offshore", "wind-onshore", "hydro-ror", "hydro-reservoir-inflow"],
        ),
        definition = "src/template/model.yaml"
    output:
        model = "build/model/{resolution}/model.yaml"
    shell:
        "cp {input.definition} {output}"
