PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

configfile: "./config/default.yaml"
include: "./rules/sync.smk"
include: "./rules/construct.smk"
include: "./rules/analyse.smk"
localrules: all, clean, copy_report_file, report, supplementary_material
onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "pushcut_secret" in config.keys():
        trigger_pushcut(event_name="snakemake_succeeded", secret=config["pushcut_secret"])
onerror:
    if "pushcut_secret" in config.keys():
        trigger_pushcut(event_name="snakemake_failed", secret=config["pushcut_secret"])
wildcard_constraints:
    resolution = "((continental)|(national))", # supported spatial resolutions
    plot_suffix = "((png)|(svg)|(tiff))"


rule all:
    message: "Run entire analysis and compile report."
    input:
        f"build/output/{config['resolution']['space']}/report.pdf",
        f"build/output/{config['resolution']['space']}/report.docx",
        f"build/output/{config['resolution']['space']}/supplementary.pdf",
        f"build/logs/{config['resolution']['space']}/test-report.html",
        f"build/output/{config['resolution']['space']}/figures/Fig1.tiff",
        f"build/output/{config['resolution']['space']}/figures/Fig2.tiff",
        f"build/output/{config['resolution']['space']}/figures/Fig3.tiff",
        f"build/output/{config['resolution']['space']}/figures/Fig4.tiff",
        f"build/output/{config['resolution']['space']}/figures/Fig5.tiff",
        f"build/output/{config['resolution']['space']}/figures/Fig6.tiff",
        f"build/output/{config['resolution']['space']}/figures/Fig7.tiff"


rule copy_report_file:
    message: "Copy file {input[0]} into dedicated report folder."
    input: "build/output/{resolution}/{filename}.{suffix}"
    wildcard_constraints: suffix = "((csv)|(png)|(svg)|(tiff))"
    output: "build/output/{resolution}/report/{filename}.{suffix}"
    shell: "ln {input} {output}"


GENERAL_DOCUMENT_DEPENDENCIES = [
    "report/literature.bib",
    "report/report.css",
    "report/plos-one.csl",
    "report/template.html",
    "report/pandoc-metadata.yaml",
    "report/fonts/KlinicSlabBook.otf",
    "report/fonts/KlinicSlabBookIt.otf",
    "report/fonts/KlinicSlabMedium.otf",
    "report/fonts/KlinicSlabMediumIt.otf",
]


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--self-contained --template template.html --to html5"
    elif suffix == "pdf":
        return "--template template.html --pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule report:
    message: "Compile report.{wildcards.suffix}."
    input:
        GENERAL_DOCUMENT_DEPENDENCIES,
        "report/report.md",
        "build/output/{resolution}/report/supply-shares.png",
        "build/output/{resolution}/report/land-use/observations.png",
        "build/output/{resolution}/report/land-use/ternary.png",
        "build/output/{resolution}/report/land-use/technology.png",
        "build/output/{resolution}/report/land-use/boxenplot-absolute.png",
        "build/output/{resolution}/report/land-use/wind.png",
        "build/output/{resolution}/report/flexibility.png",
        "build/output/{resolution}/report/overview-cost-assumptions.csv",
        "build/output/{resolution}/report/overview-uncertain-parameters.csv",
    params: options = pandoc_options
    output: "build/output/{resolution}/report.{suffix}"
    wildcard_constraints:
        suffix = "((html)|(pdf)|(docx))"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build/output/{wildcards.resolution}/report .
        {PANDOC} report.md {params.options} \
        --metadata-file=pandoc-metadata.yaml \
        -o ../build/output/{wildcards.resolution}/report.{wildcards.suffix}
        """


rule supplementary_material:
    message: "Compile the supplementary material."
    input:
        GENERAL_DOCUMENT_DEPENDENCIES,
        "report/supplementary.md",
        "build/output/{resolution}/report/land-use/map-land-requirements.png",
        "build/output/{resolution}/report/footprint-only/ternary.png",
    params: options = pandoc_options
    output: "build/output/{resolution}/supplementary.{suffix}"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build/output/{wildcards.resolution}/report .
        {PANDOC} supplementary.md {params.options} --table-of-contents --number-sections \
        --metadata-file=pandoc-metadata.yaml \
        -o ../build/output/{wildcards.resolution}/supplementary.{wildcards.suffix}
        """


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

rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """


def trigger_pushcut(event_name, secret):
    import requests
    response = requests.post(
            f'https://api.pushcut.io/{secret}/notifications/{event_name}'
    )
