PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

configfile: "./config/default.yaml"
include: "./rules/sync.smk"
include: "./rules/construct.smk"
include: "./rules/analyse.smk"


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/report.html",
        f"build/output/{config['resolution']['space']}/aggregation.nc",
        f"build/logs/{config['resolution']['space']}/test-report.html"


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--self-contained --css=report.css --template template.html --to html5"
    elif suffix == "pdf":
        return "--css=report.css --template template.html --pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule report:
    message: "Compile report.{wildcards.suffix}."
    input:
        "report/literature.bib",
        "report/report.md",
        "report/pandoc-metadata.yml",
        "report/energy-policy.csl",
        "report/template.html",
        "report/fonts/KlinicSlabBook.otf",
        "report/fonts/KlinicSlabBookIt.otf",
        "report/fonts/KlinicSlabMedium.otf",
        "report/fonts/KlinicSlabMediumIt.otf",
        "report/report.css"
    params: options = pandoc_options
    output: "build/report.{suffix}"
    wildcard_constraints:
        suffix = "((html)|(pdf)|(docx))"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} report.md pandoc-metadata.yml {params.options} \
        -o ../build/report.{wildcards.suffix}
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """

