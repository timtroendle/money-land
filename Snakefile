PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

configfile: "./config/default.yaml"
include: "./rules/sync.smk"
include: "./rules/construct.smk"
include: "./rules/analyse.smk"
localrules: all, clean, report
onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'money-land succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'money-land crashed' {config[email]}")
wildcard_constraints:
    resolution = "((continental))", # supported spatial resolutions
    plot_suffix = "((png)|(svg))"


rule all:
    message: "Run entire analysis and compile report."
    input:
        f"build/output/{config['resolution']['space']}/report.html",
        f"build/logs/{config['resolution']['space']}/test-report.html"


rule copy_report_file:
    message: "Copy file {input[0]} into dedicated report folder."
    input: "build/output/{resolution}/{filename}.{suffix}"
    wildcard_constraints: suffix = "((csv)|(png)|(svg))"
    output: "build/output/{resolution}/report/{filename}.{suffix}"
    shell: "ln {input} {output}"


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
        "report/report.css",
        "build/output/{resolution}/report/ternary.svg"
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
        {PANDOC} report.md pandoc-metadata.yml {params.options} \
        -o ../build/output/{wildcards.resolution}/report.{wildcards.suffix}
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """

