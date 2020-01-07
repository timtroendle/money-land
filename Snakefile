PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

configfile: "./config/default.yaml"
include: "./rules/sync.smk"
include: "./rules/construct.smk"
include: "./rules/analyse.smk"
localrules: all, clean, copy_report_file, report, supplementary_material
onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'money-land succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'money-land crashed' {config[email]}")
wildcard_constraints:
    resolution = "((continental)|(national))", # supported spatial resolutions
    plot_suffix = "((png)|(svg))"


rule all:
    message: "Run entire analysis and compile report."
    input:
        f"build/output/{config['resolution']['space']}/report.html",
        f"build/output/{config['resolution']['space']}/supplementary.html",
        f"build/logs/{config['resolution']['space']}/test-report.html",
        f"build/output/{config['resolution']['space']}/land-use/technology-stats.csv"


rule copy_report_file:
    message: "Copy file {input[0]} into dedicated report folder."
    input: "build/output/{resolution}/{filename}.{suffix}"
    wildcard_constraints: suffix = "((csv)|(png)|(svg))"
    output: "build/output/{resolution}/report/{filename}.{suffix}"
    shell: "ln {input} {output}"


GENERAL_DOCUMENT_DEPENDENCIES = [
    "report/literature.bib",
    "report/report.css",
    "report/plos-one.csl",
    "report/template.html",
    "report/fonts/KlinicSlabBook.otf",
    "report/fonts/KlinicSlabBookIt.otf",
    "report/fonts/KlinicSlabMedium.otf",
    "report/fonts/KlinicSlabMediumIt.otf",
]


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
        GENERAL_DOCUMENT_DEPENDENCIES,
        "report/report.md",
        "report/pandoc-metadata.yml",
        "build/output/{resolution}/report/map.png",
        "build/output/{resolution}/report/land-use/ternary-roof.svg",
        "build/output/{resolution}/report/land-use/scatter-roof.svg",
        "build/output/{resolution}/report/land-use/ternary-offshore.svg",
        "build/output/{resolution}/report/land-use/scatter-offshore.svg",
        "build/output/{resolution}/report/flexibility.svg"
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


rule supplementary_material:
    message: "Compile the supplementary material."
    input:
        GENERAL_DOCUMENT_DEPENDENCIES,
        "report/supplementary.md",
        "build/output/{resolution}/report/footprint-only/ternary-roof.svg",
        "build/output/{resolution}/report/footprint-only/scatter-roof.svg",
        "build/output/{resolution}/report/footprint-only/ternary-offshore.svg",
        "build/output/{resolution}/report/footprint-only/scatter-offshore.svg"
    params: options = pandoc_options
    output: "build/output/{resolution}/supplementary.{suffix}"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build/output/{wildcards.resolution}/report .
        {PANDOC} supplementary.md {params.options} --table-of-contents \
        -o ../build/output/{wildcards.resolution}/supplementary.{wildcards.suffix}
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """

