# Supply-side options to reduce land requirements of fully renewable electricity in Europe

Trade-off between land requirements and cost of fully renewable electricity in Europe.

This repository contains the entire scientific project, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

[![article DOI](https://img.shields.io/badge/article-10.1371/journal.pone.0236958-blue)](https://doi.org/10.1371/journal.pone.0236958)
[![data DOI](https://img.shields.io/badge/data-10.5281%2Fzenodo.3707812-blue)](https://doi.org/10.5281/zenodo.3707812)
[![code DOI](https://img.shields.io/badge/code-10.5281%2Fzenodo.3956531-blue)](https://doi.org/10.5281/zenodo.3956531)

This study uses the model of the European electricity system [euro-calliope v1.0](https://zenodo.org/record/3949794).

## Getting ready

1. Clone the repo. Because `euro-calliope` is added as a git submodule, you may want to clone using `git clone --recurse-submodules <link-to-this-repo>`.

2. Create an environment to run the analysis. You need [conda](https://conda.io/docs/index.html) to run the analysis. Using conda, you can create a conda environment from within you can run it:

    `conda env create -f environment.yaml`

3. Make sure you have a Gurobi license, or install and configure another solver.

4. You need an account at the Copernicus Climate Data Service and you need to create a `$HOME/.cdsapirc` file with your credentials, see their [How To](https://cds.climate.copernicus.eu/api-how-to) (you do not need to manually install the client).

5. Provide the input data for Euro-Calliope, as defined in "Getting Ready" in  `./euro-calliope/README.md`.

## Run the analysis

    snakemake  --use-conda

This will run all analysis steps to reproduce results and eventually build the report.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps, and if you have `dot` installed, run:

    snakemake --rulegraph | dot -Tpdf > dag.pdf

## Run on Euler cluster

To run on Euler, use the following command:

    snakemake --use-conda --profile config/euler [--config email=<you@provider.org>]

By providing an email address, you will be informed by mail when Snakemake finishes execution.

If you prefer working locally, you can sync this repository to Euler and receive build changes by running `snakemake send` and `snakemake receive`.

If you want to run on another cluster, read [snakemake's documentation on cluster execution](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) and take `config/euler` as a starting point.

## Run the tests

    snakemake test --use-conda

## Be notified of successes or fails

As the execution of this workflow may take long, you can get notified whenever the execution terminates either successfully or unsuccessfully. Notification are handled by the webservice [Pushcut](https://pushcut.io/) for which you need a free account. To activate notifications, add your Pushcut secret to the configuration using the configuration key `pushcut_secret`. For example, you may want to choose running the workflow the following way to receive notifications:

    snakemake --use-conda --config pushcut_secret=<your-secret>

This workflow will then trigger the Pushcut notifications `snakemake_succeeded` and `snakemake_failed`.

## Units and scaling

The default units for euro-calliope are `MW`, `MWh`, `EUR`, and `km2`, but you can scale all of these using the configuration values in `config/default.yaml`. Apart from convenience, this may be important to handle numerical issues with your solver.

## Repo structure

* `report`: contains all files necessary to build the report; plots and result files are not in here but generated automatically
* `rules`: contains Snakemake rule definition files
* `src`: contains the Python source code
* `envs`: contains execution environments
* `tests`: contains the test code
* `config`: configurations used in the study
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)

## License

The code in this repo is MIT licensed, see `./LICENSE.md`. This excludes the KlinicSlab font family (all files in `./report/fonts/`) which is copyright Lost Type.
