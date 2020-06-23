# Money Land

Trade-off between land requirements and cost of renewable electricity in Europe.

This repository contains the entire scientific project, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

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

As the execution of this workflow may take long, you can get notified whenever the execution terminates either successfully or unsuccessfully. Notification are handled by the webservice [IFTTT](https://ifttt.com/) for which you need a free account. To activate notifications, add your IFTTT webhooks apikey to the configuration using the configuration key `ifttt_apikey`. For example, you may want to choose running the workflow the following way to receive notifications:

    snakemake --use-conda --config ifttt_apikey=<your-apikey>

This workflow will then trigger the IFTTT events `snakemake_succeeded` and `snakemake_failed`. You can build IFTTT applets that react to these events, and for example let the notification be pushed to your phone or let your internet-connected light bulb flash red or green. [Here's an example](https://pimylifeup.com/using-ifttt-with-the-raspberry-pi/) that you can use as a starting point.

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
