import pandas as pd

DISTRIBUTION_MAP = {
    "unif": "uniform",
    "norm": "normal"
}


def uncertainty_parameters_overview(parameters, path_to_output):
    pd.DataFrame(
        data={
            "Name": [p["descriptive-name"] for p in parameters.values()],
            "Description": [p["description"] for p in parameters.values()],
            "Distribution": [DISTRIBUTION_MAP[p["distribution"]] for p in parameters.values()],
            "Min/Mean": [p["bound1"] for p in parameters.values()],
            "Max/Std": [p["bound2"] for p in parameters.values()],
            "Source": [p["source"] for p in parameters.values()],
        }
    ).to_csv(path_to_output, index=False, header=True, float_format="%.3f")


if __name__ == "__main__":
    uncertainty_parameters_overview(
        parameters=snakemake.params.parameters,
        path_to_output=snakemake.output[0]
    )
