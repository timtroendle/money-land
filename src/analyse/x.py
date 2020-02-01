import pandas as pd
from SALib.sample import saltelli
import scipy as sp
from scipy import stats
import SALib.util


def sample_x(uncertain_parameters, number_uncertainty_runs, path_to_x):
    """Create samples from input uncertainty."""
    monkey_patch_salib()
    problem = {
        'num_vars': len(uncertain_parameters.keys()),
        'names': [param_name for param_name in uncertain_parameters.keys()],
        'bounds': [(attrs["bound1"], attrs["bound2"]) for attrs in uncertain_parameters.values()],
        'dists': [param_attributes["distribution"] for param_attributes in uncertain_parameters.values()]
    }
    param_values = saltelli.sample(problem, number_uncertainty_runs, calc_second_order=False)

    (
        pd
        .DataFrame(param_values, columns=problem["names"])
        .to_csv(path_to_x, index=True, header=True)
    )


def monkey_patch_salib():
    # SALib throws an Error when using normal distribution, because scipy.stats isn't
    # properly imported. I am monkey patching this here.
    sp.stats = stats
    SALib.util.sp = sp


if __name__ == "__main__":
    sample_x(
        uncertain_parameters=snakemake.params.uncertain_parameters,
        number_uncertainty_runs=snakemake.params.runs,
        path_to_x=snakemake.output.x,
    )
