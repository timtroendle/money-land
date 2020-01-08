from dataclasses import dataclass

import numpy as np
import pandas as pd
import xarray as xr
from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.util import compute_groups_matrix


@dataclass
class CostFactors:
    wind_onshore: float
    wind_offshore: float
    roof_mounted_pv: float
    open_field_pv: float

    @classmethod
    def from_x(cls, x):
        return CostFactors(
            wind_onshore=x[0],
            wind_offshore=x[1],
            roof_mounted_pv=x[2],
            open_field_pv=x[3]
        )

    def apply_to_cost(self, cost):
        factors = xr.ones_like(cost)
        factors.loc[{"techs": ["wind_onshore_competing", "wind_onshore_monopoly"]}] = self.wind_onshore
        factors.loc[{"techs": "wind_offshore"}] = self.wind_offshore
        factors.loc[{"techs": "roof_mounted_pv"}] = self.roof_mounted_pv
        factors.loc[{"techs": "open_field_pv"}] = self.open_field_pv
        return cost * factors


def uncertainty_analysis(path_to_results, land_use_factors, path_to_xy, path_to_sobol):
    problem = { # FIXME inject the problem from config
        'num_vars': 4,
        'names': ['cost-wind-onshore',
                  'cost-wind-offshore',
                  'cost-roof-mounted-pv',
                  'cost-open-field-pv'],
        'bounds': [[800 / 1100, 1700 / 1100], # FIXME linear scaling imprecise due to discounting
                   [1790 / 2280, 3270 / 2280],
                   [760 / 880, 1000 / 880],
                   [280 / 520, 580 / 520]] # TODO link to roof mounted PV?
    }
    param_values = saltelli.sample(problem, 100)

    data = xr.open_dataset(path_to_results)
    y = np.zeros([param_values.shape[0]])
    for i, x in enumerate(param_values):
        y[i] = evaluate_model(data, land_use_factors, x)

    (
        pd
        .DataFrame(param_values, columns=problem["names"])
        .assign(land_use=y)
        .to_csv(path_to_xy, index=True, header=True)
    )

    indices = sobol.analyze(problem, y)

    with open(path_to_sobol, "w") as f_out:
        print_indices(indices, problem, True, f_out)


def evaluate_model(data, land_use_factors, x):
    cost_factors = CostFactors.from_x(x)
    data = read_data(data, land_use_factors, cost_factors)
    return data.land_use.loc[data.cost.idxmin()]


def read_data(data, land_use_factors, cost_factors):
    cost_data = (
        cost_factors
        .apply_to_cost(data.cost.sum("locs"))
        .sum("techs")
        .to_series()
    )
    land_use_data = ((
        data
        .energy_cap
        .sum("locs")
        .sel(techs=["wind_onshore_monopoly", "wind_onshore_competing", "wind_offshore",
                    "roof_mounted_pv", "open_field_pv"])
    ) * land_use_factors).sum("techs").to_series()

    both_data = pd.DataFrame({"land_use": land_use_data, "cost": cost_data})

    return both_data


def scenario_name_to_multiindex(index):
    return index.map(wxyz).rename(["util", "wind", "roof", "offshore"])


def wxyz(scenario_name):
    roof, util, wind, offshore = tuple(int(x.split("-")[1]) for x in scenario_name.split(","))
    return (util, wind, roof, offshore)


def print_indices(S, problem, calc_second_order, file):
    # taken from the SAlib source code and modified to print to file
    # https://github.com/SALib/SALib/blob/3bc2ddfb50f091e5e5dd1ed4e7fae05853f150e5/SALib/analyze/sobol.py#L240
    # Output to console
    if not problem.get('groups'):
        title = 'Parameter'
        names = problem['names']
        D = problem['num_vars']
    else:
        title = 'Group'
        _, names = compute_groups_matrix(problem['groups'])
        D = len(names)

    print('%s S1 S1_conf ST ST_conf' % title, file=file)

    for j in range(D):
        print('%s %f %f %f %f' % (names[j], S['S1'][
            j], S['S1_conf'][j], S['ST'][j], S['ST_conf'][j]), file=file)

    if calc_second_order:
        print('\n%s_1 %s_2 S2 S2_conf' % (title, title), file=file)

        for j in range(D):
            for k in range(j + 1, D):
                print("%s %s %f %f" % (names[j], names[k],
                      S['S2'][j, k], S['S2_conf'][j, k]), file=file)


if __name__ == "__main__":
    uncertainty_analysis(
        path_to_results=snakemake.input.results,
        land_use_factors=pd.Series(snakemake.params.land_factors).to_xarray().rename(index="techs"),
        path_to_xy=snakemake.output.xy,
        path_to_sobol=snakemake.output.sobol
    )
