from dataclasses import dataclass
import functools

import numpy as np
import pandas as pd
import xarray as xr
from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.util import compute_groups_matrix

ALL_TECHS = ["roof", "util", "wind", "offshore"]


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


@dataclass
class LandUseFactors:
    wind_onshore: float
    wind_offshore: float
    roof_mounted_pv: float
    open_field_pv: float

    def apply_to_energy_cap(self, energy_cap):
        factors = xr.zeros_like(energy_cap)
        factors.loc[{"techs": ["wind_onshore_competing", "wind_onshore_monopoly"]}] = self.wind_onshore
        factors.loc[{"techs": "wind_offshore"}] = self.wind_offshore
        factors.loc[{"techs": "roof_mounted_pv"}] = self.roof_mounted_pv
        factors.loc[{"techs": "open_field_pv"}] = self.open_field_pv
        return energy_cap * factors


def uncertainty_analysis(path_to_results, land_use_factors, number_uncertainty_runs, path_to_xy, path_to_sobol):
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
    param_values = saltelli.sample(problem, number_uncertainty_runs)

    data = xr.open_dataset(path_to_results)
    ys = np.zeros([param_values.shape[0], 10])
    for i, x in enumerate(param_values):
        ys[i] = evaluate_model(data, land_use_factors, x)

    (
        pd
        .DataFrame(param_values, columns=problem["names"])
        .assign(
            optimal_cost=ys.T[0],
            optimal_land_use=ys.T[1],
            optimal_util=ys.T[2],
            optimal_wind=ys.T[3],
            optimal_roof=ys.T[4],
            optimal_offshore=ys.T[5],
            r50_cost_util=ys.T[6],
            r50_cost_offshore=ys.T[7],
            r50_cost_roof=ys.T[8],
            r50_cost_wind=ys.T[9]
        )
        .to_csv(path_to_xy, index=True, header=True)
    )

    with open(path_to_sobol, "w") as f_out:
        for i, y in enumerate(ys.T):
            indices = sobol.analyze(problem, y)
            print_indices(indices, problem, True, f_out)


def evaluate_model(data, land_use_factors, x):
    cost_factors = CostFactors.from_x(x)
    data = read_data(data, land_use_factors, cost_factors)
    optimal_scenario = data.cost[data.cost == data.cost.min()].isel(scenario=0).scenario.item()
    optimal_data = data.sel(scenario=optimal_scenario)
    optimal_data = {
        key: optimal_data[key].item()
        for key in ("cost", "land_use", "roof", "util", "wind", "offshore")
    }
    return (
        optimal_data["cost"],
        optimal_data["land_use"],
        optimal_data["util"],
        optimal_data["wind"],
        optimal_data["roof"],
        optimal_data["offshore"],
        r50_cost(data, optimal_data, "util"),
        r50_cost(data, optimal_data, "offshore"),
        r50_cost(data, optimal_data, "roof"),
        r50_cost(data, optimal_data, "wind")
    )


def r50_cost(data, optimal_data, tech):
    """Calculates relative cost of reducing land use by 50% using a single technology."""
    conditions = [
        data[other_tech] <= optimal_data[other_tech]
        for other_tech in ALL_TECHS
        if other_tech != tech
    ]
    conditions.append(data.land_use <= 0.5 * optimal_data["land_use"])
    mask = functools.reduce(lambda x, y: x & y, conditions)
    if len(data.sel(scenario=mask).scenario) > 0:
        optimal_data_50 = data.sel(scenario=mask).sortby("cost").isel(scenario=0)
        delta_cost_eur = optimal_data_50.cost.item() - optimal_data["cost"]
        delta_land_m2 = (optimal_data_50.land_use.item() - optimal_data["land_use"]) * 1e6
        return delta_cost_eur / delta_land_m2
    else:
        return np.nan


def read_data(data, land_use_factors, cost_factors):
    cost_data = (
        cost_factors
        .apply_to_cost(data.cost.sum("locs"))
        .sum("techs")
    )
    land_use_data = (
        land_use_factors
        .apply_to_energy_cap(data.energy_cap.sum("locs"))
        .sum("techs")
        .rename("land_use")
    )
    return xr.Dataset({"cost": cost_data, "land_use": land_use_data})


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
        land_use_factors=LandUseFactors(
            wind_onshore=snakemake.params.land_factors["wind_onshore_monopoly"],
            wind_offshore=snakemake.params.land_factors["wind_offshore"],
            roof_mounted_pv=snakemake.params.land_factors["roof_mounted_pv"],
            open_field_pv=snakemake.params.land_factors["open_field_pv"]
        ),
        number_uncertainty_runs=snakemake.params.runs,
        path_to_xy=snakemake.output.xy,
        path_to_sobol=snakemake.output.sobol
    )
