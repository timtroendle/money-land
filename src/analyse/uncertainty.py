from dataclasses import dataclass
import functools

import numpy as np
import pandas as pd
import xarray as xr
from SALib.sample import saltelli
from SALib.analyze import sobol
import scipy as sp
from scipy import stats
import SALib.util

ALL_TECHS = ["roof", "util", "wind", "offshore"]
MAX_IRRADIATION = 1000 # Wp / m2


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

    @classmethod
    def validate_problem(cls, problem):
        # otherwise the from_x classmethod above will silently deliver wrong results
        assert problem["names"][0] == "cost-wind-onshore"
        assert problem["names"][1] == "cost-wind-offshore"
        assert problem["names"][2] == "cost-roof-mounted-pv"
        assert problem["names"][3] == "cost-open-field-pv"

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
    open_field_pv: float
    wind_offshore: float = 0
    roof_mounted_pv: float = 0

    @classmethod
    def from_x(cls, x):
        return LandUseFactors(
            wind_onshore=x[4],
            open_field_pv=LandUseFactors.pv_land_use_factor(module_efficiency=x[5], module_share=x[6])
        )

    @classmethod
    def validate_problem(cls, problem):
        # otherwise the from_x classmethod above will silently deliver wrong results
        assert problem["names"][4] == "capacity-density-wind-onshore"
        assert problem["names"][5] == "efficiency-open-field-pv"
        assert problem["names"][6] == "module-share-open-field-pv"

    @classmethod
    def pv_land_use_factor(cls, module_efficiency, module_share):
        installable_watt_per_m2 = MAX_IRRADIATION * module_share * module_efficiency
        return 1 / installable_watt_per_m2

    def apply_to_energy_cap(self, energy_cap):
        factors = xr.zeros_like(energy_cap)
        factors.loc[{"techs": ["wind_onshore_competing", "wind_onshore_monopoly"]}] = self.wind_onshore
        factors.loc[{"techs": "wind_offshore"}] = self.wind_offshore
        factors.loc[{"techs": "roof_mounted_pv"}] = self.roof_mounted_pv
        factors.loc[{"techs": "open_field_pv"}] = self.open_field_pv
        return energy_cap * factors


def uncertainty_analysis(path_to_results, uncertain_parameters, number_uncertainty_runs, path_to_xy, path_to_sobol):
    monkey_patch_salib()
    problem = {
        'num_vars': len(uncertain_parameters.keys()),
        'names': [param_name for param_name in uncertain_parameters.keys()],
        'bounds': [(attrs["bound1"], attrs["bound2"]) for attrs in uncertain_parameters.values()],
        'dists': [param_attributes["distribution"] for param_attributes in uncertain_parameters.values()]
    }
    param_values = saltelli.sample(problem, number_uncertainty_runs, calc_second_order=False)

    data = xr.open_dataset(path_to_results)
    ys = np.zeros([param_values.shape[0], 14])
    for i, x in enumerate(param_values):
        ys[i] = evaluate_model(data, x)
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
            r25_cost_util=ys.T[6],
            r25_cost_offshore=ys.T[7],
            r25_cost_roof=ys.T[8],
            r25_cost_wind=ys.T[9],
            r50_cost_util=ys.T[10],
            r50_cost_offshore=ys.T[11],
            r50_cost_roof=ys.T[12],
            r50_cost_wind=ys.T[13]
        )
        .to_csv(path_to_xy, index=True, header=True)
    )

    with open(path_to_sobol, "w") as f_out:
        for i, y in enumerate(ys.T):
            indices = sobol.analyze(problem, y, calc_second_order=False)
            print_indices(indices, problem, False, f_out)


def evaluate_model(data, x):
    cost_factors = CostFactors.from_x(x)
    land_use_factors = LandUseFactors.from_x(x)
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
        cost_of_reducing_land_use(data, optimal_data, "util", reduction_level=0.25),
        cost_of_reducing_land_use(data, optimal_data, "offshore", reduction_level=0.25),
        cost_of_reducing_land_use(data, optimal_data, "roof", reduction_level=0.25),
        cost_of_reducing_land_use(data, optimal_data, "wind", reduction_level=0.25),
        cost_of_reducing_land_use(data, optimal_data, "util", reduction_level=0.5),
        cost_of_reducing_land_use(data, optimal_data, "offshore", reduction_level=0.5),
        cost_of_reducing_land_use(data, optimal_data, "roof", reduction_level=0.5),
        cost_of_reducing_land_use(data, optimal_data, "wind", reduction_level=0.5)
    )


def cost_of_reducing_land_use(data, optimal_data, tech, reduction_level):
    """Calculates relative cost of reducing land use by given level using a single technology."""
    assert 0 < reduction_level <= 1
    conditions = [
        data[other_tech] <= optimal_data[other_tech]
        for other_tech in ALL_TECHS
        if other_tech != tech
    ]
    conditions.append(data.land_use <= (1 - reduction_level) * optimal_data["land_use"])
    mask = functools.reduce(lambda x, y: x & y, conditions)
    if len(data.sel(scenario=mask).scenario) > 0:
        selected_scenarios = data.sel(scenario=mask)
        delta_cost = selected_scenarios.cost - optimal_data["cost"]
        delta_land_m2 = (selected_scenarios.land_use - optimal_data["land_use"]) * 1e6
        return (delta_cost / delta_land_m2).max().item()
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
        _, names = SALib.util.compute_groups_matrix(problem['groups'])
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


def monkey_patch_salib():
    # SALib throws an Error when using normal distribution, because scipy.stats isn't
    # properly imported. I am monkey patching this here.
    sp.stats = stats
    SALib.util.sp = sp


if __name__ == "__main__":
    uncertainty_analysis(
        path_to_results=snakemake.input.results,
        uncertain_parameters=snakemake.params.uncertain_parameters,
        number_uncertainty_runs=snakemake.params.runs,
        path_to_xy=snakemake.output.xy,
        path_to_sobol=snakemake.output.sobol
    )
