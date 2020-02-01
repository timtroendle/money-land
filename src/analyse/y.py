from dataclasses import dataclass

import numpy as np
import pandas as pd
import xarray as xr

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


def sample_y(path_to_x, path_to_model_output, path_to_xy):
    model_output = xr.open_dataset(path_to_model_output)
    model_output = model_output[["energy_cap", "cost"]].sum("locs")

    x = pd.read_csv(path_to_x, index_col=0)

    xy = (
        xr
        .ones_like(model_output.sum("techs").expand_dims(sample_id=x.index))
        .rename({"energy_cap": "land_use"})
    ) * np.nan

    for sample_id in x.index:
        add_y(x.loc[sample_id].values, sample_id, model_output, xy)

    for col in x.columns:
        xy.coords[col] = ("sample_id", x[col].values)

    xy.to_netcdf(
        path_to_xy,
        encoding={var: dict(zlib=True, complevel=5) for var in xy.data_vars}
    )


def add_y(x, sample_id, model_output, xy):
    land_use = (
        LandUseFactors
        .from_x(x)
        .apply_to_energy_cap(model_output["energy_cap"])
        .sum("techs")
    )
    cost = (
        CostFactors
        .from_x(x)
        .apply_to_cost(model_output["cost"])
        .sum("techs")
    )
    xy["land_use"].loc[{"sample_id": sample_id}] = land_use
    xy["cost"].loc[{"sample_id": sample_id}] = cost


if __name__ == "__main__":
    sample_y(
        path_to_x=snakemake.input.x,
        path_to_model_output=snakemake.input.model_output,
        path_to_xy=snakemake.output[0]
    )
