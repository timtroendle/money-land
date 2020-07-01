from dataclasses import dataclass

import numpy as np
import geopandas as gpd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

MAX_IRRADIATION = 1000 # Wp / m2

BLUE = "#4F6DB8"
CMAP = sns.light_palette(BLUE, as_cmap=True)
EDGE_WIDTH = 0.4
EDGE_COLOR = "white"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "
MAP_MIN_X = 2200000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6300000
MAP_MAX_Y = 5500000


@dataclass
class LandUseFactors:
    wind_onshore: float
    open_field_pv: float
    wind_offshore: float = 0
    roof_mounted_pv: float = 0

    @classmethod
    def from_x(cls, x):
        return LandUseFactors(
            wind_onshore=1 / x[4],
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
        factors.loc[{"techs": ["roof_mounted_pv_n",
                               "roof_mounted_pv_e_w",
                               "roof_mounted_pv_s_flat"]}] = self.roof_mounted_pv
        factors.loc[{"techs": "open_field_pv"}] = self.open_field_pv
        return energy_cap * factors


def map_land_req(path_to_xy, path_to_aggregated_results, path_to_shapes, path_to_output):
    xy = xr.open_dataset(path_to_xy)
    expected_cost_minimal_scenario = xy.mean("sample_id").cost.to_series().idxmin()
    expected_land_use_factors = LandUseFactors.from_x([
        0, 0, 0, 0, # cost factors
        xy["capacity-density-wind-onshore"].mean().item(),
        xy["efficiency-open-field-pv"].mean().item(),
        xy["module-share-open-field-pv"].mean().item()]
    )

    det = xr.open_dataset(path_to_aggregated_results)

    land_use = (
        expected_land_use_factors
        .apply_to_energy_cap(det.sel(scenario=expected_cost_minimal_scenario).energy_cap)
        .sum("techs")
    )
    expected_land_use = xy.sel(scenario=expected_cost_minimal_scenario).mean("sample_id").land_use.item()
    assert land_use.sum() > expected_land_use * 0.9
    assert land_use.sum() < expected_land_use * 1.1

    countries = (
        gpd
        .read_file(path_to_shapes)
        .to_crs(EPSG_3035_PROJ4)
        .set_index("id")
    )
    countries["total_area_km2"] = countries.area / 1e6
    countries["land_req_km2"] = land_use.to_series().reindex(countries.index)
    countries["land_req_rel"] = (countries["land_req_km2"] / countries["total_area_km2"]) * 100

    fig = plt.figure(figsize=(8, 7))
    axes = fig.subplots(1, 2, gridspec_kw={'width_ratios': [50, 1]}).flatten()
    _plot_map(countries, axes[0])
    _plot_colorbar(fig, axes[1].inset_axes([0, 0.175, 1, 0.65]), vmin=0, vmax=countries.land_req_rel.max())
    axes[1].axis("off")
    sns.despine(fig, left=True, bottom=True)
    fig.savefig(path_to_output, dpi=600)


def _plot_map(countries, ax):
    ax.set_aspect('equal')
    countries.plot(
        ax=ax,
        column="land_req_rel",
        cmap=CMAP,
        edgecolor=EDGE_COLOR,
        linewidth=EDGE_WIDTH,
        legend=False
    )
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)
    ax.set_xticks([])
    ax.set_yticks([])


def _plot_colorbar(fig, axis, vmin, vmax):
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    s_m = matplotlib.cm.ScalarMappable(cmap=CMAP, norm=norm)
    cmap = s_m.get_cmap()
    colors = cmap(np.linspace(1 / 4, 1, cmap.N))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cut_jet', colors)
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax))
    s_m.set_array([])
    cbar = fig.colorbar(s_m, cax=axis)
    cbar.outline.set_linewidth(0)
    cbar.ax.set_ylabel('Land requirements (% total national)', rotation=90)


if __name__ == "__main__":
    map_land_req(
        path_to_xy=snakemake.input.xy,
        path_to_aggregated_results=snakemake.input.aggregated_results,
        path_to_shapes=snakemake.input.shapes,
        path_to_output=snakemake.output[0]
    )
