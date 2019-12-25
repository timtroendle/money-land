from dataclasses import dataclass

import ternary
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns


PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"


@dataclass
class PlotData:
    data: pd.Series
    right_corner_label: str
    top_corner_label: str
    left_corner_label: str
    norm: matplotlib.colors.Normalize


def plot_both_ternary(path_to_data, case, land_use_factors, path_to_plot):
    assert case in ["roof", "offshore"]
    plot_data_cost, plot_data_land_use = read_data(path_to_data, case, land_use_factors)

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 3.5))
    gs = gridspec.GridSpec(1, 5, width_ratios=[5, 0.1, 1, 5, 0.1])
    ax_1 = fig.add_subplot(gs[0])
    ax_2 = fig.add_subplot(gs[3])
    cbar_ax_1 = fig.add_subplot(gs[1])
    cbar_ax_2 = fig.add_subplot(gs[4])

    plot_ternary(plot_data_cost, ax=ax_1)
    plot_colorbar(fig, cbar_ax_1, plot_data_cost.norm, "viridis")
    ax_1.annotate('a', xy=[-0.08, 1.05], xycoords='axes fraction',
                  fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    plot_ternary(plot_data_land_use, ax=ax_2)
    plot_colorbar(fig, cbar_ax_2, plot_data_land_use.norm, "viridis")
    ax_2.annotate('b', xy=[-0.08, 1.05], xycoords='axes fraction',
                  fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    plt.subplots_adjust(
        left=0.05,
        bottom=0.05,
        right=0.95,
        top=0.90,
        wspace=0.2,
        hspace=0.2
    )
    fig.savefig(path_to_plot, dpi=600)


def read_data(path_to_data, case, land_use_factors):
    data = xr.open_dataset(path_to_data)
    cost_data = data.cost.sum(["locs", "techs"]).to_series()
    cost_data = (cost_data / cost_data.min())
    land_use_data = ((
        data
        .energy_cap
        .sum("locs")
        .sel(techs=["wind_onshore_monopoly", "wind_onshore_competing", "wind_offshore",
                    "roof_mounted_pv", "open_field_pv"])
    ) * land_use_factors).sum("techs").to_series()
    land_use_data = (land_use_data / land_use_data[cost_data[cost_data == cost_data.min()].index].values[0])

    cost_data.index = scenario_name_to_multiindex(cost_data.index)
    cost_data = filter_three_dimensions(cost_data, case)
    land_use_data.index = scenario_name_to_multiindex(land_use_data.index)
    land_use_data = filter_three_dimensions(land_use_data, case)
    if case == "roof":
        left_corner_label = "Rooftop PV"
    else:
        left_corner_label = "Offshore wind"

    plot_data_cost = PlotData(
        data=cost_data,
        right_corner_label="Utility-scale PV",
        top_corner_label="Onshore wind",
        left_corner_label=left_corner_label,
        norm=matplotlib.colors.Normalize(vmin=cost_data.min(), vmax=cost_data.max())
    )
    plot_data_land_use = PlotData(
        data=land_use_data,
        right_corner_label="Utility-scale PV",
        top_corner_label="Onshore wind",
        left_corner_label=left_corner_label,
        norm=matplotlib.colors.Normalize(vmin=land_use_data.min(), vmax=land_use_data.max())
    )
    return plot_data_cost, plot_data_land_use


def scenario_name_to_multiindex(index):
    return index.map(wxyz).rename(["util", "wind", "roof", "offshore"])


def wxyz(scenario_name):
    roof, util, wind, offshore = tuple(int(x.split("-")[1]) // 10 for x in scenario_name.split(","))
    return (util, wind, roof, offshore)


def filter_three_dimensions(data, case):
    if case == "roof":
        column = "offshore"
    else:
        column = "roof"
    return (
        data
        .reset_index()[data.reset_index()[column] == 0]
        .drop(columns=[column])
        .set_index(["util", "wind", case])
        .iloc[:, 0]
    )


def plot_ternary(plot_data, ax):
    scale = 10
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    figure, tax = ternary.figure(ax=ax, scale=scale)
    tax.boundary(linewidth=1.0)
    tax.heatmap(
        plot_data.data.to_dict(),
        scale=10,
        style="triangular",
        colorbar=False,
        vmin=plot_data.norm.vmin,
        vmax=plot_data.norm.vmax
    )
    tax.right_corner_label(plot_data.right_corner_label, ha="center", rotation=25)
    tax.top_corner_label(plot_data.top_corner_label, offset=0.2)
    tax.left_corner_label(plot_data.left_corner_label, ha="center", rotation=-25)
    tax.ticks(ticks=range(0, 110, 20), axis='lbr', linewidth=1, multiple=1, offset=0.02)
    tax.clear_matplotlib_ticks()
    tax._redraw_labels()


def plot_colorbar(fig, ax, norm, cmap):
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    cbar = fig.colorbar(s_m, ax=ax, fraction=1, aspect=35, shrink=0.65)
    cbar_ticks = np.linspace(
        start=norm.vmin,
        stop=norm.vmax,
        num=4
    )
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(["{:.1f}".format(tick)
                         for tick in cbar.get_ticks()])
    cbar.outline.set_linewidth(0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.axis('off')


if __name__ == "__main__":
    plot_both_ternary(
        path_to_data=snakemake.input.results,
        case=snakemake.wildcards.case,
        land_use_factors=pd.Series(snakemake.params.land_factors).to_xarray().rename(index="techs"),
        path_to_plot=snakemake.output[0]
    )
