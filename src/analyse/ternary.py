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
TICK_FONT_SIZE = 9
RED = "#A01914"
BLUE = "#4F6DB8"
SEQUENTIAL_PALETTE = sns.light_palette(RED, as_cmap=True)
RED_TO_BLUE = [ # from https://gka.github.io using lightness correction
    '#002d6e', '#375aa2', '#6f8ad1', '#a7bffa',
    '#f5f5f5', '#fdad97', '#e36b55', '#b23125', '#720000'
]
DIVERGING_PALETTE = matplotlib.colors.LinearSegmentedColormap.from_list("signature-BlRd", RED_TO_BLUE)


@dataclass
class PlotData:
    data: pd.Series
    norm: matplotlib.colors.Normalize
    left_axis_label: str
    bottom_axis_label: str = "Utility-scale PV (%) →"
    right_axis_label: str = "← Onshore wind (%)"


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

    plot_ternary(plot_data_cost, ax=ax_1, cmap=SEQUENTIAL_PALETTE)
    plot_colorbar(fig, cbar_ax_1, plot_data_cost.norm, cmap=SEQUENTIAL_PALETTE)
    ax_1.annotate('a - Cost', xy=[-0.08, 1.05], xycoords='axes fraction',
                  fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    plot_ternary(plot_data_land_use, ax=ax_2, cmap=DIVERGING_PALETTE)
    plot_colorbar(fig, cbar_ax_2, plot_data_land_use.norm, cmap=DIVERGING_PALETTE)
    ax_2.annotate('b - Land use', xy=[-0.08, 1.05], xycoords='axes fraction',
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
        left_axis_label = "← Rooftop PV (%)"
    else:
        left_axis_label = "← Offshore wind (%)"

    plot_data_cost = PlotData(
        data=cost_data,
        left_axis_label=left_axis_label,
        norm=matplotlib.colors.Normalize(vmin=cost_data.min(), vmax=cost_data.max())
    )
    plot_data_land_use = PlotData(
        data=land_use_data,
        left_axis_label=left_axis_label,
        norm=matplotlib.colors.Normalize(vmin=land_use_data.min(), vmax=1 + (1 - land_use_data.min()))
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


def plot_ternary(plot_data, ax, cmap):
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
        cmap=cmap,
        vmin=plot_data.norm.vmin,
        vmax=plot_data.norm.vmax
    )
    tax.bottom_axis_label(plot_data.bottom_axis_label, ha="center")
    tax.right_axis_label(plot_data.right_axis_label, offset=0.16)
    tax.left_axis_label(plot_data.left_axis_label, ha="center", offset=0.14)
    tax.ticks(ticks=range(0, 110, 20), axis='b', linewidth=1, multiple=1, offset=0.02, fontsize=TICK_FONT_SIZE)
    tax.ticks(ticks=range(0, 110, 20), axis='l', linewidth=1, multiple=1, offset=0.03, fontsize=TICK_FONT_SIZE)
    tax.ticks(ticks=range(0, 110, 20), axis='r', linewidth=1, multiple=1, offset=0.04, fontsize=TICK_FONT_SIZE)
    tax.clear_matplotlib_ticks()
    tax._redraw_labels()


def plot_colorbar(fig, ax, norm, cmap):
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    cbar = fig.colorbar(s_m, ax=ax, fraction=1, aspect=35, shrink=1.0)
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


def plot_other_colorbar(fig, ax, norm, s_m):
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
