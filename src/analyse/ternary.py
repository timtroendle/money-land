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

idx = pd.IndexSlice


@dataclass
class PlotData:
    data: pd.Series
    norm: matplotlib.colors.Normalize
    left_axis_label: str
    panel_name: str
    bottom_axis_label: str = "Utility-scale PV (%) →"
    right_axis_label: str = "← Onshore wind (%)"


def plot_both_ternary(path_to_data, land_use_factors, path_to_plot):
    plot_datas = read_data(path_to_data, land_use_factors)
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(3, 2, width_ratios=[5, 5], height_ratios=[25, 25, 1])
    ax_1 = fig.add_subplot(gs[0, 0])
    ax_2 = fig.add_subplot(gs[0, 1])
    ax_3 = fig.add_subplot(gs[1, 0])
    ax_4 = fig.add_subplot(gs[1, 1])
    cbar_ax_1 = fig.add_subplot(gs[2, 0])
    cbar_ax_2 = fig.add_subplot(gs[2, 1])

    plot_ternary(plot_datas[0], ax=ax_1, cmap=SEQUENTIAL_PALETTE)
    plot_ternary(plot_datas[1], ax=ax_2, cmap=DIVERGING_PALETTE)
    plot_ternary(plot_datas[2], ax=ax_3, cmap=SEQUENTIAL_PALETTE)
    plot_ternary(plot_datas[3], ax=ax_4, cmap=DIVERGING_PALETTE)

    plot_sequential_colorbar(fig, cbar_ax_1, plot_datas[0].norm, cmap=SEQUENTIAL_PALETTE,
                             label="Cost relative to cost minimal case")
    plot_diverging_colorbar(fig, cbar_ax_2, plot_datas[1].norm, cmap=DIVERGING_PALETTE,
                            label="Land use relative to cost minimal case",
                            land_use_data=plot_datas[1].data)

    plt.subplots_adjust(
        left=0.05,
        bottom=0.07,
        right=0.95,
        top=0.98,
        wspace=0.2,
        hspace=0.05
    )
    fig.savefig(path_to_plot, dpi=600)


def read_data(path_to_data, land_use_factors):
    data = xr.open_dataset(path_to_data).sum("locs")
    data["roof"] = data["roof"] // 10
    data["util"] = data["util"] // 10
    data["wind"] = data["wind"] // 10
    data["offshore"] = data["offshore"] // 10
    data["land_use"] = data["energy_cap"] * land_use_factors
    data = (
        data[["cost", "land_use"]]
        .sum("techs")
        .sel(scenario=(data.roof == 0) | (data.offshore == 0))
        .to_dataframe()
        .set_index(["util", "wind", "roof", "offshore"])
    )
    data = data / data.loc[data.cost.idxmin()]
    return [
        PlotData(
            data=filter_three_dimensions(data.cost, "roof"),
            left_axis_label="← Rooftop PV (%)",
            norm=matplotlib.colors.Normalize(vmin=data.cost.min(), vmax=data.cost.max()),
            panel_name="a - Cost without\n     offshore wind"
        ),
        PlotData(
            data=filter_three_dimensions(data.land_use, "roof"),
            left_axis_label="← Rooftop PV (%)",
            norm=matplotlib.colors.Normalize(vmin=data.land_use.min(), vmax=1 + (1 - data.land_use.min())),
            panel_name="b - Land use without\n     offshore wind"
        ),
        PlotData(
            data=filter_three_dimensions(data.cost, "offshore"),
            left_axis_label="← Offshore wind (%)",
            norm=matplotlib.colors.Normalize(vmin=data.cost.min(), vmax=data.cost.max()),
            panel_name="c - Cost without\n     rooftop PV"
        ),
        PlotData(
            data=filter_three_dimensions(data.land_use, "offshore"),
            left_axis_label="← Offshore wind (%)",
            norm=matplotlib.colors.Normalize(vmin=data.land_use.min(), vmax=1 + (1 - data.land_use.min())),
            panel_name="d - Land use without\n     rooftop PV"
        )
    ]


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
    ax.annotate(plot_data.panel_name, xy=[-0.08, 0.95], xycoords='axes fraction',
                fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    ax.set_aspect(1)


def plot_sequential_colorbar(fig, ax, norm, cmap, label):
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    cbar = fig.colorbar(s_m, ax=ax, fraction=1, aspect=35, shrink=1.0, orientation="horizontal")
    cbar_ticks = np.linspace(
        start=norm.vmin,
        stop=norm.vmax,
        num=4
    )
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(["{:.1f}".format(tick)
                         for tick in cbar.get_ticks()])
    cbar.outline.set_linewidth(0)
    cbar.set_label(label)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.axis('off')


def plot_diverging_colorbar(fig, ax, norm, cmap, label, land_use_data):
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    cmap = s_m.get_cmap()
    rel_max = (land_use_data.max() - land_use_data.min()) / (norm.vmax - norm.vmin)
    colors = cmap(np.linspace(0, rel_max, cmap.N))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cut_jet', colors)
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(vmin=0, vmax=land_use_data.max()))
    s_m.set_array([])
    cbar = fig.colorbar(s_m, ax=ax, fraction=1, aspect=35, shrink=1.0, orientation="horizontal")
    cbar_ticks = np.linspace(
        start=land_use_data.min(),
        stop=land_use_data.max(),
        num=6
    )
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(["{:.1f}".format(tick)
                         for tick in cbar.get_ticks()])
    cbar.outline.set_linewidth(0)
    cbar.set_label(label)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.axis('off')


if __name__ == "__main__":
    plot_both_ternary(
        path_to_data=snakemake.input.results,
        land_use_factors=pd.Series(snakemake.params.land_factors).to_xarray().rename(index="techs"),
        path_to_plot=snakemake.output[0]
    )
