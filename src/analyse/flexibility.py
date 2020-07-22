from dataclasses import dataclass

import ternary
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns


RED = "#A01914"
BLUE = "#4F6DB8"
PALETTE = sns.light_palette(RED, as_cmap=True)


@dataclass
class PlotData:
    data: pd.Series
    name: str
    left_axis_label: str
    bottom_axis_label: str = "Utility-scale PV →"
    right_axis_label: str = "← Onshore wind"


def flexibility(path_to_results, path_to_plot):
    roof_plot_datas = read_data(path_to_results, "roof")
    offshore_plot_datas = read_data(path_to_results, "offshore")

    fig = plt.figure(figsize=(8, 4))
    gs = gridspec.GridSpec(3, 4, width_ratios=[5, 5, 5, 5], height_ratios=[4, 4, 1])
    roof_axes = [fig.add_subplot(gs[0, x]) for x in [0, 1, 2, 3]]
    offshore_axes = [fig.add_subplot(gs[1, x]) for x in [0, 1, 2, 3]]
    cbar_axes = [fig.add_subplot(gs[2, x]) for x in [0, 1, 2, 3]]

    for plot_datas, axes, cbar_ax, panel_id in zip(zip(roof_plot_datas, offshore_plot_datas),
                                                   zip(roof_axes, offshore_axes),
                                                   cbar_axes,
                                                   list("abcd")):
        norm = matplotlib.colors.Normalize(
            vmin=min(plot_data.data.min() for plot_data in plot_datas),
            vmax=max(plot_data.data.max() for plot_data in plot_datas)
        )
        plot_ternary(plot_datas[0], ax=axes[0], norm=norm)
        plot_ternary(plot_datas[1], ax=axes[1], norm=norm)
        plot_colorbar(fig, cbar_ax, norm, PALETTE)
        axes[0].set_title(panel_id + " - " + plot_datas[0].name, loc="left")

    plt.subplots_adjust(
        left=0.02,
        bottom=0.03,
        right=0.98,
        top=0.90,
        wspace=0.0,
        hspace=0.2
    )

    fig.savefig(path_to_plot, dpi=600)


def plot_ternary(plot_data, ax, norm):
    scale = 10
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    figure, tax = ternary.figure(ax=ax, scale=scale)
    tax.boundary(linewidth=.5)
    tax.heatmap(
        plot_data.data.to_dict(),
        scale=10,
        style="triangular",
        colorbar=False,
        cmap=PALETTE,
        vmin=norm.vmin,
        vmax=norm.vmax
    )
    tax.bottom_axis_label(plot_data.bottom_axis_label, ha="center", offset=0.05)
    tax.right_axis_label(plot_data.right_axis_label, offset=0.10)
    tax.left_axis_label(plot_data.left_axis_label, ha="center", offset=0.10)
    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    ax.set_aspect(1)


def plot_colorbar(fig, ax, norm, cmap):
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    cbar = fig.colorbar(s_m, ax=ax, fraction=1, aspect=35, shrink=0.80, orientation="horizontal")
    cbar_ticks = np.linspace(
        start=norm.vmin,
        stop=norm.vmax,
        num=4
    )
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(["{:.0f}".format(tick)
                         for tick in cbar.get_ticks()])
    cbar.outline.set_linewidth(0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.axis('off')


def read_data(path_to_data, case):
    data = xr.open_dataset(path_to_data)
    energy_cap = data.energy_cap.sum("locs") / 1_000
    storage_cap = data.storage_cap.sum("locs") / 1_000_000

    bio = energy_cap.sel(techs="biofuel").to_series()
    stor = energy_cap.sel(techs=["battery", "hydrogen"]).sum("techs").to_series()
    e_stor = storage_cap.sel(techs=["battery", "hydrogen"]).sum("techs").to_series()
    if "ac_transmission" in energy_cap.techs:
        trans = energy_cap.sel(techs="ac_transmission").to_series()
    else:
        trans = xr.zeros_like(energy_cap.sel(techs="biofuel")).to_series()

    for da in [bio, stor, e_stor, trans]:
        da.index = scenario_name_to_multiindex(da.index)
    bio = filter_three_dimensions(bio, case)
    stor = filter_three_dimensions(stor, case)
    e_stor = filter_three_dimensions(e_stor, case)
    trans = filter_three_dimensions(trans, case)
    if case == "roof":
        left_axis_label = "← Rooftop PV"
    else:
        left_axis_label = "← Offshore wind"

    return [
        PlotData(
            data=da,
            name=name,
            left_axis_label=left_axis_label,
        )
        for da, name in ((stor, "Storage (GW)"),
                         (e_stor, "Storage (TWh)"),
                         (trans, "Transmission (GW)"),
                         (bio, "Bioenergy (GW)"))
    ]


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


if __name__ == "__main__":
    flexibility(
        path_to_results=snakemake.input.results,
        path_to_plot=snakemake.output[0]
    )
