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
FACTORS = pd.Series({ # FIXME inject
    "wind_onshore_monopoly": 1 / 8,
    "wind_onshore_competing": 1 / 8,
    "roof_mounted_pv": 0,
    "open_field_pv": 1 / 80
}).to_xarray().rename(index="techs")


def plot_both_ternary(path_to_data, path_to_plot):
    cost_data, land_use_data = read_data(path_to_data)

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 3.5))
    gs = gridspec.GridSpec(1, 5, width_ratios=[5, 0.1, 1, 5, 0.1])
    ax_1 = fig.add_subplot(gs[0])
    ax_2 = fig.add_subplot(gs[3])
    cbar_ax_1 = fig.add_subplot(gs[1])
    cbar_ax_2 = fig.add_subplot(gs[4])

    vmin = cost_data.min()
    vmax = cost_data.max()
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    plot_ternary(cost_data.to_dict(), vmin=vmin, vmax=vmax, ax=ax_1)
    plot_colorbar(fig, cbar_ax_1, norm, "viridis")
    ax_1.annotate('a', xy=[-0.08, 1.05], xycoords='axes fraction',
                  fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    vmin = land_use_data.min()
    vmax = land_use_data.max()
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    plot_ternary(land_use_data.to_dict(), vmin=vmin, vmax=vmax, ax=ax_2)
    plot_colorbar(fig, cbar_ax_2, norm, "viridis")
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


def read_data(path_to_data):
    data = xr.open_dataset(path_to_data)
    cost_data = data.cost.sum(["locs", "techs"]).to_series()
    cost_data.index = cost_data.index.map(xyz).rename(["util", "wind", "roof"])
    cost_data = (cost_data / cost_data.min())
    land_use_data = ((
        data
        .energy_cap
        .sum("locs")
        .sel(techs=["wind_onshore_monopoly", "wind_onshore_competing", "roof_mounted_pv", "open_field_pv"])
    ) * FACTORS).sum("techs").to_series()
    land_use_data.index = land_use_data.index.map(xyz).rename(["util", "wind", "roof"])
    land_use_data = (land_use_data / land_use_data[cost_data[cost_data == cost_data.min()].index].values[0])
    return cost_data, land_use_data


def xyz(scenario_name):
    roof, util, wind = tuple(int(x.split("-")[1]) // 10 for x in scenario_name.split(","))
    return (util, wind, roof)


def plot_ternary(data, vmin, vmax, ax):
    scale = 10
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    figure, tax = ternary.figure(ax=ax, scale=scale)
    tax.boundary(linewidth=1.0)
    tax.heatmap(data, scale=10, style="triangular", colorbar=False, vmin=vmin, vmax=vmax)
    tax.right_corner_label("Utility-scale PV", ha="center", rotation=25)
    tax.top_corner_label("Onshore wind", offset=0.2)
    tax.left_corner_label("Rooftop PV", ha="center", rotation=-25)
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
        path_to_plot=snakemake.output[0]
    )
