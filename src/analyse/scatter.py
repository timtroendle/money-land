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

GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"
ANTHRACITE = "#424242"
PALETTE = sns.light_palette(BLUE, n_colors=11)


def scatter(path_to_results, path_to_plot):
    cost_data, land_use_data = read_data(path_to_results)
    both_data = pd.DataFrame({"land_use": land_use_data, "cost": cost_data}).reset_index()

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 4))
    gs = gridspec.GridSpec(3, 2, width_ratios=[2, 1])
    ax_main = fig.add_subplot(gs[:, 0])
    ax_roof = fig.add_subplot(gs[0, 1], sharex=ax_main)
    ax_util = fig.add_subplot(gs[1, 1], sharex=ax_main)
    ax_wind = fig.add_subplot(gs[2, 1], sharex=ax_main)
    ax_aux = fig.add_subplot(gs[:, 1])
    ax_aux.axis('off')

    pareto_points = both_data.set_index(["util", "wind", "roof"]).loc[
        [(0, 100, 0), (10, 90, 0), (20, 80, 0), (30, 70, 0), (40, 60, 0), (50, 50, 0), (60, 40, 0),
         (70, 30, 0), (80, 20, 0), (90, 10, 0), (100, 0, 0), (0, 0, 100)], :]
    cost_optimal_points = both_data.set_index(["util", "wind", "roof"]).loc[(30, 70, 0), :]
    ax_main.plot(
        pareto_points.land_use,
        pareto_points.cost,
        label="Pareto frontier",
        color=RED
    )
    ax_main.plot(
        cost_optimal_points.land_use,
        cost_optimal_points.cost,
        label="Cost minimum",
        marker="o",
        linestyle="",
        color=RED
    )
    sns.scatterplot(
        data=both_data,
        y="cost",
        x="land_use",
        legend=False,
        ax=ax_main,
        color=BLUE
    )
    scatters = [child for child in ax_main.get_children()
                if isinstance(child, matplotlib.collections.PathCollection)][0]
    scatters.set_label("Observations")
    ax_main.legend()
    ax_main.set_xlabel("Land use relative to cost minimal case")
    ax_main.set_ylabel("Cost relative to cost minimal case")
    ax_main.annotate('a', xy=[-0.08, 1.05], xycoords='axes fraction',
                     fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    sns.scatterplot(
        data=both_data,
        y="cost",
        x="land_use",
        hue="roof",
        legend=False,
        ax=ax_roof,
        palette=PALETTE
    )
    ax_roof.set_title("Rooftop PV share")
    for tick in ax_roof.xaxis.get_major_ticks():
        tick.set_visible(False)
    for tick in ax_roof.yaxis.get_major_ticks():
        tick.set_visible(False)
    ax_roof.set_xlabel("")
    ax_roof.set_ylabel("")
    ax_aux.annotate('b', xy=[-0.08, 1.05], xycoords='axes fraction',
                    fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    sns.scatterplot(
        data=both_data,
        y="cost",
        x="land_use",
        hue="util",
        legend=False,
        ax=ax_util,
        palette=PALETTE
    )
    ax_util.set_title("Utility PV share")
    for tick in ax_util.xaxis.get_major_ticks():
        tick.set_visible(False)
    for tick in ax_util.yaxis.get_major_ticks():
        tick.set_visible(False)
    ax_util.set_xlabel("")
    ax_util.set_ylabel("")

    sns.scatterplot(
        data=both_data,
        y="cost",
        x="land_use",
        hue="wind",
        legend=False,
        ax=ax_wind,
        palette=PALETTE
    )
    ax_wind.set_title("Onshore wind share")
    ax_wind.set_ylabel("")
    for tick in ax_wind.yaxis.get_major_ticks():
        tick.set_visible(False)
    ax_wind.set_xlabel("Land use relative to cost minimal case")

    sns.despine()
    fig.tight_layout()
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
    roof, util, wind = tuple(int(x.split("-")[1]) for x in scenario_name.split(","))
    return (util, wind, roof)


if __name__ == "__main__":
    scatter(
        path_to_results=snakemake.input.results,
        path_to_plot=snakemake.output[0]
    )
