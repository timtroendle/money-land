from dataclasses import dataclass

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns

PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"

GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"
ANTHRACITE = "#424242"
PALETTE = sns.light_palette(BLUE, n_colors=11)


@dataclass
class Annotation:
    text: str
    xy: tuple
    xytext: tuple


@dataclass
class PlotData:
    data: pd.DataFrame
    pareto_points: pd.DataFrame
    case: str
    case_title: str
    annotations: tuple


ANNOTATIONS = {
    "roof": {
        "footprint-only": (),
        "land-use": (
            Annotation(text="30% utility-scale PV\n70% onshore wind", xy=(3, 7, 0), xytext=(0.4, 1.4)),
            Annotation(text="100% utility-scale PV", xy=(10, 0, 0), xytext=(0.25, 1.7)),
            Annotation(text="100% rooftop PV", xy=(0, 0, 10), xytext=(0.25, 1.8)),
            Annotation(text="100% onshore wind", xy=(0, 10, 0), xytext=(0.9, 1.3))
        )
    },
    "offshore": {
        "footprint-only": (),
        "land-use": (
            Annotation(text="30% utility-scale PV\n70% onshore wind", xy=(3, 7, 0), xytext=(1.0, 1.2)),
            Annotation(text="30% utility-scale PV\n70% offshore wind", xy=(3, 0, 7), xytext=(0.2, 1.3)),
            Annotation(text="40% utility-scale PV\n30% onshore wind\n30% offshore wind", xy=(4, 3, 3),
                       xytext=(0.5, 1.18))
        )
    }
}
PARETO_POINTS_INDEX = {
    "roof": {
        "footprint-only": (),
        "land-use": [
            (3, 7, 0), (4, 6, 0), (5, 5, 0), (6, 4, 0), (7, 3, 0), (8, 2, 0), (9, 1, 0),
            (10, 0, 0), (9, 0, 1), (8, 0, 2), (7, 0, 3), (6, 0, 4), (5, 0, 5), (4, 0, 6),
            (3, 0, 7), (2, 0, 8), (1, 0, 9), (0, 0, 10)
        ]
    },
    "offshore": {
        "footprint-only": (),
        "land-use": ()
    }
}


def scatter(path_to_results, case, land_definition, land_use_factors, path_to_plot):
    assert case in ["roof", "offshore"]
    assert land_definition in ["footprint-only", "land-use"]
    plot_data = read_data(
        path_to_results,
        case=case,
        land_definition=land_definition,
        land_use_factors=land_use_factors
    )

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 4))
    gs = gridspec.GridSpec(3, 2, width_ratios=[2, 1])
    ax_main = fig.add_subplot(gs[:, 0])
    ax_roof = fig.add_subplot(gs[0, 1], sharex=ax_main)
    ax_util = fig.add_subplot(gs[1, 1], sharex=ax_main)
    ax_wind = fig.add_subplot(gs[2, 1], sharex=ax_main)
    ax_aux = fig.add_subplot(gs[:, 1])
    ax_aux.axis('off')

    ax_main.plot(
        plot_data.pareto_points.land_use,
        plot_data.pareto_points.cost,
        label="Pareto frontier",
        color=RED
    )
    ax_main.plot(
        [1.0],
        [1.0],
        label="Cost minimum",
        marker="o",
        linestyle="",
        color=RED
    )
    sns.scatterplot(
        data=plot_data.data.reset_index(),
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
    ax_main.set_xlabel("Land requirements relative to cost minimal case")
    ax_main.set_ylabel("Cost relative to cost minimal case")
    ax_main.annotate('a', xy=[-0.08, 1.05], xycoords='axes fraction',
                     fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    for annotation in plot_data.annotations:
        ax_main.annotate(
            annotation.text,
            xy=(plot_data.data.loc[annotation.xy].values),
            xytext=annotation.xytext,
            arrowprops={"arrowstyle": "->"}
        )

    sns.scatterplot(
        data=plot_data.data.reset_index(),
        y="cost",
        x="land_use",
        hue=plot_data.case,
        legend=False,
        ax=ax_roof,
        palette=PALETTE
    )
    ax_roof.set_title(plot_data.case_title)
    for tick in ax_roof.xaxis.get_major_ticks():
        tick.set_visible(False)
    for tick in ax_roof.yaxis.get_major_ticks():
        tick.set_visible(False)
    ax_roof.set_xlabel("")
    ax_roof.set_ylabel("")
    ax_aux.annotate('b', xy=[-0.08, 1.05], xycoords='axes fraction',
                    fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    sns.scatterplot(
        data=plot_data.data.reset_index(),
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
        data=plot_data.data.reset_index(),
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
    ax_wind.set_xlabel("Land requirements relative to cost minimal case")

    sns.despine()
    fig.tight_layout()
    fig.savefig(path_to_plot, dpi=600)


def read_data(path_to_data, case, land_definition, land_use_factors):
    data = xr.open_dataset(path_to_data)
    cost_data = data.cost.sum(["locs", "techs"]).to_series()
    cost_data = (cost_data / cost_data.min())
    land_use_data = ((
        data
        .energy_cap
        .sum("locs")
        .sel(techs=["wind_onshore_monopoly", "wind_onshore_competing", "wind_offshore",
                    "roof_mounted_pv_n", "roof_mounted_pv_e_w", "roof_mounted_pv_s_flat", "open_field_pv"])
    ) * land_use_factors).sum("techs").to_series()
    land_use_data = (land_use_data / land_use_data[cost_data[cost_data == cost_data.min()].index].values[0])

    cost_data.index = scenario_name_to_multiindex(cost_data.index)
    cost_data = filter_three_dimensions(cost_data, case)
    land_use_data.index = scenario_name_to_multiindex(land_use_data.index)
    land_use_data = filter_three_dimensions(land_use_data, case)

    both_data = pd.DataFrame({"land_use": land_use_data, "cost": cost_data})

    if case == "roof":
        case_title = "Rooftop PV"
    else:
        case_title = "Offshore wind"

    pareto_index = PARETO_POINTS_INDEX[case][land_definition]
    if len(pareto_index) > 0: # manually defined
        pareto_points = both_data.loc[pareto_index]
    else: # derive automatically
        pareto_points = both_data.loc[is_pareto_efficient(both_data.values), :].sort_values("cost")

    return PlotData(
        data=both_data,
        pareto_points=pareto_points,
        case=case,
        case_title=case_title,
        annotations=ANNOTATIONS[case][land_definition]
    )


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


def is_pareto_efficient(costs):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """
    # taken from https://stackoverflow.com/a/40239615/1856079
    is_efficient = np.ones(costs.shape[0], dtype=bool)
    for i, c in enumerate(costs):
        is_efficient[i] = np.all(np.any(costs[:i] > c, axis=1)) and np.all(np.any(costs[i + 1:] > c, axis=1))
    return is_efficient


if __name__ == "__main__":
    scatter(
        path_to_results=snakemake.input.results,
        case=snakemake.wildcards.case,
        land_definition=snakemake.wildcards.land,
        land_use_factors=pd.Series(snakemake.params.land_factors).to_xarray().rename(index="techs"),
        path_to_plot=snakemake.output[0]
    )
