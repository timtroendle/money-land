from dataclasses import dataclass
import pandas as pd
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt

EPSILON = 0.00001
GREY = "#7F7F7F"
BLUE = "#4F6DB8"
BRIGHT_COLOR = sns.light_palette(BLUE, 3)[0]
DARK_COLOR = BLUE
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"

idx = pd.IndexSlice


@dataclass
class PlotData:
    name: str
    data: pd.DataFrame
    optimal_path: pd.Series
    average_slope: float
    ylabel: str = "Cost relative to\n cost minimal case"
    xlabel: str = "Land use relative to\n cost minimal case"


def technology_plot(path_to_results, path_to_plot):
    sns.set_context("paper")
    plot_datas = read_data(path_to_results)
    fig = plt.figure(figsize=(8, 3))
    axes = fig.subplots(1, 3, sharey=True, sharex=True)

    for plot_data, ax in zip(plot_datas, axes):
        sns.scatterplot(
            data=plot_data.data,
            x="land_use",
            y="cost",
            color=BRIGHT_COLOR,
            ax=ax
        )
        ax.plot(
            plot_data.data.loc[plot_data.optimal_path.index, "land_use"],
            plot_data.data.loc[plot_data.optimal_path.index, "cost"],
            color=DARK_COLOR,
            linestyle="--",
            marker="o",
            markersize=5
        )
        first_point = plot_data.data.loc[plot_data.optimal_path.index[0]]
        last_point = plot_data.data.loc[plot_data.optimal_path.index[-1]]
        ax.vlines(
            x=first_point["land_use"],
            ymin=first_point["cost"],
            ymax=last_point["cost"],
            color=GREY,
            linestyle=":",
        )
        ax.hlines(
            y=last_point["cost"],
            xmin=last_point["land_use"],
            xmax=first_point["land_use"],
            color=GREY,
            linestyle=":",
        )
        ax.plot(
            (first_point["land_use"], last_point["land_use"]),
            (first_point["cost"], last_point["cost"]),
            color=GREY,
            linestyle=":"
        )
        ax.annotate(
            s=f"{plot_data.average_slope:.2f}" + r" $ \frac{EUR}{m^2 \cdot yr} $",
            xy=(first_point["land_use"], last_point["cost"]),
            xytext=(0.35, last_point["cost"]),
            verticalalignment="bottom" if last_point["cost"] < plot_data.data.cost.max() else "top",
            color=GREY
        )
        ax.annotate(
            plot_data.name,
            xy=[-0.08, 1.05],
            xycoords='axes fraction',
            fontsize=PANEL_FONT_SIZE,
            weight=PANEL_FONT_WEIGHT
        )
        ax.set_ylabel(plot_data.ylabel)
        ax.set_xlabel(plot_data.xlabel)
    sns.despine()
    fig.tight_layout()
    fig.savefig(path_to_plot, dpi=300)


def read_data(path_to_data):
    data = (
        xr
        .open_dataset(path_to_data)
        .mean("sample_id")
        .to_dataframe()
        .set_index(["util", "wind", "roof", "offshore"])
    )
    rel_data = data / data.loc[data.cost.idxmin()]

    return [
        PlotData(
            name="a - Offshore wind",
            data=rel_data,
            optimal_path=optimal_path_series(data, "offshore"),
            average_slope=average_slope(data, optimal_path_series(data, "offshore"))
        ),
        PlotData(
            name="b - Utility-scale PV",
            data=rel_data,
            optimal_path=optimal_path_series(data, "util"),
            average_slope=average_slope(data, optimal_path_series(data, "util"))
        ),
        PlotData(
            name="c - Rooftop PV",
            data=rel_data,
            optimal_path=optimal_path_series(data, "roof"),
            average_slope=average_slope(data, optimal_path_series(data, "roof"))
        )
    ]


def slope(index1, index2, data):
    cost_delta = data.loc[index1, "cost"] - data.loc[index2, "cost"]
    land_use_delta = (data.loc[index1, "land_use"] - data.loc[index2, "land_use"]) * 1e6 # to m2
    if abs(land_use_delta) > EPSILON:
        return cost_delta / land_use_delta
    else:
        return None


def optimal_path(data, tech):
    assert tech in ["util", "wind", "roof", "offshore"]
    current_index = data["cost"].idxmin() # start at cost minimum
    yield current_index, 0
    while data.loc[[current_index]].index.get_level_values(tech)[0] < 100:
        if tech == "util":
            index = idx[current_index[0] + 10, :current_index[1], :current_index[2], :current_index[3]]
        elif tech == "wind":
            index = idx[:current_index[0], current_index[1] + 10, :current_index[2], :current_index[3]]
        elif tech == "roof":
            index = idx[:current_index[0], :current_index[1], :current_index[2] + 10, :current_index[3]]
        elif tech == "offshore":
            index = idx[:current_index[0], :current_index[1], :current_index[2], current_index[3] + 10]
        slopes = (
            data
            .sort_index(level=['util', 'wind', 'roof', 'offshore'])
            .loc[index, ]
            .apply(lambda row: slope(current_index, row.name, data), axis=1)
        )
        slopes = slopes[slopes < 0] # don't follow slopes that lead to worse points
        if len(slopes.index) == 0:
            return
        else:
            optimal_index = slopes.abs().idxmin()
            yield optimal_index, slopes.loc[optimal_index]
            current_index = optimal_index


def optimal_path_series(data, tech):
    slopes = pd.Series({index: slope for index, slope in optimal_path(data, tech)})
    slopes.index = slopes.index.set_names(["util", "wind", "roof", "offshore"])
    return slopes


def average_slope(data, optimal_path):
    cost_optimal_index = data["cost"].idxmin()
    if len(optimal_path.index) > 0:
        return slope(cost_optimal_index, optimal_path.index[-1], data)
    else:
        return None


def slope_stats_eur_per_m2(data):
    util_slopes = optimal_path_series(data, "util")
    roof_slopes = optimal_path_series(data, "roof")
    offshore_slopes = optimal_path_series(data, "offshore")
    wind_slopes = optimal_path_series(data, "wind")
    return pd.DataFrame(
        data={
            name: [slopes[slopes < 0].max(),
                   slopes[slopes < 0].min(),
                   average_slope(data, slopes[slopes < 0]),
                   slopes.iloc[0]]
            for name, slopes in [("Utility-scale PV", util_slopes),
                                 ("Rooftop PV", roof_slopes),
                                 ("Offshore Wind", offshore_slopes),
                                 ("Onshore Wind", wind_slopes)]
        },
        index=["min", "max", "average", "first"]
    )


if __name__ == "__main__":
    technology_plot(
        path_to_results=snakemake.input.results,
        path_to_plot=snakemake.output[0]
    )
