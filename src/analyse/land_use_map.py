import xarray as xr
import pandas as pd
import geopandas as gpd
import seaborn as sns
import matplotlib.pyplot as plt

BLUE = "#4F6DB8"
GREY = "#EBEBEB"
EDGE_WIDTH = 0.5
EDGE_COLOR = "white"
ANNOTATION_FONT_SIZE = 8
ANNOTATION_FONT_WEIGHT = "normal"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000

COMPARISON_COUNTRIES = ["SVN", "CHE", "AUT"]


def land_use_map(path_to_results, path_to_shapes, land_factors, path_to_plot):
    countries = (
        gpd
        .read_file(path_to_shapes)
        .to_crs(EPSG_3035_PROJ4)
        .set_index("id")
    )
    cost_optimal_land_use = read_cost_optimal_land_use(path_to_results, land_factors)
    cost_optimal_land_use_relative = cost_optimal_land_use / (countries.area.sum() / 1e6)
    comparison_countries = choose_countries(countries, cost_optimal_land_use)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    countries.plot(color=GREY, ax=ax, edgecolor=EDGE_COLOR, linewidth=EDGE_WIDTH)
    countries.loc[comparison_countries].plot(color=BLUE, ax=ax, edgecolor=EDGE_COLOR, linewidth=EDGE_WIDTH)
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)
    ax.set_xticks([])
    ax.set_yticks([])

    ax.annotate(
        f'             Total: {round_area(countries.area.sum() / 1000 / 1000)} km²             ',
        xy=[-0.20, 0.80],
        xycoords='axes fraction',
        color="black",
        backgroundcolor=GREY,
        fontsize=ANNOTATION_FONT_SIZE,
        fontweight=ANNOTATION_FONT_WEIGHT
    ).set_horizontalalignment("left")
    ax.annotate(
        f'Wind + Solar: ~{round_area(cost_optimal_land_use)} km² ({cost_optimal_land_use_relative * 100:.1f} %)',
        xy=[-0.2, 0.75],
        xycoords='axes fraction',
        color="white",
        fontsize=ANNOTATION_FONT_SIZE,
        fontweight=ANNOTATION_FONT_WEIGHT,
        backgroundcolor=BLUE
    ).set_horizontalalignment("left")

    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)
    fig.savefig(path_to_plot, dpi=600)


def read_cost_optimal_land_use(path_to_results, factors):
    data = xr.open_dataset(path_to_results)
    cost_data = data.cost.sum(["locs", "techs"]).to_series()
    land_use_data = ((
        data
        .energy_cap
        .sum("locs")
        .sel(techs=["wind_onshore_monopoly", "wind_onshore_competing", "wind_offshore",
                    "roof_mounted_pv", "open_field_pv"])
    ) * factors).sum("techs").to_series()
    return land_use_data[cost_data.sort_values(ascending=True).head(1).index][0]


def choose_countries(countries, land_area):
    comparison_countries = [COMPARISON_COUNTRIES.pop()]
    comparison_area = countries.loc[comparison_countries].area.sum() / 1e6
    while not land_area * 0.8 < comparison_area < land_area * 1.2:
        comparison_countries = comparison_countries + [COMPARISON_COUNTRIES.pop()]
        comparison_area = countries.loc[comparison_countries].area.sum() / 1e6
    return comparison_countries


def round_area(area):
    return int(round(area / 10000)) * 10000


if __name__ == "__main__":
    land_use_map(
        path_to_results=snakemake.input.results,
        path_to_shapes=snakemake.input.shapes,
        land_factors=pd.Series(snakemake.params.land_factors).to_xarray().rename(index="techs"),
        path_to_plot=snakemake.output[0]
    )
