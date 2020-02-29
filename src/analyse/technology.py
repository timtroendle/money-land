import pandas as pd
import xarray as xr

EPSILON = 0.00001

idx = pd.IndexSlice


def technology(path_to_data, land_use_factors, path_to_output):
    data = read_data(path_to_data, land_use_factors)
    data = slope_stats_eur_per_m2(data)
    data.to_csv(path_to_output, index=True, header=True)


def read_data(path_to_data, land_use_factors):
    data = xr.open_dataset(path_to_data)
    cost_data = data.cost.sum(["locs", "techs"]).to_series()
    land_use_data = ((
        data
        .energy_cap
        .sum("locs")
        .sel(techs=["wind_onshore_monopoly", "wind_onshore_competing", "wind_offshore",
                    "roof_mounted_pv_n", "roof_mounted_pv_e_w", "roof_mounted_pv_s_flat", "open_field_pv"])
    ) * land_use_factors).sum("techs").to_series()

    cost_data.index = scenario_name_to_multiindex(cost_data.index)
    land_use_data.index = scenario_name_to_multiindex(land_use_data.index)

    both_data = pd.DataFrame({"land_use": land_use_data, "cost": cost_data})
    return both_data


def scenario_name_to_multiindex(index):
    return index.map(wxyz).rename(["util", "wind", "roof", "offshore"])


def wxyz(scenario_name):
    roof, util, wind, offshore = tuple(int(x.split("-")[1]) for x in scenario_name.split(","))
    return (util, wind, roof, offshore)


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
        optimal_index = slopes.abs().idxmin()
        yield optimal_index, slopes.loc[optimal_index]
        current_index = optimal_index


def optimal_path_series(data, tech):
    slopes = pd.Series({index: slope for index, slope in optimal_path(data, tech)})
    if len(slopes.index) > 0:
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
                   slopes.iloc[0] if len(slopes.index) > 0 else None]
            for name, slopes in [("Utility-scale PV", util_slopes),
                                 ("Rooftop PV", roof_slopes),
                                 ("Offshore Wind", offshore_slopes),
                                 ("Onshore Wind", wind_slopes)]
        },
        index=["min", "max", "average", "first"]
    )


if __name__ == "__main__":
    technology(
        path_to_data=snakemake.input.results,
        land_use_factors=pd.Series(snakemake.params.land_factors).to_xarray().rename(index="techs"),
        path_to_output=snakemake.output[0]
    )
