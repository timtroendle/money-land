import functools

import numpy as np
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt

GREY = "#7F7F7F"
GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"

TOTAL_EUROPEAN_LAND_MASS_KM2 = 4920000
THRESHOLDS = [0.005, 0.01, 0.015, 0.02, 0.03]
TECHS = ["roof", "util", "offshore"]
ALL_TECHS = TECHS + ["wind"]


def boxenplot(path_to_xy_data, path_to_plot):
    data = calculate_data(path_to_xy_data)

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 4))
    ax = fig.subplots(1, 1)

    sns.boxenplot(
        data=data,
        x="Land area limit (%)",
        y="cost",
        hue="Supply technology",
        hue_order=["Offshore wind", "Utility-scale PV", "Rooftop PV"],
        palette=[BLUE, RED, GREEN],
        outlier_prop=0.01,
        scale="area",
        ax=ax
    )
    ax.set_ylabel("Cost penalty relative\nto cost minimum (%)")
    sns.despine()
    ax.legend(frameon=False)
    for rectangle in ax.get_legend().legendHandles:
        rectangle.set_linewidth(0)

    fig.tight_layout()
    fig.savefig(path_to_plot, dpi=300)


def calculate_data(path_to_xy_data):
    xy = xr.open_dataset(path_to_xy_data)
    data = (
        xr
        .ones_like(xy.cost.sum("scenario"))
        .expand_dims(tech=TECHS, threshold=THRESHOLDS)
    ) * np.nan

    cost_optimal_data = xy.isel(scenario=xy.cost.argmin("scenario"))

    for tech in TECHS:
        conditions = [
            xy[other_tech] <= cost_optimal_data[other_tech]
            for other_tech in ALL_TECHS
            if other_tech != tech
        ]
        tech_mask = functools.reduce(lambda x, y: x & y, conditions)
        for threshold in THRESHOLDS:
            absolute_threshold = threshold * TOTAL_EUROPEAN_LAND_MASS_KM2
            mask = tech_mask & (xy.land_use <= absolute_threshold)
            cost = xy.where(mask).cost.min("scenario")
            data.loc[{"tech": tech, "threshold": threshold}] = (cost - cost_optimal_data.cost) / cost_optimal_data.cost
    missing_values = data.isnull()
    assert not missing_values.sel(tech="roof").any()
    assert not missing_values.sel(tech="util").any()
    assert missing_values.sel(tech="offshore", threshold=0.01).sum() < 20 # in rare cases, 1% can not be done
    assert missing_values.sel(tech="offshore", threshold=0.005).sum() < 1000 # in 1% of cases, 0.05% can not be done
    return (
        data
        .to_series()
        .reset_index()
        .assign(
            threshold=data.to_series().reset_index().threshold * 100, # to percent
            cost=data.to_series().reset_index().cost * 100 # to percent
        )
        .rename(columns={"threshold": "Land area limit (%)", "tech": "Supply technology"})
        .replace({"offshore": "Offshore wind", "util": "Utility-scale PV", "roof": "Rooftop PV"})
    )


if __name__ == "__main__":
    boxenplot(
        path_to_xy_data=snakemake.input.xy,
        path_to_plot=snakemake.output[0]
    )
