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


REDUCTION_LEVELS = [0.1, 0.25, 0.5]
TECHS = ["roof", "util", "offshore"]
ALL_TECHS = TECHS + ["wind"]


def boxenplot(path_to_xy_data, path_to_plot):
    data = calculate_data(path_to_xy_data)

    fig = plt.figure(figsize=(7.5, 3.75))
    ax = fig.add_subplot(111)

    sns.boxenplot(
        data=data,
        y="Land reduction from cost-minimal case",
        hue="Supply technology",
        x="cost",
        palette=[BLUE, RED, GREEN],
        orient="h",
        outlier_prop=0.01,
        scale="area",
        ax=ax
    )
    sns.despine(fig=fig)

    ax.set_xlabel(r"Average cost to reduce land $\mathrm{\left(\frac{EUR}{m^2 \cdot yr}\right)} $")
    ax.legend(frameon=False, loc='center', bbox_to_anchor=(1., 0., 0.25, 1))
    for rectangle in ax.get_legend().legendHandles:
        rectangle.set_linewidth(0)

    fig.tight_layout()
    fig.savefig(path_to_plot, pil_kwargs={"compression": "tiff_lzw"})


def calculate_data(path_to_xy_data):
    xy = xr.open_dataset(path_to_xy_data)
    data = (
        xr
        .ones_like(xy.cost.sum("scenario"))
        .expand_dims(tech=TECHS, reduction_level=REDUCTION_LEVELS)
    ) * np.nan

    cost_optimal_data = xy.isel(scenario=xy.cost.argmin("scenario"))

    for tech in TECHS:
        conditions = [
            xy[other_tech] <= cost_optimal_data[other_tech]
            for other_tech in ALL_TECHS
            if other_tech != tech
        ]
        tech_mask = functools.reduce(lambda x, y: x & y, conditions)
        for reduction_level in REDUCTION_LEVELS:
            absolute_threshold = (1 - reduction_level) * cost_optimal_data["land_use"]
            mask = tech_mask & (xy.land_use <= absolute_threshold)
            delta_cost = xy.where(mask).cost - cost_optimal_data.cost
            delta_land_m2 = (xy.where(mask).land_use - cost_optimal_data.land_use) * 1e6
            rel_cost = delta_cost / delta_land_m2

            optimal_index = rel_cost.dropna("sample_id", how="all").argmax("scenario")
            optimal_cost = (
                rel_cost
                .sel(sample_id=optimal_index.sample_id)
                .isel(scenario=optimal_index)
                .reindex(sample_id=rel_cost.sample_id)
            )
            data.loc[{"tech": tech, "reduction_level": reduction_level}] = optimal_cost
    return (
        data
        .to_series()
        .mul(-1)
        .reset_index()
        .assign(reduction_level=[f'{level * 100:.0f} %' for level in data.to_series().reset_index().reduction_level])
        .rename(columns={
            "reduction_level": "Land reduction from cost-minimal case",
            "tech": "Supply technology"})
        .replace({"offshore": "Offshore wind", "util": "Utility-scale PV", "roof": "Rooftop PV"})
    )


if __name__ == "__main__":
    boxenplot(
        path_to_xy_data=snakemake.input.xy,
        path_to_plot=snakemake.output[0]
    )
