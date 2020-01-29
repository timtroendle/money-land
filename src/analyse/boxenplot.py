import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

GREY = "#7F7F7F"
GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"


def boxenplot(path_to_xy_data, path_to_plot):
    data = pd.read_csv(path_to_xy_data, index_col=0)
    clean_data = (
        data[["r25_cost_offshore", "r25_cost_util", "r25_cost_roof",
              "r50_cost_offshore", "r50_cost_util", "r50_cost_roof"]]
        .mul(-1)
        .dropna()
        .unstack()
        .reset_index()
        .rename(columns={"level_0": "scenario", 0: "cost"})
    )
    clean_data = (
        clean_data
        .assign(
            supply_technology=[tech.split("_")[-1] for tech in clean_data.scenario],
            reduction_level=[f'{tech.split("_")[0][1:]} %' for tech in clean_data.scenario])
        .rename(columns={
            "reduction_level": "Land reduction from cost-minimal case",
            "supply_technology": "Supply technology"})
        .replace({"offshore": "Offshore wind", "util": "Utility-scale PV", "roof": "Rooftop PV"})
    )

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)

    sns.boxenplot(
        data=clean_data,
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
    fig.savefig(path_to_plot, dpi=300)


if __name__ == "__main__":
    boxenplot(
        path_to_xy_data=snakemake.input.xy,
        path_to_plot=snakemake.output[0]
    )
