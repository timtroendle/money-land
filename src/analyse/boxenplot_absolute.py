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
        data
        .loc[:, [col for col in data if col.startswith('a')]]
        .mul(100) # to percent
        .unstack()
        .reset_index()
        .rename(columns={"level_0": "scenario", 0: "cost"})
    )
    clean_data = (
        clean_data
        .assign(
            supply_technology=[tech.split("_")[-1] for tech in clean_data.scenario],
            threshold=[float(tech.split("_")[0][1:]) / 10 for tech in clean_data.scenario])
        .rename(columns={"threshold": "Land use limit (%)", "supply_technology": "Supply technology"})
        .replace({"offshore": "Offshore wind", "util": "Utility-scale PV", "roof": "Rooftop PV"})
    )

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 4))
    ax = fig.subplots(1, 1)

    sns.boxenplot(
        data=clean_data,
        x="Land use limit (%)",
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


if __name__ == "__main__":
    boxenplot(
        path_to_xy_data=snakemake.input.xy,
        path_to_plot=snakemake.output[0]
    )
