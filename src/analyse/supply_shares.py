import xarray as xr

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec

GREY = "#7F7F7F"
GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"


def visualise_supply_shares(path_to_aggregated_results, path_to_plot):
    scenarios = xr.open_dataset(path_to_aggregated_results).scenario

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 3))
    gs = gridspec.GridSpec(2, 4, height_ratios=[30, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[0, 3])
    cbar_ax1 = fig.add_subplot(gs[1, 0])
    cbar_ax2 = fig.add_subplot(gs[1, 1])
    cbar_ax3 = fig.add_subplot(gs[1, 2])
    cbar_ax4 = fig.add_subplot(gs[1, 3])

    sns.heatmap(
        scenarios.wind.values.reshape(22, 13),
        cmap=sns.light_palette(GREY, 11),
        xticklabels=False,
        yticklabels=False,
        ax=ax1,
        cbar=True,
        cbar_ax=cbar_ax1,
        cbar_kws={"orientation": "horizontal"}
    )
    sns.heatmap(
        scenarios.offshore.values.reshape(22, 13),
        cmap=sns.light_palette(BLUE, 11),
        xticklabels=False,
        yticklabels=False,
        ax=ax2,
        cbar=True,
        cbar_ax=cbar_ax2,
        cbar_kws={"orientation": "horizontal"}
    )
    sns.heatmap(
        scenarios.util.values.reshape(22, 13),
        cmap=sns.light_palette(RED, 11),
        xticklabels=False,
        yticklabels=False,
        ax=ax3,
        cbar=True,
        cbar_ax=cbar_ax3,
        cbar_kws={"orientation": "horizontal"}
    )
    sns.heatmap(
        scenarios.roof.values.reshape(22, 13),
        cmap=sns.light_palette(GREEN, 11),
        xticklabels=False,
        yticklabels=False,
        ax=ax4,
        cbar=True,
        cbar_ax=cbar_ax4,
        cbar_kws={"orientation": "horizontal"}
    )

    for ax in [ax1, ax4]:
        colorbar = ax.collections[0].colorbar
        colorbar.set_ticks([(x + 1 / 2) * 100 / 11 for x in [0, 5, 10]])
        colorbar.set_ticklabels(["0%", "50%", "100%"])
    for ax in [ax2, ax3]:
        colorbar = ax.collections[0].colorbar
        colorbar.set_ticks([(x + 1 / 2) * 100 / 11 for x in [2, 8]])
        colorbar.set_ticklabels(["20%", "80%"])

    ax1.set_title("Onshore wind")
    ax2.set_title("Offshore wind")
    ax3.set_title("Utility-scale PV")
    ax4.set_title("Rooftop PV")

    cbar_ax1.set_xlabel('Supply share')
    cbar_ax2.set_xlabel('Supply share')
    cbar_ax3.set_xlabel('Supply share')
    cbar_ax4.set_xlabel('Supply share')

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.1, hspace=0.1)
    fig.savefig(path_to_plot, dpi=300)


if __name__ == "__main__":
    visualise_supply_shares(
        path_to_aggregated_results=snakemake.input.results,
        path_to_plot=snakemake.output[0]
    )
