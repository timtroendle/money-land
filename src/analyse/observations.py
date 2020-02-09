import seaborn as sns
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import gridspec

GREY = "#7F7F7F"
GREEN = "#679436"
RED = "#A01914"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"
CONTRAST_COLOR = BLUE
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"

TOTAL_EUROPEAN_LAND_MASS_KM2 = 4_920_000 # FIXME inject
TOTAL_DEMAND_KWH = 3_180_000_000_000 # FIXME inject


def plot_observations(path_to_xy, path_to_plot):
    ds = xr.open_dataset(path_to_xy)

    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 3))
    gs = gridspec.GridSpec(1, 2)
    ax_cost = fig.add_subplot(gs[0, 0])
    ax_land_use = fig.add_subplot(gs[0, 1], sharey=ax_cost)

    sns.distplot(
        ds.cost.to_series() / TOTAL_DEMAND_KWH,
        kde=False,
        color=RED,
        ax=ax_cost
    )
    ax_cost.set_xlabel("Cost (EUR / kWh)")
    ax_cost.set_ylabel("Frequency (millions)")
    ax_cost.annotate("a", xy=[-0.08, 1.05], xycoords='axes fraction',
                     fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)

    land_use = ds.land_use / TOTAL_EUROPEAN_LAND_MASS_KM2 * 100
    land_use = land_use.where(land_use < 5, drop=True) # drop upper 0.1%
    assert land_use.count() / ds.land_use.count() > 0.999
    sns.distplot(
        land_use.to_series(),
        kde=False,
        color=BLUE,
        ax=ax_land_use
    )
    ax_land_use.set_xlabel("Land requirements (% total land)")

    ax_land_use.annotate("b", xy=[-0.08, 1.05], xycoords='axes fraction',
                         fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)
    for tick in ax_land_use.yaxis.get_major_ticks():
        tick.set_visible(False)

    sns.despine(fig=fig)
    fig.tight_layout()
    ax_cost.yaxis.set_ticklabels([int(tick.get_text()) / 1e6 for tick in ax_cost.yaxis.get_ticklabels()])

    plt.subplots_adjust(
        left=0.11,
    )
    fig.savefig(path_to_plot, dpi=300)


if __name__ == "__main__":
    plot_observations(
        path_to_xy=snakemake.input.xy,
        path_to_plot=snakemake.output[0]
    )
