import pandas as pd
import seaborn as sns
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


def uncertainty_plot(path_to_xy_data, path_to_plot):
    data = pd.read_csv(path_to_xy_data, index_col=0)
    data = pd.DataFrame({
        "util": data["r50_cost_util"] * -1,
        "offshore": data["r50_cost_offshore"] * -1
    }).dropna()
    sns.set_context("paper")
    fig = plt.figure(figsize=(8, 5))
    gs = gridspec.GridSpec(2, 4, width_ratios=[9, 1, 2, 6], height_ratios=[1, 10])
    ax_joint = fig.add_subplot(gs[1, 0])
    ax_marginal_top = fig.add_subplot(gs[0, 0], sharex=ax_joint)
    ax_marginal_right = fig.add_subplot(gs[1, 1], sharey=ax_joint)
    ax_diff = fig.add_subplot(gs[1, 3])

    plot_joint(data, ax_joint, ax_marginal_top, ax_marginal_right)
    plot_diff(data, ax_diff)

    ax_joint.annotate(
        "a - Joint cost distribution",
        xy=[-0.08, 1.2],
        xycoords='axes fraction',
        fontsize=PANEL_FONT_SIZE,
        weight=PANEL_FONT_WEIGHT
    )
    ax_diff.annotate(
        "b - Cost difference distribution",
        xy=[-0.25, 1.2],
        xycoords='axes fraction',
        fontsize=PANEL_FONT_SIZE,
        weight=PANEL_FONT_WEIGHT
    )
    fig.tight_layout()
    fig.subplots_adjust(
        left=0.1,
        top=0.9,
        wspace=0.1,
        hspace=0.1
    )

    fig.savefig(path_to_plot, dpi=300)


def plot_joint(data, ax_joint, ax_marginal_top, ax_marginal_right):
    sns.kdeplot(
        data2=data["util"],
        data=data["offshore"],
        shade=True,
        shade_lowest=False,
        cmap=sns.light_palette(CONTRAST_COLOR, as_cmap=True),
        gridsize=50,
        levels=10,
        ax=ax_joint
    )
    ax_joint.plot((0, 1.1), (0, 1.1), color=GREY, linestyle=":")

    sns.distplot(
        data["offshore"],
        color=CONTRAST_COLOR,
        hist=False,
        ax=ax_marginal_top
    )
    lines = ax_marginal_top.lines[0]
    x = lines.get_xydata()[:, 0]
    y = lines.get_xydata()[:, 1]
    ax_marginal_top.fill_between(x, y, color=CONTRAST_COLOR, alpha=0.3)
    sns.distplot(
        data["util"],
        color=CONTRAST_COLOR,
        vertical=True,
        hist=False,
        ax=ax_marginal_right
    )
    lines = ax_marginal_right.lines[0]
    x = lines.get_xydata()[:, 0]
    y = lines.get_xydata()[:, 1]
    ax_marginal_right.fill_between(x, y, color=CONTRAST_COLOR, alpha=0.3)

    sns.despine(ax=ax_marginal_top, left=True)
    sns.despine(ax=ax_marginal_right, top=True, bottom=True, right=True)
    sns.despine(ax=ax_joint)
    ax_marginal_top.set_xlabel("")
    ax_marginal_top.set_yticks([])
    for tk in ax_marginal_top.get_xticklabels():
        tk.set_visible(False)

    ax_marginal_right.set_ylabel("")
    ax_marginal_right.set_xticks([])
    for tk in ax_marginal_right.get_yticklabels():
        tk.set_visible(False)
    ax_joint.set_xlim(0, 1.1)
    ax_joint.set_ylim(0, 1.1)

    ax_joint.set_ylabel(r"Average cost utility-scale PV $\mathrm{\left(\frac{EUR}{m^2 \cdot yr}\right)} $")
    ax_joint.set_xlabel(r"Average cost offshore wind $\mathrm{\left(\frac{EUR}{m^2 \cdot yr}\right)} $")


def plot_diff(data, ax):
    sns.distplot(
        data["util"] - data["offshore"],
        ax=ax,
        kde=True,
        hist=False,
        color=CONTRAST_COLOR
    )
    ax.vlines(
        x=0,
        ymin=0,
        ymax=ax.get_ylim()[1],
        color=GREY,
        linestyle=":",
    )
    lines = ax.lines[0]
    x = lines.get_xydata()[:, 0]
    y = lines.get_xydata()[:, 1]
    ax.fill_between(x, y, color=CONTRAST_COLOR, alpha=0.3)

    sns.despine(ax=ax)
    ax.set_xlabel(
        "Extra cost utility-scale PV\n"
        + r"compared to offshore wind $\mathrm{\left(\frac{EUR}{m^2 \cdot yr}\right)} $"
    )
    ax.set_ylabel("Density")


if __name__ == "__main__":
    uncertainty_plot(
        path_to_xy_data=snakemake.input.xy,
        path_to_plot=snakemake.output[0]
    )
