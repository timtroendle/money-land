import functools

import numpy as np
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt

import matplotlib as mpl
from matplotlib.collections import PatchCollection
import matplotlib.patches as Patches

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
    # monkey patch seaborn
    # * to allow smaller rhombs in boxenplot
    # * to handle zero width/height boxes
    # * to ensure smallest box is _not_ white
    sns.categorical._LVPlotter._lvplot = _lvplot
    sns.categorical._LVPlotter._width_functions = _width_functions

    data = calculate_data(path_to_xy_data)

    fig = plt.figure(figsize=(7.5, 3))
    ax = fig.subplots(1, 1)

    sns.boxenplot(
        data=data,
        x="Land area limit (% total land)",
        y="cost",
        hue="Supply technology",
        hue_order=["Offshore wind", "Utility-scale PV", "Rooftop PV"],
        palette=[BLUE, RED, GREEN],
        k_depth="proportion",
        outlier_prop=0.1,
        scale="area",
        s=1.0,
        ax=ax
    )
    ax.set_ylabel("Cost penalty relative to cost minimum (%)")
    sns.despine()
    ax.legend(frameon=False)
    for rectangle in ax.get_legend().legendHandles:
        rectangle.set_linewidth(0)

    fig.tight_layout()
    fig.savefig(path_to_plot, pil_kwargs={"compression": "tiff_lzw"})


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
        .rename(columns={"threshold": "Land area limit (% total land)", "tech": "Supply technology"})
        .replace({"offshore": "Offshore wind", "util": "Utility-scale PV", "roof": "Rooftop PV"})
    )


def _width_functions(self, width_func):
    # taken from mwaskom/seaborn, BSD-3 licensed
    # 2020-02-10: adapted to cover cases in which box height and width is 0
    # Dictionary of functions for computing the width of the boxes
    width_functions = {'linear': lambda h, i, k: (i + 1.) / k,
                       'exponential': lambda h, i, k: 2**(-k+i-1),
                       'area': lambda h, i, k: (1 - 2**(-k+i-2)) / h if h > 0 else 0}
    return width_functions[width_func]


def _lvplot(self, box_data, positions,
            color=[255. / 256., 185. / 256., 0.],
            vert=True, widths=1, k_depth='proportion',
            ax=None, outlier_prop=None, scale='exponential',
            **kws):

    x = positions[0]
    box_data = np.asarray(box_data)

    # If we only have one data point, plot a line
    if len(box_data) == 1:
        kws.update({'color': self.gray, 'linestyle': '-'})
        ys = [box_data[0], box_data[0]]
        xs = [x - widths / 2, x + widths / 2]
        if vert:
            xx, yy = xs, ys
        else:
            xx, yy = ys, xs
        ax.plot(xx, yy, **kws)
    else:
        # Get the number of data points and calculate "depth" of
        # letter-value plot
        box_ends, k = self._lv_box_ends(box_data, k_depth=k_depth,
                                        outlier_prop=outlier_prop)

        # Anonymous functions for calculating the width and height
        # of the letter value boxes
        width = self._width_functions(scale)

        # Function to find height of boxes
        def height(b):
            return b[1] - b[0]

        # Functions to construct the letter value boxes
        def vert_perc_box(x, b, i, k, w):
            rect = Patches.Rectangle((x - widths*w / 2, b[0]),
                                     widths*w,
                                     height(b), fill=True)
            return rect

        def horz_perc_box(x, b, i, k, w):
            rect = Patches.Rectangle((b[0], x - widths*w / 2),
                                     height(b), widths*w,
                                     fill=True)
            return rect

        # Scale the width of the boxes so the biggest starts at 1
        w_area = np.array([width(height(b), i, k)
                           for i, b in enumerate(box_ends)])
        # CATCH CASES WHERE BOX AREA IS ZERO
        if np.max(w_area) > 0:
            w_area = w_area / np.max(w_area)

        # Calculate the medians
        y = np.median(box_data)

        # Calculate the outliers and plot
        outliers = self._lv_outliers(box_data, k)
        hex_color = mpl.colors.rgb2hex(color)

        if vert:
            boxes = [vert_perc_box(x, b[0], i, k, b[1])
                     for i, b in enumerate(zip(box_ends, w_area))]

            # Plot the medians
            ax.plot([x - widths / 2, x + widths / 2], [y, y],
                    c='.15', alpha=.45)

            ax.scatter(np.repeat(x, len(outliers)), outliers,
                       marker='d', c=hex_color, **kws)
        else:
            boxes = [horz_perc_box(x, b[0], i, k, b[1])
                     for i, b in enumerate(zip(box_ends, w_area))]

            # Plot the medians
            ax.plot([y, y], [x - widths / 2, x + widths / 2],
                    c='.15', alpha=.45, **kws)

            ax.scatter(outliers, np.repeat(x, len(outliers)),
                       marker='d', c=hex_color, **kws)

        # Construct a color map from the input color
        rgb = sns.light_palette(hex_color, n_colors=10)[1:]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('new_map', rgb)
        collection = PatchCollection(boxes, cmap=cmap)

        # Set the color gradation
        collection.set_array(np.array(np.linspace(0, 1, len(boxes))))

        # Plot the boxes
        ax.add_collection(collection)


if __name__ == "__main__":
    boxenplot(
        path_to_xy_data=snakemake.input.xy,
        path_to_plot=snakemake.output[0]
    )
