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


def plot_wind_share(path_to_xy, path_to_plot):
    # monkey patch seaborn
    # * to allow smaller rhombs in boxenplot
    # * to ensure smallest box is _not_ white
    sns.categorical._LVPlotter._lvplot = _lvplot
    sns.categorical._LVPlotter._width_functions = _width_functions

    data = calculate_data(path_to_xy)
    fig = plt.figure(figsize=(7.5, 3))
    ax = fig.subplots(1, 1)

    sns.boxenplot(
        data=data,
        y="Max capacity share wind (%)",
        x="Land area limit (% total land)",
        hue="Technologies in addition to onshore wind",
        palette=sns.light_palette(GREY, n_colors=3, reverse=False)[1:],
        outlier_prop=0.01,
        scale="area",
        width=0.7,
        s=1.0,
        ax=ax
    )

    sns.despine()
    ax.legend(frameon=False)
    for rectangle in ax.get_legend().legendHandles:
        rectangle.set_linewidth(0)

    fig.tight_layout()
    fig.savefig(path_to_plot, pil_kwargs={"compression": "tiff_lzw"})


def calculate_data(path_to_xy_data):
    xy = xr.open_dataset(path_to_xy_data).reset_coords(["util", "wind", "roof", "offshore"])
    data = (
        xr
        .ones_like(xy.cost.sum("scenario"))
        .rename("share")
        .expand_dims(tech=["Utility-scale PV & onshore wind", "Rooftop PV, off- & onshore wind"], threshold=THRESHOLDS)
    ) * np.nan

    for tech in ["Utility-scale PV & onshore wind", "Rooftop PV, off- & onshore wind"]:
        if tech == "Utility-scale PV & onshore wind":
            tech_mask = (xy["roof"] == 0) & (xy["offshore"] == 0)
        else:
            tech_mask = xy["util"] == 0
        for threshold in THRESHOLDS:
            absolute_threshold = threshold * TOTAL_EUROPEAN_LAND_MASS_KM2
            mask = tech_mask & (xy.land_use <= absolute_threshold)
            data.loc[{"tech": tech, "threshold": threshold}] = xy.where(mask)["wind"].max("scenario")

    return (
        data
        .to_series()
        .reset_index()
        .assign(
            threshold=data.to_series().reset_index().threshold * 100, # to percent
        )
        .rename(columns={
            "threshold": "Land area limit (% total land)",
            "share": "Max capacity share wind (%)",
            "tech": "Technologies in addition to onshore wind"
        })
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
                    c='.15', alpha=.45)

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
    plot_wind_share(
        path_to_xy=snakemake.input.xy,
        path_to_plot=snakemake.output[0]
    )
