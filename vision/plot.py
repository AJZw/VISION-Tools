##############################################################################     ##    ######
#    A.J. Zwijnenburg                   2020-09-21           v1.2                 #  #      ##
#    Copyright (C) 2023 - AJ Zwijnenburg          MIT license                    ######   ##
##############################################################################  ##    ## ######

## Copyright notice ##########################################################
# Copyright 2023 AJ Zwijnenburg
#
# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal 
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in  
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE # WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
# THE SOFTWARE.
##############################################################################

"""
Provides a VISION look-a-like plotting interface into the VISION data

:class: Plot
Convenience interface using Data and API for generating VISION-like plots
.mask           -   Adds a boolean mask to 

.scatter()      -   Builds scatter plots of the VISION data
.violin()       -   Builds violin plots of the VISION data
.projection()   -   Builds scatter plots of the specified projection
.comparison()   -   Builds scatter plots comparing the x, and y property
.bar()          -   Builds a bar graph
.bar_stacked()  -   Builds a stacked bar graph of a categorical value

"""

from __future__ import annotations

from .data import Data

import pandas as pd
import numpy as np
import plotnine as p9
import copy

#p9.options.figure_size=(12.8, 9.6)
p9.options.figure_size=(6.4, 4.8)

class Plot:
    """
    Provides a plotting interface into the VISION data
        :param session: the VISION session link
    """
    # Additional discrete colorscale
    tab10 = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
    #tab20 = ["#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"]
    tab20 = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a",  "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"]

    def __init__(self, session: str) -> None:
        self.data = Data(session)

        self._mask: np.Array[bool] = None

    @property
    def mask(self) -> np.Array[bool]:
        """
        Returns the current boolean mask or, if no mask has been added, returns None
        """
        return self._mask

    @mask.setter
    def mask(self, mask: np.Array[bool]):
        """
        Set a boolean mask to the plot. Is used to filter the output data before plotting.
        If set to None, mask is disabled
        This does not change the x and y limits
        """
        if mask is None:
            self._mask = None

        if len(mask) != len(self.data.cells):
            raise ValueError(f"length of mask {len(mask)} unequal to rows in data {len(self.data.cells)}")

        self._mask = np.array(mask)

    def scatter(self, data: pd.DataFrame, x: str, y: str, c: str=None, c_map: dict=None) -> p9.ggplot:
        """
        Builds a basic scatter plot from the data
            :param data: the dataframe containing parameter x and y
            :param x: the x dimension
            :param y: the y dimension 
            :param c: the c(olor) dimension
            :param c_map: only used for factorized color parameters. Uses the color_map to map the levels
        """
        plot = p9.ggplot(
            data, 
            p9.aes(x, y, color=c)
        )

        plot = plot + p9.theme_bw() + p9.theme(
            text=p9.element_text(family="sans-serif", weight="normal"),
            plot_title=p9.element_text(ha="center", weight="bold", size=14),
            axis_text_x=p9.element_text(ha="center", va="top"),
            axis_text_y=p9.element_text(ha="right", va="center"),
            axis_ticks_major_x=p9.element_blank(),
            axis_ticks_major_y=p9.element_blank(),
            axis_ticks_minor_x=p9.element_blank(),
            axis_ticks_minor_y=p9.element_blank(),
            #panel_grid_major_x=p9.element_line(color="#DFDFDFFF"),
            #panel_grid_major_y=p9.element_line(color="#DFDFDFFF"),
            panel_grid_major_x=p9.element_blank(),
            panel_grid_major_y=p9.element_blank(),
            panel_grid_minor_x=p9.element_blank(),
            panel_grid_minor_y=p9.element_blank(),
            panel_background=p9.element_rect(fill="#F8F8F8FF", color="#FFFFFFFF"),
            legend_title=p9.element_blank(),
            legend_key=p9.element_blank(),
            legend_key_width=8,
            legend_key_height=35,
            legend_entry_spacing_x=-10,
            legend_entry_spacing_y=-20
        )

        if c:
            plot = plot + p9.ggtitle(
                c
            )
        else:
            plot = plot + p9.ggtitle(
                self.data.api.session_info.name
            )

        plot = plot + p9.labs(
            x=x, 
            y=y
        )

        # Check if x or y are factors if so, add jitter in geom_point
        width=0
        height=0
        if x in self.data.meta.names and x in self.data.meta.levels:
            width = 0.4
        if y in self.data.meta.names and y in self.data.meta.levels:
            height = 0.4
        if width != 0 or height != 0:
            position = p9.positions.position_jitter(width=width, height=height)
        else:
            position = p9.positions.position_identity()

        plot = plot + p9.geom_point(
            position=position,
            size=1
        )
        
        # Get color properties
        full_color_range: bool = True
        diverging_colormap: bool = False
        is_factor: bool = False
        if c:
            if c in self.data.expression.names:
                full_color_range = True
                diverging_colormap = False
            elif c in self.data.meta.names:
                full_color_range = False
                diverging_colormap = True
                # Check if factor
                if c in self.data.meta.levels:
                    is_factor = True
            elif c in self.data.protein.names:
                full_color_range = False
                diverging_colormap = True
            elif c in self.data.signature.names:
                full_color_range = False
                diverging_colormap = True
            else:
                raise ValueError(f"color parameter '{c}' couldnt be found in the dataset")

        if is_factor:
            if c_map:
                for level in self.data.meta.levels[c]:
                    if level not in c_map:
                        raise ValueError(f"level '{level}' of color '{c}' not found in c_map")
                plot = plot + p9.scales.scale_color_manual(
                    values = c_map
                )
            elif len(self.data.meta.levels[c]) <= 10:
                # hardcode the tab colorscales, as i dont know how else i can use them for discrete scales... scale_color_cmap doesnt work...
                plot = plot + p9.scales.scale_color_manual(
                    values = self.tab10
                )
            elif len(self.data.meta.levels[c]) <= 20:
                # hardcode the tab colorscales, as i dont know how else i can use them for discrete scales... scale_color_cmap doesnt work...
                plot = plot + p9.scales.scale_color_manual(
                    values = self.tab20
                )
            else:
                # Use default - (dont care much and) couldnt quickly see which colormap is used for >10 discrete values
                pass
            
            plot = plot + p9.theme(
                legend_background=p9.element_rect(color="#FFFFFFFF", size=5)
            )
        
        elif c:
            quantiles = data[c].quantile([0.0, 0.02, 0.98, 1.0])
            if full_color_range:
                min_color = quantiles[0.0]
                max_color = quantiles[1.0]
            else:
                min_color = quantiles[0.02]
                max_color = quantiles[0.98]

            if diverging_colormap:
                plot = plot + p9.scales.scale_color_cmap(
                    cmap_name="viridis",
                    #cmap_name="hsv", - TRIcycle
                    limits=(min_color, max_color),
                    guide=p9.guide_colorbar(
                        ticks=False
                    )
                )
            else:
                plot = plot + p9.scales.scale_color_gradientn(
                    colors=(
                        "#d8d8d8",
                        "#395252",
                        "#000000"
                    ),
                    values=(
                        0.0,
                        0.5,
                        1.0
                    ),
                    limits=(min_color, max_color),
                    guide=p9.guide_colorbar(
                        ticks=False
                    )
                )

        if self.mask is not None:
            plot.data = plot.data.loc[self.mask]

        return plot

    def violin(self, data: pd.DataFrame, x: str, y: str, x_map: dict=None) -> p9.ggplot:
        """
        Builds a basic violin/boxplot plot from the data
            :param data: the dataframe containing parameter x and y
            :param x: the x (factorized) dimension
            :param y: the y dimension 
            :param x_map: only used for x dimension. Uses the x_map to map the levels
        """
        plot = p9.ggplot(
            data, 
            p9.aes(x=x, y=y, fill=x)
        )

        plot = plot + p9.theme_bw() + p9.theme(
            text=p9.element_text(family="sans-serif", weight="normal"),
            plot_title=p9.element_text(ha="center", weight="bold", size=14),
            axis_text_x=p9.element_text(rotation=90, ha="center", va="top"),
            axis_text_y=p9.element_text(ha="right", va="center"),
            axis_ticks_major_x=p9.element_blank(),
            axis_ticks_major_y=p9.element_blank(),
            axis_ticks_minor_x=p9.element_blank(),
            axis_ticks_minor_y=p9.element_blank(),
            panel_grid_major_x=p9.element_line(color="#DFDFDFFF"),
            panel_grid_major_y=p9.element_line(color="#DFDFDFFF"),
            panel_grid_minor_x=p9.element_blank(),
            panel_grid_minor_y=p9.element_blank(),
            panel_background=p9.element_rect(fill="#EEEEEEFF", color="#FFFFFFFF"),
            legend_position="none"
        )
        plot += p9.geom_violin()
        plot += p9.geom_boxplot(p9.aes(x=x, y=y), fill=None, width=0.1, inherit_aes=False)

        if x_map:
            plot += p9.scales.scale_fill_manual(values=x_map)

        if self.mask is not None:
            plot.data = plot.data.loc[self.mask]

        return plot

    def projection(self, projection: str, x: str, y: str, c: str=None, c_map: dict=None) -> p9.ggplot:
        """
        Creates a ggplot geom_point object with the correct data and axis.
        This plots a projection x-y with color overlay according to VISION
            :param projection: the projection to use for plotting
            :param x: the x dimension
            :param y: the y dimension
            :param c: the parameter to use for color mapping. Searches entire dataset for fitting parameter
            :param c_map: only used for factorized color parameters. Uses the color_map to map the levels
        """
        plot_data = copy.deepcopy(self.data.projection[projection])
        
        if x not in plot_data or y not in plot_data:
            raise ValueError(f"x '{x}' or y '{y}' are not dimension of projection '{projection}''")
        
        plot_data[c] = self.data[c]

        plot = self.scatter(plot_data, x, y, c, c_map)

        plot = plot + p9.labs(
            x=f"{projection}: {x}", 
            y=f"{projection}: {y}"
        )
        
        return plot

    def comparison(self, x: str, y: str, c: str=None, c_map: dict=None, jitter: bool=False) -> p9.ggplot:
        """
        Creates a ggplot geom_point object with the correct data and axis.
        This plots x-y data (expression/signature) with color overlay
            :param projection: the projection to use for plotting
            :param x: the x dimension
            :param y: the y dimension
            :param c: the parameter to use for color mapping. Searches entire dataset for fitting parameter
            :param color_map: only used for factorized color parameters. Uses the color_map to map the levels
            :param jitter: whether to add a random jitter for 0 values
        """
        plot_data = copy.deepcopy(self.data[[x, y, c]])

        if jitter:
            # Check if factor
            if x in self.data.meta.names and x in self.data.meta.levels:
                # jitter is taken care off automatically .in scatter()
                pass
            else:
                is_zero = plot_data[x] == 0
                jitter = np.random.uniform(low=-0.05, high=0.05, size=sum(is_zero))
                plot_data.loc[is_zero, x] += jitter
            
            # Check if factor
            if y in self.data.meta.names and y in self.data.meta.levels:
                # jitter is taken care off automatically .in scatter()
                pass
            else:
                is_zero = plot_data[y] == 0
                jitter = np.random.uniform(low=-0.05, high=0.05, size=sum(is_zero))
                plot_data.loc[is_zero, y] += jitter

        plot = self.scatter(plot_data, x, y, c, c_map)

        return plot

    def bar(self, x: str, log: bool=False) -> p9.ggplot:
        """
        Builds a bar graph
            :param x: the x-dimension
            :param log: whether to log-transform the y-axis
        """
        plot_data = copy.deepcopy(self.data[x])

        bin_count = 40
        is_factor = False

        # Determine if factor
        if x in self.data.meta.names:
            if x in self.data.meta.levels:
                is_factor = True
                bin_range = self.data.meta.levels[x]
                bin_count = len(bin_range)

        # Do necessary binning:
        # Mainly to allow for manual log-scaling, to be able to show 0-values
        if not is_factor:
            plot_data["__bin"], bin_range = pd.cut(
                plot_data[x],
                bins=bin_count,
                right=False,
                include_lowest=True,
                retbins=True,
                labels=False
            )
            # now set the bins to the middle of the bin_range for proper x-visualization
            bin_value = {}
            bin_width = (bin_range[1] - bin_range[0])
            for i in range(0, len(bin_range)-1):
                bin_value[i] = ((i+0.5)*bin_width) + bin_range[0]

            plot_data["__bin"] = plot_data["__bin"].apply(lambda x: bin_value[x])

            bin_data = pd.DataFrame(pd.value_counts(plot_data["__bin"]))
            bin_data.columns = ["__count"]
            bin_data["__bin"] = bin_data.index
        else:
            bin_data = pd.DataFrame(pd.value_counts(plot_data[x]))
            bin_data.columns = ["__count"]
            bin_data["__bin"] = bin_data.index
        
        plot = p9.ggplot(
            bin_data, 
            p9.aes(
                x="__bin",
                y="__count"
            )
        )

        # hides axis titles
        plot = plot + p9.theme_bw() + p9.theme(
            text=p9.element_text(family="sans-serif", weight="normal"),
            plot_title=p9.element_text(ha="center", weight="bold", size=14),
            axis_title_x=p9.element_blank(),
            axis_title_y=p9.element_blank(),
            axis_text_x=p9.element_text(ha="center", va="top"),
            axis_text_y=p9.element_text(ha="right", va="center"),
            axis_line_x=p9.element_line(color="#000000FF"),
            axis_line_y=p9.element_blank(),
            axis_ticks_major_x=p9.element_blank(),
            axis_ticks_major_y=p9.element_blank(),
            axis_ticks_minor_x=p9.element_blank(),
            axis_ticks_minor_y=p9.element_blank(),
            panel_grid_major_x=p9.element_blank(),
            panel_grid_major_y=p9.element_line(color="#DFDFDFFF"),
            panel_grid_minor_x=p9.element_blank(),
            panel_grid_minor_y=p9.element_blank(),
            panel_background=p9.element_rect(fill="#FFFFFFFF", color="#FFFFFFFF"),
            panel_border=p9.element_blank(),
            legend_title=p9.element_blank(),
            legend_key=p9.element_blank(),
            legend_key_width=8,
            legend_key_height=35,
            legend_entry_spacing_x=-10,
            legend_entry_spacing_y=-20
        )

        plot = plot + p9.ggtitle(
            x
        )

        # Set x & y scales
        if not is_factor:
            min_range = min(plot_data[x])
            max_range = max(plot_data[x])
            plot = plot + p9.scale_x_continuous(
                limits=(min_range, max_range),
                expand=(0,0)
            )

        if log:
            # Build a fake log-scale
            plot.data["__count"] = plot.data["__count"].apply(lambda x: np.log10(x + 1))

            # Label mapping
            labels = pd.Series(["0", "1", "2", "5", "10", "20", "50", "100", "200", "500", "1 000", "2 000", "5 000", "10 000", "20 000", "50 000", "100 000", "200 000", "500 000", "1 000 000"])
            values = pd.Series([0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000, 500_000, 1_000_000])
            values = values.apply(lambda x: np.log10(x + 1))

            # Remove unused labels
            cutoff = max(plot.data["__count"]) * 1.05
            for i in range(0, len(values)):
                if values[i] > cutoff:
                    values = values[:i]
                    labels = labels[:i]
                    break

            plot = plot + p9.scale_y_continuous(
                breaks=values,
                labels=labels,
                expand=(0, 0, 0.05, 0)
            ) 
        else:
            plot = plot + p9.scale_y_continuous(
                expand=(0, 0, 0.05, 0)
            )

        # Set binning
        if is_factor:
            plot = plot + p9.geom_bar(
                stat="identity",
                fill="#1F77B4"
            )
        else:
            plot = plot + p9.geom_bar(
                stat="identity",
                fill="#1F77B4",
                na_rm=True
            )

        return 
        
    def bar_stacked(self, x: str, y: str, y_map: dict=None, fraction: bool=True) -> p9.ggplot:
        """
        Builds a bar graph
            :param x: the x-dimension
            :param y: the y-dimension (categorical)
            :param y_map: (optional) the color_map to use for the y-dimension
            :param fraction: whether to report the fraction (instead of absolute counts) / x dimension
        """
        if y not in self.data.meta.names or y not in self.data.meta.levels:
            raise ValueError(f"the y-dimension '{y}' has to be a categorical parameter")

        plot_data = copy.deepcopy(self.data[[x, y]])

        bin_count = 40
        is_factor = False

        # Determine if factor
        if x in self.data.meta.names:
            if x in self.data.meta.levels:
                is_factor = True
                bin_range = self.data.meta.levels[x]
                bin_count = len(bin_range)

        # Do necessary binning:
        if not is_factor:
            plot_data["__bin"], bin_range = pd.cut(
                plot_data[x],
                bins=bin_count,
                right=False,
                include_lowest=True,
                retbins=True,
                labels=False
            )
            # now set the bins to the middle of the bin_range for proper x-visualization
            bin_value = {}
            bin_width = (bin_range[1] - bin_range[0])
            for i in range(0, len(bin_range)-1):
                bin_value[i] = ((i+0.5)*bin_width) + bin_range[0]

            plot_data["__bin"] = plot_data["__bin"].apply(lambda x: bin_value[x])
        else:
            plot_data["__bin"] = plot_data[x]

        # Generate x&y bins
        plot_data["__xy"] = plot_data.loc[:,["__bin", y]].apply(lambda x: '__'.join(x), axis=1)

        # Count bins
        bin_data = pd.value_counts(plot_data["__xy"])

        bin_data = pd.concat([
            bin_data,
            pd.DataFrame(bin_data.index.str.split("__").to_list(), index=bin_data.index, columns=["__x", "__y"]),
        ], axis=1)

        # Fractionize
        if fraction:
            output = []
            for key, group in bin_data.groupby("__x"):
                total = sum(group["__xy"])
                group["__xy"] = group["__xy"].apply(lambda x: x/total)
                output.append(group)
            bin_data = pd.concat(output)

        plot = p9.ggplot(
            bin_data
        )

        # hides axis titles
        plot = plot + p9.theme_bw() + p9.theme(
            text=p9.element_text(family="sans-serif", weight="normal"),
            plot_title=p9.element_text(ha="center", weight="bold", size=14),
            axis_title_x=p9.element_blank(),
            axis_title_y=p9.element_blank(),
            axis_text_x=p9.element_text(ha="center", va="top"),
            axis_text_y=p9.element_text(ha="right", va="center"),
            axis_line_x=p9.element_line(color="#000000FF"),
            axis_line_y=p9.element_blank(),
            axis_ticks_major_x=p9.element_blank(),
            axis_ticks_major_y=p9.element_blank(),
            axis_ticks_minor_x=p9.element_blank(),
            axis_ticks_minor_y=p9.element_blank(),
            panel_grid_major_x=p9.element_blank(),
            panel_grid_major_y=p9.element_line(color="#DFDFDFFF"),
            panel_grid_minor_x=p9.element_blank(),
            panel_grid_minor_y=p9.element_blank(),
            panel_background=p9.element_rect(fill="#FFFFFFFF", color="#FFFFFFFF"),
            panel_border=p9.element_blank(),
            legend_title=p9.element_blank(),
            legend_key=p9.element_blank(),
            legend_key_width=8,
            legend_key_height=35,
            legend_entry_spacing_x=-10,
            legend_entry_spacing_y=-20
        )

        plot = plot + p9.ggtitle(
            x
        )

        # Set x scale
        if not is_factor:
            min_range = min(plot_data[x])
            max_range = max(plot_data[x])
            plot = plot + p9.scale_x_continuous(
                limits=(min_range, max_range),
                expand=(0,0)
            )

        # Set y scale 
        if fraction:
            max_y = 1.0
        else:
            max_y = 0
            for key, group in bin_data.groupby("__x"):
                total = sum(group["__xy"])
                if total > max_y:
                    max_y = total

        plot = plot + p9.scale_y_continuous(
            limits=(0, max_y),
            expand=(0, 0, 0.05, 0)
        )

        # Set fill-scaling
        if y_map:
            for level in self.data.meta.levels[y]:
                if level not in y_map:
                    raise ValueError(f"level '{level}' of y-dimension '{y}' not found in y_map")
            plot = plot + p9.scales.scale_fill_manual(
                values = y_map
            )
        elif len(self.data.meta.levels[y]) <= 10:
            # hardcode the tab colorscales, as i dont know how else i can use them for discrete scales... scale_color_cmap doesnt work...
            plot = plot + p9.scales.scale_fill_manual(
                values = self.tab10
            )
        else:
            # Use default - (dont care much and) couldnt quickly see which colormap is used for >10 discrete values
            pass
        
        # Add bar graphs
        if is_factor:
            plot = plot + p9.geom_bar(
                p9.aes(
                    x="__x",
                    y="__xy",
                    fill="__y"
                ),
                stat="stat_identity",
                inherit_aes=False
            )
        else:
            plot = plot + p9.geom_bar(
                p9.aes(
                    x="__x",
                    y="__xy",
                    fill="__y"
                ),
                stat="stat_identity",
                inherit_aes=False,
                na_rm=True
            )

        return plot
