##############################################################################     ##    ######
#    A.J. Zwijnenburg                   2020-09-03           v1.0                 #  #      ##
#    Copyright (C) 2020 - AJ Zwijnenburg          MIT license                    ######   ##
##############################################################################  ##    ## ######

## Copyright notice ##########################################################
# Copyright 2020 AJ Zwijnenburg
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
.scatter()  -   Builds scatter plots of the VISION data

"""

from __future__ import annotations

from .data import Data

import pandas as pd
import plotnine as p9

p9.options.figure_size=(12.8, 9.6)

class Plot:
    """
    Provides a plotting interface into the VISION data
        :param session: the VISION session link
    """
    def __init__(self, session: str) -> None:
        self.data = Data(session)

    def scatter(self, projection:str, x: str, y: str, color: str, color_map: dict=None) -> p9.ggplot:
        """
        Creates a ggplot dotplot object with the correct data and axis
            :param projection: the projection to use for plotting
            :param x: the x dimension
            :param y: the y dimension
            :param color: the parameter to use for color mapping. Searches entire dataset for fitting parameter
            :param color_map: only used for factorized color parameters. Uses the color_map to map the levels
        """
        plot_data = self.data.projection[projection]

        if x not in plot_data or y not in plot_data:
            raise ValueError(f"x '{x}' or y '{y}' are not dimension of projection '{projection}''")

        if color:
            # Vision has two color options
            full_color_range: bool = True
            diverging_colormap: bool = False

            # factors need to be handled differently
            is_factor: bool = False

            if color in self.data.expression.names:
                plot_data[color] = self.data.expression[color]
                full_color_range = True
                diverging_colormap = False
            elif color in self.data.meta.names:
                plot_data[color] = self.data.meta[color]
                full_color_range = False
                diverging_colormap = True

                # Check if factor
                if color in self.data.meta.levels:
                    is_factor = True

            elif color in self.data.protein.names:
                plot_data[color] = self.data.protein[color]
                full_color_range = False
                diverging_colormap = True
            elif color in self.data.signature.names:
                plot_data[color] = self.data.signature[color]
                full_color_range = False
                diverging_colormap = True
            else:
                raise ValueError(f"color parameter '{color}' couldnt be found in the dataset")
                
        plot = p9.ggplot(
            plot_data, 
            p9.aes(x, y, color=color)
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
            panel_grid_major_x=p9.element_line(color="#DFDFDFFF"),
            panel_grid_major_y=p9.element_line(color="#DFDFDFFF"),
            panel_grid_minor_x=p9.element_blank(),
            panel_grid_minor_y=p9.element_blank(),
            panel_background=p9.element_rect(fill="#EEEEEEFF", color="#FFFFFFFF"),
            legend_title=p9.element_blank(),
            legend_key=p9.element_blank(),
            legend_key_width=8,
            legend_key_height=35,
            legend_entry_spacing_x=-10,
            legend_entry_spacing_y=-20
        )
        
        if color:
            plot = plot + p9.ggtitle(
                color
            )
        else:
            plot = plot + p9.ggtitle(
                self.data.api.session_info.name
            )
        
        plot = plot + p9.labs(
            x=f"{projection}: {x}", 
            y=f"{projection}: {y}"
        )
        plot = plot + p9.geom_point(
            size=1
        )

        if is_factor:
            if color_map:
                for level in self.data.meta.levels[color]:
                    if level not in color_map:
                        raise ValueError(f"level '{level}' of color '{color}' not found in color_map")

                plot = plot + p9.scales.scale_color_manual(
                    values = color_map
                )

            elif len(self.data.meta.levels[color]) <= 10:
                # hardcode the tab colorscales, as i dont know how else i can use them for discrete scales... scale_color_cmap doesnt work...
                tab10 = ["#1f77b4","#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
                tab20 = ["#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"]

                plot = plot + p9.scales.scale_color_manual(
                    values = tab10
                )
           
            else:
                pass
            
            plot = plot + p9.theme(
                legend_background=p9.element_rect(color="#FFFFFFFF", size=5)
            )
        
        elif color:
            quantiles = plot_data[color].quantile([0.0, 0.02, 0.98, 1.0])
            if full_color_range:
                min_color = quantiles[0.0]
                max_color = quantiles[1.0]
            else:
                min_color = quantiles[0.02]
                max_color = quantiles[0.98]

            if diverging_colormap:
                plot = plot + p9.scales.scale_color_cmap(
                    cmap_name="viridis",
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

        return plot
