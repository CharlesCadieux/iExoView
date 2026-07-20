"use strict";

var data_available = source_available.data;
var data_visible = source_visible.data;

Object.keys(data_available).forEach(function (prop) {
    data_visible[prop] = data_available[prop];
});

x_axis.value = "Period (days)";
y_axis.value = "Radius (R⊕)";
color_select.value = default_color_label;
search.value = "";
ra_min_slider.value = 0;
ra_max_slider.value = 24;
dec_slider.value = [-90, 90];
teff_slider.value = [2500, 15000];
rade_slider.value = [0.1, 40];
checkbox.active = [0, 1, 2];

// Re-point the active colormap field/legend back to the default parameter.
data_visible.bokeh_colors = data_visible['bokeh_colors_' + default_color_key];
color_bar.color_mapper.low = color_ranges[default_color_key][0];
color_bar.color_mapper.high = color_ranges[default_color_key][1];
color_bar.title = default_color_label;

source_visible.data = data_visible;
source_visible.change.emit();
p.reset.emit();
title.text = titlestr;

inverted.active = [];