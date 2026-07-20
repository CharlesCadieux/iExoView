"use strict";

// Map a raw value to a Viridis hex color using the same log-normalization
// as the Python-side `logscale()` function.
function valueToColor(rawValue, low, high) {
    const value = parseFloat(rawValue);
    if (isNaN(value) || value <= 0 || low <= 0 || high <= 0) {
        return '#707070';
    }
    let t = (Math.log10(value) - Math.log10(low)) / (Math.log10(high) - Math.log10(low));
    t = Math.min(1, Math.max(0, t));
    const palette = color_bar.color_mapper.palette;
    const idx = Math.min(palette.length - 1, Math.floor(t * palette.length));
    return palette[idx];
}

const selected_label = cb_obj.value;
const key = axis_map[selected_label];
const data_visible = source_visible.data;

let low = color_ranges[key][0];
let high = color_ranges[key][1];

if (key === 'teff') {
    low = teff_slider.value[0];
    high = teff_slider.value[1];
    data_visible.bokeh_colors = data_visible.teff.map(v => valueToColor(v, low, high));
} else if (key === 'rade') {
    low = rade_slider.value[0];
    high = rade_slider.value[1];
    data_visible.bokeh_colors = data_visible.rade.map(v => valueToColor(v, low, high));
} else {
    data_visible.bokeh_colors = data_visible['bokeh_colors_' + key];
}

source_visible.change.emit();

color_bar.color_mapper.low = low;
color_bar.color_mapper.high = high;
color_bar.title = selected_label;