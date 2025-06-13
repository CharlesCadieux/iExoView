"use strict";

const selected_x_axis = cb_obj.value;
const data_visible = source_visible.data;
const inverted_axis = inverted.active;

data_visible.x = data_visible[axis_map[selected_x_axis]] || [];
source_visible.change.emit();
xaxis.axis_label = selected_x_axis;
xaxis.change.emit();

p.reset.emit();
inverted.active = inverted_axis;
