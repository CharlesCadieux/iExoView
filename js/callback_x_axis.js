"use strict";

const selected_x_axis = cb_obj.value;
const data_visible = source_visible.data;
const inverted_axis = inverted.active;

data_visible.x = data_visible[axis_map[selected_x_axis]] || [];
source_visible.change.emit();
xaxis.axis_label = selected_x_axis;
xaxis.change.emit();

p.reset.emit();

// p.reset.emit() restores the default (non-inverted) range direction,
// so re-apply any active axis inversion manually.
if (inverted_axis.includes(0)) {
    const tmp = p.x_range.start;
    p.x_range.start = p.x_range.end;
    p.x_range.end = tmp;
}
if (inverted_axis.includes(1)) {
    const tmp = p.y_range.start;
    p.y_range.start = p.y_range.end;
    p.y_range.end = tmp;
}
inverted.active = inverted_axis;