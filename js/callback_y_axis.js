var selected_y_axis = cb_obj.value
var data_visible = source_visible.data
const inverted_axis = inverted.active

data_visible.y = data_visible[axis_map[selected_y_axis]]
source_visible.change.emit()
yaxis.attributes.axis_label = selected_y_axis
yaxis.change.emit()

p.reset.emit()
inverted.active = []
inverted.active = inverted_axis
