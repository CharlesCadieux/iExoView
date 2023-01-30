var data_visible = source_visible.data
var len = data_visible.name.length
const inverted_axis = inverted.active

search.value = ""
data_visible.edge_alpha = []
data_visible.fill_alpha = []
for (var i = 0; i < len; i++) {
data_visible.edge_alpha[i] = 1
data_visible.fill_alpha[i] = 0.8
}
source_visible.change.emit()

p.reset.emit()
inverted.active = []
inverted.active = inverted_axis
