var text = cb_obj.value
var data_visible = source_visible.data
var len = data_visible.name.length
var name = data_visible.name
const inverted_axis = inverted.active

index = []
for (var i = 0; i < len; i++) {
str = name[i].toLowerCase()
str = str.trim()
str = str.replace(/\s/g, '')
str = str.replace('-', '')
text = text.toLowerCase()
text = text.trim()
text = text.replace(/\s/g, '')
text = text.replace('-', '')
if (str.search(text) > -1){
   index.push(i)
}
}


if (index.length > 0) {
data_visible.edge_alpha = []
data_visible.fill_alpha = []
for (var i = 0; i < len; i++) {
   data_visible.edge_alpha[i] = 0.10
   data_visible.fill_alpha[i] = 0.10
}
for (var i = 0; i < index.length; i++) {
data_visible.edge_alpha[index[i]] = 1
data_visible.fill_alpha[index[i]] = 1
}
}

else {
data_visible.edge_alpha = []
data_visible.fill_alpha = []
for (var i = 0; i < len; i++) {
   data_visible.edge_alpha[i] = 1
   data_visible.fill_alpha[i] = 0.8
}

}


source_visible.change.emit()
p.reset.emit()
inverted.active = []
inverted.active = inverted_axis
