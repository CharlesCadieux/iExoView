"use strict";

var data_visible = source_visible.data;
var len = data_visible.name.length;

search.value = "";
data_visible.edge_alpha = new Array(len).fill(1);
data_visible.fill_alpha = new Array(len).fill(0.8);

source_visible.change.emit();