"use strict";

const text = cb_obj.value;
const data_visible = source_visible.data;
const len = data_visible.name.length;
const name = data_visible.name;

let index = [];
for (let i = 0; i < len; i++) {
    let str = name[i].toLowerCase().trim().replace(/\s/g, '').replace(/-/g, '');
    let searchText = text.toLowerCase().trim().replace(/\s/g, '').replace(/-/g, '');
    if (str.includes(searchText)) {
        index.push(i);
    }
}

if (index.length > 0) {
    data_visible.edge_alpha = new Array(len).fill(0.10);
    data_visible.fill_alpha = new Array(len).fill(0.10);
    for (let i = 0; i < index.length; i++) {
        data_visible.edge_alpha[index[i]] = 1;
        data_visible.fill_alpha[index[i]] = 1;
    }
} else {
    data_visible.edge_alpha = new Array(len).fill(1);
    data_visible.fill_alpha = new Array(len).fill(0.8);
}

source_visible.change.emit();