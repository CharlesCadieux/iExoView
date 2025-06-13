"use strict";

const xi = xr.start;
const xf = xr.end;
const yi = yr.start;
const yf = yr.end;

let previous = [];

if (xi < xf && yi < yf) {
    previous = [];
} else if (xi > xf && yi < yf) {
    previous = [0];
} else if (xi < xf && yi > yf) {
    previous = [1];
} else if (xi > xf && yi > yf) {
    previous = [0, 1];
}

if (previous.length === 0 && cb_obj.active.length === 0) {
    // No action needed
} else if (previous.length === 0 && cb_obj.active.includes(0)) {
    xr.start = xf;
    xr.end = xi;
} else if (previous.length === 0 && cb_obj.active.includes(1)) {
    yr.start = yf;
    yr.end = yi;
} else if (previous.length === 0 && cb_obj.active.includes(0) && cb_obj.active.includes(1)) {
    xr.start = xf;
    xr.end = xi;
    yr.start = yf;
    yr.end = yi;
} else if (previous.includes(0) && cb_obj.active.includes(0) && cb_obj.active.includes(1)) {
    yr.start = yf;
    yr.end = yi;
} else if (previous.includes(0) && cb_obj.active.length === 0) {
    xr.start = xf;
    xr.end = xi;
} else if (previous.includes(1) && cb_obj.active.includes(0) && cb_obj.active.includes(1)) {
    xr.start = xf;
    xr.end = xi;
} else if (previous.includes(1) && cb_obj.active.length === 0) {
    yr.start = yf;
    yr.end = yi;
} else if (previous.includes(0) && previous.includes(1) && cb_obj.active.includes(0)) {
    yr.start = yf;
    yr.end = yi;
} else if (previous.includes(0) && previous.includes(1) && cb_obj.active.includes(1)) {
    xr.start = xf;
    xr.end = xi;
}