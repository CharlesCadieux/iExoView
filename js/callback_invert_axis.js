const xi = xr.start
const xf = xr.end
const yi = yr.start
const yf = yr.end

if (xi < xf && yi < yf){
previous = []

}

else if (xi > xf && yi < yf){
previous = 0
}

else if (xi < xf && yi > yf){
previous = 1
}

else if (xi > xf && yi > yf){
previous = [0,1]
}


if (previous.length == 0 && cb_obj.active.length == 0){


}

else if (previous.length == 0 && cb_obj.active == 0){
xr.start = xf;
xr.end = xi;

}

else if (previous.length == 0 && cb_obj.active == 1){
yr.start = yf;
yr.end = yi;

}

else if (previous.length == 0 && cb_obj.active == "0,1"){
xr.start = xf;
xr.end = xi;
yr.start = yf;
yr.end = yi;

}

else if (previous == 0 && cb_obj.active == "0,1"){
yr.start = yf;
yr.end = yi;

}

else if (previous == 0 && cb_obj.active.length == 0){
xr.start = xf;
xr.end = xi;

}

else if (previous == 1 && cb_obj.active == "0,1"){
xr.start = xf;
xr.end = xi;

}

else if (previous == 1 && cb_obj.active.length == 0){
yr.start = yf;
yr.end = yi;

}

else if (previous == "0,1" && cb_obj.active == 0){
yr.start = yf;
yr.end = yi;

}

else if (previous == "0,1" && cb_obj.active == 1){
xr.start = xf;
xr.end = xi;

}
