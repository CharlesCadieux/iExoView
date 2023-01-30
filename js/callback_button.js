var data_available = source_available.data
var data_visible = source_visible.data
var ra_min = ra_min_slider.value
var ra_max = ra_max_slider.value
var dec_values = dec_slider.value
var teff_values = teff_slider.value
var rade_values = rade_slider.value
var len_available = data_available.name.length
var selected_x_axis = x_axis.value
var selected_y_axis = y_axis.value
const inverted_axis = inverted.active


if (checkbox.active == 0){
targettype = "Exoplanet"
}
else if (checkbox.active == 0){
targettype = "Exoplanet"
}

else if (checkbox.active == 1){
targettype = "TOI"
}

else if (checkbox.active == 2){
targettype = "K2"
}

else if (checkbox.active == "0,1"){
targettype = ["Exoplanet", "TOI"]
}

else if (checkbox.active == "0,2"){
targettype = ["Exoplanet", "K2"]
}

else if (checkbox.active == "1,2"){
targettype = ["TOI", "K2"]
}

else if (checkbox.active == "0,1,2"){
targettype = ["Exoplanet", "TOI","K2"]
}

else if (checkbox.active.length == 0){
targettype = ""
}

data_visible.x = []
data_visible.y = []
data_visible.hostname = []
data_visible.name = []
data_visible.targettype = []
data_visible.ra_str = []
data_visible.dec_str = []
data_visible.stmass = []
data_visible.stradius = []
data_visible.stlum = []
data_visible.teff = []
data_visible.spt = []
data_visible.dist = []
data_visible.mass = []
data_visible.rade = []
data_visible.rho = []
data_visible.period = []
data_visible.depth = []
data_visible.dur = []
data_visible.semima = []
data_visible.teq = []
data_visible.insol = []
data_visible.atmosignal = []
data_visible.tsm = []
data_visible.esm = []
data_visible.K = []
data_visible.vmag = []
data_visible.gmag = []
data_visible.jmag = []
data_visible.hmag = []
data_visible.kmag = []
data_visible.size = []
data_visible.bokeh_colors = []
data_visible.edge_color = []
data_visible.edge_thick = []
data_visible.edge_alpha = []
data_visible.fill_alpha = []
data_visible.semima_colors = []
data_visible.teq_colors = []
data_visible.rho_colors = []
data_visible.insol_colors = []
data_visible.spt_colors = []
data_visible.K_colors = []
data_visible.mass_colors = []
data_visible.red_colors = []
data_visible.lum_colors = []
data_visible.depth_colors = []

index = []

if (ra_min <= ra_max){

for (var i = 0; i < len_available; i++) {
   if  (data_available.ra[i] >= ra_min * 15 && data_available.ra[i] <= ra_max * 15
   && data_available.dec[i] >= dec_values[0] && data_available.dec[i] <= dec_values[1]
   && data_available.teff[i] >= teff_values[0] && data_available.teff[i] <= teff_values[1]
   && data_available.rade[i] >= rade_values[0] && data_available.rade[i] <= rade_values[1]
   && targettype.includes(data_available.targettype[i])){
       index.push(i)
   }
}

data_visible.x = index.map((item) => data_available.x[item])
data_visible.y = index.map((item) => data_available.y[item])
data_visible.hostname = index.map((item) => data_available.hostname[item])
data_visible.name = index.map((item) => data_available.name[item])
data_visible.targettype = index.map((item) => data_available.targettype[item])
data_visible.ra_str = index.map((item) => data_available.ra_str[item])
data_visible.dec_str = index.map((item) => data_available.dec_str[item])
data_visible.stmass = index.map((item) => data_available.stmass[item])
data_visible.stradius = index.map((item) => data_available.stradius[item])
data_visible.stlum = index.map((item) => data_available.stlum[item])
data_visible.teff = index.map((item) => data_available.teff[item])
data_visible.spt = index.map((item) => data_available.spt[item])
data_visible.dist = index.map((item) => data_available.dist[item])
data_visible.mass = index.map((item) => data_available.mass[item])
data_visible.rade = index.map((item) => data_available.rade[item])
data_visible.rho = index.map((item) => data_available.rho[item])
data_visible.period = index.map((item) => data_available.period[item])
data_visible.depth = index.map((item) => data_available.depth[item])
data_visible.dur = index.map((item) => data_available.dur[item])
data_visible.semima = index.map((item) => data_available.semima[item])
data_visible.teq = index.map((item) => data_available.teq[item])
data_visible.insol = index.map((item) => data_available.insol[item])
data_visible.atmosignal = index.map((item) => data_available.atmosignal[item])
data_visible.tsm = index.map((item) => data_available.tsm[item])
data_visible.esm = index.map((item) => data_available.esm[item])
data_visible.K = index.map((item) => data_available.K[item])
data_visible.vmag = index.map((item) => data_available.vmag[item])
data_visible.gmag = index.map((item) => data_available.gmag[item])
data_visible.jmag = index.map((item) => data_available.jmag[item])
data_visible.hmag = index.map((item) => data_available.hmag[item])
data_visible.kmag = index.map((item) => data_available.kmag[item])
data_visible.size = index.map((item) => data_available.size[item])
data_visible.bokeh_colors = index.map((item) => data_available.bokeh_colors[item])
data_visible.edge_color = index.map((item) => data_available.edge_color[item])
data_visible.edge_thick = index.map((item) => data_available.edge_thick[item])
data_visible.edge_alpha = index.map((item) => data_available.edge_alpha[item])
data_visible.fill_alpha = index.map((item) => data_available.fill_alpha[item])
data_visible.semima_colors = index.map((item) => data_available.semima_colors[item])
data_visible.teq_colors = index.map((item) => data_available.teq_colors[item])
data_visible.rho_colors = index.map((item) => data_available.rho_colors[item])
data_visible.insol_colors = index.map((item) => data_available.insol_colors[item])
data_visible.spt_colors = index.map((item) => data_available.spt_colors[item])
data_visible.K_colors = index.map((item) => data_available.K_colors[item])
data_visible.mass_colors = index.map((item) => data_available.mass_colors[item])
data_visible.red_colors = index.map((item) => data_available.red_colors[item])
data_visible.lum_colors = index.map((item) => data_available.lum_colors[item])
data_visible.depth_colors = index.map((item) => data_available.depth_colors[item])
}

if (ra_min > ra_max){

for (var i = 0; i < len_available; i++) {
   if  (data_available.dec[i] >= dec_values[0] && data_available.dec[i] <= dec_values[1]
   && data_available.teff[i] >= teff_values[0] && data_available.teff[i] <= teff_values[1]
   && data_available.rade[i] >= rade_values[0] && data_available.rade[i] <= rade_values[1]
   && targettype.includes(data_available.targettype[i])
   && (data_available.ra[i] >= ra_min * 15 || data_available.ra[i] <= ra_max * 15)){
       index.push(i)
   }
}

data_visible.x = index.map((item) => data_available.x[item])
data_visible.y = index.map((item) => data_available.y[item])
data_visible.hostname = index.map((item) => data_available.hostname[item])
data_visible.name = index.map((item) => data_available.name[item])
data_visible.targettype = index.map((item) => data_available.targettype[item])
data_visible.ra_str = index.map((item) => data_available.ra_str[item])
data_visible.dec_str = index.map((item) => data_available.dec_str[item])
data_visible.stmass = index.map((item) => data_available.stmass[item])
data_visible.stradius = index.map((item) => data_available.stradius[item])
data_visible.stlum = index.map((item) => data_available.stlum[item])
data_visible.teff = index.map((item) => data_available.teff[item])
data_visible.spt = index.map((item) => data_available.spt[item])
data_visible.dist = index.map((item) => data_available.dist[item])
data_visible.mass = index.map((item) => data_available.mass[item])
data_visible.rade = index.map((item) => data_available.rade[item])
data_visible.rho = index.map((item) => data_available.rho[item])
data_visible.period = index.map((item) => data_available.period[item])
data_visible.depth = index.map((item) => data_available.depth[item])
data_visible.dur = index.map((item) => data_available.dur[item])
data_visible.semima = index.map((item) => data_available.semima[item])
data_visible.teq = index.map((item) => data_available.teq[item])
data_visible.insol = index.map((item) => data_available.insol[item])
data_visible.atmosignal = index.map((item) => data_available.atmosignal[item])
data_visible.tsm = index.map((item) => data_available.tsm[item])
data_visible.esm = index.map((item) => data_available.esm[item])
data_visible.K = index.map((item) => data_available.K[item])
data_visible.vmag = index.map((item) => data_available.vmag[item])
data_visible.gmag = index.map((item) => data_available.gmag[item])
data_visible.jmag = index.map((item) => data_available.jmag[item])
data_visible.hmag = index.map((item) => data_available.hmag[item])
data_visible.kmag = index.map((item) => data_available.kmag[item])
data_visible.size = index.map((item) => data_available.size[item])
data_visible.bokeh_colors = index.map((item) => data_available.bokeh_colors[item])
data_visible.edge_color = index.map((item) => data_available.edge_color[item])
data_visible.edge_thick = index.map((item) => data_available.edge_thick[item])
data_visible.edge_alpha = index.map((item) => data_available.edge_alpha[item])
data_visible.fill_alpha = index.map((item) => data_available.fill_alpha[item])
data_visible.semima_colors = index.map((item) => data_available.semima_colors[item])
data_visible.teq_colors = index.map((item) => data_available.teq_colors[item])
data_visible.rho_colors = index.map((item) => data_available.rho_colors[item])
data_visible.insol_colors = index.map((item) => data_available.insol_colors[item])
data_visible.spt_colors = index.map((item) => data_available.spt_colors[item])
data_visible.K_colors = index.map((item) => data_available.K_colors[item])
data_visible.mass_colors = index.map((item) => data_available.mass_colors[item])
data_visible.red_colors = index.map((item) => data_available.red_colors[item])
data_visible.lum_colors = index.map((item) => data_available.lum_colors[item])
data_visible.depth_colors = index.map((item) => data_available.depth_colors[item])
}
data_visible.x = data_visible[axis_map[selected_x_axis]]
data_visible.y = data_visible[axis_map[selected_y_axis]]


count_exo = 0
for(var i = 0; i < data_visible.targettype.length; ++i){
if(data_visible.targettype[i] == "Exoplanet")
   count_exo++;
}

count_toi = 0
for(var i = 0; i < data_visible.targettype.length; ++i){
if(data_visible.targettype[i] == "TOI")
   count_toi++;
}

count_k2 = 0
for(var i = 0; i < data_visible.targettype.length; ++i){
if(data_visible.targettype[i] == "K2")
   count_k2++;
}

templateString = "{name1} Transiting Exoplanets with Known Mass + {name2} TOI + {name3} K2 Planets (as of {name4})";

templateString = templateString.replace('{name1}', count_exo);
templateString = templateString.replace('{name2}', count_toi);
templateString = templateString.replace('{name3}', count_k2);
templateString = templateString.replace('{name4}', date);

title.text = templateString
source_visible.change.emit()
search.value = ""
p.reset.emit()
inverted.active = []
inverted.active = inverted_axis
