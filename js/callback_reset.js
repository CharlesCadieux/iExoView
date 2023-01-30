var data_available = source_available.data
var data_visible = source_visible.data


data_visible.x = data_available.x
data_visible.y = data_available.y
data_visible.hostname = data_available.hostname
data_visible.name = data_available.name
data_visible.targettype = data_available.targettype
data_visible.ra_str = data_available.ra_str
data_visible.dec_str = data_available.dec_str
data_visible.stmass = data_available.stmass
data_visible.stradius = data_available.stradius
data_visible.stlum = data_available.stlum
data_visible.teff = data_available.teff
data_visible.mass = data_available.mass
data_visible.spt = data_available.spt
data_visible.dist = data_available.dist
data_visible.rade = data_available.rade
data_visible.rho = data_available.rho
data_visible.period = data_available.period
data_visible.depth = data_available.depth
data_visible.dur = data_available.dur
data_visible.semima = data_available.semima
data_visible.teq = data_available.teq
data_visible.insol = data_available.insol
data_visible.atmosignal = data_available.atmosignal
data_visible.tsm = data_available.tsm
data_visible.esm = data_available.esm
data_visible.K = data_available.K
data_visible.vmag = data_available.vmag
data_visible.gmag = data_available.gmag
data_visible.jmag = data_available.jmag
data_visible.hmag = data_available.hmag
data_visible.kmag = data_available.kmag
data_visible.size = data_available.size
data_visible.bokeh_colors = data_available.bokeh_colors
data_visible.edge_color = data_available.edge_color
data_visible.edge_thick = data_available.edge_thick
data_visible.edge_alpha = data_available.edge_alpha
data_visible.fill_alpha = data_available.fill_alpha
data_visible.semima_colors = data_available.semima_colors
data_visible.teq_colors = data_available.teq_colors
data_visible.rho_colors = data_available.rho_colors
data_visible.insol_colors = data_available.insol_colors
data_visible.spt_colors = data_available.spt_colors
data_visible.K_colors = data_available.K_colors
data_visible.mass_colors = data_available.mass_colors
data_visible.red_colors = data_available.red_colors
data_visible.lum_colors = data_available.lum_colors
data_visible.depth_colors = data_available.depth_colors

x_axis.value = "Period (days)"
y_axis.value = "Radius (R⊕)"
search.value = ""
ra_min_slider.value = 0
ra_max_slider.value = 24
dec_slider.value = [-90,90]
teff_slider.value = [2500,15000]
rade_slider.value = [0.1,40]
checkbox.active = [0,1,2]


source_visible.change.emit()
p.reset.emit()
title.text = titlestr

inverted.active = []
