var selected_index = cb_data.source.selected.indices
var data_visible = source_visible.data
var tabular = TOOLTIPS

if (selected_index.length == 1){
   tabular = tabular.replace('@hostname', data_visible.hostname[selected_index])
   tabular = tabular.replace('@name', data_visible.name[selected_index])
   tabular = tabular.replace('@ra_str', data_visible.ra_str[selected_index])
   tabular = tabular.replace('@mass_colors', data_visible.mass_colors[selected_index])
   if (data_visible.mass_colors[selected_index] == '#ff0000'){
       tabular = tabular.replace('@mass', data_visible.mass[selected_index].toFixed(2))
   }
   else tabular = tabular.replace('@mass', data_visible.mass[selected_index])
   tabular = tabular.replace('@dec_str', data_visible.dec_str[selected_index])
   tabular = tabular.replace('@rade', data_visible.rade[selected_index].toFixed(3))
   tabular = tabular.replace('@stmass', data_visible.stmass[selected_index].toFixed(3))
   tabular = tabular.replace('@rho_colors', data_visible.rho_colors[selected_index])
   if (data_visible.rho_colors[selected_index] == '#ff0000'){
       tabular = tabular.replace('@rho', data_visible.rho[selected_index].toFixed(2))
   }
   else tabular = tabular.replace('@rho', data_visible.rho[selected_index].toFixed(2))
   tabular = tabular.replace('@stradius', data_visible.stradius[selected_index].toFixed(3))
   tabular = tabular.replace('@lum_colors', data_visible.lum_colors[selected_index])
   if (data_visible.lum_colors[selected_index] == '#ff0000'){
       tabular = tabular.replace('@stlum', data_visible.stlum[selected_index].toFixed(3))
   }
   else tabular = tabular.replace('@stlum', data_visible.stlum[selected_index].toFixed(3))
   tabular = tabular.replace('@depth_colors', data_visible.depth_colors[selected_index])
   if (data_visible.depth_colors[selected_index] == '#ff0000'){
       tabular = tabular.replace('@depth', data_visible.depth[selected_index].toFixed(1))
   }
   else tabular = tabular.replace('@depth', data_visible.depth[selected_index].toFixed(1))
   tabular = tabular.replace('@period', data_visible.period[selected_index].toFixed(5))
   tabular = tabular.replace('@teff', data_visible.teff[selected_index].toFixed(0))
   tabular = tabular.replace('@dur', data_visible.dur[selected_index].toFixed(1))
   tabular = tabular.replace('@spt_colors', data_visible.spt_colors[selected_index])
   tabular = tabular.replace('@spt', data_visible.spt[selected_index])
   tabular = tabular.replace('@semima_colors', data_visible.semima_colors[selected_index])
   if (data_visible.semima_colors[selected_index] == '#ff0000'){
       tabular = tabular.replace('@semima', data_visible.semima[selected_index].toFixed(3))
   }
   else tabular = tabular.replace('@semima', data_visible.semima[selected_index])
   tabular = tabular.replace('@dist', data_visible.dist[selected_index])
   tabular = tabular.replace('@teq_colors', data_visible.teq_colors[selected_index])
   if (data_visible.teq_colors[selected_index] == '#ff0000'){
       tabular = tabular.replace('@teq', data_visible.teq[selected_index].toFixed(0))
   }
   else tabular = tabular.replace('@teq', data_visible.teq[selected_index])
   tabular = tabular.replace('@vmag', data_visible.vmag[selected_index])
   tabular = tabular.replace('@gmag', data_visible.gmag[selected_index])
   tabular = tabular.replace('@insol_colors', data_visible.insol_colors[selected_index])
   if (data_visible.insol_colors[selected_index] == '#ff0000'){
       tabular = tabular.replace('@insol', data_visible.insol[selected_index].toFixed(1))
   }
   else tabular = tabular.replace('@insol', data_visible.insol[selected_index])
   tabular = tabular.replace('@jmag', data_visible.jmag[selected_index])
   tabular = tabular.replace('@K_colors', data_visible.K_colors[selected_index])
   if (data_visible.K_colors[selected_index] == '#ff0000'){
       tabular = tabular.replace('@K', data_visible.K[selected_index].toFixed(2))
   }
   else tabular = tabular.replace('@K', data_visible.K[selected_index])
   tabular = tabular.replace('@hmag', data_visible.hmag[selected_index])
   tabular = tabular.replace('@red_colors', data_visible.red_colors[selected_index])
   tabular = tabular.replace('@red_colors', data_visible.red_colors[selected_index])
   tabular = tabular.replace('@red_colors', data_visible.red_colors[selected_index])
   tabular = tabular.replace('@atmosignal', data_visible.atmosignal[selected_index].toFixed(0))
   tabular = tabular.replace('@kmag', data_visible.kmag[selected_index])
   tabular = tabular.replace('@tsm', data_visible.tsm[selected_index].toFixed(1))
   tabular = tabular.replace('@esm', data_visible.esm[selected_index].toFixed(1))

   var x = window.open("", 'windowName','height=375,width=500')
   x.document.open();
   x.document.write(tabular);
   x.document.close();

}

else {
   ;
}
