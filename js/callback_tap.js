"use strict";

const selected_index = cb_data.source.selected.indices;
const data_visible = source_visible.data;
const safeFixed = (val, digits) =>
  Number.isFinite(parseFloat(val)) ? parseFloat(val).toFixed(digits) : '';
let tabular = TOOLTIPS;

if (selected_index.length === 1) {
    const index = selected_index[0];
    tabular = tabular.replace('@hostname', data_visible.hostname[index] ?? '');
    tabular = tabular.replace('@name', data_visible.name[index] ?? '');
    tabular = tabular.replace('@ra_str', data_visible.ra_str[index] ?? '');
    tabular = tabular.replace('@mass_colors', data_visible.mass_colors[index] ?? '');
    tabular = tabular.replace('@mass', data_visible.mass[index] != null ? 
        (data_visible.mass_colors[index] === '#ff0000' ? 
         safeFixed(data_visible.mass[index], 2) : 
         data_visible.mass[index]) : '');
    tabular = tabular.replace('@dec_str', data_visible.dec_str[index] ?? '');
    tabular = tabular.replace('@rade', data_visible.rade[index] != null ? 
        safeFixed(data_visible.rade[index], 3) : '');
    tabular = tabular.replace('@stmass', data_visible.stmass[index] != null ? 
        safeFixed(data_visible.stmass[index], 3) : '');
    tabular = tabular.replace('@rho_colors', data_visible.rho_colors[index] ?? '');
    tabular = tabular.replace('@rho', data_visible.rho[index] != null ? 
        safeFixed(data_visible.rho[index], 2) : '');
    tabular = tabular.replace('@stradius', data_visible.stradius[index] != null ? 
        safeFixed(data_visible.stradius[index], 3) : '');
    tabular = tabular.replace('@lum_colors', data_visible.lum_colors[index] ?? '');
    tabular = tabular.replace('@stlum', data_visible.stlum[index] != null ? 
        safeFixed(data_visible.stlum[index], 3) : '');
    tabular = tabular.replace('@depth_colors', data_visible.depth_colors[index] ?? '');
    tabular = tabular.replace('@depth', data_visible.depth[index] != null ? 
        safeFixed(data_visible.depth[index], 1) : '');
    tabular = tabular.replace('@period', data_visible.period[index] != null ? 
        safeFixed(data_visible.period[index], 5) : '');
    tabular = tabular.replace('@teff', data_visible.teff[index] != null ? 
        safeFixed(data_visible.teff[index], 0) : '');
    tabular = tabular.replace('@dur', data_visible.dur[index] != null ? 
        safeFixed(data_visible.dur[index], 1) : '');
    tabular = tabular.replace('@spt_colors', data_visible.spt_colors[index] ?? '');
    tabular = tabular.replace('@spt', data_visible.spt[index] ?? '');
    tabular = tabular.replace('@semima_colors', data_visible.semima_colors[index] ?? '');
    tabular = tabular.replace('@semima', data_visible.semima[index] != null ? 
        (data_visible.semima_colors[index] === '#ff0000' ? 
         safeFixed(data_visible.semima[index], 3) : 
         data_visible.semima[index]) : '');
    tabular = tabular.replace('@dist', data_visible.dist[index] ?? '');
    tabular = tabular.replace('@teq_colors', data_visible.teq_colors[index] ?? '');
    tabular = tabular.replace('@teq', data_visible.teq[index] != null ? 
        (data_visible.teq_colors[index] === '#ff0000' ? 
         safeFixed(data_visible.teq[index], 0) : 
         data_visible.teq[index]) : '');
    tabular = tabular.replace('@vmag', data_visible.vmag[index] ?? '');
    tabular = tabular.replace('@gmag', data_visible.gmag[index] ?? '');
    tabular = tabular.replace('@insol_colors', data_visible.insol_colors[index] ?? '');
    tabular = tabular.replace('@insol', data_visible.insol[index] != null ? 
        (data_visible.insol_colors[index] === '#ff0000' ? 
         safeFixed(data_visible.insol[index], 1) : 
         data_visible.insol[index]) : '');
    tabular = tabular.replace('@jmag', data_visible.jmag[index] ?? '');
    tabular = tabular.replace('@K_colors', data_visible.K_colors[index] ?? '');
    tabular = tabular.replace('@K', data_visible.K[index] != null ? 
        (data_visible.K_colors[index] === '#ff0000' ? 
         safeFixed(data_visible.K[index], 2) : 
         data_visible.K[index]) : '');
    tabular = tabular.replace('@hmag', data_visible.hmag[index] ?? '');
    tabular = tabular.replace(/@red_colors/g, data_visible.red_colors[index] ?? '');
    tabular = tabular.replace('@atmosignal', data_visible.atmosignal[index] != null ? 
        safeFixed(data_visible.atmosignal[index], 0) : '');
    tabular = tabular.replace('@kmag', data_visible.kmag[index] ?? '');
    tabular = tabular.replace('@tsm', data_visible.tsm[index] != null ? 
        safeFixed(data_visible.tsm[index], 1) : '');
    tabular = tabular.replace('@esm', data_visible.esm[index] != null ? 
        safeFixed(data_visible.esm[index], 1) : '');

    const x = window.open("", 'windowName', 'height=375,width=500');
    x.document.open();
    x.document.write(tabular);
    x.document.close();
}