"use strict";
try {
    const data_available = source_available.data;
    const data_visible = source_visible.data;
    const ra_min = ra_min_slider.value;
    const ra_max = ra_max_slider.value;
    const dec_values = dec_slider.value;
    const teff_values = teff_slider.value;
    const rade_values = rade_slider.value;
    const len_available = data_available.name.length;
    const selected_x_axis = x_axis.value;
    const selected_y_axis = y_axis.value;
    const inverted_axis = inverted.active;

    let targettype = [];
    if (checkbox.active.includes(0)) targettype.push("Exoplanet");
    if (checkbox.active.includes(1)) targettype.push("TOI");
    if (checkbox.active.includes(2)) targettype.push("K2");
    if (checkbox.active.length === 0) targettype = "";

    const properties = [
        'x', 'y', 'hostname', 'name', 'targettype', 'ra_str', 'dec_str', 'stmass', 
        'stradius', 'stlum', 'teff', 'spt', 'dist', 'mass', 'rade', 'rho', 'period', 
        'depth', 'dur', 'semima', 'teq', 'insol', 'atmosignal', 'tsm', 'esm', 'K', 
        'vmag', 'gmag', 'jmag', 'hmag', 'kmag', 'size', 'bokeh_colors', 'edge_color', 
        'edge_thick', 'edge_alpha', 'fill_alpha', 'semima_colors', 'teq_colors', 
        'rho_colors', 'insol_colors', 'spt_colors', 'K_colors', 'mass_colors', 
        'red_colors', 'lum_colors', 'depth_colors'
    ];

    let index = [];
    console.log("Filtering data...");
    if (ra_min <= ra_max) {
        for (let i = 0; i < len_available; i++) {
            const teff_val = parseFloat(data_available.teff[i]);
            const rade_val = parseFloat(data_available.rade[i]);
            if (data_available.ra[i] >= ra_min * 15 && data_available.ra[i] <= ra_max * 15
                && data_available.dec[i] >= dec_values[0] && data_available.dec[i] <= dec_values[1]
                && (!isNaN(teff_val) && teff_val >= teff_values[0] && teff_val <= teff_values[1])
                && (!isNaN(rade_val) && rade_val >= rade_values[0] && rade_val <= rade_values[1])
                && (typeof targettype === "string" || targettype.includes(data_available.targettype[i]))) {
                index.push(i);
            }
        }
    } else {
        for (let i = 0; i < len_available; i++) {
            const teff_val = parseFloat(data_available.teff[i]);
            const rade_val = parseFloat(data_available.rade[i]);
            if (data_available.dec[i] >= dec_values[0] && data_available.dec[i] <= dec_values[1]
                && (!isNaN(teff_val) && teff_val >= teff_values[0] && teff_val <= teff_values[1])
                && (!isNaN(rade_val) && rade_val >= rade_values[0] && rade_val <= rade_values[1])
                && (typeof targettype === "string" || targettype.includes(data_available.targettype[i]))
                && (data_available.ra[i] >= ra_min * 15 || data_available.ra[i] <= ra_max * 15)) {
                index.push(i);
            }
        }
    }

    Object.keys(data_available).forEach(prop => {
        data_visible[prop] = index.map(i => data_available[prop][i]);
    });

    // Compute x/y based on selected axes
    data_visible.x = data_visible[axis_map[selected_x_axis]];
    data_visible.y = data_visible[axis_map[selected_y_axis]];

    // Update data
    source_visible.data = data_visible;


    let count_exo = 0, count_toi = 0, count_k2 = 0;
    for (let i = 0; i < data_visible.targettype.length; i++) {
        if (data_visible.targettype[i] === "Exoplanet") count_exo++;
        if (data_visible.targettype[i] === "TOI") count_toi++;
        if (data_visible.targettype[i] === "K2") count_k2++;
    }

    const date = new Date().toLocaleDateString();
    let templateString = "{name1} Transiting Exoplanets with Known Mass + {name2} TOI + {name3} K2 Planets (as of {name4})";
    templateString = templateString.replace('{name1}', count_exo);
    templateString = templateString.replace('{name2}', count_toi);
    templateString = templateString.replace('{name3}', count_k2);
    templateString = templateString.replace('{name4}', date);

    title.text = templateString;
    source_visible.change.emit();
    search.value = "";
    p.reset.emit();
    inverted.active = inverted_axis;
} catch (e) {
    console.log("Error in callback_button.js:", e);
}