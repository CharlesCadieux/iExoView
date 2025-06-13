"use strict";

function table_to_csv(source) {
    const columns = ['hostname', 'name', 'targettype', 'ra_str', 'dec_str', 'stmass', 'stradius', 'vmag', 'jmag', 'hmag', 'kmag', 'teff', 'spt', 'mass', 'rade', 'period', 'K', 'teq', 'tsm', 'esm'];
    const nrows = source.get_length();
    const lines = [columns.join(',')];

    for (let i = 0; i < nrows; i++) {
        let row = [];
        for (let j = 0; j < columns.length; j++) {
            const column = columns[j];
            const value = source.data[column][i];
            // Handle null/undefined values and ensure proper CSV escaping
            row.push(value == null ? '' : `"${value.toString().replace(/"/g, '""')}"`);
        }
        lines.push(row.join(','));
    }
    return lines.join('\n') + '\n';
}

const filename = 'table.csv';
const filetext = table_to_csv(source);
const blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' });

// Addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename);
} else {
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = filename;
    link.target = '_blank';
    link.style.visibility = 'hidden';
    link.dispatchEvent(new MouseEvent('click'));
    URL.revokeObjectURL(link.href); // Clean up the URL object
}