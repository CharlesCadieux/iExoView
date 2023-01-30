function table_to_csv(source) {
    const columns = ['hostname', 'name', 'targettype', 'ra_str', 'dec_str', 'stmass', 'stradius', 'vmag', 'jmag', 'hmag', 'kmag', 'teff', 'spt', 'mass', 'rade', 'period', 'K', 'teq', 'tsm', 'esm']
    const nrows = source.get_length()
    const lines = [columns.join(',')]

    for (let i = 0; i < nrows; i++) {
    let row = [];
    for (let j = 0; j < columns.length; j++) {
        const column = columns[j]
        row.push(source.data[column][i].toString())
    }
    lines.push(row.join(','))
    }
    return lines.join('\n').concat('\n')
}


const filename = 'table.csv'
filetext = table_to_csv(source)
const blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' })

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename)
} else {
    const link = document.createElement('a')
    link.href = URL.createObjectURL(blob)
    link.download = filename
    link.target = '_blank'
    link.style.visibility = 'hidden'
    link.dispatchEvent(new MouseEvent('click'))
}
