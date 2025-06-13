# %%
"""
iExoView, the Interactive Exoplanet Viewer
astro.umontreal.ca/~charles/iExoView.html

@author: Charles Cadieux 2019-2022
@updated: 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import requests
from datetime import datetime
from astropy.table import Table, Column
from astropy import units as u
from astropy.coordinates import Angle
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import (Circle, CustomJS, Label, ColorBar, LogColorMapper,
                         LogAxis, LogTicker, Slider, RangeSlider, Select, Div,
                         TextInput, Button, CheckboxGroup, CheckboxButtonGroup,
                         Scatter)
from bokeh.models.tools import TapTool
from bokeh.layouts import layout, row, column
from bokeh import events
from exofile.archive import ExoFile
from utils import bulkrho, semiamp, mass_est, semima, Teq, TSM, ESM, atmosig, find_nearest

# %%
# Import data
tbl = ExoFile.load()
tbl['pl_trandep'] *= 10  # Transit depth: percent to ppt
tbl['st_lum'] = 10**tbl['st_lum']  # Stellar luminosity in solar units

# TOI data from ExoFOP
toi_response = requests.get(
    "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
)

# Parse TOI data into Table
tbl_toi = Table.read(toi_response.text, format="ascii.csv", data_start=1, delimiter=',')
tbl_toi['TFOPWG Disposition'].mask = False
valid_toi_idx = np.where(
    (tbl_toi['TFOPWG Disposition'] == '0') |
    (tbl_toi['TFOPWG Disposition'] == 'PC')
)[0]
tbl_toi = tbl_toi[valid_toi_idx]
number_toi = len(tbl_toi)

# Spectral type data (Pecaut & Mamajek 2013)
spectral_response = requests.get(
    "http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt"
)
tbl_spt = Table.read(
    spectral_response.text, format='ascii.csv', delimiter=' ', comment='&',
    header_start=22, data_start=23, data_end=141
)
del tbl_spt['#SpT_1']
tbl_spt.rename_column('#SpT', 'SpT')

# %%
# Populate TOIs
# DEBUG: Limit to 1 TOI for testing
#number_toi = 1
col_starts_stellar = (
    0, 25, 43, 65, 87, 109, 131, 153, 175, 197, 219, 241, 263, 285, 307,
    329, 351, 373, 395, 417, 439, 461, 483, 505, 530, 552, 574, 599, 621,
    643, 665, 687, 709, 731, 749, 767
)
col_starts_mag = (0, 18, 36, 54, 79, 99, 117, 135)

new_row_number = 0
for i in range(number_toi):
    tic_id = tbl_toi['TIC ID'][i]
    try:
        exofop_response = requests.get(
            f"https://exofop.ipac.caltech.edu/tess/download_target.php?id={tic_id}"
        )
        exofop_text = exofop_response.text

        # Parse stellar parameters
        idx1 = exofop_text.find("STELLAR PARMAMETERS")
        idx2 = exofop_text.find("MAGNITUDES")
        stellar_params = Table.read(
            exofop_text[idx1:idx2], format='ascii.fixed_width',
            col_starts=col_starts_stellar, header_start=1, data_start=2
        )

        # Parse magnitudes
        idx1 = exofop_text.find("MAGNITUDES")
        idx2 = exofop_text.find("IMAGING OBSERVATIONS")
        magnitudes = Table.read(
            exofop_text[idx1:idx2], format='ascii.fixed_width',
            col_starts=col_starts_mag, header_start=1, data_start=2
        )
        idx_v = np.where(magnitudes['Band'] == 'V')[0]
        idx_g = np.where(magnitudes['Band'] == 'Gaia')[0]
        idx_j = np.where(magnitudes['Band'] == 'J')[0]
        idx_h = np.where(magnitudes['Band'] == 'H')[0]
        idx_k = np.where(magnitudes['Band'] == 'K')[0]

        # Estimate planetary parameters
        masse = tbl_toi['Predicted Mass (M_Earth)'][i]
        density = bulkrho(masse, tbl_toi['Planet Radius (R_Earth)'][i])
        rv_amplitude = semiamp(
            tbl_toi['Period (days)'][i], masse,
            tbl_toi['Planet Radius (R_Earth)'][i],
            float(stellar_params['Mass (M_Sun)'][0]), 90, 0
        )

        # Add row to table
        new_row_number += 1
        print(f"Adding TOI {tbl_toi['TOI'][i]} ({new_row_number}/{number_toi})")
        tbl.add_row()
        tbl[-1]['hostname'] = f'TIC {tbl_toi["TIC ID"][i]}'
        tbl[-1]['rastr'] = str(tbl_toi['RA'][i])
        tbl[-1]['decstr'] = str(tbl_toi['Dec'][i])
        tbl[-1]['ra'] = Angle(tbl_toi['RA'][i], unit=u.hour).deg
        tbl[-1]['dec'] = Angle(tbl_toi['Dec'][i], unit=u.deg).deg
        tbl[-1]['st_rad'] = tbl_toi['Stellar Radius (R_Sun)'][i]
        tbl[-1]['st_mass'] = stellar_params['Mass (M_Sun)'][0]
        tbl[-1]['st_lum'] = stellar_params['Luminosity'][0]
        tbl[-1]['st_teff'] = tbl_toi['Stellar Eff Temp (K)'][i]
        tbl[-1]['st_spectype'] = '0'
        tbl[-1]['sy_dist'] = tbl_toi['Stellar Distance (pc)'][i]
        if len(idx_v) == 1:
            tbl[-1]['sy_vmag'] = magnitudes['Value'][idx_v]
        if len(idx_g) == 1:
            tbl[-1]['sy_gaiamag'] = magnitudes['Value'][idx_g]
        if len(idx_j) == 1:
            tbl[-1]['sy_jmag'] = magnitudes['Value'][idx_j]
        if len(idx_h) == 1:
            tbl[-1]['sy_hmag'] = magnitudes['Value'][idx_h]
        if len(idx_k) == 1:
            tbl[-1]['sy_kmag'] = magnitudes['Value'][idx_k]
        tbl[-1]['pl_name'] = f'TOI {tbl_toi["TOI"][i]}'
        tbl[-1]['pl_bmassprov'] = 'Mass Estimate'
        tbl[-1]['pl_bmasse'] = np.round(masse, 2)
        tbl[-1]['pl_rade'] = tbl_toi['Planet Radius (R_Earth)'][i]
        tbl[-1]['pl_dens'] = np.round(density, 2)
        tbl[-1]['pl_orbper'] = tbl_toi['Period (days)'][i]
        tbl[-1]['pl_trandur'] = tbl_toi['Duration (hours)'][i]
        tbl[-1]['pl_orbsmax'] = 0  # Estimated later
        tbl[-1]['pl_eqt'] = tbl_toi['Planet Equil Temp (K)'][i]
        tbl[-1]['pl_insol'] = tbl_toi['Planet Insolation (Earth Flux)'][i]
        tbl[-1]['pl_rvamp'] = np.round(rv_amplitude, 2)
        tbl[-1]['pl_bmasselim'] = 0
        tbl[-1]['tran_flag'] = 1
        tbl[-1]['pl_trandep'] = tbl_toi['Depth (ppm)'][i] / 1000  # in ppt
        tbl[-1]['disc_facility'] = 'Transiting Exoplanet Survey Satellite (TESS)'
    except Exception as e:
        print(f"Error processing TIC {tic_id}: {e}")
        continue

# Making sure we have the right number of TOIs added
number_toi = new_row_number

# Populate K2 planets with unknown mass
k2_idx, = np.where((tbl['disc_facility'] == 'K2') &
                    (tbl['tran_flag'] == 1))
mass_idx, = np.where((tbl['pl_bmassprov'] == 'Mass') &
                      (tbl['pl_bmasselim'] == 0))
k2_idx = np.setdiff1d(k2_idx, mass_idx)
number_k2 = len(k2_idx)

for i in k2_idx:
    masse = mass_est(tbl['pl_rade'][i])
    density = bulkrho(masse, tbl['pl_rade'][i])
    rv_amplitude = semiamp(
        tbl['pl_orbper'][i], masse, tbl['pl_rade'][i], tbl['st_mass'][i], 90, 0
    )
    tbl[i]['pl_bmassprov'] = 'Mass Estimate'
    tbl[i]['pl_bmasse'] = masse
    tbl[i]['pl_dens'] = density
    tbl[i]['pl_rvamp'] = rv_amplitude
    tbl[i]['pl_bmasselim'] = 0


# Clean table
valid_idx = np.where((tbl['tran_flag'] == 1) & (tbl['pl_rade'] > 0))[0]
tbl = tbl[valid_idx]
valid_idx = np.where(
    (tbl['pl_bmassprov'] == 'Mass') |
    (tbl['pl_bmassprov'] == 'Msini') |
    (tbl['pl_bmassprov'] == 'Msin(i)/sin(i)') |
    (tbl['pl_bmassprov'] == 'Mass Estimate')
)[0]
tbl = tbl[valid_idx]
valid_idx = np.where(tbl['pl_bmasselim'] == 0)[0]
tbl = tbl[valid_idx]
k2_idx = np.where(
    (tbl['disc_facility'] == 'K2') & (tbl['pl_bmassprov'] == 'Mass Estimate')
)[0]

# %%
# Fill missing parameters
c_dict = {}
c_dict['mass'] = np.array(len(tbl['pl_name']) * ['default'])
c_dict['mass'][k2_idx] = '#ff0000'
c_dict['mass'][-number_toi:] = '#ff0000'

# Semi-major axis
tbl['pl_orbsmax'].mask = False
sma = semima(tbl['pl_orbper'].data, tbl['st_mass'].data)
sma_idx = np.where(tbl['pl_orbsmax'] == 0)[0]
tbl['pl_orbsmax'][sma_idx] = sma[sma_idx]
c_dict['semima'] = np.array(len(tbl['pl_name']) * ['default'])
c_dict['semima'][sma_idx] = '#ff0000'
mask_idx = np.where(tbl['pl_orbsmax'].mask)[0]
tbl['pl_orbsmax'][mask_idx] = 0

# Equilibrium temperature
tbl['pl_eqt'].mask = False
teq = Teq(tbl['st_rad'].data, tbl['st_teff'].data, tbl['pl_orbsmax'].data)
teq_idx = np.where(tbl['pl_eqt'] == 0)[0]
tbl['pl_eqt'][teq_idx] = teq[teq_idx]
c_dict['teq'] = np.array(len(tbl['pl_name']) * ['default'])
c_dict['teq'][teq_idx] = '#ff0000'
mask_idx = np.where(tbl['pl_eqt'].mask)[0]
tbl['pl_eqt'][mask_idx] = 0

# Bulk density
tbl['pl_dens'].mask = False
dens_idx = np.where(tbl['pl_dens'].data == 0)[0]
density = bulkrho(tbl['pl_bmasse'].data, tbl['pl_rade'].data)
tbl['pl_dens'][dens_idx] = density[dens_idx]
c_dict['rho'] = np.array(len(tbl['pl_name']) * ['default'])
c_dict['rho'][dens_idx] = '#ff0000'
c_dict['rho'][k2_idx] = '#ff0000'
c_dict['rho'][-number_toi:] = '#ff0000'
mask_idx = np.where(tbl['pl_dens'].mask)[0]
tbl['pl_dens'][mask_idx] = 0

# Insolation
tbl['pl_insol'].mask = False
insol_idx = np.where(tbl['pl_insol'].data == 0)[0]
insol = (
    tbl['st_rad'].data**2 * (tbl['st_teff'].data / 5777)**4 /
    tbl['pl_orbsmax'].data**2
)
tbl['pl_insol'][insol_idx] = insol[insol_idx]
c_dict['insol'] = np.array(len(tbl['pl_name']) * ['default'])
c_dict['insol'][insol_idx] = '#ff0000'
mask_idx = np.where(tbl['pl_insol'].mask)[0]
tbl['pl_insol'][mask_idx] = 0

# RV semi-amplitude
tbl['pl_rvamp'].mask = False
rv_idx = np.where(tbl['pl_rvamp'].data == 0)[0]
rv_amplitude = semiamp(
    tbl['pl_orbper'].data, tbl['pl_bmasse'].data, tbl['pl_rade'].data,
    tbl['st_mass'].data, 90, 0
)
tbl['pl_rvamp'][rv_idx] = rv_amplitude[rv_idx]
c_dict['K'] = np.array(len(tbl['pl_name']) * ['default'])
c_dict['K'][rv_idx] = '#ff0000'
c_dict['K'][k2_idx] = '#ff0000'
c_dict['K'][-number_toi:] = '#ff0000'
mask_idx = np.where(tbl['pl_rvamp'].mask)[0]
tbl['pl_rvamp'][mask_idx] = 0

# Atmospheric signal
atmosig_col = Column(name='pl_atmosig', data=np.zeros(len(tbl['pl_name'])))
tbl.add_column(atmosig_col, 66)
mu_array = np.zeros(len(tbl['pl_name']))
for i, rade in enumerate(tbl['pl_rade']):
    mu_array[i] = (
        28.97 if rade <= 2 else 2.61 if 2 < rade <= 6 else 2.22
    )
tbl['pl_atmosig'] = atmosig(
    tbl['pl_eqt'].data, tbl['st_rad'].data, tbl['pl_dens'].data, mu_array
)
c_dict['red'] = np.array(len(tbl['pl_name']) * ['#ff0000'])
mask_idx = np.where(tbl['pl_atmosig'].mask)[0]
tbl['pl_atmosig'][mask_idx] = 0

# TSM
tsm_col = Column(name='pl_tsm', data=np.zeros(len(tbl['pl_name'])))
tbl.add_column(tsm_col, 67)
tsm = np.zeros(len(tbl['pl_name']))
for i in range(len(tsm)):
    tsm[i] = TSM(
        tbl['pl_rade'][i], tbl['pl_bmasse'][i], tbl['st_rad'][i],
        tbl['st_teff'][i], tbl['pl_orbsmax'][i], tbl['sy_jmag'][i]
    )
tsm[np.isnan(tsm) | np.isinf(tsm)] = 0
tbl['pl_tsm'] = tsm
mask_idx = np.where(tbl['pl_tsm'].mask)[0]
tbl['pl_tsm'][mask_idx] = 0

# ESM
esm_col = Column(name='pl_esm', data=np.zeros(len(tbl['pl_name'])))
tbl.add_column(esm_col, 68)
esm = np.zeros(len(tbl['pl_name']))
for i in range(len(esm)):
    esm[i] = ESM(
        tbl['pl_rade'][i], tbl['st_rad'][i], tbl['st_teff'][i],
        tbl['pl_orbsmax'][i], tbl['sy_kmag'][i]
    )
esm[np.isnan(esm) | np.isinf(esm)] = 0
tbl['pl_esm'] = esm
mask_idx = np.where(tbl['pl_esm'].mask)[0]
tbl['pl_esm'][mask_idx] = 0

# Spectral type
tbl['st_spectype'].mask = False
tbl['st_teff'].mask = False
spt_idx = np.where(tbl['st_spectype'].data == '0')[0]
c_dict['spt'] = np.array(len(tbl['pl_name']) * ['default'])
for i in spt_idx:
    if tbl['st_teff'][i] > 0:
        tbl['st_spectype'][i] = tbl_spt['SpT'][
            find_nearest(tbl_spt['Teff'], tbl['st_teff'][i])
        ]
        c_dict['spt'][i] = '#ff0000'
mask_idx = np.where(tbl['st_spectype'].mask)[0]
tbl['st_spectype'][mask_idx] = 'N/A'

# Transit depth
tbl['pl_trandep'].mask = False
depth_idx = np.where(tbl['pl_trandep'].data == 0)[0]
depth = (
    (6.371e6 / 6.956e8)**2 * (tbl['pl_rade'].data / tbl['st_rad'].data)**2 * 1000
)
tbl['pl_trandep'][depth_idx] = depth[depth_idx]
c_dict['depth'] = np.array(len(tbl['pl_name']) * ['default'])
c_dict['depth'][depth_idx] = '#ff0000'
mask_idx = np.where(tbl['pl_trandep'].mask)[0]
tbl['pl_trandep'][mask_idx] = 0

# Luminosity
tbl['st_lum'].mask = False
lum_idx = np.where(tbl['st_lum'].data == 10)[0]
lum = (tbl['st_rad'].data)**2 * (tbl['st_teff'].data / 5777)**4
tbl['st_lum'][lum_idx] = lum[lum_idx]
c_dict['lum'] = np.array(len(tbl['pl_name']) * ['default'])
c_dict['lum'][lum_idx] = '#ff0000'
mask_idx = np.where(tbl['st_lum'].mask)[0]
tbl['st_lum'][mask_idx] = 0

# %%
def min_val(x):
    """
    Return the minimum value in array (skipping masked elements)
    """
    return np.min(x[np.nonzero(x)])

def max_val(x):
    """
    Return the maximum value in array (skipping masked elements)
    """
    return np.max(x[np.nonzero(x)])

def logscale(x):
    """
    Return the normalized value in array (between 0 and 1) in log scale.
    """
    min_value = min_val(x)
    max_value = max_val(x)

    return (np.log10(x)-np.log10(min_value))/np.log10(max_value/min_value)

def scaleMag(x):
    """
    Arbitrary scale to plot the size of a data point according to the
    magnitude of the host star.
    """
    scale = 5*((20-x)/8)**3

    if type(scale) == float or type(scale) == int:
        return scale

    index, = np.where(x == 0)
    scale[index] = 5*((20-12)/8)**3 # Small size for system with unknown mag.
    index, = np.where(x > 14)
    scale[index] = 5*((20-14)/8)**3 # Minimum size for mag. > 14

    return scale

def find_na(x):
    """
    Finds element that are masked or NaN and return an empty string to display
    """
    x.mask = False
    index = np.array([])
    try:
        index, = np.where(np.isnan(x))
    except:
        pass

    index = np.append(index, np.where(x == 0)[0])
    listx = list(x)
    for i in index:
        listx[i] = ' '
    return listx

from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import CustomJS, Label, ColorBar, LogColorMapper, LogAxis, LogTicker
from bokeh.models.widgets import Slider, RangeSlider, Select, Div, TextInput, Button, CheckboxGroup, CheckboxButtonGroup
from bokeh.models.tools import TapTool
from bokeh.layouts import layout, row, column, Spacer
from bokeh.models import Scatter
from bokeh import events

# Current date
date = datetime.today().strftime('%Y-%m-%d')

# Available parameters
hostname = find_na(tbl['hostname'])
name = find_na(tbl['pl_name'])
ra_str = find_na(tbl['rastr'])
dec_str = find_na(tbl['decstr'])
ra = find_na(tbl['ra'])
dec = find_na(tbl['dec'])
stmass = find_na(tbl['st_mass'])
stradius = find_na(tbl['st_rad'])
stlum = find_na(tbl['st_lum'])
teff = find_na(tbl['st_teff'])
spt = find_na(tbl['st_spectype'])
dist = find_na(tbl['sy_dist'])
mass = find_na(tbl['pl_bmasse'])
rade = find_na(tbl['pl_rade'])
rho = find_na(tbl['pl_dens'])
period = find_na(tbl['pl_orbper'])
depth = find_na(tbl['pl_trandep'])
dur = find_na(tbl['pl_trandur'])
semima = find_na(tbl['pl_orbsmax'])
teq = find_na(tbl['pl_eqt'])
insol = find_na(tbl['pl_insol'])
K = find_na(tbl['pl_rvamp'])
atmosignal = find_na(tbl['pl_atmosig'])
tsm = find_na(tbl['pl_tsm'])
esm = find_na(tbl['pl_esm'])
vmag = find_na(tbl["sy_vmag"])
gmag = find_na(tbl["sy_gaiamag"])
jmag = find_na(tbl["sy_jmag"])
hmag = find_na(tbl["sy_hmag"])
kmag = find_na(tbl["sy_kmag"])

# Size of data points
size = scaleMag(tbl["sy_jmag"])

# Color of data points
c = logscale(tbl['pl_eqt'])
colors = plt.cm.viridis(c, 1, True)
c_dict['bokeh'] = np.array(["#%02x%02x%02x" % (r, g, b) for r, g, b in colors[:,0:3]])
index, = np.where(tbl['pl_eqt'] == 0)
c_dict['bokeh'][index] = '#707070'
c_dict['edge'] = np.array(len(name)*['black'])
c_dict['edge'][k2_idx] = 'blue'
c_dict['edge'][-number_toi:] = 'red'
edge_thick = 1 * np.ones(len(name))
edge_thick[k2_idx] = 1
edge_thick[-number_toi:] = 1
edge_alpha = 1 * np.ones(len(name))
fill_alpha = 0.8 * np.ones(len(name))

# Target type
targettype = np.array(len(name) * ["Exoplanet"])
targettype[k2_idx] = "K2"
targettype[-number_toi:] = "TOI"

# Changing x and y axis mapping
axis_map = {
    "Stellar Mass (M☉)": "stmass",
    "Stellar Radius (R☉)": "stradius",
    "Stellar Lum. (L☉)": "stlum",
    "Stellar Teff. (K)": "teff",
    "Stellar Distance (pc)": "dist",
    "Mass (M⊕)": "mass",
    "Radius (R⊕)": "rade",
    "Density (g/cm³)": "rho",
    "Period (days)": "period",
    "Orbital Separation (au)": "semima",
    "Planetary Teq (K)": "teq",
    "Insolation (S⊕)": "insol",
    "Transit Depth (ppt)": "depth",
    "Transit Duration (hr)": "dur",
    "Atmo. signal (ppm)": "atmosignal",
    "TSM": "tsm",
    "ESM": "esm",
    "RV Semi-amplitude (m/s)": "K",
    "V mag.": "vmag",
    "Gaia mag.": "gmag",
    "J mag.": "jmag",
    "H mag.": "hmag",
    "K mag.": "kmag"}
axis_map_list = [i for i in axis_map.keys()]

# Create Widgets
x_axis = Select(title = "X Axis",
                options = axis_map_list,
                value = "Period (days)",
                width = 260)
y_axis = Select(title="Y Axis",
                options = axis_map_list,
                value = "Radius (R⊕)",
                width = 260)
inverted = CheckboxButtonGroup(labels=["Inverse X", "Inverse Y"], width = 260)
color_select = Select(title = "Colorbar",
                      options = ["Planetary Teq (K)", "In a future release"],
                      value = "Planetary Teq (K)", width = 260)
ra_min_slider = Slider(title = "RA min (hr)",
                       start = 0,
                       end = 24,
                       value = 0,
                       step = 1,
                       width = 125)
ra_max_slider = Slider(title = "RA max (hr)",
                       start = 0,
                       end = 24,
                       value = 24,
                       step = 1,
                       width = 125)
dec_slider = RangeSlider(title = "Declination (deg.)",
                         start = -90,
                         end = 90,
                         value = (-90, 90),
                         step = 1,
                         width = 260)
teff_slider = RangeSlider(title = "Stellar Teff (K)",
                          start = 2500,
                          end = 15000,
                          value = (2500, 15000),
                          step = 100,
                          width = 260)
rade_slider = RangeSlider(title = "Radius (R⊕)",
                          start = 0.1,
                          end = 20,
                          value = (0.1, 20),
                          step = .1,
                          width = 260)
checkbox = CheckboxGroup(labels=["Exoplanets", "TOI", "K2"],
                         active = [0, 1, 2],
                         width = 260)
button = Button(label = "Submit", button_type = "success", width = 80)
reset = Button(label = "Reset", width = 80)
help_button = Button(label = "Help", button_type = "warning", width = 80)
search = TextInput(value = "", title = "Search:", width = 170)
search_button = Button(label = "Search", button_type = "primary", width = 80)
clear = Button(label = "Clear", width = 80)
download = Button(label = "Download to CSV", button_type = "success", width = 80)
text = Div(text="""
<p class="big">
Data from <a href="https://exoplanetarchive.ipac.caltech.edu/">Exoplanet Archive</a>
and <a href="https://exofop.ipac.caltech.edu/">ExoFop TESS</a>
</p>

<p class="big" style="margin-top: 1em;">
Code available on <a href="https://github.com/CharlesCadieux/iExoView" target="_blank">GitHub</a><br>
Contact: <a href="mailto:charles.cadieux.1@umontreal.ca">charles.cadieux.1@umontreal.ca</a>
</p>
""", height=60)
figure_title = Div(text = """<font size="+1"> <b> iExoView: The Interactive Exoplanet Viewer </b> </font>""",
                   height = 30)

source_available = ColumnDataSource(data = dict(x = period,
                                                y = rade,
                                                hostname = hostname,
                                                name = name,
                                                targettype = targettype,
                                                ra = ra,
                                                dec = dec,
                                                ra_str = ra_str,
                                                dec_str = dec_str,
                                                stmass = stmass,
                                                stradius = stradius,
                                                stlum = stlum,
                                                teff = teff,
                                                spt = spt,
                                                dist = dist,
                                                mass = mass,
                                                rade = rade,
                                                rho = rho,
                                                period = period,
                                                depth = depth,
                                                semima = semima,
                                                teq = teq,
                                                insol = insol,
                                                atmosignal = atmosignal,
                                                tsm = tsm,
                                                esm = esm,
                                                K = K,
                                                dur = dur,
                                                vmag = vmag,
                                                gmag = gmag,
                                                jmag = jmag,
                                                hmag = hmag,
                                                kmag = kmag,
                                                size = size,
                                                bokeh_colors = c_dict['bokeh'],
                                                edge_color = c_dict['edge'],
                                                edge_alpha = edge_alpha,
                                                edge_thick = edge_thick,
                                                fill_alpha = fill_alpha,
                                                semima_colors = c_dict['semima'],
                                                teq_colors = c_dict['teq'],
                                                rho_colors = c_dict['rho'],
                                                insol_colors = c_dict['insol'],
                                                spt_colors = c_dict['spt'],
                                                mass_colors = c_dict['mass'],
                                                K_colors = c_dict['K'],
                                                red_colors = c_dict['red'],
                                                lum_colors = c_dict['lum'],
                                                depth_colors = c_dict['depth']))

source_visible = ColumnDataSource(data = dict(x = source_available.data[axis_map[x_axis.value]],
                                              y = source_available.data[axis_map[y_axis.value]],
                                              hostname = hostname,
                                              name = name,
                                              targettype = targettype,
                                              ra = ra,
                                              dec = dec,
                                              ra_str = ra_str,
                                              dec_str = dec_str,
                                              stmass = stmass,
                                              stradius = stradius,
                                              stlum = stlum,
                                              teff = teff,
                                              spt = spt,
                                              dist = dist,
                                              mass = mass,
                                              rade = rade,
                                              rho = rho,
                                              period = period,
                                              depth = depth,
                                              semima = semima,
                                              teq = teq,
                                              insol = insol,
                                              atmosignal = atmosignal,
                                              tsm = tsm,
                                              esm = esm,
                                              K = K,
                                              dur = dur,
                                              vmag = vmag,
                                              gmag = gmag,
                                              jmag = jmag,
                                              hmag = hmag,
                                              kmag = kmag,
                                              size = size,
                                              bokeh_colors = c_dict['bokeh'],
                                              edge_color = c_dict['edge'],
                                              edge_alpha = edge_alpha,
                                              edge_thick = edge_thick,
                                              fill_alpha = fill_alpha,
                                              semima_colors = c_dict['semima'],
                                              teq_colors = c_dict['teq'],
                                              rho_colors = c_dict['rho'],
                                              insol_colors = c_dict['insol'],
                                              spt_colors = c_dict['spt'],
                                              mass_colors = c_dict['mass'],
                                              K_colors = c_dict['K'],
                                              red_colors = c_dict['red'],
                                              lum_colors = c_dict['lum'],
                                              depth_colors = c_dict['depth']))

# Tools to interact with the figure (documentation available at
# https://docs.bokeh.org/en/latest/index.html)
tap = TapTool()
TOOLS = ["crosshair",
         "hover",
         "pan",
         "box_zoom",
         "undo",
         "redo",
         "reset",
         "save",
         tap]

TOOLTIPS = """
<table>
  <tr>
    <th colspan="2">Stellar Parameters</th>
    <th colspan="2">Planetary Parameters</th>
  </tr>
  <tr>
    <td><span style="color: #2874a6;">Host name</span></td>
    <td>@hostname</td>
    <td><span style="color: #2874a6;">Name</span></td>
    <td>@name</td>
  </tr>
  <tr>
    <td><span style="color: #2874a6;">RA</span></td>
    <td>@ra_str</td>
    <td><span style="color: #2874a6;">Mass (M⊕)</span></td>
    <td><span style="color: @mass_colors;">@mass</span></td>
  </tr>
  <tr>
    <td><span style="color: #2874a6;">Dec.</span></td>
    <td>@dec_str</td>
    <td><span style="color: #2874a6;">Radius (R⊕)</span></td>
    <td>@rade</td>
  </tr>
  <tr>
    <td><span style="color: #2874a6;">Mass (M⊙)</span></td>
    <td>@stmass</td>
    <td><span style="color: #2874a6;">Density (g/cm³)</span></td>
    <td><span style="color: @rho_colors;">@rho</span></td>
  </tr>
    <tr>
    <td><span style="color: #2874a6;">Radius (R⊙)</span></td>
    <td>@stradius</td>
    <td><span style="color: #2874a6;">Period (days)</span></td>
    <td>@period</td>
  </tr>
  <tr>
    <td><span style="color: #2874a6;">Lum. (L☉)</span></td>
    <td><span style="color: @lum_colors;">@stlum</span></td>
    <td><span style="color: #2874a6;">Transit depth (ppt)</span></td>
    <td><span style="color: @depth_colors;">@depth</span></td>
  </tr>
  <tr>
    <td><span style="color: #2874a6;">Teff (K)</span></td>
    <td>@teff</td>
    <td><span style="color: #2874a6;">Transit dur. (hr)</span></td>
    <td>@dur</td>
  </tr>
  <tr>
    <td><span style="color: #2874a6;">Spect. type</span></td>
    <td><span style="color: @spt_colors;">@spt</span></td>
    <td><span style="color: #2874a6;">Orb. sep. (au)</span></td>
    <td><span style="color: @semima_colors;">@semima</span></td>
  </tr>
  <tr>
    <td><span style="color: #2874a6;">Dist. (pc)</span></td>
    <td>@dist</td>
    <td><span style="color: #2874a6;">Teq (K)</span></td>
    <td><span style="color: @teq_colors;">@teq</span></td>
  </tr>
  <tr>
    <td><span style="color: #2874a6;">V</span></td>
    <td>@vmag</td>
    <td><span style="color: #2874a6;">Insolation (S⊕)</span></td>
    <td><span style="color: @insol_colors;">@insol</span></td>
  </tr>
    <tr>
    <td><span style="color: #2874a6;">Gaia</span></td>
    <td>@gmag</td>
    <td><span style="color: #2874a6;">RV K (m/s)</span></td>
    <td><span style="color: @K_colors;">@K</span></td>

  </tr>
  <tr>
    <td><span style="color: #2874a6;">J</span></td>
    <td>@jmag</td>
    <td><span style="color: #2874a6;">Atmo. sig. (ppm)</span></td>
    <td><span style="color: @red_colors;">@atmosignal</span></td>
  </tr>
    <tr>
    <td><span style="color: #2874a6;">H</span></td>
    <td>@hmag</td>
    <td><span style="color: #2874a6;">TSM</span></td>
    <td><span style="color: @red_colors;">@tsm</span></td>

  </tr>
    <tr>
    <td><span style="color: #2874a6;">K</span></td>
    <td>@kmag</td>
    <td><span style="color: #2874a6;">ESM</span></td>
    <td><span style="color: @red_colors;">@esm</span></td>
  </tr>
</table>
"""

# Generating the figure
titlestr = ("{0} Transiting Exoplanets with Known Mass + {1} TOI + {2} K2 Planets (as of {3})"
           ).format(len(tbl) - number_toi - number_k2, number_toi, number_k2, date)
p = figure(width = 1120,
           height=630,
           tools = TOOLS,
           x_axis_type = "log",
           y_axis_type = "log",
           x_axis_location = 'above',
           y_axis_location = 'right',
           tooltips = TOOLTIPS,
           title = titlestr,
           title_location = 'below',
           toolbar_location = "left")
p.title.text_font_size = '14pt'
p.axis.axis_label_text_font_size = '14pt'
p.axis.major_label_text_font_size = '12pt'
circle = p.scatter('x',
                   'y',
                   source=source_visible,
                   size='size',
                   fill_color='bokeh_colors',
                   fill_alpha='fill_alpha',
                   line_color='edge_color',
                   line_alpha='edge_alpha',
                   line_width='edge_thick')

selected_circle = Scatter(
    marker="circle",
    size="size",
    fill_color='bokeh_colors',
    fill_alpha='fill_alpha',
    line_color='edge_color',
    line_alpha='edge_alpha',
    line_width='edge_thick'
)

nonselected_circle = Scatter(
    marker="circle",
    size="size",
    fill_color='bokeh_colors',
    fill_alpha='fill_alpha',
    line_color='edge_color',
    line_alpha='edge_alpha',
    line_width='edge_thick'
)

circle.selection_glyph = selected_circle
circle.nonselection_glyph = nonselected_circle
xaxis = LogAxis(axis_label = x_axis.value,
               ticker = LogTicker(num_minor_ticks = 10),
               axis_label_text_font_size = '14pt')
xaxis.major_label_text_font_size = '12pt'
yaxis = LogAxis(axis_label = y_axis.value,
                ticker = LogTicker(num_minor_ticks = 10),
                axis_label_text_font_size = '14pt')
yaxis.major_label_text_font_size = '12pt'
p.add_layout(xaxis, 'below')
p.add_layout(yaxis, 'left')

#### User interaction ####

# Opens a new window when clicking on a data point
with open('js/callback_tap.js') as f:
    script = f.read()
callback_tap = CustomJS(args = dict(source_visible = source_visible,
                                    TOOLTIPS = TOOLTIPS),
                        code = script)

# Submit button. Update the figure based on user inputs
with open('js/callback_button.js') as f:
    script = f.read()
callback_button = CustomJS(args = dict(source_available = source_available,
                                       source_visible = source_visible,
                                       axis_map = axis_map,
                                       ra_min_slider = ra_min_slider,
                                       ra_max_slider = ra_max_slider,
                                       dec_slider = dec_slider,
                                       teff_slider = teff_slider,
                                       rade_slider = rade_slider,
                                       x_axis = x_axis,
                                       y_axis = y_axis,
                                       search = search,
                                       checkbox = checkbox,
                                       button = button,
                                       title = p.title,
                                       p = p,
                                       inverted = inverted,
                                       date = date),
                           code = script)

# Change x axis
with open('js/callback_x_axis.js') as f:
    script = f.read()
callback_x_axis = CustomJS(args = dict(source_visible = source_visible,
                                       axis_map = axis_map,
                                       xaxis = xaxis,
                                       p = p,
                                       inverted = inverted),
                           code = script)

# Change y axis
with open('js/callback_y_axis.js') as f:
    script = f.read()
callback_y_axis = CustomJS(args = dict(source_visible = source_visible,
                                       axis_map = axis_map,
                                       yaxis = yaxis,
                                       p = p,
                                       inverted = inverted),
                           code = script)

# Invert x and y axis
with open('js/callback_invert_axis.js') as f:
    script = f.read()
callback_invert_axis = CustomJS(args = dict(xr = p.x_range,
                                            yr = p.y_range),
                                code = script)

# Search engine
with open('js/callback_search.js') as f:
    script = f.read()
callback_search = CustomJS(args = dict(source_visible = source_visible,
                                       p = p,
                                       inverted = inverted),
                           code = script)

# Reset figure
with open('js/callback_reset.js') as f:
    script = f.read()
callback_reset = CustomJS(args = dict(p = p,
                                      source_available = source_available,
                                      source_visible = source_visible,
                                      ra_min_slider = ra_min_slider,
                                      ra_max_slider = ra_max_slider,
                                      dec_slider = dec_slider,
                                      teff_slider = teff_slider,
                                      rade_slider = rade_slider,
                                      x_axis = x_axis,
                                      y_axis = y_axis,
                                      search = search,
                                      button = button,
                                      checkbox = checkbox,
                                      inverted = inverted,
                                      title = p.title,
                                      titlestr = titlestr),
                           code = script)

# Clear search
with open('js/callback_clear.js') as f:
    script = f.read()
callback_clear = CustomJS(args = dict(source_visible = source_visible,
                                      search = search,
                                      p = p,
                                      inverted = inverted),
                           code = script)

# Help button
with open('js/callback_help.js') as f:
    script = f.read()
callback_help = CustomJS(code = script)

# Download to .CSV button
with open('js/callback_download.js') as f:
    script = f.read()
callback_download = CustomJS(args=dict(source = source_visible),
                             code = script)
tap.callback = callback_tap
x_axis.js_on_change('value', callback_x_axis)
y_axis.js_on_change('value', callback_y_axis)
search.js_on_change('value', callback_search)
button.js_on_event(events.ButtonClick, callback_button)
inverted.js_on_change('active', callback_invert_axis)
reset.js_on_event(events.ButtonClick, callback_reset)
help_button.js_on_event(events.ButtonClick, callback_help)
clear.js_on_event(events.ButtonClick, callback_clear)
download.js_on_event(events.ButtonClick, callback_download)



# %%

# Figure legend
x0, y0 = 0, 1
legend = figure(title = "Legend",
                align = 'end',
                x_range = (-0.1, 0.5),
                y_range = (0.1, 1.5),
                y_axis_type = 'log',
                width = 100,
                height = 110,
                tools = '',
                toolbar_location = None,
                min_border = 0,
                outline_line_color = None)
legend.title.text_font_size = '10pt'
legend.title.text_font_style = "normal"
legend.axis.visible = False
legend.xgrid.grid_line_color = None
legend.ygrid.grid_line_color = None
legend.scatter(x0,
              y0,
              size = 15,
              line_color = "black",
              fill_color = "white",
              line_width = 1)
legend.scatter(x0,
              y0 - 0.5,
              size = 15,
              line_color = "red" ,
              fill_color = "white",
              line_width = 1)
legend.scatter(x0,
              y0 - 0.75,
              size = 15,
              line_color = "blue",
              fill_color = "white",
              line_width=1)
label1 = Label(x = x0,
               y = y0,
               text = "Exoplanets",
               x_offset = 10,
               y_offset = -8,
               text_font_size = '10pt')
label2 = Label(x = x0,
               y = y0 - 0.5,
               text = "TOI",
               text_color = "red",
               x_offset = 10,
               y_offset = -8,
               text_font_size = '10pt')
label3 = Label(x = x0,
               y = y0 - 0.75,
               text = "K2",
               text_color = "blue",
               x_offset = 10,
               y_offset = -8,
               text_font_size = '10pt')
legend.add_layout(label1)
legend.add_layout(label2)
legend.add_layout(label3)

# J mag legend
mags = range(6, 13, 2)
j_mag_legend = figure(title = "J mag.",
                      align = "end",
                      x_range = (-0.2, 1),
                      y_range = (0.2, 1.3),
                      y_axis_type = 'log',
                      width = 100,
                      height = 180,
                      tools = '',
                      toolbar_location = None,
                      min_border = 0,
                      outline_line_color = None)
j_mag_legend.title.text_font_size = '10pt'
j_mag_legend.title.text_font_style = "normal"
j_mag_legend.axis.visible = False
j_mag_legend.xgrid.grid_line_color = None
j_mag_legend.ygrid.grid_line_color = None
for i, m in enumerate(mags):
    j_mag_legend.scatter(x0,
                        y0/1.4**i,
                        size = scaleMag(m),
                        line_color = "black",
                        fill_color="white",
                        line_width=1)
    label = Label(x = x0,
                  y = y0/1.4**i,
                  text = "{0}".format(m),
                  x_offset=15,
                  y_offset = -10,
                  text_font_size = '10pt')
    j_mag_legend.add_layout(label)

j_mag_legend.scatter(x0,
                    y0/1.4**4,
                    size = scaleMag(14),
                    line_color = "black",
                    fill_color="white",
                    line_width=1)
label = Label(x = x0,
              y = y0/1.4**4,
              text = "≥14",
              x_offset=15,
              y_offset = -10,
              text_font_size = '10pt')
j_mag_legend.add_layout(label)

# Color map
color_bar = ColorBar(color_mapper = LogColorMapper(palette = "Viridis256",
                                            low = min_val(tbl['pl_eqt']),
                                            high = max_val(tbl['pl_eqt'])),
                     ticker = LogTicker(),
                     label_standoff = 10,
                     border_line_color = None,
                     location = (0, 0),
                     title = 'Planetary Teq (K)',
                     title_standoff = 10,
                     title_text_align = "center",
                     title_text_font_style = 'italic',
                     title_text_font_size='12pt')
p.add_layout(color_bar, 'right')

# Layout
ra_sliders = row(ra_min_slider, ra_max_slider)
buttons = row(button, reset, help_button)
legends = column(legend, j_mag_legend, align = "end")
search_buttons = row(search_button, clear)
widgets = column(x_axis,
                 y_axis,
                 inverted,
                 color_select,
                 ra_sliders,
                 dec_slider,
                 teff_slider,
                 rade_slider,
                 checkbox,
                 buttons,
                 search,
                 search_buttons,
                 download,
                 text,
                 width=300)
# Layout
legends = column(legend, j_mag_legend, align="end")
layout = row(widgets, Spacer(width=10), p, legends)
layout = column(figure_title, layout)

# %%
# Save as .html
output_file("iExoView.html", title = 'iExoView')
show(layout)

# %%



