"""
iExoView, the Interactive Exoplanet Viewer
astro.umontreal.ca/~charles/iExoView.html

@author: CharlesCadieux 2022
"""

import numpy as np
import matplotlib.pyplot as plt

import requests

from astropy.table import Table, Row, Column
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import Angle

from datetime import datetime

from exofile.archive import ExoFile

from utils import *

##### Import data #####

# NASA Exoplanet Archive data using exofile (github.com/AntoineDarveau/exofile)
tbl = ExoFile.load()
# Unit conversion
tbl['pl_trandep'] *= 10 # Transit depth from percents to part-per-thousands (ppt)
tbl['st_lum'] = 10**tbl['st_lum'] # Stellar luminosity in solar lum.

# TOI data from ExoFOP
TOI = requests.get(
    "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv")
# From txt to Table
colnames = TOI.text[0:965]
colnames = list(colnames.split(","))
tblTOI = Table.read(TOI.text, format = "ascii", names = colnames, data_start = 1)
# Remove Known Planets and False Positive TOIs
tblTOI['TFOPWG Disposition'].mask = False
indexTOI, = np.where((tblTOI['TFOPWG Disposition'] == '0') |
                     (tblTOI['TFOPWG Disposition'] == 'PC') |
                     (tblTOI['TFOPWG Disposition'] == 'CP'))
tblTOI = tblTOI[indexTOI]
numberTOI = len(tblTOI)

# Table 5 from Pecaut & Mamajek 2013
spectral = requests.get(
    "http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt")
tblSpT = Table.read(spectral.text,
                    format = 'ascii.csv',
                    delimiter = ' ',
                    comment = '&',
                    header_start = 22,
                    data_start = 23,
                    data_end = 141)
del tblSpT['#SpT_1']
tblSpT.rename_column('#SpT', 'SpT')


##### Create master table with Known exoplanets, TOIs and K2 planets #####

#DEBUG
numberTOI = 1

# Populate TOIs
for i in range(numberTOI):
    TICid = tblTOI['TIC ID'][i]
    try:
        # Open the TOI page on ExoFOP to obtain additional parameters
        ExoFop = requests.get(
            "https://exofop.ipac.caltech.edu/tess/download_target.php?id={0}"
            .format(TICid))

        # Host Star Parameters
        index1 = ExoFop.text.find("STELLAR PARMAMETERS")
        index2 = ExoFop.text.find("MAGNITUDES")
        st = Table.read(ExoFop.text[index1:index2:],
               format = 'ascii.fixed_width',
               col_starts=(0, 25, 43, 65, 87, 109, 131, 153, 175, 197, 219, 241,
               263, 285, 307, 329, 351, 373, 395, 417, 439, 461, 483, 505, 530,
               552, 574, 599, 621, 643, 665, 687, 709, 731, 749, 767),
               header_start = 1,
               data_start = 2)

        # Magnitudes
        index1 = ExoFop.text.find("MAGNITUDES")
        index2 = ExoFop.text.find("IMAGING OBSERVATIONS")
        mag = Table.read(ExoFop.text[index1:index2:],
            format = 'ascii.fixed_width',
            col_starts=(0, 18, 36, 54, 79, 99, 117, 135),
            header_start = 1,
            data_start = 2)
        indexV, = np.where(mag['Band'] == 'V')
        indexG, = np.where(mag['Band'] == 'Gaia')
        indexJ, = np.where(mag['Band'] == 'J')
        indexH, = np.where(mag['Band'] == 'H')
        indexK, = np.where(mag['Band'] == 'K')

        # TOI planetary parameter estimation
        masse = mass_est(tblTOI['Planet Radius (R_Earth)'][i])
        rho = bulkrho(masse, tblTOI['Planet Radius (R_Earth)'][i])
        K = semiamp(tblTOI['Period (days)'][i],
                    masse,
                    tblTOI['Planet Radius (R_Earth)'][i],
                    float(st['Mass (M_Sun)'][0]),
                    90,
                    0)

    except:
        print('Error with TIC {0}'.format(TICid))
        numberTOI -= 1
        continue

    # Add row to master tbl
    tbl.add_row()
    tbl[-1]['hostname'] = 'TIC {0}'.format(tblTOI['TIC ID'][i])
    tbl[-1]['rastr'] = str(tblTOI["RA"][i])
    tbl[-1]['decstr'] = str(tblTOI["Dec"][i])
    tbl[-1]['ra'] = Angle(tblTOI["RA"][i], unit=u.hour).deg
    tbl[-1]['dec'] = Angle(tblTOI["Dec"][i], unit = u.deg).deg
    tbl[-1]['st_rad'] = tblTOI['Stellar Radius (R_Sun)'][i]
    tbl[-1]['st_mass'] = st['Mass (M_Sun)'][0]
    tbl[-1]['st_lum'] = st['Luminosity'][0]
    tbl[-1]['st_teff'] = tblTOI['Stellar Eff Temp (K)'][i]
    tbl[-1]['st_spectype'] = 0
    tbl[-1]['sy_dist'] = tblTOI['Stellar Distance (pc)'][i]
    if len(indexV) == 1 : tbl[-1]['sy_vmag'] = mag['Value'][indexV]
    if len(indexG) == 1 : tbl[-1]['sy_gaiamag'] = mag['Value'][indexG]
    if len(indexJ) == 1 : tbl[-1]['sy_jmag'] = mag['Value'][indexJ]
    if len(indexH) == 1 : tbl[-1]['sy_hmag'] = mag['Value'][indexH]
    if len(indexK) == 1 : tbl[-1]['sy_kmag'] = mag['Value'][indexK]
    tbl[-1]['pl_name'] = 'TOI {0}'.format(tblTOI['TOI'][i])
    tbl[-1]['pl_bmassprov'] = 'Mass Estimate'
    tbl[-1]['pl_bmasse'] = np.round(masse, 2)
    tbl[-1]['pl_rade'] = tblTOI['Planet Radius (R_Earth)'][i]
    tbl[-1]['pl_dens'] = np.round(rho, 2)
    tbl[-1]['pl_orbper'] = tblTOI['Period (days)'][i]
    tbl[-1]['pl_trandur'] = tblTOI['Duration (hours)'][i]
    tbl[-1]['pl_orbsmax'] = 0 # To be estimated later
    tbl[-1]['pl_eqt'] = tblTOI['Planet Equil Temp (K)'][i]
    tbl[-1]['pl_insol'] = tblTOI['Planet Insolation (Earth Flux)'][i]
    tbl[-1]['pl_rvamp'] = np.round(K, 2)
    tbl[-1]['pl_bmasselim'] = 0
    tbl[-1]['tran_flag'] = 1
    tbl[-1]['pl_trandep'] = tblTOI['Depth (ppm)'][i] / 1000 # in ppt
    tbl[-1]['disc_facility'] = 'Transiting Exoplanet Survey Satellite (TESS)'


# Populate K2 planets with unknown mass
indexK2, = np.where((tbl['disc_facility'] == 'K2') &
                    (tbl['tran_flag'] == 1))
massindex, = np.where((tbl['pl_bmassprov'] == 'Mass') &
                      (tbl['pl_bmasselim'] == 0))
indexK2 = np.setdiff1d(indexK2, massindex)
numberK2 = len(indexK2)

for i in indexK2:
    # K2 planets parameter estimation
    masse = mass_est(tbl['pl_rade'][i])
    rho = bulkrho(masse, tbl['pl_rade'][i])
    K = semiamp(tbl['pl_orbper'][i],
                masse,
                tbl['pl_rade'][i],
                tbl['st_mass'][i],
                90,
                0)
    tbl[i]['pl_bmassprov'] = 'Mass Estimate'
    tbl[i]['pl_bmasse'] = masse
    tbl[i]['pl_dens'] = rho
    tbl[i]['pl_rvamp'] = K
    tbl[i]['pl_bmasselim'] = 0


##### Clean master table #####

# Keep only transiting planets with known radius
good, = np.where((tbl['tran_flag'] == 1) &
                 (tbl['pl_rade'] > 0))
tbl = tbl[good]
good, = np.where((tbl['pl_bmassprov'] == 'Mass') |
                 (tbl['pl_bmassprov'] == 'Msini') |
                 (tbl['pl_bmassprov'] == 'Msin(i)/sin(i)') |
                 (tbl['pl_bmassprov'] == 'Mass Estimate'))
tbl = tbl[good]
# Keep only planets with known mass (no upper/lower limit)
good, = np.where(tbl['pl_bmasselim'] == 0)
tbl = tbl[good]
# Update indexK2
indexK2, = np.where((tbl['disc_facility'] == 'K2') &
                    (tbl['pl_bmassprov'] == 'Mass Estimate'))


##### Fill in missing parameters #####

# Color dictionary for print output (black: observed value, red: estimation)
c_dict = {}
c_dict['mass'] = np.array(len(tbl['pl_name'])*['default'])
c_dict['mass'][indexK2] = '#ff0000'
c_dict['mass'][-1 * numberTOI:] = '#ff0000'

# Semi-major axis
tbl['pl_orbsmax'].mask = False
sma = semima(tbl['pl_orbper'].data, tbl['st_mass'].data)
index, = np.where(tbl['pl_orbsmax'] == 0)
tbl['pl_orbsmax'][index] = sma[index]
c_dict['semima'] = np.array(len(tbl['pl_name'])*['default'])
c_dict['semima'][index] = '#ff0000'
index, = np.where(tbl['pl_orbsmax'].mask == 1)
tbl['pl_orbsmax'][index] = 0

# Equilibrium temperature
tbl['pl_eqt'].mask = False
teq = Teq(tbl['st_rad'].data, tbl['st_teff'].data, tbl['pl_orbsmax'].data)
index, = np.where(tbl['pl_eqt'] == 0)
tbl['pl_eqt'][index] = teq[index]
c_dict['teq'] = np.array(len(tbl['pl_name'])*['default'])
c_dict['teq'][index] = '#ff0000'
index, = np.where(tbl['pl_eqt'].mask == 1)
tbl['pl_eqt'][index] = 0

# Bulk density
tbl['pl_dens'].mask = False
index, = np.where(tbl['pl_dens'].data == 0)
density = bulkrho(tbl['pl_bmasse'].data, tbl['pl_rade'].data)
tbl['pl_dens'][index] = density[index]
c_dict['rho'] = np.array(len(tbl['pl_name'])*['default'])
c_dict['rho'][index] = '#ff0000'
c_dict['rho'][indexK2] = '#ff0000'
c_dict['rho'][-1 * numberTOI:] = '#ff0000'
index, = np.where(tbl['pl_dens'].mask == 1)
tbl['pl_dens'][index] = 0

# Insolation
tbl['pl_insol'].mask = False
index, = np.where(tbl['pl_insol'].data == 0)
insol = (tbl['st_rad'].data**2 * (tbl['st_teff'].data / 5777)**4
        / tbl['pl_orbsmax'].data**2)
tbl['pl_insol'][index] = insol[index]
c_dict['insol'] = np.array(len(tbl['pl_name'])*['default'])
c_dict['insol'][index] = '#ff0000'
index, = np.where(tbl['pl_insol'].mask == 1)
tbl['pl_insol'][index] = 0

# RV semi-amplitude
tbl['pl_rvamp'].mask = False
index, = np.where(tbl['pl_rvamp'].data == 0)
K = semiamp(tbl['pl_orbper'].data,
            tbl['pl_bmasse'].data,
            tbl['pl_rade'].data,
            tbl['st_mass'].data,
            90,
            0)
tbl['pl_rvamp'][index] = K[index]
c_dict['K'] = np.array(len(tbl['pl_name'])*['default'])
c_dict['K'][index] = '#ff0000'
c_dict['K'][indexK2] = '#ff0000'
c_dict['K'][-1 * numberTOI:] = '#ff0000'
index, = np.where(tbl['pl_rvamp'].mask == 1)
tbl['pl_rvamp'][index] = 0

# Add atmospheric signal in master table
atmosigcol = Column(name = 'pl_atmosig', data = np.zeros(len(tbl['pl_name'])))
tbl.add_column(atmosigcol, 66)
mu_array = np.zeros(len(tbl['pl_name']))
for i in range(len(tbl['pl_name'])):
    if tbl['pl_rade'][i] <= 2 : mu = 28.97 # Earth's atmosphere mu
    elif 2 < tbl['pl_rade'][i] <= 6 : mu = 2.61  # Neptune's atmosphere mu
    elif tbl['pl_rade'][i] > 6 : mu = 2.22  # Jupiter's atmosphere mu
    else: print('Planetary radius unknown for {0}'.format(tbl['pl_name'][i]))
    mu_array[i] = mu

tbl['pl_atmosig'] = atmosig(tbl['pl_eqt'].data,
                            tbl['st_rad'].data,
                            tbl['pl_dens'].data,
                            mu_array)
c_dict['red'] = np.array(len(tbl['pl_name'])*['#ff0000'])
index, = np.where(tbl['pl_atmosig'].mask == 1)
tbl['pl_atmosig'][index] = 0

# Add TSM column in master table
tsmcol = Column(name = 'pl_tsm', data = np.zeros(len(tbl['pl_name'])))
tbl.add_column(tsmcol, 67)
tsm = np.zeros(len(tbl['pl_name']))
for i in range(len(tsm)):
    tsm[i] = TSM(tbl['pl_rade'][i], tbl['pl_bmasse'][i], tbl['st_rad'][i],
             tbl['st_teff'][i], tbl['pl_orbsmax'][i], tbl['sy_jmag'][i])
index, = np.where((np.isnan(tsm)) | (np.isinf(tsm)))
tsm[index] = 0
tbl['pl_tsm'] = tsm
index, = np.where(tbl['pl_tsm'].mask == 1)
tbl['pl_tsm'][index] = 0

# Add ESM column in master table
esmcol = Column(name = 'pl_esm', data = np.zeros(len(tbl['pl_name'])))
tbl.add_column(esmcol, 68)
esm = np.zeros(len(tbl['pl_name']))
for i in range(len(esm)):
    esm[i] = ESM(tbl['pl_rade'][i], tbl['st_rad'][i], tbl['st_teff'][i],
                 tbl['pl_orbsmax'][i], tbl["sy_kmag"][i])
index, = np.where((np.isnan(esm)) | (np.isinf(esm)))
esm[index] = 0
tbl['pl_esm'] = esm
index, = np.where(tbl['pl_esm'].mask == 1)
tbl['pl_esm'][index] = 0

# Spectral type
tbl['st_spectype'].mask = False
tbl['st_teff'].mask = False
index, = np.where(tbl['st_spectype'].data == '0')
c_dict['spt'] = np.array(len(tbl['pl_name'])*['default'])
for i in index:
    if tbl['st_teff'][i] > 0:
        tbl['st_spectype'][i] = tblSpT['SpT'][find_nearest(tblSpT['Teff'], tbl['st_teff'][i])]
        c_dict['spt'][i] = '#ff0000'
    else : continue
index, = np.where(tbl['st_spectype'].mask == 1)
tbl['st_spectype'][index] = 0

# Transit depth
tbl['pl_trandep'].mask = False
index, = np.where(tbl['pl_trandep'].data == 0)
depth = ((6.371e6 / 6.956e8)**2
        * (tbl['pl_rade'].data / tbl['st_rad'].data)**2 * 1000) # in ppt
tbl['pl_trandep'][index] = depth[index]
c_dict['depth'] = np.array(len(tbl['pl_name'])*['default'])
c_dict['depth'][index] = '#ff0000'
index, = np.where(tbl['pl_trandep'].mask == 1)
tbl['pl_trandep'][index] = 0

# Luminosity
tbl['st_lum'].mask = False
index, = np.where(tbl['st_lum'].data == 10)
lum = (tbl['st_rad'].data)**2 * (tbl['st_teff'].data / 5777)**4
tbl['st_lum'][index] = lum[index]
c_dict['lum'] = np.array(len(tbl['pl_name'])*['default'])
c_dict['lum'][index] = '#ff0000'
index, = np.where(tbl['st_lum'].mask == 1)
tbl['st_lum'][index] = 0


#### Bokeh interactive plot #####

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

def findNA(x):
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
    return np.array(listx) # Array works better with bokeh

from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models.markers import Circle
from bokeh.models import CustomJS, Label, ColorBar, LogColorMapper, LogAxis, LogTicker
from bokeh.models.widgets import Slider, RangeSlider, Select, Div, TextInput, Button, CheckboxGroup, CheckboxButtonGroup
from bokeh.models.tools import TapTool
from bokeh.layouts import layout, row, column
from bokeh.models.annotations import Title
from bokeh import events

# Current date
date = datetime.today().strftime('%Y-%m-%d')

# Available parameters
hostname = findNA(tbl['hostname'].data)
name = findNA(tbl['pl_name'].data)
ra_str = findNA(tbl['rastr'].data)
dec_str = findNA(tbl['decstr'].data)
ra = findNA(tbl['ra'].data)
dec = findNA(tbl['dec'].data)
stmass = findNA(tbl['st_mass'].data)
stradius = findNA(tbl['st_rad'].data)
stlum = findNA(tbl['st_lum'].data)
teff = findNA(tbl['st_teff'].data)
spt = findNA(tbl['st_spectype'].data)
dist = findNA(tbl['sy_dist'].data)
mass = findNA(tbl['pl_bmasse'].data)
rade = findNA(tbl['pl_rade'].data)
rho = findNA(tbl['pl_dens'].data)
period = findNA(tbl['pl_orbper'].data)
depth = findNA(tbl['pl_trandep'].data)
dur = findNA(tbl['pl_trandur'].data)
semima = findNA(tbl['pl_orbsmax'].data)
teq = findNA(tbl['pl_eqt'].data)
insol = findNA(tbl['pl_insol'].data)
K = findNA(tbl['pl_rvamp'].data)
atmosignal = findNA(tbl['pl_atmosig'].data)
tsm = findNA(tbl['pl_tsm'].data)
esm = findNA(tbl['pl_esm'].data)
vmag = findNA(tbl["sy_vmag"].data)
gmag = findNA(tbl["sy_gaiamag"].data)
jmag = findNA(tbl["sy_jmag"].data)
hmag = findNA(tbl["sy_hmag"].data)
kmag = findNA(tbl["sy_kmag"].data)

# Size of data points
size = scaleMag(tbl["sy_jmag"].data)

# Color of data points
c = logscale(tbl['pl_eqt'].data)
colors = plt.cm.viridis(c, 1, True)
c_dict['bokeh'] = np.array(["#%02x%02x%02x" % (r, g, b) for r, g, b in colors[:,0:3]])
index, = np.where(tbl['pl_eqt'] == 0)
c_dict['bokeh'][index] = '#707070'
c_dict['edge'] = np.array(len(name)*['black'])
c_dict['edge'][indexK2] = 'blue'
c_dict['edge'][-1 * numberTOI:] = 'red'
edge_thick = 1 * np.ones(len(name))
edge_thick[indexK2] = 1
edge_thick[-1 * numberTOI:] = 1
edge_alpha = 1 * np.ones(len(name))
fill_alpha = 0.8 * np.ones(len(name))

# Target type
targettype = np.array(len(name) * ["Exoplanet"])
targettype[indexK2] = "K2"
targettype[-1 * numberTOI:] = "TOI"

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
                          end = 40,
                          value = (0.1, 40),
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
Update every Monday<br>
Exoplanets data from <a href="https://exoplanetarchive.ipac.caltech.edu/">Exoplanet Archive</a> <br>
TOI data from <a href="https://exofop.ipac.caltech.edu/">ExoFop TESS</a><br>
</p>""", height = 40)
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
           ).format(len(tbl) - numberTOI - numberK2, numberTOI, numberK2, date)
p = figure(plot_width = 1200,
           plot_height=700,
           tools = TOOLS,
           x_axis_type = "log",
           y_axis_type = "log",
           x_axis_location = 'above',
           y_axis_location = 'right',
           tooltips = TOOLTIPS,
           title = titlestr,
           title_location = 'below',
           toolbar_location = "left")
p.title.text_font_size = '12pt'
circle = p.circle('x',
                  'y',
                  source = source_visible,
                  size = 'size',
                  fill_color = 'bokeh_colors',
                  fill_alpha = 'fill_alpha',
                  line_color = 'edge_color',
                  line_alpha = 'edge_alpha',
                  line_width = 'edge_thick')
selected_circle = Circle(fill_color = 'bokeh_colors',
                         fill_alpha = 'fill_alpha',
                         line_color = 'edge_color',
                         line_alpha = 'edge_alpha',
                         line_width = 'edge_thick')
nonselected_circle = Circle(fill_color = 'bokeh_colors',
                            fill_alpha = 'fill_alpha',
                            line_color = 'edge_color',
                            line_alpha = 'edge_alpha',
                            line_width = 'edge_thick')
circle.selection_glyph = selected_circle
circle.nonselection_glyph = nonselected_circle
xaxis = LogAxis(axis_label = x_axis.value,
               ticker = LogTicker(num_minor_ticks = 10),
               axis_label_text_font_size = '12pt')
yaxis = LogAxis(axis_label = y_axis.value,
                ticker = LogTicker(num_minor_ticks = 10),
                axis_label_text_font_size = '12pt')
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
legend.circle(x0,
              y0,
              size = 15,
              line_color = "black",
              fill_color = "white",
              line_width = 1)
legend.circle(x0,
              y0 - 0.5,
              size = 15,
              line_color = "red" ,
              fill_color = "white",
              line_width = 1)
legend.circle(x0,
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
    j_mag_legend.circle(x0,
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

j_mag_legend.circle(x0,
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
                                            low = min_val(tbl['pl_eqt'].data),
                                            high = max_val(tbl['pl_eqt'].data)),
                     ticker = LogTicker(),
                     label_standoff = 10,
                     border_line_color = None,
                     location = (5, 0),
                     title = 'Teq (K)',
                     title_standoff = 10,
                     title_text_align = "center",
                     title_text_font_style = 'italic')
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
                 width = 300)
layout = row(widgets, p, legends)
layout = column(figure_title, layout)

# Save as .html
output_file("iExoView.html", title = 'iExoView')
show(layout)
print("Done")

##### Insert Google Analytics #####
from bs4 import BeautifulSoup

with open("iExoView.html", "r") as f:
    contents = f.read()
    soup = BeautifulSoup(contents, 'html.parser')

google_str = BeautifulSoup("""
<!-- Global site tag (gtag.js) - Google Analytics -->
 <script async src="https://www.googletagmanager.com/gtag/js?id=UA-166761240-1"></script>
 <script>
   window.dataLayer = window.dataLayer || [];
   function gtag(){dataLayer.push(arguments);}
   gtag('js', new Date());

   gtag('config', 'UA-166761240-1');
 </script>
""", 'html.parser')
soup.head.insert(0, google_str)

with open("iExoView.html", "w") as file:
    file.write(str(soup))

# import os
# os.system('rsync -avz -e "ssh -oport=5822" iExoView.html charles@venus.astro.umontreal.ca:/home/charles/www/')
