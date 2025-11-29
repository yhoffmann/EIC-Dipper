import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colormaps
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.patches import Circle
from matplotlib.gridspec import GridSpec
import numpy as np
from math import *
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import scipy.special as sfs
from scipy.optimize import curve_fit
from matplotlib.ticker import (
    MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import pandas as pd
from matplotlib.lines import Line2D
import os

YBY2 = ["#00274cff", "#132472ff", "#261895ff", "#4d1aa8ff", "#711bbbff", "#a71cc9ff", "#d81db6ff",
        "#e22277ff", "#e62d46ff", "#e74d37ff", "#eb6e3dff", "#ec9646ff", "#edbd4fff", "#eed558ff", "#f4e471ff"]
scheme = YBY2

EXTDATA_COLOR = '#333333'

sdec10 = pd.read_csv('data/sdec10.csv')
sdic10 = pd.read_csv('data/sdic10.csv')

t = sdic10['t']

DEMIRCI = pd.read_csv('data/dsdt_demirci.dat', sep=' ')

H1_co = pd.read_csv('data/dsigmadt_H1_co.dat', sep=' ')
H1_inco = pd.read_csv('data/dsigmadt_H1_inco.dat', sep=' ')
H1_inco_ht = pd.read_csv('data/dsigmadt_H1_inco_high_t.dat', sep=' ')


def err(x):
    return np.abs(x)


fig, sub0 = plt.subplots(1, 1, figsize=(7, 5))

sub0.set_title('Incoherent cross section')

sub0.set_xlabel(r'$|t|~[{\rm GeV}^2]$')
sub0.set_ylabel(r'${\rm d}\sigma/{\rm d}t~[{\rm nb~GeV}^{-2}]$')

plt.yscale('log')
sub0.minorticks_on()

sub0.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
sub0.grid(which='major', linestyle='dashed', linewidth='0.5')

sub0.set_xlim(0.0, 8.0)
sub0.set_ylim(2.0e-1, 1.0e2)

SCHEME = YBY2

INCO_COLOR = SCHEME[11]
COLOR_COLOR = SCHEME[8]
HOTSPOT_COLOR = SCHEME[1]
DIDE_COLOR = '#aaaaaa'

labels = [
    Line2D([0], [0], color=DIDE_COLOR, linestyle=':', label='Dilute'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='-', label='Dense'),
    Line2D([0], [0], label='Analytical', ls='none',
           marker='+', ms=7, color=DIDE_COLOR),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='none',
           marker='none', alpha=0.0, ms=7, label=''),

    Line2D([0], [0], color=INCO_COLOR, linestyle='-',
           linewidth=5, label=r'$3.0 \times$'+'Incoherent'),
    Line2D([0], [0], color=COLOR_COLOR, linestyle='-',
           linewidth=5, label=r'$3.0 \times$'+'Color fluc.'),
    Line2D([0], [0], color=HOTSPOT_COLOR, linestyle='-',
           linewidth=5, label=r'$3.0 \times$'+'Hotspot fluc.'),
    Line2D([0], [0], color=EXTDATA_COLOR, linestyle='none', marker='s', fillstyle='none',
           label='H1, low '+r'$|t|$'),
    Line2D([0], [0], color=EXTDATA_COLOR, linestyle='none', marker='o', fillstyle='none',
           label='H1, high '+r'$|t|$'),
]
sub0.legend(handles=labels)

SAMP_DATA_FACTOR = 3.0
sub0.plot(t, SAMP_DATA_FACTOR*sdic10['Inco'],
          label='Dilute', linestyle=':', color=INCO_COLOR)
sub0.fill_between(t, SAMP_DATA_FACTOR*err(sdic10['Inco']-sdic10['IncoErr']), SAMP_DATA_FACTOR*err(
    sdic10['Inco']+sdic10['IncoErr']), color=INCO_COLOR, alpha=0.4)

sub0.plot(t, SAMP_DATA_FACTOR*sdec10['Inco'], label='Dense', color=INCO_COLOR)
sub0.fill_between(t, SAMP_DATA_FACTOR*err(sdec10['Inco']-sdec10['IncoErr']), SAMP_DATA_FACTOR*err(
    sdec10['Inco']+sdec10['IncoErr']), color=INCO_COLOR, alpha=0.4)

sub0.plot(t, SAMP_DATA_FACTOR *
          sdic10['Color'], label='Dilute, color fluc.', linestyle=':', color=COLOR_COLOR)
sub0.fill_between(t, SAMP_DATA_FACTOR*err(sdic10['Color']-sdic10['ColorErr']), SAMP_DATA_FACTOR*err(
    sdic10['Color']+sdic10['ColorErr']), color=COLOR_COLOR, alpha=0.4)

sub0.plot(t, SAMP_DATA_FACTOR*sdec10['Color'],
          label='Dense, color fluc.', color=COLOR_COLOR)
sub0.fill_between(t, SAMP_DATA_FACTOR*err(sdec10['Color']-sdec10['ColorErr']), SAMP_DATA_FACTOR*err(
    sdec10['Color']+sdec10['ColorErr']), color=COLOR_COLOR, alpha=0.4)

sub0.plot(t, SAMP_DATA_FACTOR*sdic10['Hotspot'],
          label='Dilute, hotspot fluc.', linestyle=':', color=HOTSPOT_COLOR)
sub0.fill_between(t, SAMP_DATA_FACTOR*err(sdic10['Hotspot']-sdic10['HotspotErr']), SAMP_DATA_FACTOR*err(
    sdic10['Hotspot']+sdic10['HotspotErr']), color=HOTSPOT_COLOR, alpha=0.4)

sub0.plot(t, SAMP_DATA_FACTOR*sdec10['Hotspot'],
          label='Dense, hotspot fluc.', color=HOTSPOT_COLOR)
sub0.fill_between(t, SAMP_DATA_FACTOR*err(sdec10['Hotspot']-sdec10['HotspotErr']), SAMP_DATA_FACTOR*err(
    sdec10['Hotspot']+sdec10['HotspotErr']), color=HOTSPOT_COLOR, alpha=0.4)

sub0.plot(DEMIRCI['t'], SAMP_DATA_FACTOR*DEMIRCI['Inco'],
          label='Analytical', ls='none', marker='+', ms=7, color=INCO_COLOR)
sub0.plot(DEMIRCI['t'], SAMP_DATA_FACTOR*DEMIRCI['Color'],
          label='Analytical, color fluc.', ls='none', marker='+', ms=7, color=COLOR_COLOR)
sub0.plot(DEMIRCI['t'], SAMP_DATA_FACTOR*DEMIRCI['Hotspot'],
          label='Analytical, hotspot fluc.', ls='none', marker='+', ms=7, color=HOTSPOT_COLOR)

H1_inco = pd.read_csv('data/dsigmadt_H1_inco.dat', sep=' ')
H1_inco_ht = pd.read_csv('data/dsigmadt_H1_inco_high_t.dat', sep=' ')
H1_FACTOR = 1.0
H1_inco['Inco'] = H1_inco['Inco']*H1_FACTOR
H1_inco['IncoErrStat'] = H1_inco['IncoErrStat']*H1_FACTOR
H1_inco['IncoErrSys'] = H1_inco['IncoErrSys']*H1_FACTOR
H1_inco_ht['Inco'] = H1_inco_ht['Inco']*H1_FACTOR
H1_inco_ht['IncoErrStat'] = H1_inco_ht['IncoErrStat']*H1_FACTOR
H1_inco_ht['IncoErrSys'] = H1_inco_ht['IncoErrSys']*H1_FACTOR

sub0.errorbar(H1_inco['t'], H1_inco['Inco'], np.sqrt(np.square(H1_inco['IncoErrStat']) + np.square(H1_inco['IncoErrSys'])),
              marker='s', fillstyle='none', color=EXTDATA_COLOR, linestyle='none', label='H1'+r'$\cdot \frac{1}{3}$')
sub0.errorbar(H1_inco_ht['t'], H1_inco_ht['Inco'], np.sqrt(np.square(H1_inco_ht['IncoErrStat']) + np.square(H1_inco_ht['IncoErrSys'])),
              marker='o', fillstyle='none', color=EXTDATA_COLOR, linestyle='none', label='H1'+r'$\cdot \frac{1}{3}$'+', high t')

fig.savefig('figure4.pdf', bbox_inches='tight')
