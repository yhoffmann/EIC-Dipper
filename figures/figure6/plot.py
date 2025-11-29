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

sdec05 = pd.read_csv('data/sdec05.csv')
sdec20 = pd.read_csv('data/sdec20.csv')

t = sdic10['t']


def err(x):
    return np.abs(x)


fig1 = plt.figure(figsize=(7, 7))
gs1 = GridSpec(2, 1, height_ratios=[1, 0.4])
fig1.subplots_adjust(hspace=0.0)

sub1 = fig1.add_subplot(gs1[0])

sub1.set_title(
    'Color and hotspot fluctuations with different ' + r'$g^2\mu_0^2$')

# sub1.set_xlabel(r'$|t|~[{\rm GeV}^2]$')
sub1.set_ylabel(r'${\rm d}\sigma/{\rm d}t~[{\rm nb~GeV}^{-2}]$')

sub1.set_yscale('log')
sub1.minorticks_on()

sub1.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
sub1.grid(which='major', linestyle='dashed', linewidth='0.5')
plt.setp(sub1.get_xticklabels(), visible=False)

sub1.set_xlim(0.0, 8.0)
sub1.set_ylim(3.0e-2, 1.0e2)

SCHEME = YBY2
COLOR1 = SCHEME[11]
COLOR2 = SCHEME[8]
COLOR3 = SCHEME[1]
DIDE_COLOR = '#aaaaaa'

labels = [
    Line2D([0], [0], color=DIDE_COLOR, linestyle=':',
           label='Dilute, color fluc.'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='-',
           label='Dense, color fluc.'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle=(
        0, (3, 2, 1, 2, 1, 2)), label='Dilute, hotspot fluc.'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle=(
        0, (10, 4)), label='Dense, hotspot fluc.'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='none',
           marker='none', alpha=0.0, ms=7, label=''),
    Line2D([0], [0], color=COLOR1, lw=5,
           label=r'$g^2\mu_0^2 = 2.0 \times (g^2\mu_0^2)_0$'),
    Line2D([0], [0], color=COLOR2, lw=5,
           label=r'$g^2\mu_0^2 = 1.0 \times (g^2\mu_0^2)_0$'),
    Line2D([0], [0], color=COLOR3, lw=5,
           label=r'$g^2\mu_0^2 = 0.5 \times (g^2\mu_0^2)_0$'),
]
sub1.legend(handles=labels)

sub1.plot(t, 4.0*sdic10['Color'], label='Dilute, color fluc., ' +
          r'$g^2\mu_0^2+100\%$', color=COLOR1, linestyle=':')
sub1.fill_between(t, 4.0*err(sdic10['Color']-sdic10['ColorErr']), 4.0*err(
    sdic10['Color']+sdic10['ColorErr']), color=COLOR1, alpha=0.4)

sub1.plot(t, sdec20['Color'], label='Dense, color fluc., ' +
          r'$g^2\mu_0^2+100\%$', color=COLOR1)
sub1.fill_between(t, err(sdec20['Color']-sdec20['ColorErr']),
                  err(sdec20['Color']+sdec20['ColorErr']), color=COLOR1, alpha=0.4)

sub1.plot(t, 4.0*sdic10['Hotspot'], label='Dilute, hotspot fluc., ' +
          r'$g^2\mu_0^2+100\%$', color=COLOR1, linestyle=(0, (3, 2, 1, 2, 1, 2)))
sub1.fill_between(t, 4.0*err(sdic10['Hotspot']-sdic10['HotspotErr']), 4.0*err(
    sdic10['Hotspot']+sdic10['HotspotErr']), color=COLOR1, alpha=0.4)

sub1.plot(t, sdec20['Hotspot'], label='Dense, hotspot fluc., ' +
          r'$g^2\mu_0^2+100\%$', color=COLOR1, ls=(0, (10, 4)))
sub1.fill_between(t, err(sdec20['Hotspot']-sdec20['HotspotErr']), err(
    sdec20['Hotspot']+sdec20['HotspotErr']), color=COLOR1, alpha=0.4)


sub1.plot(t, sdic10['Color'], label='Dilute, color fluc., ' +
          r'$g^2\mu_0^2$', color=COLOR2, linestyle=':')
sub1.fill_between(t, err(sdic10['Color']-sdic10['ColorErr']),
                  err(sdic10['Color']+sdic10['ColorErr']), color=COLOR2, alpha=0.4)

sub1.plot(t, sdec10['Color'], label='Dense, color fluc., ' +
          r'$g^2\mu_0^2$', color=COLOR2)
sub1.fill_between(t, err(sdec10['Color']-sdec10['ColorErr']),
                  err(sdec10['Color']+sdec10['ColorErr']), color=COLOR2, alpha=0.4)

sub1.plot(t, sdic10['Hotspot'], label='Dilute, hotspot fluc., ' +
          r'$g^2\mu_0^2$', color=COLOR2, linestyle=(0, (3, 2, 1, 2, 1, 2)))
sub1.fill_between(t, err(sdic10['Hotspot']-sdic10['HotspotErr']), err(
    sdic10['Hotspot']+sdic10['HotspotErr']), color=COLOR2, alpha=0.4)

sub1.plot(t, sdec10['Hotspot'], label='Dense, hotspot fluc., ' +
          r'$g^2\mu_0^2$', color=COLOR2, ls=(0, (10, 4)))
sub1.fill_between(t, err(sdec10['Hotspot']-sdec10['HotspotErr']), err(
    sdec10['Hotspot']+sdec10['HotspotErr']), color=COLOR2, alpha=0.4)


sub1.plot(t, 0.25*sdic10['Color'], label='Dilute, color fluc., ' +
          r'$g^2\mu_0^2-50\%$', color=COLOR3, linestyle=':')
sub1.fill_between(t, 0.25*err(sdic10['Color']-sdic10['ColorErr']), 0.25*err(
    sdic10['Color']+sdic10['ColorErr']), color=COLOR3, alpha=0.4)

sub1.plot(t, sdec05['Color'], label='Dense, color fluc., ' +
          r'$g^2\mu_0^2-50\%$', color=COLOR3)
sub1.fill_between(t, err(sdec05['Color']-sdec05['ColorErr']),
                  err(sdec05['Color']+sdec05['ColorErr']), color=COLOR3, alpha=0.4)

sub1.plot(t, 0.25*sdic10['Hotspot'], label='Dilute, hotspot fluc., ' +
          r'$g^2\mu_0^2-50\%$', color=COLOR3, linestyle=(0, (3, 2, 1, 2, 1, 2)))
sub1.fill_between(t, 0.25*err(sdic10['Hotspot']-sdic10['HotspotErr']), 0.25*err(
    sdic10['Hotspot']+sdic10['HotspotErr']), color=COLOR3, alpha=0.4)

sub1.plot(t, sdec05['Hotspot'], label='Dense, hotspot fluc., ' +
          r'$g^2\mu_0^2-50\%$', color=COLOR3, ls=(0, (10, 4)))
sub1.fill_between(t, err(sdec05['Hotspot']-sdec05['HotspotErr']), err(
    sdec05['Hotspot']+sdec05['HotspotErr']), color=COLOR3, alpha=0.4)


sub1ratio = fig1.add_subplot(gs1[1], sharex=sub1)

sub1ratio.set_xlabel(r'$|t|~[{\rm GeV}^2]$')
sub1ratio.set_ylabel('Dense/Dilute')

sub1ratio.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
sub1ratio.grid(which='major', linestyle='dashed', linewidth='0.5')

sub1ratio.set_ylim(0.5, 0.999)

labels = [
    Line2D([0], [0], color=DIDE_COLOR, linestyle='-.',
           label='Dense/Dilute, color fluc.'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='--',
           label='Dense/Dilute, hotspot fluc.'),
]
sub1ratio.legend(handles=labels, loc='lower right')

sub1ratio.plot(t, 1.0/(4.0*sdic10['Color']/sdec20['Color']),
               label=r'Color fluc., $g^2\mu_0^2$ + 100%', color=COLOR1, linestyle='dashdot')
sub1ratio.plot(t, 1.0/(4.0*sdic10['Hotspot']/sdec20['Hotspot']),
               label=r'Hotspot fluc., $g^2\mu_0^2$ + 100%', color=COLOR1, linestyle='--')


sub1ratio.plot(t, 1.0/(sdic10['Color']/sdec10['Color']),
               label=r'Color fluc., $g^2\mu_0^2$', color=COLOR2, linestyle='dashdot')
sub1ratio.plot(t, 1.0/(sdic10['Hotspot']/sdec10['Hotspot']),
               label=r'Hotspot fluc., $g^2\mu_0^2$', color=COLOR2, linestyle='--')

sub1ratio.plot(t, 1.0/(0.25*sdic10['Color']/sdec05['Color']),
               label=r'Color fluc., $g^2\mu_0^2$ - 50%', color=COLOR3, linestyle='dashdot')
sub1ratio.plot(t, 1.0/(0.25*sdic10['Hotspot']/sdec05['Hotspot']),
               label=r'Hotspot fluc., $g^2\mu_0^2$ - 50%', color=COLOR3, linestyle='--')

fig1.savefig('figure6.pdf', bbox_inches='tight')
