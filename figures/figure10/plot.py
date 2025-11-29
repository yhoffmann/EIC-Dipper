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

YBY2 = ["#00274cff", "#132472ff", "#261895ff", "#4d1aa8ff", "#711bbbff", "#a71cc9ff", "#d81db6ff",
        "#e22277ff", "#e62d46ff", "#e74d37ff", "#eb6e3dff", "#ec9646ff", "#edbd4fff", "#eed558ff", "#f4e471ff"]

EXTDATA_COLOR = '#333333'

sdib10 = pd.read_csv('data/sdib10.csv')
sdeb10 = pd.read_csv('data/sdeb10.csv')

sdec10 = pd.read_csv('data/sdec10.csv')
sdic10 = pd.read_csv('data/sdic10.csv')


t = sdic10['t']


def err(x):
    return np.abs(x)


M_C = 1.275
M_B = 4.18

C_FACTOR = 1.0
# factor of 4 comes from (e_c/e_b)^2 = ((2/3) / (1/3))^2 = 4
B_FACTOR = 1.0*np.power(M_B/M_C, 7)*4.0


fig1 = plt.figure(figsize=(7, 7))
gs1 = GridSpec(2, 1, height_ratios=[1, 0.4])
fig1.subplots_adjust(hspace=0.0)

sub1 = fig1.add_subplot(gs1[0])

sub1.set_title('Color and hotspot fluctuations for different quarks')

sub1.set_ylabel(
    r'${\rm d}\sigma/{\rm d}t \times (\Gamma_{{\rm J}/\Psi} / \Gamma_{\rm V}) \times (m_{\rm q}/m_{\rm c})^5 ~[{\rm nb~GeV}^{-2}]$')

sub1.set_yscale('log')
sub1.minorticks_on()

sub1.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
sub1.grid(which='major', linestyle='dashed', linewidth='0.5')
plt.setp(sub1.get_xticklabels(), visible=False)

sub1.set_xlim(0.0, 8.0)
sub1.set_ylim(1.0e-1, 2.0e2)

SCHEME = YBY2

CHARM_COLOR = SCHEME[8]
BOTTOM_COLOR = SCHEME[1]
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
    Line2D([0], [0], color=CHARM_COLOR, lw=5, label='Charm'),
    Line2D([0], [0], color=BOTTOM_COLOR, lw=5, label='Bottom'),
]
sub1.legend(handles=labels)

sub1.plot(t, C_FACTOR*sdic10['Color'],
          label='Dilute, color fluc., charm', color=CHARM_COLOR, linestyle=':')
sub1.fill_between(t, C_FACTOR*err(sdic10['Color']-sdic10['ColorErr']), C_FACTOR*err(
    sdic10['Color']+sdic10['ColorErr']), color=CHARM_COLOR, alpha=0.4)

sub1.plot(t, C_FACTOR*sdec10['Color'],
          label='Dense, color fluc., charm', color=CHARM_COLOR, ls='-')
sub1.fill_between(t, C_FACTOR*err(sdec10['Color']-sdec10['ColorErr']), C_FACTOR*err(
    sdec10['Color']+sdec10['ColorErr']), color=CHARM_COLOR, alpha=0.4)

sub1.plot(t, C_FACTOR*sdic10['Hotspot'], label='Dilute, hotspot fluc., charm',
          color=CHARM_COLOR, linestyle=(0, (3, 2, 1, 2, 1, 2)))
sub1.fill_between(t, C_FACTOR*err(sdic10['Hotspot']-sdic10['HotspotErr']), C_FACTOR*err(
    sdic10['Hotspot']+sdic10['HotspotErr']), color=CHARM_COLOR, alpha=0.4)

sub1.plot(t, C_FACTOR*sdec10['Hotspot'],
          label='Dense, hotspot fluc., charm', color=CHARM_COLOR, ls=(0, (10, 4)))
sub1.fill_between(t, C_FACTOR*err(sdec10['Hotspot']-sdec10['HotspotErr']), C_FACTOR*err(
    sdec10['Hotspot']+sdec10['HotspotErr']), color=CHARM_COLOR, alpha=0.4)


sub1.plot(t, B_FACTOR*sdib10['Color'],
          label='Dilute, color fluc., bottom', color=BOTTOM_COLOR, linestyle=':')
sub1.fill_between(t, B_FACTOR*err(sdib10['Color']-sdib10['ColorErr']), B_FACTOR*err(
    sdib10['Color']+sdib10['ColorErr']), color=BOTTOM_COLOR, alpha=0.4)

sub1.plot(t, B_FACTOR*sdeb10['Color'],
          label='Dense, color fluc., bottom', color=BOTTOM_COLOR, ls='-')
sub1.fill_between(t, B_FACTOR*err(sdeb10['Color']-sdeb10['ColorErr']), B_FACTOR*err(
    sdeb10['Color']+sdeb10['ColorErr']), color=BOTTOM_COLOR, alpha=0.4)

sub1.plot(t, B_FACTOR*sdib10['Hotspot'], label='Dilute, hotspot fluc., bottom',
          color=BOTTOM_COLOR, linestyle=(0, (3, 2, 1, 2, 1, 2)))
sub1.fill_between(t, B_FACTOR*err(sdib10['Hotspot']-sdib10['HotspotErr']), B_FACTOR*err(
    sdib10['Hotspot']+sdib10['HotspotErr']), color=BOTTOM_COLOR, alpha=0.4)

sub1.plot(t, B_FACTOR*sdeb10['Hotspot'],
          label='Dense, hotspot fluc., bottom', color=BOTTOM_COLOR, ls=(0, (10, 4)))
sub1.fill_between(t, B_FACTOR*err(sdeb10['Hotspot']-sdeb10['HotspotErr']), B_FACTOR*err(
    sdeb10['Hotspot']+sdeb10['HotspotErr']), color=BOTTOM_COLOR, alpha=0.4)


sub1ratio = fig1.add_subplot(gs1[1], sharex=sub1)
sub1ratio.set_xlabel(r'$|t|~[{\rm GeV}^2]$')
sub1ratio.set_ylabel('Dense/Dilute')

sub1ratio.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
sub1ratio.grid(which='major', linestyle='dashed', linewidth='0.5')

sub1ratio.set_ylim(0.65, 0.999)

labels = [
    Line2D([0], [0], color=DIDE_COLOR, linestyle='-.',
           label='Dense/Dilute, color fluc.'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='--',
           label='Dense/Dilute, hotspot fluc.'),
]
sub1ratio.legend(handles=labels, loc='lower right')

sub1ratio.plot(t, 1.0/(sdic10['Color']/sdec10['Color']),
               color=CHARM_COLOR, linestyle='dashdot')
sub1ratio.plot(t, 1.0/(sdib10['Color']/sdeb10['Color']),
               color=BOTTOM_COLOR, linestyle='dashdot')

sub1ratio.plot(
    t, 1.0/(sdic10['Hotspot']/sdec10['Hotspot']), color=CHARM_COLOR, linestyle='--')
sub1ratio.plot(t, 1.0/(sdib10['Hotspot']/sdeb10['Hotspot']),
               color=BOTTOM_COLOR, linestyle='--')

fig1.savefig('figure10.pdf', bbox_inches='tight')
