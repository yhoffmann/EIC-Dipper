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


fig = plt.figure(figsize=(7, 7))
gs = GridSpec(2, 1, height_ratios=[1, 0.4])
plt.subplots_adjust(hspace=0.00)

sub0 = fig.add_subplot(gs[0])

sub0.set_title('Coherent cross section with different ' + r'$g^2\mu_0^2$')

sub0.set_ylabel(r'${\rm d}\sigma_{\rm coh}/{\rm d}t~[{\rm nb~GeV}^{-2}]$')

plt.yscale('log')
sub0.minorticks_on()

sub0.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
sub0.grid(which='major', linestyle='dashed', linewidth='0.5')

sub0.set_xlim(0.0, 2.0)
sub0.set_ylim(1.0e-2, 3.0e3)

SCHEME = YBY2
COLOR1 = SCHEME[11]
COLOR2 = SCHEME[8]
COLOR3 = SCHEME[1]
DIDE_COLOR = '#aaaaaa'

INCO_COLOR = SCHEME[11]
COLOR_COLOR = SCHEME[8]
HOTSPOT_COLOR = SCHEME[1]

labels = [
    Line2D([0], [0], color=DIDE_COLOR, linestyle=':', label='Dilute'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='-', label='Dense'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='none',
           marker='none', alpha=0.0, ms=7, label=''),
    Line2D([0], [0], color=COLOR1, lw=5,
           label=r'$g^2\mu_0^2 = 2.0 \times (g^2\mu_0^2)_0$'),
    Line2D([0], [0], color=COLOR2, lw=5,
           label=r'$g^2\mu_0^2 = 1.0 \times (g^2\mu_0^2)_0$'),
    Line2D([0], [0], color=COLOR3, lw=5,
           label=r'$g^2\mu_0^2 = 0.5 \times (g^2\mu_0^2)_0$'),
]
sub0.legend(handles=labels)

sub0.plot(t, 4.0*sdic10['Co'], label='Dilute, ' +
          r'$g^2\mu_0^2+100\%$', color=SCHEME[11], linestyle=':')
sub0.fill_between(t, 4.0*err(sdic10['Co']-sdic10['CoErr']), 4.0*err(
    sdic10['Co']+sdic10['CoErr']), color=SCHEME[11], alpha=0.4)

sub0.plot(t, sdec20['Co'], label='Dense, ' +
          r'$g^2\mu_0^2+100\%$', color=SCHEME[11])
sub0.fill_between(t, err(sdec20['Co']-sdec20['CoErr']),
                  err(sdec20['Co']+sdec20['CoErr']), color=SCHEME[11], alpha=0.4)


sub0.plot(t, sdic10['Co'], label='Dilute, ' +
          r'$g^2\mu_0^2$', color=SCHEME[8], linestyle=':')
sub0.fill_between(t, err(sdic10['Co']-sdic10['CoErr']),
                  err(sdic10['Co']+sdic10['CoErr']), color=SCHEME[8], alpha=0.4)

sub0.plot(t, sdec10['Co'], label='Dense, ' + r'$g^2\mu_0^2$', color=SCHEME[8])
sub0.fill_between(t, err(sdec10['Co']-sdec10['CoErr']),
                  err(sdec10['Co']+sdec10['CoErr']), color=SCHEME[8], alpha=0.4)


sub0.plot(t, 0.25*sdic10['Co'], label='Dilute, ' +
          r'$g^2\mu_0^2-50\%$', color=SCHEME[1], linestyle=':')
sub0.fill_between(t, 0.25*err(sdic10['Co']-sdic10['CoErr']), 0.25*err(
    sdic10['Co']+sdic10['CoErr']), color=SCHEME[1], alpha=0.4)

sub0.plot(t, sdec05['Co'], label='Dense, ' +
          r'$g^2\mu_0^2-50\%$', color=SCHEME[1])
sub0.fill_between(t, err(sdec05['Co']-sdec05['CoErr']),
                  err(sdec05['Co']+sdec05['CoErr']), color=SCHEME[1], alpha=0.4)

ratio = fig.add_subplot(gs[1], sharex=sub0)


ratio.set_xlabel(r'$|t|~[{\rm GeV}^2]$')
ratio.set_ylabel('Dense/Dilute')

ratio.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
ratio.grid(which='major', linestyle='dashed', linewidth='0.5')

ratio.set_ylim(0.0, 0.999)

labels = [
    Line2D([0], [0], color=DIDE_COLOR, linestyle='--', label='Dense/Dilute'),
]
ratio.legend(handles=labels)

ratio.plot(t, 1.0/(4.0*sdic10['Co']/sdec20['Co']), ls='--', color=SCHEME[11])
ratio.plot(t, 1.0/(sdic10['Co']/sdec10['Co']), ls='--', color=SCHEME[8])
ratio.plot(t, 1.0/(0.25*sdic10['Co']/sdec05['Co']), ls='--', color=SCHEME[1])

plt.setp(sub0.get_xticklabels(), visible=False)

fig.savefig('figure5.pdf', bbox_inches='tight')
