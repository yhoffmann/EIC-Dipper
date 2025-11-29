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


fig0 = plt.figure(figsize=(7, 7))

gs0 = GridSpec(2, 1, height_ratios=[1, 0.4])

fig0.subplots_adjust(hspace=0.0)

sub0 = fig0.add_subplot(gs0[0])

sub0.set_title('Incoherent cross section with different ' + r'$g^2\mu_0^2$')

sub0.set_ylabel(r'${\rm d}\sigma_{\rm incoh}/{\rm d}t~[{\rm nb~GeV}^{-2}]$')

sub0.set_yscale('log')
sub0.minorticks_on()

sub0.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
sub0.grid(which='major', linestyle='dashed', linewidth='0.5')
plt.setp(sub0.get_xticklabels(), visible=False)

sub0.set_xlim(0.0, 8.0)
sub0.set_ylim(3.0e-2, 1.0e2)

SCHEME = YBY2
COLOR1 = SCHEME[11]
COLOR2 = SCHEME[8]
COLOR3 = SCHEME[1]
DIDE_COLOR = '#aaaaaa'

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

sub0.plot(t, 4.0*sdic10['Inco'], label='Dilute, ' +
          r'$g^2\mu_0^2+100\%$', color=SCHEME[11], linestyle=':')
sub0.fill_between(t, 4.0*err(sdic10['Inco']-sdic10['IncoErr']), 4.0*err(
    sdic10['Inco']+sdic10['IncoErr']), color=SCHEME[11], alpha=0.4)

sub0.plot(t, sdec20['Inco'], label='Dense, ' +
          r'$g^2\mu_0^2+100\%$', color=SCHEME[11])
sub0.fill_between(t, err(sdec20['Inco']-sdec20['IncoErr']), err(
    sdec20['Inco']+sdec20['IncoErr']), color=SCHEME[11], alpha=0.4)


sub0.plot(t, sdic10['Inco'], label='Dilute, ' +
          r'$g^2\mu_0^2$', color=SCHEME[8], linestyle=':')
sub0.fill_between(t, err(sdic10['Inco']-sdic10['IncoErr']),
                  err(sdic10['Inco']+sdic10['IncoErr']), color=SCHEME[8], alpha=0.4)

sub0.plot(t, sdec10['Inco'], label='Dense, ' +
          r'$g^2\mu_0^2$', color=SCHEME[8])
sub0.fill_between(t, err(sdec10['Inco']-sdec10['IncoErr']),
                  err(sdec10['Inco']+sdec10['IncoErr']), color=SCHEME[8], alpha=0.4)


sub0.plot(t, 0.25*sdic10['Inco'], label='Dilute, ' +
          r'$g^2\mu_0^2-50\%$', color=SCHEME[1], linestyle=':')
sub0.fill_between(t, 0.25*err(sdic10['Inco']-sdic10['IncoErr']), 0.25*err(
    sdic10['Inco']+sdic10['IncoErr']), color=SCHEME[1], alpha=0.4)

sub0.plot(t, sdec05['Inco'], label='Dense, ' +
          r'$g^2\mu_0^2-50\%$', color=SCHEME[1])
sub0.fill_between(t, err(sdec05['Inco']-sdec05['IncoErr']),
                  err(sdec05['Inco']+sdec05['IncoErr']), color=SCHEME[1], alpha=0.4)


sub0ratio = fig0.add_subplot(gs0[1], sharex=sub0)

sub0ratio.set_xlabel(r'$|t|~[{\rm GeV}^2]$')
sub0ratio.set_ylabel('Dense/Dilute')

sub0ratio.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
sub0ratio.grid(which='major', linestyle='dashed', linewidth='0.5')

sub0ratio.set_ylim(0.5, 0.999)

labels = [
    Line2D([0], [0], color=DIDE_COLOR, linestyle='--', label='Dense/Dilute'),
]
sub0ratio.legend(handles=labels)

sub0ratio.plot(t, 1.0/(4.0*sdic10['Inco']/sdec20['Inco']),
               label=r'$g^2\mu_0^2$ + 100%', color=SCHEME[11], ls='--')
sub0ratio.plot(t, 1.0/(sdic10['Inco']/sdec10['Inco']),
               label=r'$g^2\mu_0^2$', color=SCHEME[8], ls='--')
sub0ratio.plot(t, 1.0/(0.25*sdic10['Inco']/sdec05['Inco']),
               label=r'$g^2\mu_0^2$ - 50%', color=SCHEME[1], ls='--')

fig0.savefig('figure7.pdf', bbox_inches='tight')
