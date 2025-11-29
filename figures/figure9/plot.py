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


fig = plt.figure(figsize=(7, 7))
gs = GridSpec(2, 1, height_ratios=[1, 0.4])
fig.subplots_adjust(hspace=0.0)

sub0 = fig.add_subplot(gs[0])
sub0.set_title('Coherent cross section for different quarks')
sub0.set_ylabel(
    r'${\rm d}\sigma_{\rm coh}/{\rm d}t \times (\Gamma_{{\rm J}/\Psi} / \Gamma_{\rm V}) \times (m_{\rm q}/m_{\rm c})^5 ~[{\rm nb~GeV}^{-2}]$')

sub0.set_yscale('log')
sub0.minorticks_on()

sub0.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
sub0.grid(which='major', linestyle='dashed', linewidth='0.5')
plt.setp(sub0.get_xticklabels(), visible=False)

sub0.set_xlim(0.0, 2.0)
sub0.set_ylim(1.0e-2, 1.0e4)

SCHEME = YBY2

CHARM_COLOR = SCHEME[8]
BOTTOM_COLOR = SCHEME[1]
DIDE_COLOR = '#aaaaaa'

labels = [
    Line2D([0], [0], color=DIDE_COLOR, linestyle=':', label='Dilute'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='-', label='Dense'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='none',
           marker='none', alpha=0.0, ms=7, label=''),
    Line2D([0], [0], color=CHARM_COLOR, linestyle='-', lw=5, label='Charm'),
    Line2D([0], [0], color=BOTTOM_COLOR, linestyle='-', lw=5, label='Bottom'),
]
sub0.legend(handles=labels)

sub0.plot(t, C_FACTOR*sdic10['Co'],
          label='Dilute, charm', color=SCHEME[8], linestyle=':')
sub0.fill_between(t, C_FACTOR*err(sdic10['Co']-sdic10['CoErr']), C_FACTOR*err(
    sdic10['Co']+sdic10['CoErr']), color=SCHEME[8], alpha=0.4)

sub0.plot(t, C_FACTOR*sdec10['Co'], label='Dense, charm', color=SCHEME[8])
sub0.fill_between(t, C_FACTOR*err(sdec10['Co']-sdec10['CoErr']), C_FACTOR*err(
    sdec10['Co']+sdec10['CoErr']), color=SCHEME[8], alpha=0.4)

sub0.plot(t, B_FACTOR*sdib10['Co'],
          label='Dilute, bottom', color=SCHEME[1], linestyle=':')
sub0.fill_between(t, B_FACTOR*err(sdib10['Co']-sdib10['CoErr']), B_FACTOR*err(
    sdib10['Co']+sdib10['CoErr']), color=SCHEME[1], alpha=0.4)

sub0.plot(t, B_FACTOR*sdeb10['Co'], label='Dense, bottom', color=SCHEME[1])
sub0.fill_between(t, B_FACTOR*err(sdeb10['Co']-sdeb10['CoErr']), B_FACTOR*err(
    sdeb10['Co']+sdeb10['CoErr']), color=SCHEME[1], alpha=0.4)

sub0ratio = fig.add_subplot(gs[1], sharex=sub0)
sub0ratio.set_xlabel(r'$|t|~[{\rm GeV}^2]$')
sub0ratio.set_ylabel('Dense/Dilute')
sub0ratio.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
sub0ratio.grid(which='major', linestyle='dashed', linewidth='0.5')


sub0ratio.set_ylim(0.5, 0.999)

labels = [
    Line2D([0], [0], color=DIDE_COLOR, linestyle='--', label='Dense/Dilute')
]
sub0ratio.legend(handles=labels, loc='lower left')

sub0ratio.plot(t, 1.0/(sdic10['Co']/sdec10['Co']),
               color=CHARM_COLOR, linestyle='--')
sub0ratio.plot(t, 1.0/(sdib10['Co']/sdeb10['Co']),
               color=BOTTOM_COLOR, linestyle='--')

fig.savefig('figure9.pdf', bbox_inches='tight')
