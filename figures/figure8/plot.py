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
SCHEME = YBY2

EXTDATA_COLOR = '#333333'

sdic10 = pd.read_csv('data/sdic10.csv')

sgc = pd.read_csv('data/sgc.csv')


def err(x):
    return np.abs(x)


g2mu02 = np.sqrt(sgc['t'])
INDEX = 13
CO = sdic10['Co'][INDEX]
INCO = sdic10['Inco'][INDEX]
COLOR = sdic10['Color'][INDEX]
HOTSPOT = sdic10['Hotspot'][INDEX]
DIDE_COLOR = '#aaaaaa'

fig = plt.figure(figsize=(7, 7))

fig.text(0.35, 0.815, r'$|t| = 1.0{\rm GeV}^2$' +
         '\n'+r'$(g^2\mu_0^2)_0 = \sqrt{43.22}$')

gs = GridSpec(2, 1, height_ratios=[1, 0.4])
fig.subplots_adjust(hspace=0.0)
sub0 = fig.add_subplot(gs[0])
sub0.set_title(
    r'Normalized cross sections at single momentum transfer for different $g^2\mu_0^2$')

sub0.grid(which='minor', axis='both', linestyle='dashed', linewidth='0.2')
sub0.grid(which='major', linestyle='dashed', linewidth='0.5')
plt.setp(sub0.get_xticklabels(), visible=False)
sub0.set_ylabel(
    r'${\rm d}\sigma/{\rm d}t / (g^2\mu_0^2) ~[{\rm nb~ GeV}^{-2}]$')

sub0.set_xlim(0.0, 2.5)
sub0.set_ylim(0.0, 15.0)

CO_COLOR = SCHEME[1]
INCO_COLOR = SCHEME[8]

labels = [
    Line2D([0], [0], color=DIDE_COLOR, linestyle=':', label='Dilute'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='-', label='Dense'),
    Line2D([0], [0], color=DIDE_COLOR, linestyle='none',
           marker='none', alpha=0.0, ms=7, label=''),
    Line2D([0], [0], color=CO_COLOR, linestyle='-', lw=5, label='Coherent'),
    Line2D([0], [0], color=INCO_COLOR,
           linestyle='-', lw=5, label='Incoherent'),
]
sub0.legend(handles=labels)

sub0.plot(g2mu02, g2mu02*CO, color=CO_COLOR, linestyle=':', label='Dilute')
sub0.plot(g2mu02, sgc['Co']/g2mu02, color=CO_COLOR, label='Dense')

sub0.plot(g2mu02, g2mu02*INCO, ls=':', color=INCO_COLOR)
sub0.plot(g2mu02, sgc['Inco']/g2mu02, color=INCO_COLOR)

sub0ratio = fig.add_subplot(gs[1], sharex=sub0)
sub0ratio.grid(which='major', linestyle='dashed', linewidth='0.5')
sub0ratio.set_xlabel(r'$g^2 \mu_0^2 / (g^2\mu_0^2)_0$')
sub0ratio.set_ylabel('Dense/Dilute')

sub0ratio.set_ylim(0.5, 0.999)

labels = [
    Line2D([0], [0], color=DIDE_COLOR, linestyle='--', label='Dense/Dilute')
]
sub0ratio.legend(handles=labels)

sub0ratio.plot(g2mu02, 1.0/(g2mu02*g2mu02*CO /
               sgc['Co']), color=CO_COLOR, linestyle='--')
sub0ratio.plot(g2mu02, 1.0/(g2mu02*g2mu02*INCO /
               sgc['Inco']), color=INCO_COLOR, linestyle='--')

fig.savefig('figure8.pdf', bbox_inches='tight')
