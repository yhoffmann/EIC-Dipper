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
scheme = YBY2

EXTDATA_COLOR = '#333333'

sdec10 = pd.read_csv('data/sdec10.csv')
sdic10 = pd.read_csv('data/sdic10.csv')

t = sdic10['t']

DEMIRCI = pd.read_csv('data/dsdt_demirci.dat', sep=' ')

H1_co = pd.read_csv('data/dsigmadt_H1_co.dat', sep=' ')


def err(x):
    return np.abs(x)


fig, sub0 = plt.subplots(1, 1, figsize=(7, 5))

sub0.set_title('Coherent cross section')

sub0.set_xlabel(r'$|t|~[{\rm GeV}^2]$')
sub0.set_ylabel(r'${\rm d}\sigma/{\rm d}t~[{\rm nb~GeV}^{-2}]$')

plt.yscale('log')
sub0.minorticks_on()

sub0.grid(which='minor', axis='both',
          linestyle='dashed', linewidth='0.3', alpha=0.6)
sub0.grid(which='major', linestyle='dashed', linewidth='0.5')

sub0.set_xlim(0.0, 2.0)
sub0.set_ylim(1.0e-2, 1.0e3)


SCHEME = YBY2

DI_COLOR = SCHEME[1]
DE_COLOR = SCHEME[1]

labels = [
    Line2D([0], [0], color=DI_COLOR, linestyle=':', label='Dilute'),
    Line2D([0], [0], color=DE_COLOR, linestyle='-', label='Dense'),
    Line2D([0], [0], label='Analytical', ls='none',
           marker='+', ms=7, color=DI_COLOR),
    Line2D([0], [0], color=EXTDATA_COLOR, linestyle='none',
           marker='s', fillstyle='none', label='H1'),
]

sub0.legend(handles=labels)

sub0.plot(sdic10['t'], sdic10['Co'], label='Dilute',
          color=DI_COLOR, linestyle=':')
sub0.fill_between(sdic10['t'], err(sdic10['Co']-sdic10['CoErr']),
                  err(sdic10['Co']+sdic10['CoErr']), color=DI_COLOR, alpha=0.4)

sub0.plot(sdec10['t'], sdec10['Co'], label='Dense', color=DE_COLOR)
sub0.fill_between(sdec10['t'], err(sdec10['Co']-sdec10['CoErr']),
                  err(sdec10['Co']+sdec10['CoErr']), color=DE_COLOR, alpha=0.4)

sub0.plot(DEMIRCI['t'], DEMIRCI['Co'], label='Analytical',
          ls='none', marker='+', ms=7, color=DI_COLOR)

sub0.errorbar(H1_co['t'], H1_co['Co'], np.sqrt(np.square(H1_co['CoErrStat']) + np.square(H1_co['CoErrSys'])),
              marker='s', fillstyle='none', color=EXTDATA_COLOR, linestyle='none', label='H1')

fig.savefig('figure3.pdf', bbox_inches='tight')
