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

G_DATA = pd.read_csv('data/G_m_022.dat', sep=' ')

fig, sub0 = plt.subplots(1, 1, figsize=(7, 5))

sub0.set_title('Dipole correlation for a set of dipole configurations')
sub0.set_ylabel(r'$G_{\bf x y}$')
sub0.set_xlabel(r'$|{\bf b}|~[{\rm GeV}^{-1}]$')
sub0.grid(which='major', linestyle='dashed', linewidth=0.5)

PARAMS_TEXT = r'$b_1=0 b_2=b r_1=r\cos(\varphi) r_2=r\sin(\varphi) r=1$'

fig.text(0.6, 0.6, r'$b_1=0$'
         '\n'
         r'$b_2=|{\bf b}|$'
         '\n'
         r'$|{\bf r}|=1{\rm GeV}^{-1}$'
         '\n'
         r'$r_1=|{\bf r}|\cos(\varphi)$'
         '\n'
         r'$r_2=|{\bf r}|\sin(\varphi)$'
         '\n'
         r'${\bf x}={\bf b} + {\bf r}/2$'
         '\n'
         r'${\bf y}={\bf b} - {\bf r}/2$')

# plt.xscale('log')

sub0.set_xlim(0.0, 5.0)
sub0.set_ylim(0.0, -0.07)

# sub0.set_yscale('log')
# sub0.set_ylim(1.0e-6, 1.0e-1)

for i in range(5):
    sub0.plot(G_DATA['x'], G_DATA['angle'+str(i)+'integration'],
              color=SCHEME[i*2+3],  label=r'$\varphi=$'+str(i)+r'$/8\pi$')

plt.legend()

fig.savefig("figure13.pdf", bbox_inches='tight')
