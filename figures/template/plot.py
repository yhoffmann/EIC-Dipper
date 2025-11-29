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

sdib10 = pd.read_csv('data/sdib10.csv')
sdeb10 = pd.read_csv('data/sdeb10.csv')

sdec10 = pd.read_csv('data/sdec10.csv')
sdic10 = pd.read_csv('data/sdic10.csv')

t = sdic10['t']

DEMIRCI = pd.read_csv('data/dsdt_demirci.dat', sep=' ')

H1_co = pd.read_csv('data/dsigmadt_H1_co.dat', sep=' ')
H1_inco = pd.read_csv('data/dsigmadt_H1_inco.dat', sep=' ')
H1_inco_ht = pd.read_csv('data/dsigmadt_H1_inco_high_t.dat', sep=' ')


def err(x):
    return np.abs(x)


M_C = 1.275
M_B = 4.18
