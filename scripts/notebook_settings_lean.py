#Initial configuration, probably overkill in imports.
import sys, os, re
import numpy as np
import pandas as pd
import gc

from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina', 'png')
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
#from horizonplot import horizonplot
import seaborn as sns
sns.set()
sns.set_theme()
sns.set_style("white")
sns.set_context("notebook")