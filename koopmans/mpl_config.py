import matplotlib
import os

if 'DISPLAY' not in os.environ:
    matplotlib.use('Agg')
