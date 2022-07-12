"""
Import this matplotlib config file prior to importing matplotlib.pyplot in order to appropriately configure the
matplotlib backend
"""

import os

import matplotlib

if 'DISPLAY' not in os.environ:
    matplotlib.use('Agg')
