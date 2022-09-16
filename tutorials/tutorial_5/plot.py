from pathlib import Path
import matplotlib

from koopmans import io

matplotlib.use('Agg')  # nopep8
import matplotlib.pyplot as plt  # nopep8

if __name__ == '__main__':
    # Read in the workflow
    wf = io.read(Path('tutorial_5a') / 'h2o_conv.kwf')
