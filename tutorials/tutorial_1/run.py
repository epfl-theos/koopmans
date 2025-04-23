"""Run the ozone workflow with a python script."""

from koopmans.io import read

wf = read('ozone.json')
wf.run()
