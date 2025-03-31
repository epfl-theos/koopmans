"""Update the version in the CFF file."""
import sys

import yaml

if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    import importlib_metadata as metadata


# Load the CITATION.cff file
with open('CITATION.cff', 'r') as fd:
    data = yaml.load(fd, Loader=yaml.Loader)

# Update version
data['version'] = metadata.version('koopmans')

# Overwrite the CITATION.cff file
with open('CITATION.cff', 'w') as fd:
    yaml.dump(data, fd, allow_unicode=True)
