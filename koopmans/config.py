"""

Contains the global variables for python_KI

Written by Edward Linscott Nov 2020

"""

# Contain all globals in an init function so that globals don't get reset upon every "import config"


def init(**kwargs):

    # if False, we try and pick up from a calculation that is partially complete
    global from_scratch

    for k, v in kwargs.items():
        if k == 'from_scratch':
            from_scratch = v
        else:
            raise ValueError('Could not set the unrecognised global {k}')
