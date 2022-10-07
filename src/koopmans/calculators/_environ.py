import copy
import os

from ._pw import PWCalculator

_default_settings = {
    'ENVIRON': {
        'verbose': 0,
        'environ_type': 'input',
        'environ_type': 'input',
        'env_surface_tension': 0,
        'env_pressure': 0},
    'BOUNDARY': {},
    'ELECTROSTATIC': {}}


class EnvironCalculator(PWCalculator):
    # Create an environ calculator that inherits from the vanilla pw.x calculator

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Ensure we're using an environ-enabled version of pw.x
        self.command.flags = '--environ'

        # Add dictionary of environ settings
        self.environ_settings = copy.deepcopy(_default_settings)

    def calculate(self):
        # Generic function for running a calculation
        self.write_environ_in()
        super().calculate()

    def set_environ_settings(self, settings, use_defaults=True):
        self.environ_settings = settings

        if use_defaults:
            # cycle through blocks in default settings
            for block_name, block in _default_settings.items():
                # if an entire block is missing, add it
                if block_name not in settings.keys():
                    self.environ_settings[block_name] = {}
                # if a particular keyword is missing, add it
                for key, value in block.items():
                    if key not in settings[block_name].keys():
                        self.environ_settings[block_name][key] = value

    def write_environ_in(self):
        # Write an environ.in file
        with open(f'{self.directory}/environ.in', 'w') as f:
            # cycle through blocks
            for block_name, block in self.environ_settings.items():
                # add header
                f.write(f'&{block_name}\n')
                # add key-value pairs
                for key, value in block.items():
                    if isinstance(value, str):
                        value = f"'{value}'"
                    f.write(f'   {key} = {value}\n')
                # add footer
                f.write('/\n')

    def check_code_is_installed(self):
        super().check_code_is_installed()
        if not environ_addon_is_installed(self.command.path.parent):
            raise OSError('The pw add-on "environ" is not installed')


def environ_addon_is_installed(qe_directory):
    # This is written in this way so it can be externally imported by the test suite
    return (qe_directory / 'Environ_PATCH').is_file()
