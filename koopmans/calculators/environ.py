import os
from koopmans.calculators.pw import PW_calc

_default_settings = {
   'ENVIRON': {
   'verbose': 0,
   'environ_type': 'input',
   'environ_type': 'input',
   'env_surface_tension': 0,
   'env_pressure': 0},
   'BOUNDARY': {},
   'ELECTROSTATIC': {}}

class Environ_calc(PW_calc):
   # Create an environ calculator that inherits from the vanilla pw.x calculator

   environ_settings = _default_settings

   def calculate(self):
      # Generic function for running a calculation

      # If pseudo_dir is a relative path then make sure it accounts for self.directory
      if self.pseudo_dir is not None and self.pseudo_dir[0] != '/':
          directory_depth = self.directory.strip('./').count('/') + 1
          self.pseudo_dir = '../'*directory_depth + self.pseudo_dir

      self.write_environ_in()

      # Ensure we're using an environ-enabled version of pw.x
      if '--environ' not in self._ase_calc.command:
          self._ase_calc.command = self._ase_calc.command.replace('pw.x', 'pw.x --environ')

      self._ase_calculate()

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
             
