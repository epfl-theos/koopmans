'''

utils module for koopmans

Written by Edward Linscott May 2020

'''

from ._figures import savefig
from ._io import (construct_atomic_positions_block,
                  construct_atomic_species_block,
                  construct_cell_parameters_block,
                  generate_wannier_centers_file_contents,
                  generate_wannier_hr_file_contents,
                  generate_wannier_u_file_contents, indented_print, parse_dict,
                  parse_wannier_centers_file_contents,
                  parse_wannier_hr_file_contents,
                  parse_wannier_u_file_contents, print_alert, read_alpha_file,
                  read_atomic_positions, read_cell_parameters,
                  read_wannier_hr_file, write_alpha_file)
from ._misc import flatten, update_nested_dict
from ._os import (HasDirectory, chdir, copy, find_executable,
                  get_binary_content, get_content, remove, set_env, symlink,
                  symlink_tree, system_call, write_binary_content,
                  write_content)
from ._units import units
from ._warnings import CalculatorNotConvergedWarning, warn
from ._xml import read_xml_array, read_xml_nr
from ._xsf import write_xsf
