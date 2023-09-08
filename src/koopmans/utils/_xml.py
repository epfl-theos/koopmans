import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Tuple

import numpy as np


def _load_xml_element(xml_file: Path, string: str) -> ET.Element:
    with open(xml_file, 'r') as fd:
        tree = ET.parse(fd)
    root = tree.getroot()
    branch = root.find(string)
    assert isinstance(branch, ET.Element)
    return branch


def _get_nr(element: ET.Element) -> List[int]:
    # Extract the grid dimension
    info = element.find('INFO')
    assert isinstance(info, ET.Element)
    nr = []
    for i in range(3):
        nr_str = info.get(f'nr{i+1}')
        assert isinstance(nr_str, str)
        nr.append(int(nr_str) + 1)
    return nr


def read_xml_nr(xml_file: Path, string: str = 'EFFECTIVE-POTENTIAL') -> List[int]:
    """
    Loads the dimensions of an array from an xml file
    """

    branch = _load_xml_element(xml_file, string)
    return _get_nr(branch)


def read_xml_array(
    xml_file: Path, norm_const: float, string: str = 'EFFECTIVE-POTENTIAL', retain_final_element: bool=False
    ) -> np.ndarray:
    """
    Loads an array from an xml file.
    
    :param xml_file: The xml file to read from
    :param norm_const: The normalization constant to multiply the array with (in our case 1/((Bohr radii)^3)
    :param string: The name of the field in the xml file that contains the array, in our case either 
    'EFFECTIVE-POTENTIAL' or 'CHARGE-DENSITY'
    :param retain_final_element: If True, the array is returned in with periodic boundary conditions, i.e. the last 
    element in each dimension is equal to the first element in each dimension. This is required for the xsf format.
    
    :return: The array
    """

    # Load the branch of the xml tree
    branch = _load_xml_element(xml_file, string)

    # Extract the nr grid
    nr_xml = _get_nr(branch)

    # Extract the array
    array_xml = np.zeros((nr_xml[2], nr_xml[1], nr_xml[0]), dtype=float)
    for k in range(nr_xml[2]):
        current_name = 'z.' + str(k % (nr_xml[2]-1)+1)
        entry = branch.find(current_name)
        assert isinstance(entry, ET.Element)
        text = entry.text
        assert isinstance(text, str)
        rho_tmp = np.array(text.split('\n')[1:-1], dtype=float)
        for j in range(nr_xml[1]):
            for i in range(nr_xml[0]):
                array_xml[k, j, i] = rho_tmp[(j % (nr_xml[1]-1))*(nr_xml[0]-1)+(i % (nr_xml[0]-1))]
    array_xml *= norm_const

    if retain_final_element: 
        # the xsf format requires an array where the last element is equal to the first element in each dimension
        return array_xml
    else:
        return array_xml[:-1, :-1, :-1]
