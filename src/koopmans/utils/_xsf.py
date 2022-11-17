from pathlib import Path
from typing import List, Tuple

import numpy as np
from ase import Atoms


def write_xsf(filename: Path, atoms: Atoms, arrays: List[np.ndarray], nr_xml: Tuple[int, int, int]):
    """
    Writes a quantity defined on the real space grid to a xsf file, which can be plotted with xcrysden
    """

    # Convert the arrays to xsf-format, where the last row at for each axis is identical to the first
    arrays_xsf: List[np.ndarray] = []
    for array in arrays:
        array_xsf = np.zeros((nr_xml[2], nr_xml[1], nr_xml[0]))
        for k in range(nr_xml[2]):
            for j in range(nr_xml[1]):
                for i in range(nr_xml[0]):
                    array_xsf[k, j, i] = array[k % (nr_xml[2]-1), j % (nr_xml[1]-1), i % (nr_xml[0]-1)]
        arrays_xsf.append(array_xsf)

    cell_parameters = atoms.get_cell()
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    with open(filename, 'w') as out:
        out.write('# xsf file \n')
        out.write('CRYSTAL\n\n')
        out.write('PRIMVEC\n\n')
        for vec in cell_parameters:
            out.write("\t" + " ".join([f"{x:13.10f}" for x in vec]))
        out.write('PRIMCOORD\n')
        out.write(f"\t{len(symbols)}\t1\n")
        for symbol, pos in zip(symbols, positions):
            out.write("\t" + symbol + " " + " ".join([f"{x:13.10f}" for x in pos]) + " \n")
        out.write('BEGIN_BLOCK_DATAGRID_3D\n')
        out.write("\t" + 'my_first_example_of_3D_datagrid\n')
        for i, array_xsf in enumerate(arrays_xsf):
            out.write("\t" + 'BEGIN_DATAGRID_3D_this_is_3Dgrid#' + str(i+1) + '\n')
            out.write("\t" + "\t" + str(nr_xml[0]) + '\t' + str(nr_xml[1]) + '\t' + str(nr_xml[2]) + '\t\n')
            out.write("\t" + "\t" + str(0.0) + '\t' + str(0.0) + '\t' + str(0.0) + '\t\n')  # origin of the data grid
            # third spanning vector of the data grid
            out.write("\t" + "\t" + str(cell_parameters[0][0]) + '\t' + str(0.0) + '\t' + str(0.0) + '\t\n')
            # second spanning vector of the data grid
            out.write("\t" + "\t" + str(0.0) + '\t' + str(cell_parameters[1][1]) + '\t' + str(0.0) + '\t\n')
            out.write("\t" + "\t" + str(0.0) + '\t' + str(0.0) + '\t' +
                      str(cell_parameters[2][2]) + '\t\n')  # first spanning vector of the data grid
            for k in range(nr_xml[2]):
                for j in range(nr_xml[1]):
                    out.write("\t\t")
                    for i in range(nr_xml[0]):
                        out.write("{:.15E}\t".format(array[k, j, i]))
                    out.write('\n')
                out.write("\n\n")
            out.write("\n\t" + 'END_DATAGRID_3D\n')
        out.write('END_BLOCK_DATAGRID_3D')
