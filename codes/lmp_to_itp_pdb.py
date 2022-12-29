import pandas as pd
import read_lmp_data as relmp


class Doc:
    """read data file in LAMMPS format and convert it to pdb and itp
    files for GROMACS
    The input is the LAMMPS data file:
        Name of the atoms in Masses section should be the ones in the
        force filed file in GROMACS
    """
