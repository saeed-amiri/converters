import sys
import pandas as pd
import read_lmp_data as relmp
import lmp_to_pdb as lmpdb


class Doc:
    """read data file in LAMMPS format and convert it to pdb and itp
    files for GROMACS
    The input is the LAMMPS data file:
        Name of the atoms in Masses section should be the ones in the
        force filed file in GROMACS
    """


class ReadLmp:
    """read lammps from read_lmp_data.py"""
    def __init__(self,
                 fname: str  # Input file name
                 ) -> None:
        self.get_lmp(fname)

    def get_lmp(self,
                fname: str  # Input file name
                ) -> None:
        """get the data"""
        self.lmp_data = relmp.ReadData(fname)


if __name__ == '__main__':
    lmp: relmp.ReadData = ReadLmp(sys.argv[1])  # All data in input file
