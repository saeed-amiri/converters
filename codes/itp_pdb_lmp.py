import typing
import pandas as pd
from itp_to_df import Itp
from pdb_to_df import Pdb

class Doc:
    """call the classes to read itp and pdb files and return data
    in format that can be write into LAMMPS files"""


class ItpPdb(Pdb,  # class which give dataframe for the pdb file
             Itp  # class which give dataframe for the itp file
             ):
    """call all the classes and update the Atoms_df and all other dfs
    so can be writeed by write_lmp.py
    """
    def __init__(self,
                 fname: str  # itp file name, for pdb it shoud write
                 ) -> None:
        Pdb.__init__(self, f'{fname}.pdb')
        Itp.__init__(self, f'{fname}.itp')
        print(self.atoms_extra)

if __name__ == '__main__':
    data = ItpPdb('UND')