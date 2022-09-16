import typing
import pandas as pd
from itp_to_df import Itp, get_type
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
        self.Atoms_df = self.mk_atoms()
        print(self.Atoms_df)

    def mk_atoms(self) -> pd.DataFrame:
        """make DataFrame form the data files"""
        columns: list[str]  # Header of the LAMMPS atoms full style
        columns = ['atom_id', 'mol', 'typ', 'charge', 'x', 'y', 'z',
                   'nx', 'ny', 'nz', 'cmt', 'name']
        Atoms_df: pd.DataFrame  # The main dataframe for the LAMMPS
        Atoms_df = pd.DataFrame(columns=columns)
        Atoms_df['atom_id'] = self.atoms['atom_id']
        Atoms_df['mol'] = self.atoms['residue_number']
        Atoms_df['typ'] = get_type(self.atoms_extra['name'])
        Atoms_df['charge'] = self.atoms_extra['charge']
        Atoms_df['x'] = self.atoms['x']
        Atoms_df['y'] = self.atoms['y']
        Atoms_df['z'] = self.atoms['z']
        Atoms_df['nx'] = [0 for _ in self.atoms.index]
        Atoms_df['ny'] = Atoms_df['nx']
        Atoms_df['nz'] = Atoms_df['nx']
        Atoms_df['cmt'] = ['#' for _ in self.atoms.index]
        Atoms_df['name'] = self.atoms_extra['name']
        return Atoms_df


if __name__ == '__main__':
    data = ItpPdb('UND')
