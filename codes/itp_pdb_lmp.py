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
        self.Bonds_df = self.bonds
        self.Angles_df = self.angles
        self.Dihedrals_df = self.dihedrals
        self.Masses_df = self.mk_masses()

    def mk_atoms(self) -> pd.DataFrame:
        """make atom DataFrame form the data files"""
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

    def mk_bonds(self) -> None:
        """make bonds DataFrame from data files"""
        columns: list[str]  # Header of the bonds in LAMMPS
        columns = ['typ', 'ai', 'aj', 'cmt', 'name']
        # It made by Itp class

    def mk_angles(self) -> None:
        """make angles DataFrame from data files"""
        columns: list[str]  # Header of the angles in LAMMPS
        columns = ['typ', 'ai', 'aj', 'ak', 'cmt', 'name']
        # It made by Itp class

    def mk_masses(self) -> pd.DataFrame:
        """masses and type informatin about atoms"""
        columns: list[str] = ['typ', 'mass', 'cmt', 'name']
        df: pd.DataFrame  # temporary for the getting data
        df = self.atoms_extra.copy()
        df = df.astype({'mass': float})
        df = df.loc[df.groupby(by=['name'])['mass'].idxmin(), :]
        num_list: list[str] = list(df['atomnr'])  # Atoms index of grouped df
        typ: list[int]  # type pf atoms from main Atoms_df
        typ = [self.Atoms_df.loc[self.Atoms_df['atom_id'] ==
               int(item)]['typ'][int(item)-1] for item in num_list]
        Masses_df = pd.DataFrame(columns=columns)
        Masses_df['name'] = df['name']
        Masses_df['mass'] = df['mass']
        Masses_df['cmt'] = ['#' for _ in Masses_df.index]
        Masses_df['typ'] = typ
        Masses_df.reset_index(inplace=True)
        Masses_df.drop('index', axis=1, inplace=True)
        return Masses_df


if __name__ == '__main__':
    data = ItpPdb('UND')
