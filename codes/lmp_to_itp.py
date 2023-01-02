import pandas as pd
import read_lmp_data as relmp
import lmp_to_pdb as lmpdb
from colors_text import TextColor as bcolors


class Doc:
    """read the LAMMPS data file and get force field parameters and write
        `itp` formats for GROMACS conversion"""


class ItpStyleInfos:
    """ informations in the itp file:
        [ moleculetype ] : defines the name of your molecule in this top
    and nrexcl = 3 stands for excluding non-bonded interactions between
    atoms that are no further than 3 bonds away.

        [ atoms ] : defines the molecule, where nr and type are fixed,
    the rest is user defined. So atom can be named as you like,
    cgnr made larger or smaller (if possible, the total charge of
    a charge group should be zero), and charges can be changed here
    too.

        [ bonds ] : no comment.

        [ pairs ] : LJ and Coulomb 1-4 interactions

        [ angles ] : no comment

        [ dihedrals ] : in this case there are 9 proper dihedrals
    (funct = 1), 3 improper (funct = 4) and no Ryckaert-Bellemans type
    dihedrals."""


class Itp:
    """get data from main"""
    def __init__(self,
                 lmp: relmp.ReadData,  # LAMMPS data file
                 pdb_df: pd.DataFrame  # Final df for pdb file
                 ) -> None:
        self.__mk_itp(lmp, pdb_df)

    def __mk_itp(self,
                 lmp: relmp.ReadData,  # LAMMPS data file
                 pdb_df: pd.DataFrame  # Final df for pdb file
                 ) -> None:
        """call functions"""
        self.atoms = self.__mk_atoms(lmp, pdb_df)
        self.bonds = self.__mk_bonds(lmp)
        self.angles = self.__mk_angles(lmp)
        self.dihedrals = self.__mk_dihedrals(lmp)
        self.__mk_pairs(lmp)

    def __mk_atoms(self,
                   lmp: relmp.ReadData,  # LAMMPS data file
                   pdb_df: pd.DataFrame  # Final df for pdb file
                   ) -> pd.DataFrame:
        df: pd.DataFrame  # df in the itp format
        columns: list[str]  # Columns of the DataFrame
        columns = ['atomnr',  # Index of the atoms in lmp file
                   'atomtype',  # FF type of the atom, e.g. for c, opls_135
                   'resnr',  # Index of the residue which atom is belonged
                   'resname',  # Name of the residue which atom is belonged
                   'atomname',  # Name of the atom as in PDB file
                   'chargegrp',  # Charge group as in PDB file
                   'charge',  # Charge as in the lmp file
                   'mass',  # Mass odf the atom as in the lmp file
                   ' ',  # Comment column for ;
                   '  '  # Second column for the coments
                   ]
        df = pd.DataFrame(columns=columns)
        df['atomnr'] = pdb_df['atom_id']
        df['atomtype'] = pdb_df['ff_type']
        df['resnr'] = pdb_df['residue_id']
        df['resname'] = pdb_df['residue_name']
        df['atomname'] = pdb_df['atom_name']
        df['chargegrp'] = pdb_df['charge']
        df['charge'] = pdb_df['q']
        df['mass'] = pdb_df['mass']
        df['chargegrp'] = [1 for _ in df['atomnr']]
        df[' '] = [';' for _ in df['atomnr']]
        df['  '] = pdb_df["element"]
        return df

    def __mk_bonds(self,
                   lmp: relmp.ReadData  # LAMMPS data file
                   ) -> pd.DataFrame:
        df: pd.DataFrame  # df in the itp format
        columns: list[str]  # Columns of the DataFrame
        columns = ['ai',  # 1st atom in bond
                   'aj',  # 2nd atom in bond
                   'funct',  # not sure what is this, just set to 1, or empty!
                   'r',  # Distance parameter in harmonic bond interactions
                   'k',  # Harmonic constant in the harmonic interactions
                   ' '  # Comment: ;
                   '  '  # Comment: name of the bond
                   ]
        df = pd.DataFrame(columns=columns)
        Bonds_df: pd.DataFrame = lmp.Bonds_df.sort_values(by='ai')
        df['ai'] = Bonds_df['ai']
        df['aj'] = Bonds_df['aj']
        df['funct'] = [1 for _ in df['ai']]
        df['r'] = ['r' for _ in df['ai']]
        df['k'] = ['K' for _ in df['ai']]
        try:
            df[' '] = [';' for _ in df['ai']]
            df['  '] = lmp.Bonds_df['name']
        except KeyError:
            df.drop(columns=[' '])
            print(f'{bcolors.WARNING}{self.__class__.__name__}:\n'
                  f'\t There is no bonds` names in LAMMPS read data\n'
                  f'{bcolors.ENDC}')
        return df

    def __mk_angles(self,
                    lmp: relmp.ReadData  # LAMMPS data file
                    ) -> pd.DataFrame:
        df: pd.DataFrame  # df in the itp format
        columns: list[str]  # Columns of the DataFrame
        columns = ['ai',  # 1st atom in angle
                   'aj',  # 2nd atom in angle
                   'ak',  # 3rd atom in angle
                   'funct',  # not sure what is this, just set to 1, or empty!
                   'theta',  # The angle between bonds
                   'cth',  # Strength of the bonds
                   ' '  # Comment: name of the angle
                   ]
        df = pd.DataFrame(columns=columns)
        Angles_df: pd.DataFrame = lmp.Angles_df.sort_values(by='ai')
        df['ai'] = Angles_df['ai']
        df['aj'] = Angles_df['aj']
        df['ak'] = Angles_df['ak']
        df['funct'] = [1 for _ in df['ai']]
        df['theta'] = ['theta' for _ in df['ai']]
        df['cth'] = ['cth' for _ in df['ai']]
        try:
            df[' '] = '; ' + lmp.Angles_df['name']
        except KeyError:
            df.drop(columns=[' '], inplace=True)
            print(f'{bcolors.WARNING}{self.__class__.__name__}:\n'
                  f'\t There is no angles` names in LAMMPS read data\n'
                  f'{bcolors.ENDC}')
        return df

    def __mk_dihedrals(self,
                       lmp: relmp.ReadData  # LAMMPS data file
                       ) -> pd.DataFrame:
        df: pd.DataFrame  # df in the itp format
        columns: list[str]  # Columns of the DataFrame
        columns = ['ai',  # 1st atom in dihedrals
                   'aj',  # 2nd atom in dihedrals
                   'ak',  # 3rd atom in dihedrals
                   'ah',  # 4th atom in dihedrals
                   'funct',  # not sure what is this, just set to 1, or empty!
                   'C0',  # Dihedrals parameters
                   'C1',  # Dihedrals parameters
                   'C2',  # Dihedrals parameters
                   'C3',  # Dihedrals parameters
                   'C4',  # Dihedrals parameters
                   ' '  # Comment: name of the dihedrals
                   ]
        df = pd.DataFrame(columns=columns)
        Dihedrals_df: pd.DataFrame = lmp.Dihedrals_df.sort_values(by='ai')
        df['ai'] = Dihedrals_df['ai']
        df['aj'] = Dihedrals_df['aj']
        df['ak'] = Dihedrals_df['ak']
        df['ah'] = Dihedrals_df['ah']
        df['funct'] = [1 for _ in df['ai']]
        df['C0'] = ['C0' for _ in df['ai']]
        df['C1'] = ['C1' for _ in df['ai']]
        df['C2'] = ['C2' for _ in df['ai']]
        df['C3'] = ['C3' for _ in df['ai']]
        df['C4'] = ['C4' for _ in df['ai']]
        try:
            df[' '] = '; ' + lmp.Dihedrals_df['name']
        except KeyError:
            df.drop(columns=[' '], inplace=True)
            print(f'{bcolors.WARNING}{self.__class__.__name__}:\n'
                  f'\t There is no dihedralss` names in LAMMPS read data\n'
                  f'{bcolors.ENDC}')
        return df

    def __mk_pairs(self,
                   lmp: relmp.ReadData  # LAMMPS data file
                   ) -> pd.DataFrame:
        df: pd.DataFrame  # df in the itp format
        columns: list[str]  # Columns of the DataFrame
