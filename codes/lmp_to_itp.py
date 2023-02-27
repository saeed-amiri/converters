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
        self.atoms = self.__mk_atoms(pdb_df)
        self.bonds = self.__mk_bonds(lmp)
        self.angles = self.__mk_angles(lmp)
        self.dihedrals = self.__mk_dihedrals(lmp)
        self.__mk_pairs(lmp)

    def __mk_atoms(self,
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
                   'element'  # Second column for the coments
                   ]
        df = pd.DataFrame(columns=columns)
        df['atomnr'] = [int(item) for item in pdb_df['atom_id']]
        df['atomtype'] = [str(item) for item in pdb_df['ff_type']]
        df['resnr'] = [int(item) for item in pdb_df['residue_id']]
        df['resname'] = [str(item) for item in pdb_df['residue_name']]
        df['atomname'] = [str(item) for item in pdb_df['atom_name']]
        df['chargegrp'] = [int(1) for _ in pdb_df['charge']]
        df['charge'] = [float(item) for item in pdb_df['q']]
        df['mass'] = [float(item) for item in pdb_df['mass']]
        df[' '] = [';' for _ in df['atomnr']]
        df['element'] = [str(item) for item in pdb_df["element"]]
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
                   ' ',  # Comment: ;
                   '  ',  # Comment: name of the bond
                   'resname',  # Name of the residue which atoms belonged to
                   'resnr',  # Nr. of the residue which atoms belonged to
                   ]
        df = pd.DataFrame(columns=columns)
        Bonds_df: pd.DataFrame = lmp.Bonds_df.sort_values(by='ai')
        df['ai'] = Bonds_df['ai']
        df['aj'] = Bonds_df['aj']
        df['funct'] = [1 for _ in df['ai']]
        # df['r'] = [' ' for _ in df['ai']]
        # df['k'] = [' ' for _ in df['ai']]
        try:
            df[' '] = [';' for _ in df['ai']]
            df['  '] = lmp.Bonds_df['name']
        except KeyError:
            df.drop(columns=[' '])
            print(f'{bcolors.WARNING}{self.__class__.__name__}:\n'
                  f'\t There is no bonds` names in LAMMPS read data\n'
                  f'{bcolors.ENDC}')
        df['resname'], df['resnr'] = self.__get_bonds_res(lmp, df)
        return df

    def __get_bonds_res(self,
                        lmp: relmp.ReadData,  # LAMMPS data file
                        df: pd.DataFrame  # df contain itp info
                        ) -> tuple[list]:
        """return residues name and index"""
        resname: list[str] = []  # Name of the residues
        resnr: list[int] = []  # index of the residues
        flag_war: bool = False  # To print the warning
        for ai, aj in zip(df['ai'], df['aj']):
            ai_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == ai]['typ'][ai]
            aj_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == aj]['typ'][aj]
            mol_i: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ai_type]['residues'][ai_type]
            mol_j: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == aj_type]['residues'][aj_type]
            mol_iid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == ai]['mol'][ai]
            mol_jid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == aj]['mol'][aj]
            if mol_i != mol_j or mol_iid != mol_jid:
                if not flag_war:
                    print(f'{bcolors.WARNING}{self.__class__.__name__}:\n'
                         f'\tBond between atoms with diffrents residues '
                         f'type\n{bcolors.ENDC}')
                    flag_war = True
            resnr.append(mol_iid)
            resname.append(mol_i)
        return resname, resnr

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
                   ' ',  # Comment: name of the angle
                   'angle_name',  # Comment: name of the angle
                   'resname',  # Name of the residue which atoms belonged to
                   'resnr',  # Nr. of the residue which atoms belonged to
                   ]
        df = pd.DataFrame(columns=columns)
        Angles_df: pd.DataFrame = lmp.Angles_df.sort_values(by='ai')
        df['ai'] = Angles_df['ai']
        df['aj'] = Angles_df['aj']
        df['ak'] = Angles_df['ak']
        df['funct'] = [1 for _ in df['ai']]
        # df['theta'] = [' ' for _ in df['ai']]
        # df['cth'] = [' ' for _ in df['ai']]
        df[' '] = [';' for _ in df['ai']]
        try:
            df['angle_name'] = lmp.Angles_df['name']
        except KeyError:
            print(f'{bcolors.WARNING}{self.__class__.__name__}:\n'
                  f'\t There is no angles` names in LAMMPS read data\n'
                  f'{bcolors.ENDC}')
        df['resname'], df['resnr'] = self.__get_angles_res(lmp, df)
        return df

    def __get_angles_res(self,
                         lmp: relmp.ReadData,  # LAMMPS data file
                         df: pd.DataFrame  # df contain itp info
                         ) -> tuple[list]:
        """return residues name and index"""
        resname: list[str] = []  # Name of the residues
        resnr: list[int] = []  # index of the residues
        flag_war: bool = False  # To print the warning
        for ai, aj, ak in zip(df['ai'], df['aj'], df['ak']):
            ai_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == ai]['typ'][ai]
            aj_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == aj]['typ'][aj]
            ak_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == ak]['typ'][ak]
            mol_i: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ai_type]['residues'][ai_type]
            mol_j: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == aj_type]['residues'][aj_type]
            mol_k: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ak_type]['residues'][ak_type]
            mol_iid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == ai]['mol'][ai]
            mol_jid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == aj]['mol'][aj]
            mol_kid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == ak]['mol'][ak]
            check_list_name = set([mol_i, mol_j, mol_k])
            check_list_id = set([mol_iid, mol_jid, mol_kid])
            if len(check_list_name) != 1 or len(check_list_id) != 1:
                if not flag_war:
                    print(f'{bcolors.WARNING}{self.__class__.__name__}:\n'
                         f'\tangles between atoms with diffrents residues '
                         f'type\n{bcolors.ENDC}')
                flag_war = True
            resnr.append(mol_iid)
            resname.append(mol_i)
        return resname, resnr

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
                   ' ',  # Comment: name of the dihedrals
                   'dihedral_name',  # names
                   'resname',  # Name of the residue which atoms belonged to
                   'resnr',  # Nr. of the residue which atoms belonged to
                   ]
        df = pd.DataFrame(columns=columns)
        Dihedrals_df: pd.DataFrame = lmp.Dihedrals_df.sort_values(by='ai')
        df['ai'] = Dihedrals_df['ai']
        df['aj'] = Dihedrals_df['aj']
        df['ak'] = Dihedrals_df['ak']
        df['ah'] = Dihedrals_df['ah']
        df['funct'] = [3 for _ in df['ai']]
        # df['C0'] = [' ' for _ in df['ai']]
        # df['C1'] = [' ' for _ in df['ai']]
        # df['C2'] = [' ' for _ in df['ai']]
        # df['C3'] = [' ' for _ in df['ai']]
        # df['C4'] = [' ' for _ in df['ai']]
        df[' '] = [';' for _ in df['ai']]
        try:
            df['dihedral_name'] = lmp.Dihedrals_df['name']
        except KeyError:
            print(f'{bcolors.WARNING}{self.__class__.__name__}:\n'
                  f'\t There is no dihedralss` names in LAMMPS read data\n'
                  f'{bcolors.ENDC}')
        df['resname'], df['resnr'] = self.__get_dihedrals_res(lmp, df)
        return df

    def __get_dihedrals_res(self,
                            lmp: relmp.ReadData,  # LAMMPS data file
                            df: pd.DataFrame  # df contain itp info
                            ) -> tuple[list]:
        """return residues name and index"""
        resname: list[str] = []  # Name of the residues
        resnr: list[int] = []  # index of the residues
        flag_war: bool = False  # To print the warning
        for ai, aj, ak, ah in zip(df['ai'], df['aj'], df['ak'], df['ah']):
            ai_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == ai]['typ'][ai]
            aj_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == aj]['typ'][aj]
            ak_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == ak]['typ'][ak]
            ah_type: int = lmp.Atoms_df.loc[
                           lmp.Atoms_df['atom_id'] == ah]['typ'][ah]
            mol_i: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ai_type]['residues'][ai_type]
            mol_j: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == aj_type]['residues'][aj_type]
            mol_k: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ak_type]['residues'][ak_type]
            mol_h: str = lmp.Masses_df[
                         lmp.Masses_df['typ'] == ah_type]['residues'][ah_type]
            mol_iid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == ai]['mol'][ai]
            mol_jid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == aj]['mol'][aj]
            mol_kid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == ak]['mol'][ak]
            mol_hid: int = lmp.Atoms_df[
                           lmp.Atoms_df['atom_id'] == ah]['mol'][ah]
            check_list_name = set([mol_i, mol_j, mol_k, mol_h])
            check_list_id = set([mol_iid, mol_jid, mol_kid, mol_hid])
            if len(check_list_name) != 1 or len(check_list_id) != 1:
                if not flag_war:
                    print(f'{bcolors.WARNING}{self.__class__.__name__}:\n'
                         f'\tdihedrals between atoms with diffrents residues '
                         f'type\n{bcolors.ENDC}')
                    flag_war = True

            resnr.append(mol_iid)
            resname.append(mol_i)
        return resname, resnr

    def __mk_pairs(self,
                   lmp: relmp.ReadData  # LAMMPS data file
                   ) -> pd.DataFrame:
        df: pd.DataFrame  # df in the itp format
        columns: list[str]  # Columns of the DataFrame
