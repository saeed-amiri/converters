import itertools
from pprint import pprint
import re
import typing
import pandas as pd
import sys
import read_lmp_data as rdlmp
from colors_text import TextColor as bcolors


class Doc:
    """read the lammps data file and clean it:
        read all the data
        add names from the mass section to all atoms, bonds, ...
        remove empty sections
    Input:
        A LAMMPS data file
    Restriction:
        Mass section must contain atoms name right after a hash '#'
    """


# A helper function needed by most of the classes to get types for LAMMPS
def get_type(lst: list[str]  # list to get the number of distenguished ones
             ) -> list[int]:  # types' index
    """make type based on the unique items in the lst"""
    # eleminate the repeated names
    seen: set[str] = set()
    seen_add = seen.add
    type_set = [x for x in lst if not (x in seen or seen_add(x))]
    type_dict: dict[str, int]  # to make a list with type index
    type_dict = {item: i+1 for i, item in enumerate(type_set)}
    types: list[int]  # types to return
    types = [type_dict[item] for item in lst]
    return types


class CleanData:
    """read data and clean it"""
    def __init__(self,
                 fname: str  # Name of the input file
                 ) -> None:
        print(f'{bcolors.OKCYAN}{self.__class__.__name__}:\n'
              f'\tCleaning: `{fname}`\n')
        self.raw_data = rdlmp.ReadData(fname)
        self.clean_data()

    def clean_data(self) -> None:
        """call all the methods"""
        self.Atoms_df: pd.DataFrame = self.raw_data.Atoms_df
        self.Bonds_df: pd.DataFrame = self.get_bonds()  # Bonds df to write
        self.Angles_df: pd.DataFrame = self.get_angles()  # Angles df to write
        self.Dihedrals_df: pd.DataFrame = self.get_dihedrals()  # Dihedrals df
        self.Masses_df: pd.DataFrame = self.raw_data.Masses_df

    def get_bonds(self) -> pd.DataFrame:  # Bonds DataFrame for writing
        """correct the name and type of the bonds"""
        names: list[str]  # Bonds names
        types: list[int]  # Bonds type from names
        names = self.bonds_name()
        types = get_type(names)
        columns: list[str]  # DataFrame columns for bonds in LAMMPS
        columns = ['typ', 'ai', 'aj', 'cmt', 'name']
        df: pd.DataFrame  # Bonds df to write out
        df = pd.DataFrame(columns=columns)
        df['typ'] = types
        df.index += 1  # Since the raw data increased one
        df['ai'] = self.raw_data.Bonds_df['ai']
        df['aj'] = self.raw_data.Bonds_df['aj']
        df['name'] = names
        df['cmt'] = ['#' for _ in df.index]
        return df

    def bonds_name(self) -> list[str]:  # Name of the bonds
        """return name of the bonds by making from Atoms_df name"""
        ai_name: list[str]  # 1st atoms names
        aj_name: list[str]  # 2nd atoms names
        bonds: list[str]  # Bonds names
        ai_name = [
                self.raw_data.Atoms_df.loc[
                    self.raw_data.Atoms_df['atom_id'] == item
                    ]['name'][item]
                for item in self.raw_data.Bonds_df['ai']
                ]
        aj_name = [
                self.raw_data.Atoms_df.loc[
                    self.raw_data.Atoms_df['atom_id'] == item
                    ]['name'][item]
                for item in self.raw_data.Bonds_df['aj']
                ]
        bonds = [f'{i}_{j}' for i, j in zip(ai_name, aj_name)]
        return bonds

    def get_angles(self) -> pd.DataFrame:  # Angles DataFrame for writing
        """correct the name and type of the angles"""
        names: list[str]  # Angles names
        types: list[int]  # Angles type from names
        names = self.angles_name()
        types = get_type(names)
        angle_type, type_name = self.angles_type(names)
        columns: list[str]  # DataFrame columns for angles in LAMMPS
        columns = ['typ', 'ai', 'aj', 'ak', 'cmt', 'name']
        df: pd.DataFrame  # Angles df to write out
        df = pd.DataFrame(columns=columns)
        df['typ'] = types
        df.index += 1   # Since the raw data increased one
        df['ai'] = self.raw_data.Angles_df['ai']
        df['aj'] = self.raw_data.Angles_df['aj']
        df['ak'] = self.raw_data.Angles_df['ak']
        df['name'] = names
        df['cmt'] = ['#' for _ in df.index]
        return df

    def angles_name(self) -> list[str]:  # Name of the angles
        """return name of the angles by making from Atoms_df name"""
        ai_name: list[str]  # 1st atoms names
        aj_name: list[str]  # 2nd atoms names
        ak_name: list[str]  # 3rd atoms names
        angles: list[str]  # Angles names
        ai_name = [
                self.raw_data.Atoms_df.loc[
                    self.raw_data.Atoms_df['atom_id'] == item
                    ]['name'][item]
                for item in self.raw_data.Angles_df['ai']
                  ]
        aj_name = [
                self.raw_data.Atoms_df.loc[
                    self.raw_data.Atoms_df['atom_id'] == item
                    ]['name'][item]
                for item in self.raw_data.Angles_df['aj']
                  ]
        ak_name = [
                self.raw_data.Atoms_df.loc[
                    self.raw_data.Atoms_df['atom_id'] == item
                    ]['name'][item]
                for item in self.raw_data.Angles_df['ak']
                  ]
        angles = [f'{i}_{j}_{k}' for i, j, k in zip(ai_name, aj_name, ak_name)]
        return angles

    def angles_type(self,
                    names: list[str]  # Names of the bonds
                    ) -> list[int]:  # Type of the angles
        """make a correct type for angles
        Angle ABC is same as CBA
        Fixing this is complicated since it entirely depends on the
        names of the type of atoms in the Mass section.
        """
        # Remove digits from the name of the atoms:
        tmp_name: list[list[str]]  # name of the bonds without digits
        tmp_name = [re.sub('\d+', '', item) for item in names]
        tmp_lst: list[list[str]] = [item.split('_') for item in tmp_name]
        # uniqe names after removing digits:
        name_set: set[str]
        name_set = self.seen_set(tmp_name)  # Remove the duplicates
        name_set = [item.split('_') for item in name_set]
        name_set = [item for sublist in name_set for item in sublist]
        type_set = self.seen_set(name_set)  # Final set of uniqes atoms
        # Set an index to each atom:
        type_dict: dict[str, int]  # dict for name and type
        type_dict = {item: v+1 for v, item in enumerate(type_set)}
        # new type for all sets of angles
        name_lst: list[list[str]]  # type of each set of angles
        name_lst = [[str(type_dict[i]) for i in item] for item in tmp_lst]
        # Make one list of all the angles with thier atom types
        name_int: list[str]  # make a str list of int of types
        name_int = ['_'.join(item) for item in name_lst]
        # Get the uniqe type of all angles
        seen: set[str] = set()
        name_tup = [tuple([item]) for item in name_int]
        # name_tup = [item[0] for item in name_tup]
        # Make a list of each set of angles
        t: list[str] = [
            x for x in name_tup if frozenset(x) not in seen and
            not seen.add(frozenset(x))
            ]
        t = [list(item)[0] for item in t]  # make list of items
        # remove duplicates and reverse duplicated
        t = list({i[::-1] if i[-1] < i[0] else i: i for i in t}.values())
        # Give a name to each uniqe set of angle
        angle_dict: dict[str, int]
        angle_dict = {item: v+1 for v, item in enumerate(t)}
        # List of types for each angle
        angle_type: list[int] = []  # final list for each angle
        for item in name_int:
            typ: int  # index for each angle
            try:
                typ = angle_dict[item]
            except KeyError:
                typ = angle_dict[item[::-1]]
            angle_type.append(typ)
        name_dict: dict[str, int] = {}  # name angles based on the types
        for k, v in angle_dict.items():
            n = []
            for i in k.split('_'):
                value = self.get_key(type_dict, int(i))
                n.append(value)
            name_dict[v] = f'_'.join(n)
        type_name: list[str] = []  # main name of each angle
        for item in angle_type:
            type_name.append(name_dict[item])
        return angle_type, type_name

    # function to return key for any value
    def get_key(self,
                dic,
                val):
        for key, value in dic.items():
            if val == value:
                return key

    def seen_set(self,
                 lst: list[typing.Any]  # to drop duplicate with keping order
                 ) -> list[typing.Any]:
        """remove duplicated item with keeping order of them in the
        main list"""
        seen: set[str] = set()
        seen_add = seen.add
        return [x for x in lst if not (x in seen or seen_add(x))]

    def get_dihedrals(self) -> pd.DataFrame:  # Dihedrals DataFrame for writing
        """correct the name and type of the dihedrals"""
        names: list[str]  # Dihedrals names
        types: list[int]  # Dihedrals type from names
        names = self.dihedrals_name()
        types = get_type(names)
        columns: list[str]  # DataFrame columns for dihedrals in LAMMPS
        columns = ['typ', 'ai', 'aj', 'ak', 'ah', 'cmt', 'name']
        df: pd.DataFrame  # Dihedrals df to write out
        df = pd.DataFrame(columns=columns)
        df['typ'] = types
        df.index += 1   # Since the raw data increased one
        df['ai'] = self.raw_data.Dihedrals_df['ai']
        df['aj'] = self.raw_data.Dihedrals_df['aj']
        df['ak'] = self.raw_data.Dihedrals_df['ak']
        df['ah'] = self.raw_data.Dihedrals_df['ah']
        df['name'] = names
        df['cmt'] = ['#' for _ in df.index]
        return df

    def dihedrals_name(self) -> list[str]:  # Name of the dihedrals
        """return name of the dihedrals by making from Atoms_df name"""
        ai_name: list[str]  # 1st atoms names
        aj_name: list[str]  # 2nd atoms names
        ak_name: list[str]  # 3rd atoms names
        ah_name: list[str]  # 4th atoms names
        dihedrals: list[str]  # Dihedrals names
        ai_name = [
                self.raw_data.Atoms_df.loc[
                    self.raw_data.Atoms_df['atom_id'] == item
                    ]['name'][item]
                for item in self.raw_data.Dihedrals_df['ai']
                  ]
        aj_name = [
                self.raw_data.Atoms_df.loc[
                    self.raw_data.Atoms_df['atom_id'] == item
                    ]['name'][item]
                for item in self.raw_data.Dihedrals_df['aj']
                  ]
        ak_name = [
                self.raw_data.Atoms_df.loc[
                    self.raw_data.Atoms_df['atom_id'] == item
                    ]['name'][item]
                for item in self.raw_data.Dihedrals_df['ak']
                  ]
        ah_name = [
                self.raw_data.Atoms_df.loc[
                    self.raw_data.Atoms_df['atom_id'] == item
                    ]['name'][item]
                for item in self.raw_data.Dihedrals_df['ah']
                  ]
        dihedrals = [
            f'{i}_{j}_{k}_{h}' for i, j, k, h in
            zip(ai_name, aj_name, ak_name, ah_name)]
        return dihedrals


if __name__ == '__main__':
    fname = sys.argv[1]
    data = CleanData(fname)
