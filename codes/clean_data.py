import re
import sys
import typing
import pandas as pd
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


class Bonds:
    """clean the bonds section with the right names and types"""
    def __init__(self,
                 raw_data: rdlmp.ReadData  # raw data to clean the bonds
                 ) -> None:
        self.Bonds_df: pd.DataFrame = self.get_bonds(raw_data)  # Bonds df

    def get_bonds(self,
                  raw_data: rdlmp.ReadData  # raw data to clean the bonds
                  ) -> pd.DataFrame:  # Bonds DataFrame for writing
        """correct the name and type of the bonds"""
        names: list[str]  # Bonds names
        types: list[int]  # Bonds type from names
        types_name: list[str]  # Bonds types' names
        bonds_atom: list[str]  # Atoms name which share a bond
        names, bonds_atom = self.bonds_name(raw_data)
        typ = Types(names)
        types = typ.types
        types_name = typ.types_name
        columns: list[str]  # DataFrame columns for bonds in LAMMPS
        columns = ['typ', 'ai', 'aj', 'cmt', 'name', 'type_name']
        df: pd.DataFrame  # Bonds df to write out
        df = pd.DataFrame(columns=columns)
        df['typ'] = types
        df.index += 1  # Since the raw data increased one
        df['ai'] = raw_data.Bonds_df['ai']
        df['aj'] = raw_data.Bonds_df['aj']
        df['cmt'] = ['#' for _ in df.index]
        df['name'] = bonds_atom
        df['type_name'] = types_name
        return df

    def bonds_name(self,
                   raw_data: rdlmp.ReadData  # raw data to clean the bonds
                   ) -> tuple[list[str], list[str]]:  # Name of the bonds
        """return name of the bonds by making from Atoms_df name"""
        ai_name: list[str]  # 1st atoms names, type of atom for bond
        aj_name: list[str]  # 2nd atoms names, type of atom for bond
        ai_b_name: list[str]  # 1st atoms names which share bond
        aj_b_name: list[str]  # 2nd atoms names which share bond
        bonds: list[str]  # Bonds names as type of the bonds
        bonds_atom: list[str]  # Atoms name which share a bond
        ai_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['name'][item]
                for item in raw_data.Bonds_df['ai']
                ]
        aj_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['name'][item]
                for item in raw_data.Bonds_df['aj']
                ]
        ai_b_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['b_name'][item]
                for item in raw_data.Bonds_df['ai']
                ]
        aj_b_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['b_name'][item]
                for item in raw_data.Bonds_df['aj']
                ]
        bonds = [f'{i}_{j}' for i, j in zip(ai_b_name, aj_b_name)]
        bonds_atom = [f'{i}_{j}' for i, j in zip(ai_name, aj_name)]
        return bonds, bonds_atom


class Angles:
    """Clean the angles section by giving the right names abd types"""
    def __init__(self,
                 raw_data: rdlmp.ReadData  # raw data to clean the angles
                 ) -> None:
        try:
            self.Angles_df: pd.DataFrame = self.get_angles(raw_data)  # Angles
        except KeyError:
            pass

    def get_angles(self,
                   raw_data: rdlmp.ReadData  # raw data to clean the angles
                   ) -> pd.DataFrame:  # Angles DataFrame for writing
        """correct the name and type of the angles"""
        names: list[str]  # Angles names
        angles_atom: list[str]  # Name of the atoms in the angle
        names, angles_atom = self.angles_name(raw_data)
        angle_type: tuple[list[int]]  # type of each angle with ABC=CBA
        type_name: list[str]  # name of each type of each angle Angle
        typ = Types(names)
        angle_type, type_name = typ.types, typ.types_name
        columns: list[str]  # DataFrame columns for angles in LAMMPS
        columns = ['typ', 'ai', 'aj', 'ak', 'cmt', 'name', 'type_name']
        df: pd.DataFrame  # Angles df to write out
        df = pd.DataFrame(columns=columns)
        df['typ'] = angle_type
        df.index += 1   # Since the raw data increased one
        df['ai'] = raw_data.Angles_df['ai']
        df['aj'] = raw_data.Angles_df['aj']
        df['ak'] = raw_data.Angles_df['ak']
        df['name'] = angles_atom
        df['type_name'] = type_name
        df['cmt'] = ['#' for _ in df.index]
        return df

    def angles_name(self,
                    raw_data: rdlmp.ReadData  # raw data to clean the angles
                    ) -> tuple[list[str], list[str]]:  # Name of the angles
        """return name of the angles by making from Atoms_df name"""
        ai_a_name: list[str]  # 1st atoms names, type of atoms in the angle
        aj_a_name: list[str]  # 2nd atoms names, type of atoms in the angle
        ak_a_name: list[str]  # 3rd atoms names, type of atoms in the angle
        ai_name: list[str]  # 1st atoms names which share angles
        aj_name: list[str]  # 2nd atoms names which share angles
        ak_name: list[str]  # 3rd atoms names which share angles
        angles: list[str]  # Angles names
        angles_atom: list[str]  # Angles names based on the atom share them
        ai_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['name'][item]
                for item in raw_data.Angles_df['ai']
                  ]
        aj_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['name'][item]
                for item in raw_data.Angles_df['aj']
                  ]
        ak_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['name'][item]
                for item in raw_data.Angles_df['ak']
                  ]
        ai_a_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['b_name'][item]
                for item in raw_data.Angles_df['ai']
                  ]
        aj_a_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['b_name'][item]
                for item in raw_data.Angles_df['aj']
                  ]
        ak_a_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['b_name'][item]
                for item in raw_data.Angles_df['ak']
                  ]
        angles = [f'{i}_{j}_{k}' for i, j, k in zip(ai_a_name,
                                                    aj_a_name,
                                                    ak_a_name)]
        angles_atom = [f'{i}_{j}_{k}' for i, j, k in zip(ai_name,
                                                         aj_name,
                                                         ak_name)]
        return angles, angles_atom


class Dihedrals:
    """set the right names and types for dihedrals section"""
    def __init__(self,
                 raw_data: rdlmp.ReadData  # raw data to clean the angles
                 ) -> None:
        try:
            self.Dihedrals_df: pd.DataFrame = self.get_dihedrals(raw_data)
        except KeyError:
            pass

    def get_dihedrals(self,
                      raw_data: rdlmp.ReadData  # raw data to clean the dihedra
                      ) -> pd.DataFrame:  # Dihedrals DataFrame for writing
        """correct the name and type of the dihedrals"""
        names: list[str]  # Dihedrals names
        types: list[int]  # Dihedrals type from names
        type_names: list[str]  # Dihedrals types' names
        dihedrals_atoms: list[str]  # Names by atoms in the dihedrals
        names, dihedrals_atoms = self.dihedrals_name(raw_data)
        typ = Types(names, dihedrals=True)
        types = typ.types
        type_names = typ.types_name
        columns: list[str]  # DataFrame columns for dihedrals in LAMMPS
        columns = ['typ', 'ai', 'aj', 'ak', 'ah', 'cmt', 'name', 'type_name']
        df: pd.DataFrame  # Dihedrals df to write out
        df = pd.DataFrame(columns=columns)
        df['typ'] = types
        df.index += 1   # Since the raw data increased one
        df['ai'] = raw_data.Dihedrals_df['ai']
        df['aj'] = raw_data.Dihedrals_df['aj']
        df['ak'] = raw_data.Dihedrals_df['ak']
        df['ah'] = raw_data.Dihedrals_df['ah']
        df['cmt'] = ['#' for _ in df.index]
        df['name'] = dihedrals_atoms
        df['type_name'] = type_names
        return df

    def dihedrals_name(self,
                       raw_data: rdlmp.ReadData  # raw data to clean dihedrals
                       ) -> tuple[list[str], list[str]]:  # Names of dihedrals
        """return name of the dihedrals by making from Atoms_df name"""
        ai_d_name: list[str]  # 1st atoms names as thier type in the dihedrals
        aj_d_name: list[str]  # 2nd atoms names as thier type in the dihedrals
        ak_d_name: list[str]  # 3rd atoms names as thier type in the dihedrals
        ah_d_name: list[str]  # 4th atoms names as thier type in the dihedrals
        ai_name: list[str]  # 1st atoms names which involve in a dihedrals
        aj_name: list[str]  # 2nd atoms names which involve in a dihedrals
        ak_name: list[str]  # 3rd atoms names which involve in a dihedrals
        ah_name: list[str]  # 4th atoms names which involve in a dihedrals
        dihedrals: list[str]  # Dihedrals names
        dihedrals_atoms: list[str]  # Dihedrals names by atoms names
        ai_d_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['b_name'][item]
                for item in raw_data.Dihedrals_df['ai']
                  ]
        aj_d_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['b_name'][item]
                for item in raw_data.Dihedrals_df['aj']
                  ]
        ak_d_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['b_name'][item]
                for item in raw_data.Dihedrals_df['ak']
                  ]
        ah_d_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['b_name'][item]
                for item in raw_data.Dihedrals_df['ah']
                  ]
        ai_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['name'][item]
                for item in raw_data.Dihedrals_df['ai']
                  ]
        aj_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['name'][item]
                for item in raw_data.Dihedrals_df['aj']
                  ]
        ak_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['name'][item]
                for item in raw_data.Dihedrals_df['ak']
                  ]
        ah_name = [
                raw_data.Atoms_df.loc[
                    raw_data.Atoms_df['atom_id'] == item]['name'][item]
                for item in raw_data.Dihedrals_df['ah']
                  ]
        dihedrals = [f'{i}_{j}_{k}_{h}' for i, j, k, h in zip(ai_d_name,
                                                              aj_d_name,
                                                              ak_d_name,
                                                              ah_d_name)]
        dihedrals_atoms = [f'{i}_{j}_{k}_{h}' for i, j, k, h in zip(ai_name,
                                                                    aj_name,
                                                                    ak_name,
                                                                    ah_name)]
        return dihedrals, dihedrals_atoms


class Types:
    """class to return right names and types"""
    def __init__(self,
                 names: list[str],    # Names to get the types based on them
                 dihedrals: bool = False  # Flag to correct the dihedrals duple
                 ) -> None:
        self.types: list[int]  # types to set to the system
        self.types_name: list[str]  # Names of the types NOT names
        self.types, self.types_name = self.get_types(names, dihedrals)

    def get_types(self,
                  names: list[str],  # Names of the bonds
                  dihedrals: bool = False  # Flag to correct the dihedrals dups
                  ) -> tuple[list[str], list[str]]:  # Type of the angles
        """make a correct type for angles
        Angle ABC is same as CBA
        Fixing this is complicated since it entirely depends on the
        names of the type of atoms in the Mass section.
        """
        # Remove digits from the name of the atoms:
        tmp_name: list[str]  # name of the bonds without digits
        tmp_name = [re.sub('\d+', '', item) for item in names]
        tmp_lst: list[list[str]] = [item.split('_') for item in tmp_name]
        # uniqe names after removing digits:
        name_set: list[typing.Any]  # set[str]
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
        if dihedrals:
            name_lst = [sorted(item) for item in name_lst]
        # Make one list of all the angles with thier atom types
        name_int: list[str]  # make a str list of int of types
        name_int = ['_'.join(item) for item in name_lst]
        # Get the uniqe type of all angles
        seen: set[frozenset[str]] = set()
        name_tup = [tuple([item]) for item in name_int]
        # name_tup = [item[0] for item in name_tup]
        # Make a list of each set of angles
        t: list[typing.Any] = [
            x for x in name_tup if frozenset(x) not in seen and
            not seen.add(frozenset(x))
            ]
        t = [list(item)[0] for item in t]  # make list of items
        # remove duplicates and reverse duplicated
        t = list({i[::-1] if i[-1] < i[0] else i: i for i in t}.values())
        # Give a name to each uniqe set of names
        angle_dict: dict[str, int]
        angle_dict = {item: v+1 for v, item in enumerate(t)}
        # List of types for each angle
        final_types: list[str] = []  # final list for each types
        for item in name_int:
            typ: int  # index for each angle
            try:
                typ = angle_dict[item]
            except KeyError:
                typ = angle_dict[item[::-1]]
            final_types.append(typ)
        name_dict: dict[int, str] = {}  # names based on the types
        for k, v in angle_dict.items():
            n = []
            for i in k.split('_'):
                value = self.get_key(type_dict, int(i))
                n.append(value)
            name_dict[v] = f'({"_".join(n)})'
        type_name: list[str] = []  # main name of each set in names
        for item in final_types:
            type_name.append(str(name_dict[item]))
        return final_types, type_name

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


class CleanData(Bonds,  # To get Bonds_df to write into files
                Angles,  # To get Angles_df to write into files
                Dihedrals  # To get Dihedrals_df to write into files
                ):
    """read data and clean it"""
    def __init__(self,
                 fname: str  # Name of the input file
                 ) -> None:
        print(f'{bcolors.OKCYAN}{self.__class__.__name__}:\n'
              f'\tCleaning: `{fname}`\n')
        raw_data = rdlmp.ReadData(fname)
        Bonds.__init__(self, raw_data)
        Angles.__init__(self, raw_data)
        Dihedrals.__init__(self, raw_data)
        self.clean_data(raw_data)

    def clean_data(self,
                   raw_data  # Data read by read_lmp_data scripts
                   ) -> None:
        """call all the methods, the other DataFrame are set with
        childer classes"""
        self.Atoms_df: pd.DataFrame = raw_data.Atoms_df
        self.Masses_df: pd.DataFrame = raw_data.Masses_df


if __name__ == '__main__':
    fname = sys.argv[1]
    data = CleanData(fname)
