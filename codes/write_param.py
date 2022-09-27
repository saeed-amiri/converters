import re
import typing
import pandas as pd
import read_param as rdprm
import clean_data as cllmp


class Doc:
    """update the DataFrame from read_param by the type of each set
    from clean_data
    input:
        DataFrames from clean_data
        DataFrames from read_param
    Output:
        File in JSON with all the interactions that can be read by
        combination codes
    """


# Rteturn unique atoms
def seen_set(
             lst: list[typing.Any]  # to drop duplicate with keping order
             ) -> list[typing.Any]:
    """remove duplicated item with keeping order of them in the
    main list"""
    seen: set[str] = set()
    seen_add = seen.add
    return [x for x in lst if not (x in seen or seen_add(x))]


class GetType:
    """set type for bonds, angles, dihedrals"""
    def __init__(self,
                 data: cllmp.CleanData  # Cleaned data to write
                 ) -> None:
        param: rdprm.ReadParam  # All the data from param file
        param = rdprm.ReadParam(fname='param')
        self.set_types(param, data)

    def set_types(self,
                  param: rdprm.ReadParam,  # ForceField from param
                  data: cllmp.CleanData  # Cleaned data to write
                  ) -> None:
        """call all the methods to set types and write them"""
        self.atoms = self.set_atoms(param.atoms, data.Atoms_df)  # Atoms prm
        self.bonds = self.set_abd_types(
            param.bonds, data.Bonds_df, 'bond_name')  # Bonds prm
        self.angles = self.set_abd_types(
            param.angles, data.Angles_df, 'angle_name')  # Angles prm
        self.dihedrals = self.set_abd_types(
            param.dihedrals, data.Dihedrals_df, 'dihedral_name')  # Dihedralprm

    def set_atoms(self,
                  atoms_prm: pd.DataFrame,  # LJ parameters for atoms
                  df: pd.DataFrame  # Atoms information
                  ) -> pd.DataFrame:  # Updated atoms parameters
        """update parameters of atoms with thier types"""
        atoms_name: list[str]  # Names of the atoms != Names in bonds, etc
        atoms_name = list(df['name'])
        atoms_name = seen_set(atoms_name)
        atoms_type: dict[str, int] = {}  # Atoms name with thier types
        for item in atoms_name:
            typ: int = df.loc[df['name'] == item]['typ']  # Type of the atom
            atoms_type[item] = seen_set(typ)[0]
        i_type: list[int] = []  # Type of each atom
        for item in atoms_prm['atom_name']:
            i_type.append(atoms_type[item])
        atoms_prm['type'] = i_type
        return atoms_prm

    def set_abd_types(self,
                      abd_prm: pd.DataFrame,  # LJ parameters for bonds, etc
                      df: pd.DataFrame,  # Bonds, or Angles or Dihedrlas infos
                      char: str  # Column's in df
                      ) -> pd.DataFrame:  # Updated atoms prm
        """updated parameters of bonds, angles, dihedrlas with thier types"""
        abd_name: list[str]  # All the names
        abd_name = list(df['type_name'])
        abd_name = seen_set(abd_name)
        abd_type: dict[str, int] = {}  # Type of the each name
        for item in abd_name:
            typ: int = df.loc[df['type_name'] == item]['typ']  # Type of bonds
            abd_type[item] = seen_set(typ)[0]
        i_type: list[int] = []  # Type of each name
        for item in abd_prm[char]:
            try:
                i_type.append(abd_type[item])
            except KeyError:
                i_type.append(0)
        abd_prm['type'] = i_type
        return abd_prm

    def drop_parentheses(self,
                         lst: list[str]  # drop ()
                         ) -> list[str]:
        """drop parentheses from the names"""
        return [re.sub(r'[()]', '', item) for item in lst]


class WriteParam(GetType):
    """write parameters in json format"""
    def __init__(self,
                 data: cllmp.CleanData  # Read data to set types
                 ) -> None:
        super().__init__(data)
        self.write_params()

    def write_params(self) -> None:
        """write the parameters in json by calling methods"""
        with open('param.json', 'w') as f:
            self.write_atoms(f)
            self.write_bonds(f)
            self.write_angles(f)

    def write_atoms(self,
                    f: typing.TextIO  # To write into
                    ) -> None:
        """write atoms in the LJ section into the output file"""
        #   atom_name    mass   sigma epsilom charge  style  type
        f.write(f'"atoms": [\n')
        for i, row in self.atoms.iterrows():
            f.write(f'\t{{\n'
                    f'\t"type": {row["type"]}, '
                    f'"name": "{row["atom_name"]}", '
                    f'"sigma": {row["sigma"]},\t'
                    f'"epsilon": {row["epsilon"]}, '
                    f'"mass": {row["mass"]}, '
                    f'"charge": {row["charge"]}\n'
                    f'\t}},\n')
        f.write(f'  ],\n')  # 2 spaces before closing bracket

    def write_bonds(self,
                    f: typing.TextIO  # To write into
                    ) -> None:
        """write bonds section into the output file"""
        f.write(f'"bonds": [\n')
        for i, row in self.bonds.iterrows():
            f.write(f'\t{{\n'
                    f'\t"type": {row["type"]}, '
                    f'"name": "{row["bond_name"]}", '
                    f'"style": "{row["style"]}", '
                    f'"kbond": {row["kbond"]}, '
                    f'"r": {row["r"]}\n'
                    f'\t}},\n')

    def write_angles(self,
                    f: typing.TextIO  # To write into
                    ) -> None:
        """write angles section into the output file"""
        f.write(f'"angles": [\n')
        for i, row in self.angles.iterrows():
            f.write(f'\t{{\n'
                    f'\t"type": {row["type"]}, '
                    f'"name": "{row["angle_name"]}", '
                    f'"style": "{row["style"]}", '
                    f'"kangle": {row["kangle"]}, '
                    f'"angle": {row["angle"]}\n'
                    f'\t}},\n')
