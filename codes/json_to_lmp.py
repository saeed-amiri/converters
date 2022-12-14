from pprint import pprint
import sys
import json
import typing
import pandas as pd
import numpy as np
from itertools import combinations
import periodictabel_df as periodf
from colors_text import TextColor as bcolors


class Doc:
    """transfer the data file of atom structures in JSON format
    to lammps full atom format
    This script suppose to read data file form:
    https://pubchem.ncbi.nlm.nih.gov/compound/5974#section=Structures
    The strcture is in the following format:
    {
        "PC_Compounds": [
          {
          "id": {
            "id": {
              "cid": 2681
            }
          "atoms": {
            }
            "aid": [
            ],
            "element": [
            ],
            "charge": [
            ]
          },
          "bonds": {
          },
          "coords": [
          ],
          "props": [
          ],
          "count": {
          }
        }
      ]
    }
    """


class ReadJson:
    """get the json file to a dictionary"""
    def __init__(self,
                 fname: str  # Json file
                 ) -> None:
        print(f'\t{bcolors.OKCYAN}Get the Json file: '
              f'`{fname}`{bcolors.ENDC}\n')
        self.get_param(fname)

    def get_param(self, fname: str) -> None:
        with open(fname, 'r') as f:
            data = json.load(f)
        self.param = data
        del data


class Angle:
    """to guess angles for the input file"""
    def __init__(self,
                 bonds_df: pd.DataFrame,  # Bonds DF from ConvertJson
                 atoms_info: pd.DataFrame,  # All the atoms information
                 atoms_df: pd.DataFrame  # Coordinates of atoms
                 ) -> None:
        """make angles_df for LAMMPS data file"""
        self.angles_df: pd.DataFrame = self.mk_angels_df(bonds_df,
                                                         atoms_info,
                                                         atoms_df)

    def mk_angels_df(self,
                     bonds_df: pd.DataFrame,  # Bonds DF from ConvertJson
                     atoms_info: pd.DataFrame,  # All the atoms information
                     atoms_df: pd.DataFrame  # All the atoms coordinates
                     ) -> pd.DataFrame:
        """the dataframe for LAMMPS to write"""
        columns: list[str] = ['typ', 'ai', 'aj', 'ak', 'cmt', 'name']
        angle_df: pd.DataFrame = self.get_angle(bonds_df)  # df of just angles
        type_df: pd.DataFrame = self.get_types(angle_df, atoms_info)  # types
        types: list[int] = self.set_angle_type(type_df)  # Type of each angle
        angle_name: list[str] = self.set_angle_name(angle_df, atoms_info)
        df = pd.DataFrame(columns=columns)
        df['ai'] = angle_df['ai']
        df['aj'] = angle_df['aj']
        df['ak'] = angle_df['ak']
        df['typ'] = types
        df['cmt'] = ['#' for _ in df.index]
        df['name'] = angle_name
        df.index += 1
        self.get_radian(angle_df, atoms_df)
        return df

    def get_angle(self,
                  bonds_df: pd.DataFrame  # Bonds DF from ConvertJson
                  ) -> pd.DataFrame:
        """guess angles between atoms based on the bonds"""
        ai_set: set[int]  # Unrepeted atom index
        i_df: pd.DataFrame  # For each atom bonds
        i_angle: pd.DataFrame  # Angles for each atom
        angle_df_list: list[pd.DataFrame] = []  # list of the DataFrames
        columns: list[str] = ['ai', 'aj', 'ak']  # columns of angle DF
        i_angle = pd.DataFrame(columns=columns)
        ai_set = set(bonds_df['ai'])
        for ai in ai_set:
            i_df = bonds_df.loc[bonds_df['ai'] == ai]
            angles = self.mk_angles(ai, i_df)
            i_angle = pd.DataFrame(angles, columns=columns)
            angle_df_list.append(i_angle)
        angle_df = pd.concat(angle_df_list, ignore_index=True)
        return angle_df

    def mk_angles(self,
                  ai: int,  # The center atom to make angles
                  df: pd.DataFrame  # Slice of the bonds_df for each atom
                  ) -> list[list[int]]:
        """make a angles for each sub dataframe"""
        aj: list[int] = df['aj']  # The other atoms in the angle
        aj_comb: list[tuple[int, int]]  # All the combinations of the aj
        aj_comb = [[aj1, aj2] for aj1, aj2 in combinations(aj, 2)]
        angles: list[list[int]]  # list of angles
        angles = [[aj[0], ai, aj[1]] for aj in aj_comb]
        return angles

    def set_angle_type(self,
                       df: pd.DataFrame  # DataFrame of the types
                       ) -> list[int]:
        """make a type to distenguise the angle types"""
        types: list[tuple[int, int, int]]  # Types of atom in each set
        types = [(i, j, k) for i, j, k in
                 zip(df['ai_type'], df['aj_type'], df['ak_type'])]
        types_unique = list(set(types))
        type_dict: dict[tuple[int, int, int], int]  # Indexing the types
        type_dict = {k: v+1 for v, k in enumerate(types_unique)}
        angle_type: list[int] = []  # Type of each angle
        angle_type = [type_dict[item] for item in types]
        return angle_type

    def set_angle_name(self,
                       angles: pd.DataFrame,  # DataFrame of the angles
                       atoms_info: pd.DataFrame  # All atoms informations
                       ) -> list[str]:
        """make angels name based on the atoms"""
        angle_name: list[str] = []  # Name for each angle
        for _, row in angles.iterrows():
            ai: int = row['ai']  # index of the atom
            aj: int = row['aj']  # index of the atom
            ak: int = row['ak']  # index of the atom
            i_name: int = atoms_info.loc[atoms_info['aid'] == ai]['name'][ai-1]
            j_name: int = atoms_info.loc[atoms_info['aid'] == aj]['name'][aj-1]
            k_name: int = atoms_info.loc[atoms_info['aid'] == ak]['name'][ak-1]
            a_name: str = f'{i_name}_{j_name}_{k_name}'  # Name of the angle
            angle_name.append(a_name)
        return angle_name

    def get_types(self,
                  angels: pd.DataFrame,  # Index of atoms share an angles
                  atoms_info: pd.DataFrame  # All atoms information
                  ) -> pd.DataFrame:
        """retrun the type of atoms"""
        ai_type: list[int] = []  # type of the ai atom
        aj_type: list[int] = []  # type of the aj atom
        ak_type: list[int] = []  # type of the ak atom
        for _, row in angels.iterrows():
            ai: int = row['ai']  # index of the atom
            aj: int = row['aj']  # index of the atom
            ak: int = row['ak']  # index of the atom
            i_type: int = atoms_info.loc[atoms_info['aid'] == ai]['type'][ai-1]
            j_type: int = atoms_info.loc[atoms_info['aid'] == aj]['type'][aj-1]
            k_type: int = atoms_info.loc[atoms_info['aid'] == ak]['type'][ak-1]
            ai_type.append(i_type)
            aj_type.append(j_type)
            ak_type.append(k_type)
        columns: list[str]  # Columns of the types DataFrame
        type_df: pd.DataFrame  # To return types
        columns = ['ai_type', 'aj_type', 'ak_type']
        type_df = pd.DataFrame(columns=columns)
        type_df['ai_type'] = ai_type
        type_df['aj_type'] = aj_type
        type_df['ak_type'] = ak_type
        return type_df

    def get_radian(self,
                   angles: pd.DataFrame,  # Index of atoms share an angles
                   atoms: pd.DataFrame  # Atoms' coordinates
                   ) -> list[float]:
        """calculate the angle of each set"""
        radians: list[float] = []  # Angles between atoms to return
        for _, row in angles.iterrows():
            ai: int = row['ai']  # index 1st atoms in the angles
            aj: int = row['aj']  # index 2nd atoms in the angles
            ak: int = row['ak']  # index 3rd atoms in the angles
            ai_coords = atoms.loc[  # Coord of the 1st atom
                atoms['atom_id'] == ai][['x', 'y', 'z']].to_numpy()[0]
            aj_coords = atoms.loc[  # Coord of the 2nd atom
                atoms['atom_id'] == aj][['x', 'y', 'z']].to_numpy()[0]
            ak_coords = atoms.loc[  # Coord of the 3rd atom
                atoms['atom_id'] == ak][['x', 'y', 'z']].to_numpy()[0]
            radians.append(
                self.calculate_radian(ai_coords, aj_coords, ak_coords)
                )
        return radians

    def calculate_radian(self,
                         a: np.array,  # Coords of 1st atom in the angle
                         b: np.array,  # Coords of 2nd atom in the angle
                         c: np.array,  # Coords of 3rd atom in the angle
                         ):
        """return the angle of each set in radian"""
        ba = a - b
        bc = c - b
        cosine_angle =\
            np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        return angle


class ConvertJson(ReadJson,  # Read the main data file for atoms and bonds
                  periodf.PeriodicTable  # Periodic table for all elements
                  ):
    """read the json files"""
    def __init__(self,
                 fname: str  # Name of the input files
                 ) -> None:
        ReadJson.__init__(self, fname)
        periodf.PeriodicTable.__init__(self)
        self.compounds: dict[str, list[typing.Any]]  # Needed values
        self.compounds = self.param['PC_Compounds'][0]
        self.atom_info: pd.DataFrame = self.get_atom_info()
        self.get_atoms()
        self.get_bonds()
        self.get_angles()
        self.Masses_df: pd.DataFrame = self.mk_masses()  # Masses info

    def get_atoms(self) -> None:
        """get all the atoms coords and return a lammps version
        of full atom style"""
        coords_df: pd.DataFrame = self.get_atoms_coords()  # xyz of atoms
        element_df: pd.DataFrame = self.get_element()  # Atomic numbers
        self.Atoms_df: pd.DataFrame = self.mk_atom_df(coords_df)

    def get_bonds(self) -> None:
        """get all the bonds types and atoms ids"""
        self.Bonds_df: pd.DataFrame = self.get_bonds_df()

    def mk_masses(self) -> pd.DataFrame:
        """make masses dataframe for the datafile"""
        columns: list[str] = ['typ', 'mass', 'cmt', 'name']
        df = self.atom_info.copy()
        df = df.loc[df.groupby(by=['element'])['mass'].idxmin(), :]
        df.reset_index(inplace=True)
        df.drop(['aid', 'index'], axis=1, inplace=True)
        masses_df = pd.DataFrame(columns=columns)
        masses_df['typ'] = df['type']
        masses_df['mass'] = df['mass']
        masses_df['cmt'] = ['#' for _ in df.index]
        masses_df['name'] = df['name']
        return masses_df

    def get_bonds_df(self) -> pd.DataFrame:
        """make bonds dataframe"""
        bonds: dict[str, list[int]] = self.compounds['bonds']  # Bonds ids
        columns = ['typ', 'ai', 'aj', 'cmt', 'name']
        bonds_df = pd.DataFrame(columns=columns)
        bonds_df['typ'] = self.get_bonds_type(bonds['aid1'], bonds['aid2'])
        bonds_df['ai'] = bonds['aid1']
        bonds_df['aj'] = bonds['aid2']
        bonds_df['cmt'] = ['#' for _ in bonds_df.index]
        bonds_df['name'] = self.get_bonds_name(bonds['aid1'], bonds['aid2'])
        bonds_df.index += 1
        return bonds_df

    def get_bonds_type(self,
                       ai: list[int],  # Id of the 1st atoms in bonds
                       aj: list[int],  # Id of the 2nd atoms in bonds
                       ) -> list[int]:
        """get bonds type, since there is no infos about it in data
        file"""
        ai_type: list[int]  # type of 1st atoms in bonds
        aj_type: list[int]  # type of 2nd atoms in bonds
        bond_couple: list[tuple[int, int]]  # Bond couples
        bond_dict: dict[tuple[int, int], int]  # type of each couple
        bonds_type: list[int]  # type of each bond
        ai_type = [self.atom_info.loc[self.atom_info['aid'] == item]
                   ['type'][item-1] for item in ai]
        aj_type = [self.atom_info.loc[self.atom_info['aid'] == item]
                   ['type'][item-1] for item in aj]
        bond_couple = [(i, j) for i, j in zip(ai_type, aj_type)]
        bond_dict = {t: i+1 for i, t in enumerate(set(bond_couple))}
        bonds_type = [bond_dict[item] for item in bond_couple]
        return bonds_type

    def get_bonds_name(self,
                       ai: list[int],  # Id of the 1st atoms in bonds
                       aj: list[int],  # Id of the 2nd atoms in bonds
                       ) -> list[str]:
        """return a list of atoms atomics number since there is no
        name yet"""
        ai_name: list[str]  # Name of the atoms
        aj_name: list[str]  # Name of the atoms
        ai_name = [self.atom_info.loc[self.atom_info['aid'] == item]
                   ['name'][item-1] for item in ai]
        aj_name = [self.atom_info.loc[self.atom_info['aid'] == item]
                   ['name'][item-1] for item in aj]
        name: list[str] = [f'{i}_{j}' for i, j in zip(ai_name, aj_name)]
        return name

    def mk_atom_df(self,
                   coords_df: pd.DataFrame,  # xyz of atoms
                   ) -> pd.DataFrame:
        """return atoms coordinates in LAMMPS style"""
        mol: list[int]  # index for the moelcule, HERE 1!
        nxyz: list[int]  # nx, ny, nz flags of the full style
        charges: list[float]  # Charges for each atom
        cmt: list[str]  # Comment for the name of the atoms
        mol = [int(1) for _ in coords_df.index]
        nxyz = [int(0) for _ in coords_df.index]
        cmt = ['#' for _ in coords_df.index]
        charges = [float(0) for _ in coords_df.index]
        columns = ['atom_id', 'mol', 'typ', 'charge', 'x', 'y', 'z',
                   'nx', 'ny', 'nz', 'cmt', 'name']
        atoms_df: pd.DataFrame  # LAMMPS dataframe
        atoms_df = pd.DataFrame(columns=columns)
        atoms_df['atom_id'] = coords_df['aid']
        atoms_df['mol'] = mol
        atoms_df['typ'] = self.atom_info['type']
        atoms_df['charge'] = charges
        atoms_df['x'] = coords_df['x']
        atoms_df['y'] = coords_df['y']
        atoms_df['z'] = coords_df['z']
        atoms_df['nx'] = nxyz
        atoms_df['ny'] = nxyz
        atoms_df['nz'] = nxyz
        atoms_df['cmt'] = cmt
        atoms_df['name'] = self.atom_info['name']
        return atoms_df

    def get_element(self) -> pd.DataFrame:
        """return atoms chemical elemnt numbers"""
        atoms: dict[str, list[typing.Any]]  # contain the elemnt information
        atoms = self.compounds['atoms']
        aid: list[int]  # id of all the atoms
        a_element: list[int]  # atomic number of each elemcnt
        aid = atoms['aid']
        a_element = atoms['element']
        type_list: list[int] = self.mk_types(a_element)
        element_df = pd.DataFrame(list(zip(aid, type_list, a_element)),
                                  columns=['aid', 'typ', 'element'])
        return element_df

    def mk_types(self,
                 a_element: list[int]  # Atomic number of atoms
                 ) -> list[int]:
        """return a list for atomic types"""
        set_element: set[int] = set(a_element)
        type_dict: dict[int, int] = {k: i+1 for i, k in enumerate(set_element)}
        type_list: list[int] = [type_dict[k] for k in a_element]
        return type_list

    def get_atoms_coords(self) -> pd.DataFrame:
        """get the id of all the atoms"""
        aid: list[int]  # id of all the atoms
        x: list[float]  # x component for all the atoms
        y: list[float]  # y component for all the atoms
        z: list[float]  # z component for all the atoms
        coords: dict[str, list[float]]
        coords = self.compounds['coords'][0]['conformers'][0]
        x = coords['x']
        y = coords['y']
        z = coords['z']
        aid = self.compounds['coords'][0]['aid']
        coords_df = pd.DataFrame(list(zip(aid, x, y, z)),
                                 columns=['aid', 'x', 'y', 'z'])
        return coords_df  # xyz of atoms

    def get_atom_info(self) -> pd.DataFrame:
        """return the name, atomic nnumber and mass of each atom"""
        df: pd.DataFrame  # A dataframe for index and element id
        columns: list[str]  # Name of the columns
        name: list[str] = []  # Name of each element
        mass: list[float] = []  # Mass of each element
        atoms: list[int]  # Index of atoms (atomic number)
        columns = ['aid', 'type', 'element', 'name', 'mass']  # Columns name
        df = pd.DataFrame(columns=columns)
        df['aid'] = self.compounds['atoms']['aid']
        df['element'] = self.compounds['atoms']['element']
        for i, aid in enumerate(df['aid']):
            element = df.iloc[i]['element']
            iloc = self.peridic_table.loc[self.peridic_table['number'] ==
                                          element]
            i_name = iloc['symbol'][element]
            i_mass = iloc['atomic_mass'][element]
            name.append(i_name)
            mass.append(i_mass)
        atoms = self.compounds['atoms']['element']
        type_list: list[int] = self.mk_types(atoms)
        df['name'] = name
        df['mass'] = mass
        df['type'] = type_list
        return df

    def get_angles(self) -> None:
        """call class Angles to find the angles between particles"""
        angle = Angle(self.Bonds_df, self.atom_info, self.Atoms_df)
        self.Angles_df = angle.angles_df


if __name__ == '__main__':
    fname = sys.argv[1]
    fjson = ConvertJson(fname)
