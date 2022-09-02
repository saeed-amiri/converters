from pprint import pprint
import sys
import json
import typing
import pandas as pd
from colors_text import TextColor as bcolors
import periodictabel_df as pridf


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


class ConvertJson(ReadJson,
                  pridf.PeriodicTable):
    """read the json files"""
    def __init__(self,
                 fname: str  # Name of the input files
                 ) -> None:
        ReadJson.__init__(self, fname)
        pridf.PeriodicTable.__init__(self)
        self.compounds: dict[str, list[typing.Any]]  # Needed values
        self.compounds = self.param['PC_Compounds'][0]
        self.atom_info: pd.DataFrame = self.get_atom_info()
        self.get_atoms()
        self.get_bonds()
        self.Masses_df: pd.DataFrame = self.mk_masses()  # Masses info

    def get_atoms(self) -> None:
        """get all the atoms coords and return a lammps version
        of full atom style"""
        coords_df: pd.DataFrame = self.get_atoms_coords()  # xyz of atoms
        element_df: pd.DataFrame = self.get_element()  # Atomic numbers
        self.Atoms_df: pd.DataFrame = self.mk_atom_df(coords_df, element_df)

    def get_bonds(self) -> None:
        """get all the bonds types and atoms ids"""
        self.Bonds_df: pd.DataFrame = self.get_bonds_df()

    def mk_masses(self) -> None:
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
        bonds_df['typ'] = bonds['order']
        bonds_df['ai'] = bonds['aid1']
        bonds_df['aj'] = bonds['aid2']
        bonds_df['cmt'] = ['#' for _ in bonds_df.index]
        bonds_df['name'] = self.get_bonds_name(bonds['aid1'], bonds['aid2'])
        return bonds_df

    def get_bonds_name(self,
                       ai: list[int],  # Id of the 1st atoms in bonds
                       aj: list[int],  # Id of the 2nd atoms in bonds
                       ) -> list[str]:
        """return a list of atoms atomics number since there is no
        name yet"""
        ai_name: list[str]  # Name of the atoms
        aj_name: list[str]  # Name of the atoms
        ai_name = [self.atom_info.loc[self.atom_info['aid']==item]
                   ['name'][item-1] for item in ai]
        aj_name = [self.atom_info.loc[self.atom_info['aid']==item]
                   ['name'][item-1] for item in aj]
        name: list[str] = [f'{i}_{j}' for i, j in zip(ai_name, aj_name)]
        return name

    def mk_atom_df(self,
                   coords_df: pd.DataFrame,  # xyz of atoms
                   element_df: pd.DataFrame  # Atomic numbers
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
        columns = ['aid', 'type', 'element', 'name', 'mass']
        df = pd.DataFrame(columns=columns)
        df['aid'] = self.compounds['atoms']['aid']
        df['element'] = self.compounds['atoms']['element']
        for i, aid in enumerate(df['aid']):
            element = df.iloc[i]['element']
            iloc = self.peridic_table.loc[self.peridic_table['number']==
                                          element]
            i_name = iloc['name'][element]
            i_mass = iloc['atomic_mass'][element]
            name.append(i_name)
            mass.append(i_mass)
        atoms = self.compounds['atoms']['element']
        type_list: list[int] = self.mk_types(atoms)
        df['name'] = name
        df['mass'] = mass
        df['type'] = type_list
        return df


if __name__ == '__main__':
    fname = sys.argv[1]
    json = ConvertJson(fname)
