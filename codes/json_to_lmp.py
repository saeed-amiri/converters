import sys
import json
import typing
import pandas as pd
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


class ConvertJson(ReadJson):
    """read the json files"""
    def __init__(self,
                 fname: str  # Name of the input files
                 ) -> None:
        super().__init__(fname)
        self.compounds: dict[str, list[typing.Any]]
        self.compounds = self.param['PC_Compounds'][0]
        self.get_atoms()

    def get_atoms(self) -> None:
        """get all the atoms coords and return a lammps version
        of full atom style"""
        coords_df: pd.DataFrame = self.get_atoms_coords()  # xyz of atoms

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


if __name__ == '__main__':
    fname=sys.argv[1]
    json = ConvertJson(fname)
