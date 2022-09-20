import pandas as pd
import sys
import read_lmp_data as rdlmp


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


class CleanData:
    """read data and clean it"""
    def __init__(self,
                 fname: str  # Name of the input file
                 ) -> None:
        raw_data: rdlmp.ReadData  # Raw data read by rdlmp
        self.raw_data = rdlmp.ReadData(fname)
        self.clean_data()

    def clean_data(self) -> None:
        """call all the methods"""
        self.get_bonds()

    def get_bonds(self) -> None:
        """correct the name and type of the bonds"""
        name: list[str]  # Bonds names
        name = self.bonds_name()
        print(name)

    def bonds_name(self) -> list[str]:  # Name of the bonds
        """return name of the bonds by making from Atoms_df name"""
        ai_name: list[str]  # 1st atoms names
        aj_name: list[str]  # 2nd atoms names
        names: list[str]  # Bonds names
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


if __name__ == '__main__':
    fname = sys.argv[1]
    data = CleanData(fname)
