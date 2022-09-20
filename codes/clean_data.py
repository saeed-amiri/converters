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


# A helper function needed by most of the classes to get types for LAMMPS
def get_type(lst: list[str]  # list to get the number of distenguished ones
             ) -> list[int]:  # types' index
    """make type based on the unique items in the lst"""
    type_set: set[str] = set(lst)  # eleminate the repeated names
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
        raw_data: rdlmp.ReadData  # Raw data read by rdlmp
        self.raw_data = rdlmp.ReadData(fname)
        self.clean_data()

    def clean_data(self) -> None:
        """call all the methods"""
        self.Bonds_df: pd.DataFrame = self.get_bonds()  # Bonds df to write

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
        df['ai'] = self.raw_data.Bonds_df['ai']
        df['aj'] = self.raw_data.Bonds_df['aj']
        df['name'] = names
        df['cmt'] = ['#' for _ in df.index]
        df.index += 1
        return df

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
