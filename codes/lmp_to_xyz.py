import re
import sys
import csv
import pandas as pd
import read_lmp_data as rdlmp
from colors_text import TextColor as bcolors


class Doc:
    """convert LAMMPS data file to xyz file
    Input:
        LAMMPS data file (fname)
    Output:
        xyz file
    """


class XYZ:
    """read the data and convert it to xyz"""
    def __init__(self,
                 fname: str  # Name of the input file
                 ) -> None:
        print(f'{bcolors.OKCYAN}{self.__class__.__name__}:\n'
              f'\tConverting `{fname}` to XYZ file{bcolors.ENDC}\n')
        data = rdlmp.ReadData(fname)
        df: pd.DataFrame = self.to_xyz(data.Atoms_df)  # Atoms in xyz format
        self.write_xyz(fname, df)

    def to_xyz(self,
               atoms: pd.DataFrame  # Atoms df
               ) -> pd.DataFrame:  # In xyz format
        """convert the LAMMPS atoms to xyz format"""
        columns: list[str]  # Columns for xyz file
        columns = ['name', 'x', 'y', 'z']
        df = pd.DataFrame(columns=columns)
        df['x'] = atoms['x']
        df['y'] = atoms['y']
        df['z'] = atoms['z']
        df['name'] = self.clean_names(atoms['b_name'])
        df = df.astype({'name': str, 'x': float, 'y': float, 'z': float})
        return df

    def write_xyz(self,
                  fname: str,  # Input name to make output name
                  df: pd.DataFrame  # To write
                  ) -> None:
        """write the xyz file"""
        fout: str = f'{fname.split(".")[0]}.xyz'  # Output file name
        with open(fout, 'w') as f:
            f.write(f'{len(df)}\n')
            f.write(f'\n')
            df.to_csv(f, sep=' ', index=False, header=None,
                      quoting=csv.QUOTE_NONE)

    def clean_names(self,
                    names: list[str]  # Names of the atoms
                    ) -> list[str]:  # Cleaned names
        """remove special chars and return in capitalize first letter
        format"""
        lst = [item.capitalize() for item in names]
        lst = [re.sub(r'[^a-zA-Z]', '', s) for s in names]
        lst = [item.capitalize() for item in names]
        names = lst
        return names


if __name__ == '__main__':
    xyz = XYZ(sys.argv[1])
