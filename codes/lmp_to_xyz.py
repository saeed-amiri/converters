from ast import Index
import sys
import pandas as pd
import read_lmp_data as rdlmp


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
        data = rdlmp.ReadData(fname)
        df: pd.DataFrame = self.to_xyz(data.Atoms_df)  # Atoms in xyz format
        print(df)
        self.write_xyz(fname, df)

    def to_xyz(self,
               atoms: pd.DataFrame  # Atoms df
               ) -> pd.DataFrame:  # In xyz format
        """convert the LAMMPS atoms to xyz format"""
        columns: list[str]  # Columns for xyz file
        columns = ['name', 'x', 'y', 'z']
        df = pd.DataFrame(columns=columns)
        df['name'] = atoms['b_name']
        df['x'] = atoms['x']
        df['y'] = atoms['y']
        df['z'] = atoms['z']
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
            df.to_csv(f, sep=' ', index=False, header=None)


if __name__ == '__main__':
    xyz = XYZ(sys.argv[1])
