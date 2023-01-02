import sys
import pandas as pd
import typing
import read_lmp_data as relmp
import lmp_to_pdb as lmpdb
import lmp_to_itp as lmpitp
from colors_text import TextColor as bcolors


class Doc:
    """read data file in LAMMPS format and convert it to pdb and itp
    files for GROMACS
    The input is the LAMMPS data file:
        Name of the atoms in Masses section should be the ones in the
        force filed file in GROMACS
    """


class ReadLmp:
    """read lammps from read_lmp_data.py"""
    def __init__(self,
                 fname: str  # Input file name
                 ) -> None:
        self.get_lmp(fname)

    def get_lmp(self,
                fname: str  # Input file name
                ) -> None:
        """get the data"""
        self.lmp_data = relmp.ReadData(fname)


class WritePdb:
    """write pdb file from Pdb class in lmpdb"""
    def __init__(self,
                 pdb_df: pd.DataFrame,  # df in pdb format
                 fname: str  # Input file name of LAMMPS data
                 ) -> None:
        self.write_pdb(pdb_df, fname)

    def write_pdb(self,
                  pdb_df: pd.DataFrame,  # df in pdb format
                  fname: str  # Input file name of LAMMPS data
                  ) -> None:
        """write the dataframe into a file"""
        fout: str  # Name of the output file
        fout = rename_file(fname, extension='pdb')
        print(f'{bcolors.OKBLUE}{self.__class__.__name__}:\n'
              f'\tPDB file is `{fout}`{bcolors.ENDC}\n')
        with open(fout, 'w') as f:
            f.write(f'HEADER\n')
            for row in pdb_df.iterrows():
                line: list[str]  # line with length of pdb line fill by spaces
                line = [' '*79]
                line[0:6] = f'{row[1]["records"]:<6s}'
                line[6:11] = f'{row[1]["atom_id"]:>5d}'
                line[11:12] = f' '
                line[12:16] = f'{row[1]["atom_name"]:<4s}'
                line[16:17] = f' '
                line[17:20] = f'{row[1]["residue_name"]:<3s}'
                line[20:22] = f'{" "*2}'
                line[22:26] = f'{row[1]["residue_id"]:>4d}'
                line[26:27] = f' '
                line[27:30] = f'{" "*3}'
                line[30:38] = f'{row[1]["x"]:>8.3f}'
                line[38:46] = f'{row[1]["y"]:>8.3f}'
                line[46:54] = f'{row[1]["z"]:>8.3f}'
                line[54:60] = f'{row[1]["occupancy"]:>6.2f}'
                line[60:66] = f'{row[1]["temperature"]:>6s}'
                line[66:72] = f'{" "*6}'
                line[72:76] = f'{row[1]["Segment_id"]:<4s}'
                line[76:78] = f'{row[1]["element"]:>2s}'
                line[78:] = f'{row[1]["charge"]:2s}'
                f.write(''.join(line))
                f.write(f'\n')
            f.write(f'END\n')


def rename_file(fname: str,  # Input file name
                extension: str  # The extension of the output file
                ) -> str:  # Out put file name
    """rename file name, same name with pdb extension"""
    fout: str  # Output file name
    fout = f'{fname.strip().split(".")[0]}.{extension}'
    return fout


class WriteItp:
    """write itp file
    There is no uniqe structure one should follow
    The columns will be seperated by single space"""
    def __init__(self,
                 itp: lmpitp.Itp,  # Data frames restructerd from LAMMPS
                 fname: str  # Name of the input files
                 ) -> None:
        """call functions"""
        self.write_itp(itp, fname)

    def write_itp(self,
                  itp: lmpitp.Itp,  # Data frames restructerd from LAMMPS
                  fname: str  # Name of the input data file
                  ) -> None:
        fout: str  # Name of the input file
        fout = rename_file(fname, extension='itp')
        print(f'{bcolors.OKBLUE}{self.__class__.__name__}:\n'
              f'\tITP file is `{fout}`{bcolors.ENDC}\n')
        with open(fout, 'w') as f:
            f.write(f'; input pdb SMILES:\n')
            f.write(f'\n')
            self.write_molecule(f)
            self.write_atoms(f, itp.atoms)
            self.write_bonds(f, itp.bonds)
            self.write_angles(f, itp.angles)
            self.write_dihedrals(f, itp.dihedrals)
            self.write_pairs(f)

    def write_molecule(self,
                       f: typing.Any  # The out put file
                       ) -> None:
        """write section of the itp file"""
        f.write(f'[ moleculetype ]\n')
        f.write(f'; Name\t\tnrexcl\n')
        f.write(f'ITP\t\t3\n')
        f.write(f'\n')

    def write_atoms(self,
                    f: typing.Any,  # The out put file
                    atoms: pd.DataFrame  # Atoms information
                    ) -> None:
        """write section of the itp file"""
        header: list[str] = [item for item in atoms.columns]
        f.write(f'[ atoms ]\n')
        f.write(f'; {" ".join(header)}\n')
        atoms.to_csv(f, header=None, sep=' ', index=False)
        f.write(f'\n')

    def write_bonds(self,
                    f: typing.Any,  # The out put file
                    bonds: pd.DataFrame  # bonds information
                    ) -> None:
        """write section of the itp file"""

    def write_angles(self,
                     f: typing.Any,  # The out put file
                     bonds: pd.DataFrame  # Bonds inoformation
                     ) -> None:
        """write section of the itp file"""

    def write_dihedrals(self,
                        f: typing.Any,  # The out put file
                        angles: pd.DataFrame  # Angles inoformation
                        ) -> None:
        """write section of the itp file"""

    def write_pairs(self,
                    f: typing.Any  # The out put file
                    ) -> None:
        """write section of the itp file"""


if __name__ == '__main__':
    fname: str = sys.argv[1]  # Input file name
    lmp: relmp.ReadData = relmp.ReadData(fname)  # All data in input file
    pdb = lmpdb.Pdb(lmp.Masses_df, lmp.Atoms_df)
    pdb_w = WritePdb(pdb.pdb_df, fname)
    itp = lmpitp.Itp(lmp, pdb.pdb_df)
    itp_w = WriteItp(itp, fname)
