import sys
import typing
import numpy as np
import pandas as pd
import lmp_to_pdb as lmpdb
import lmp_to_itp as lmpitp
import read_lmp_data as relmp
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
        moles: list[str]  # Names of each mol to make files
        moles = set(itp.atoms['resname'])
        for mol in moles:
            fout = rename_file(mol, 'itp')
            fout: str  # Name of the input file
            # fout = rename_file(fname, extension='itp')
            print(f'{bcolors.OKBLUE}{self.__class__.__name__}:\n'
                  f'\tITP file is `{fout}`{bcolors.ENDC}\n')
            with open(fout, 'w') as f:
                f.write(f'; input pdb SMILES:\n')
                f.write(f'\n')
                self.write_molecule(f, mol)
                self.write_atoms(f, itp.atoms, mol)
                self.write_bonds(f, itp.bonds, mol)
                self.write_angles(f, itp.angles, mol)
                self.write_dihedrals(f, itp.dihedrals, mol)
                self.write_pairs(f)

    def write_molecule(self,
                       f: typing.Any,  # The out put file
                       mol: str  # Name of the molecule to write into file
                       ) -> None:
        """write section of the itp file"""
        f.write(f'[ moleculetype ]\n')
        f.write(f'; Name\t\tnrexcl\n')
        f.write(f'{mol:>6s}\t\t3\n')
        f.write(f'\n')

    def write_atoms(self,
                    f: typing.Any,  # The out put file
                    atoms: pd.DataFrame,  # Atoms information
                    mol: str  # Name of the molecule to write into file
                    ) -> None:
        """write section of the itp file"""
        header: list[str] = [item for item in atoms.columns]
        f.write(f'[ atoms ]\n')
        f.write(f'; {"  ".join(header)}\n')

        df_raw: pd.DataFrame  # Copy of the df with mol selected info
        # df_raw = atoms[atoms['resname'] == mol]
        df_raw = atoms.copy()
        resides_ids = set(df_raw['resnr'])
        df1: pd.DataFrame  # Copy of the df with mol_id selected info
        df1 = df_raw[df_raw['resnr'] == list(resides_ids)[0]]
        df1 = atoms.copy()
        self.__atoms_one: dict[int, int]  # Atoms index from one
        self.__atoms_one = {nr: nr-np.min(df1['atomnr']) + 1
                            for nr in df1['atomnr']}
        atomnr: list[int] = [self.__atoms_one[item] for item in df1['atomnr']]
        resnr: list[int] = [1 for _ in self.__atoms_one]
        df = pd.DataFrame({'atomnr': atomnr,
                           'atomtype': df1['atomtype'],
                           'resnr': resnr,
                           'resname': df1['resname'],
                           'atomname': df1['atomname'],
                           'chargegrp': df1['chargegrp'],
                           'charge': df1['charge'],
                           'mass': df1['mass'],
                           ' ': df1[' '],
                           'element': df1['element']
                           })
        resides_ids: set[int]  # all of the residues ids
        for row in df.iterrows():
            line: list[str]  # line with length of 85 spaces to fix output
            line = [' '*85]
            line[0:7] = f'{row[1]["atomnr"]:>7d}'
            line[7:11] = f'{" "*2}'
            line[11:19] = f'{row[1]["atomtype"]:>7s}'
            line[19:21] = f'{" "*2}'
            line[21:26] = f'{row[1]["resnr"]:5d}'
            line[26:28] = f'{" "*2}'
            line[28:35] = f'{row[1]["resname"]:>7s}'
            line[35:37] = f'{" "*2}'
            line[37:45] = f'{row[1]["atomname"]:>8s}'
            line[45:47] = f'{" "*2}'
            line[47:56] = f'{row[1]["chargegrp"]:>9d}'
            line[56:58] = f'{" "*2}'
            line[58:64] = f'{row[1]["charge"]:>6.3f}'
            line[64:66] = f'{" "*2}'
            line[66:73] = f'{row[1]["mass"]:>6.3f}'
            line[73:74] = f'{" "*1}'
            line[75:77] = f'{row[1][" "]:>2s}'
            line[77:78] = f'{" "*1}'
            line[78:] = f'{row[1]["element"]:>6s}'
            f.write(''.join(line))
            f.write(f'\n')
        f.write(f'; Total charge : {df["charge"].sum()}\n')
        f.write(f'\n')

    def write_bonds(self,
                    f: typing.Any,  # The out put file
                    bonds: pd.DataFrame,  # bonds information
                    mol: str  # Name of the molecule to write into file
                    ) -> None:
        """write section of the itp file"""
        df_raw: pd.DataFrame  # Copy of the df with mol selected info
        df_raw = bonds.copy()
        resides_ids = set(df_raw['resnr'])
        df_list: list[pd.DataFrame] = []  # For keeping all bonds
        if resides_ids:
            for id in resides_ids:
                df1: pd.DataFrame  # Copy of the df with mol_id selected info
                df1 = df_raw[df_raw['resnr'] == id]
                ai = [self.__atoms_one[item] for item in df1['ai']]
                aj = [self.__atoms_one[item] for item in df1['aj']]
                df_i = pd.DataFrame({'ai': ai,
                                     'aj': aj,
                                     'funct': df1['funct'],
                                     'r': df1['r'],
                                     'k': df1['k'],
                                     ' ': df1[' '],
                                     '  ': df1['  ']})
                df_list.append(df_i)
                del df_i
            df = pd.concat(df_list)
            header: list[str] = [item for item in df.columns]
            f.write(f'[ bonds ]\n')
            f.write(f'; {" ".join(header)}\n')
            df.to_csv(f, header=None, sep='\t', index=False)
            f.write(f'\n')

    def write_angles(self,
                     f: typing.Any,  # The out put file
                     angles: pd.DataFrame,  # Angles inoformation
                     mol: str  # Name of the molecule to write into file
                     ) -> None:
        """write section of the itp file"""
        df_raw: pd.DataFrame  # Copy of the df with mol selected info
        df_raw = angles.copy()
        resides_ids = set(df_raw['resnr'])
        df_list: list[pd.DataFrame] = []  # Keeping all the angles
        if resides_ids:
            for id in resides_ids:
                df1: pd.DataFrame  # Copy of the df with mol_id selected info
                df1 = df_raw[df_raw['resnr'] == id]
                ai = [self.__atoms_one[item] for item in df1['ai']]
                aj = [self.__atoms_one[item] for item in df1['aj']]
                ak = [self.__atoms_one[item] for item in df1['ak']]

                df_i = pd.DataFrame({'ai': ai,
                                     'aj': aj,
                                     'ak': ak,
                                     'funct': df1['funct'],
                                     'theta': df1['theta'],
                                     'cth': df1['cth'],
                                     ' ': df1[' '],
                                     'angle_name': df1['angle_name']})
                df_list.append(df_i)
                del df_i
            df: pd.DataFrame = pd.concat(df_list)
            header: list[str] = [item for item in df.columns]
            f.write(f'[ angles ]\n')
            f.write(f'; {" ".join(header)}\n')
            df.to_csv(f, header=None, sep='\t', index=False)
            f.write(f'\n')

    def write_dihedrals(self,
                        f: typing.Any,  # The out put file
                        dihedrals: pd.DataFrame,  # Dihedrals inoformation
                        mol: str  # Name of the molecule tto write into file
                        ) -> None:
        """write section of the itp file"""
        df_raw: pd.DataFrame  # Copy of the df with mol selected info
        df_raw = dihedrals.copy()
        resides_ids = set(df_raw['resnr'])
        df_list: list[pd.DataFrame] = []  # Keeping all the angles
        if resides_ids:
            for id in resides_ids:
                df1: pd.DataFrame  # Copy of the df with mol_id selected info
                df1 = df_raw[df_raw['resnr'] == id]
                ai = [self.__atoms_one[item] for item in df1['ai']]
                aj = [self.__atoms_one[item] for item in df1['aj']]
                ak = [self.__atoms_one[item] for item in df1['ak']]
                ah = [self.__atoms_one[item] for item in df1['ah']]
                df_i = pd.DataFrame({'ai': ai,
                                     'aj': aj,
                                     'ak': ak,
                                     'ah': ah,
                                     'funct': df1['funct'],
                                     'C0': df1['C0'],
                                     'C1': df1['C1'],
                                     'C2': df1['C2'],
                                     'C3': df1['C3'],
                                     'C4': df1['C4'],
                                     ' ': df1[' '],
                                     'dihedral_name': df1['dihedral_name']
                                     })
                df_list.append(df_i)
                del df_i
            df: pd.DataFrame = pd.concat(df_list)
            header: list[str] = [item for item in df.columns]
            f.write(f'[ dihedrals ]\n')
            f.write(f'; {" ".join(header)}\n')
            df.to_csv(f, header=None, sep='\t', index=False)
            f.write(f'\n')

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
