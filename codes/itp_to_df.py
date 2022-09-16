from pprint import pprint
import sys
import typing
import pandas as pd
from colors_text import TextColor as bcolors


class Doc:
    """reading itp files and return data several data frames
    Usually, the itp file contains information about atoms (not
    coordinates), bonds, angles, dihedrals, and improper dihedrals.

    informations in the itp file:
        [ moleculetype ] : defines the name of your molecule in this top
    and nrexcl = 3 stands for excluding non-bonded interactions between
    atoms that are no further than 3 bonds away.

        [ atoms ] : defines the molecule, where nr and type are fixed,
    the rest is user defined. So atom can be named as you like,
    cgnr made larger or smaller (if possible, the total charge of
    a charge group should be zero), and charges can be changed here
    too.

        [ bonds ] : no comment.

        [ pairs ] : LJ and Coulomb 1-4 interactions

        [ angles ] : no comment

        [ dihedrals ] : in this case there are 9 proper dihedrals
    (funct = 1), 3 improper (funct = 4) and no Ryckaert-Bellemans type
    dihedrals.
    """


# A helper function needed by most of the classes
def free_char_line(line: str  # line of the itp file
                   ) -> list[str]:  # Free from chars
    """cheack the lines and return the line free special chars"""
    char_list: list[str] = [';', '#', ':']  # chars to eliminate from lines
    l_line: list[str]  # Breaking the line cahrs
    l_line = line.strip().split(' ')
    l_line = [item for item in l_line if item]
    l_line = [item for item in l_line if item not in char_list]
    return l_line


class Itp:
    """read itp file and return a DataFrame of the information
    within the file"""
    def __init__(self,
                 fname: str  # Name of the itp file
                 ) -> None:
        self.get_itp(fname)

    def get_itp(self,
                fname: str  # Name of the itp file
                ) -> None:
        """read the file line by line and call related methods"""
        atoms: bool = False  # Flag of 'atoms' occurrence
        bonds: bool = False  # Flag of 'bonds' occurrence
        angles: bool = False  # Flag of 'angles' occurrence
        dihedrals: bool = False  # Flag of 'dihedrals' occurrence
        imporopers: bool = False  # Flag of 'imporopers' occurrence
        moleculetype: bool = False  # Flag of 'moleculetype' occurrence
        atoms_info: list[str] = []  # to append atoms lines
        bonds_info: list[str] = []  # to append bonds lines
        with open(fname, 'r') as f:
            while True:
                line: str = f.readline()
                if line.strip():
                    if line.strip().startswith('['):
                        wilds: list[str]  # Parts of the line
                        wilds = line.strip().split()
                        if wilds[1] == 'atoms':
                            atoms, bonds, angles, dihedrals, imporopers,\
                                moleculetype = True, False, False, False,\
                                False, False
                        elif wilds[1] == 'bonds':
                            atoms, bonds, angles, dihedrals, imporopers,\
                                moleculetype = False, True, False, False,\
                                False, False
                        elif wilds[1] == 'angles':
                            atoms, bonds, angles, dihedrals, imporopers,\
                                moleculetype = False, False, True, False,\
                                False, False
                        elif wilds[1] == 'dihedrals':
                            atoms, bonds, angles, dihedrals, imporopers,\
                                moleculetype = False, False, True, False,\
                                False, False
                        elif wilds[1] == 'moleculestype':
                            atoms, bonds, angles, dihedrals, imporopers,\
                                moleculetype = False, False, False, False,\
                                False, True
                        else:
                            atoms, bonds, angles, dihedrals, imporopers,\
                                moleculetype = False, False, False, False,\
                                False, False
                    else:
                        if atoms:
                            atoms_info.append(line)
                        if bonds:
                            bonds_info.append(line)
                if not line:
                    break
        atom = AtomsInfo(atoms_info)
        bond = BondsInfo(atoms=atom.df, bonds=bonds_info)


class AtomsInfo:
    """get atoms wild information and retrun a DataFrame"""
    def __init__(self,
                 atoms: list[str]  # lines read by Itp class
                 ) -> None:
        self.df = self.get_atoms_info(atoms)

    def get_atoms_info(self,
                       atoms: list[str]  # Lines of the atoms' section
                       ) -> pd.DataFrame:
        """get atoms info from the file"""
        l_line: list[str]  # Breaking the line cahrs
        # Check if header of the atoms section is same as the defeined one
        columns: list[str]   # columns for the atoms dict, name of each column
        columns = ['atomnr', 'atomtype', 'resnr', 'resname', 'atomname',
                   'chargegrp', 'charge', 'mass']
        atomnr: list[typing.Any] = []  # list to append info: atoms id
        atomtype: list[typing.Any] = []  # list to append info: forcefield type
        resnr: list[typing.Any] = []  # list to append info: res infos
        resname: list[typing.Any] = []  # list to append info: res number
        atomname: list[typing.Any] = []  # list to append info: atom name
        chargegrp: list[typing.Any] = []  # list to append info: charge group
        charge: list[typing.Any] = []  # list to append info: charge value
        mass: list[typing.Any] = []  # list to append info: mass value
        # The 3 columns below do not have header, I added to get useful name
        atomsty: list[typing.Any] = []  # list to append info: name with style
        chemi: list[typing.Any] = []  # list to append info: name with alkane
        name: list[typing.Any] = []  # list to append info: real name
        for line in atoms:
            if line.startswith(';'):  # line start with ';' are commets&header
                l_line = free_char_line(line)
                if 'Total' not in l_line:  # Not header!
                    if l_line == columns:
                        pass
                    else:
                        exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                             f'\tError in the [ atoms ] header of the '
                             f'itp file\n{bcolors.ENDC}')
            else:
                l_line = free_char_line(line)
                atomnr.append(l_line[0])
                atomtype.append(l_line[1])
                resnr.append(l_line[2])
                resname.append(l_line[3])
                atomname.append(l_line[4])
                chargegrp.append(l_line[5])
                charge.append(l_line[6])
                mass.append(l_line[7])
                atomsty.append(l_line[8])
                chemi.append(l_line[9])
                name.append(l_line[10])
        df: pd.DataFrame  # DataFrame from the infos
        df = pd.DataFrame(columns=columns)
        df['atomnr'] = atomnr
        df['atomtype'] = atomtype
        df['resnr'] = resnr
        df['resname'] = resname
        df['atomname'] = atomname
        df['chargegrp'] = chargegrp
        df['charge'] = charge
        df['mass'] = mass
        df['atomsty'] = atomsty
        df['chemi'] = chemi
        df['name'] = name
        return df


class BondsInfo:
    """get the bonds list from Itp class and return a dataframe"""
    def __init__(self,
                 bonds: list[str],  # lines of bonds section read by Itp class
                 atoms: pd.DataFrame  # atoms df from AtomsInfo to get names
                 ) -> None:
        """get the bonds infos"""
        self.mk_bonds_df(bonds, atoms)

    def mk_bonds_df(self,
                    bonds: list[str],  # lines of bonds section
                    atoms: pd.DataFrame  # atoms df from AtomInfo
                    ) -> None:
        """call all the methods to make the bonds DataFrame"""
        ai: list[int]  # index of the 1st atoms in the bonds
        aj: list[int]  # index of the 2nd atoms in the bonds
        names: list[str]  # name of the bonds
        ai, aj, names = self.get_bonds(bonds, atoms)
        self.mk_df(ai, aj, names, atoms)

    def get_bonds(self,
                  bonds: list[str],  # lines of bonds section read by Itp class
                  atoms: pd.DataFrame  # atoms df from AtomsInfo to get names
                  ) -> pd.DataFrame:  # DataFrame of the bonds
        """return bonds dataframe to make bonds dataframe"""
        columns: list[str]  # Columns of the bonds wild
        columns = ['ai', 'aj', 'funct', 'r', 'k']
        ai: list[int] = []  # index of the 1st atoms in the bonds
        aj: list[int] = []  # index of the 2nd atoms in the bonds
        names: list[str] = []  # name of the bonds
        for line in bonds:
            if line.startswith(';'):  # line start with ';' are commets&header
                l_line = free_char_line(line)
                if 'Total' not in l_line:  # Not header!
                    if l_line == columns:
                        pass
                    else:
                        exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                             f'\tError in the [ bonds ] header of the '
                             f'itp file\n{bcolors.ENDC}')
            else:
                l_line = free_char_line(line)
                ai.append(int(l_line[0]))
                aj.append(int(l_line[1]))
                names.append(l_line[3])
        return ai, aj, names

    def mk_df(self,
              ai: list[int],  # index of the 1st atom in the bonds
              aj: list[int],  # index of the 2nd atom in the bonds
              names: list[str],  # names of the bonds form bonds section
              atoms: pd.DataFrame  # atoms df from AtomsInfo to cehck the name
              ) -> None:
        """make DataFrame and check if they are same as atoms name"""
        df: pd.DataFrame  # to save the bonds_df
        df = pd.DataFrame(columns=['typ', 'ai', 'aj', 'cmt', 'name'])
        df['ai'] = ai
        df['aj'] = aj
        df['name'] = names
        df['cmt'] = ['#' for _ in ai]
        df['typ'] = self.get_type(names)
        self.check_names(df, atoms)
        # print(df)

    def get_type(self,
                 lst: list[str]  # list to get the number of distenguished ones
                 ) -> list[int]:  # types' index
        """make type based on the unique items in the lst"""
        type_set: set[str] = set(lst)  # eleminate the repeated names
        type_dict: dict[str, int]  # to make a list with type index
        type_dict = {item: i+1 for i, item in enumerate(type_set)}
        types: list[int]  # types to return
        types = [type_dict[item] for item in lst]
        return types

    def check_names(self,
                    df: pd.DataFrame,  # Df to check its names column
                    atoms: pd.DataFrame  # Df source (mostly atoms df)
                    ) -> None:
        """ checks the names column in the source file with comparing it
        with names from the atoms dataframe name column for each atom"""
        ai_name: list[str]  # name of the 1st atom from the list
        aj_name: list[str]  # name of the 2nd atom from the list
        ai_name = [atoms.loc[atoms['atomnr'] == str(item)]['atomsty'][item-1]
                   for item in df['ai']]
        aj_name = [atoms.loc[atoms['atomnr'] == str(item)]['atomsty'][item-1]
                   for item in df['aj']]
        names: list[str] = [f'{i}-{j}' for i, j in zip(ai_name, aj_name)]
        if list(df['name']) != names:
            exit(f'{bcolors.FAIL}\tError! in the bonds and atoms name!'
                 f'{bcolors.ENDC}')


if __name__ == '__main__':
    itp = Itp(sys.argv[1])
