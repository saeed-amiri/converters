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
        self.atoms_info: dict[int, list[str]] = {}  # Lines of atoms section
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
                            self.get_atoms_info(line.strip())
                if not line:
                    break

    def get_atoms_info(self,
                       line: str  # Line of the atoms' section)
                       ) -> None:
        """get atoms info from the file"""
        l_line: list[str]  # Breaking the line cahrs
        keys: list[str]   # Keys for the atoms dict
        keys = ['atomnr', 'atomtype', 'resnr', 'resname', 'atomname',
                'chargegrp', 'charge', 'mass']

        if line.startswith(';'):
            l_line = self.free_char_line(line)
            if l_line[0] != 'Total':
                if l_line != keys:
                    exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                         f'\tError in the [ atoms ] header of the itp file\n'
                         f'{bcolors.ENDC}')
        else:
            l_line = self.free_char_line(line)

    def free_char_line(self,
                       line: str  # line of the itp file
                       ) -> list[str]:  # Free from chars
        """cheack the lines and return the line free special chars"""
        char_list: list[str] = [';', '#', ':']  # chars to eliminate from lines
        l_line: list[str]  # Breaking the line cahrs
        l_line = line.strip().split(' ')
        l_line = [item for item in l_line if item]
        l_line = [item for item in l_line if item not in char_list]
        return l_line


if __name__ == '__main__':
    itp = Itp(sys.argv[1])
