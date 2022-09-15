import sys
import typing
import pandas as pd


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
        with open(fname, 'r') as f:
            while True:
                line: str = f.readline()
                l_line: list[str]  # Breaking the line to cahrs
                if line:
                    l_line = line.strip().split()
                    print(l_line)
                if not line:
                    break


if __name__ == '__main__':
    itp = Itp(sys.argv[1])
