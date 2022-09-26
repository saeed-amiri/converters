import readline
import sys
import pandas as pd


class Doc:
    """read param file
    The file contain Force Field informations for the interactions
    inside the system runing.
    input:
        file: param
    output:
        DataFrame(s) if there were information for any of the parts
            atoms, bonds, angles, dihedrals
    """


class ReadParam:
    """read param file"""
    def __init__(self,
                 fname: str = 'param'  # Name of the parameter file
                 ) -> None:
        self.get_params(fname)

    def get_params(self,
                   fname: str  # Name of the input file
                   ) -> None:  # set attributs to self
        """read and call all the methods"""
        with open(fname, 'r') as f:
            while True:
                line: str = f.readline()
                if line.strip():
                    line = line.strip()
                    if line.startswith('#'):  # Comments lines
                        pass
                    else:
                        if line.startswith('atom'):
                            pass
                        elif line.startswith('bond'):
                            pass
                        elif line.startswith('angle'):
                            pass
                        elif line.startswith('dihedral'):
                            pass
                else:
                    pass
                if not line:
                    break


if __name__ == '__main__':
    fname: str = sys.argv[1]  # Input file
    ReadParam(fname)
