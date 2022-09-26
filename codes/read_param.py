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
        atoms: list[str] = []  # Atoms informations
        bonds: list[str] = []  # Bonds informations
        angles: list[str] = []  # Angles informations
        dihedrals: list[str] = []  # Dihedrals informations
        with open(fname, 'r') as f:
            while True:
                line: str = f.readline()
                if line.strip():
                    line = line.strip()
                    if line.startswith('#'):  # Comments lines
                        pass
                    else:
                        if line.startswith('atom'):
                            atoms.append(line)
                        elif line.startswith('bond'):
                            bonds.append(line)
                        elif line.startswith('angle'):
                            angles.append(line)
                        elif line.startswith('dihedral'):
                            dihedrals.append(line)
                else:
                    pass
                if not line:
                    break


if __name__ == '__main__':
    fname: str = sys.argv[1]  # Input file
    ReadParam(fname)
