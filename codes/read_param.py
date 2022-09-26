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
        self.set_fields(fname)

    def set_fields(self,
                   fname: str  # Name of the input file
                   ) -> None:  # Set attributes to self
        """set force filed for each section to self"""
        atoms: list[str]  # Atoms informations
        bonds: list[str]  # Bonds informations
        angles: list[str]  # Angles informations
        dihedrals: list[str]  # Dihedrals informations
        atoms, bonds, angles, dihedrals = self.get_params(fname)
        if atoms:
            self.atoms: pd.DataFrame = self.get_atoms(atoms)  # LJ information

    def get_atoms(self,
                  atoms: list[str]  # atoms line in param file
                  ) -> pd.DataFrame:  # Atoms LJ information
        """convert info to dataframe"""
        columns: list[str]  # Columns for the df
        columns =['atom_name', 'mass', 'sigma', 'epsilom', 'charge']
        df = pd.DataFrame(columns=columns)
        atom_name: list[str] = []  # To save from each line
        mass: list[str] = []  # To save from each line
        sigma: list[str] = []  # To save from each line
        epsilom: list[str] = []  # To save from each line
        charge: list[str] = []  # To save from each line
        for item in atoms:
            l_line: list[str]  # breaking the line
            l_line = item.strip().split(' ')
            atom_name.append(l_line[1])
            mass.append(l_line[2])
            sigma.append(l_line[3])
            epsilom.append(l_line[4])
            charge.append(l_line[5])
        df['atom_name'] = atom_name
        df['mass'] = mass
        df['sigma'] = sigma
        df['epsilom'] = epsilom
        df['charge'] = charge
        return df

    def get_params(self,
                   fname: str  # Name of the input file
                   ) -> tuple[list[str], list[str], list[str], list[str]]:
                   # return lines about atoms, bonds, angles, dihedrals
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
        return atoms, bonds, angles, dihedrals
        

if __name__ == '__main__':
    fname: str = sys.argv[1]  # Input file
    ReadParam(fname)
