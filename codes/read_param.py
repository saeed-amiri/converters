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
        if bonds:
            self.bonds: pd.DataFrame = self.get_bonds(bonds)  # Bonds infos
        if angles:
            self.angles: pd.DataFrame = self.get_angles(angles)  # Angles infos
        if dihedrals:
            self.dihedrals: pd.DataFrame = self.get_dihedrals(dihedrals)

    def get_atoms(self,
                  atoms: list[str]  # atoms line in param file
                  ) -> pd.DataFrame:  # Atoms LJ information
        """convert info to dataframe"""
        columns: list[str]  # Columns for the df
        columns = ['atom_name', 'mass', 'sigma', 'epsilom', 'charge']
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

    def get_bonds(self,
                  bonds: list[str]  # List of bonds lines
                  ) -> pd.DataFrame:  # Bonds interaction, in harmonic format
        """return data for harmonic interactions"""
        columns: list[str]  # Columns for the df
        columns = ['bond_name', 'r', 'kbond']
        df = pd.DataFrame(columns=columns)
        bond_name: list[str] = []  # Name of the bonds (based on type)
        r: list[str] = []  # To save bonds length
        kbond: list[str] = []  # To save bonds strength
        for item in bonds:
            l_line: list[str]  # breacking the line
            l_line = item.strip().split(' ')
            bond_name.append(l_line[1])
            r.append(l_line[2])
            try:
                kbond.append(l_line[3])
            except IndexError:
                kbond.append('0.0')
        df['bond_name'] = bond_name
        df['r'] = r
        df['kbond'] = kbond
        return df

    def get_angles(self,
                   angles: list[str]  # List of angles lines
                   ) -> pd.DataFrame:  # Angles interaction, in harmonic format
        """return data for harmonic interactions"""
        columns: list[str]  # Columns for the df
        columns = ['angle_name', 'angle', 'kangle']
        df = pd.DataFrame(columns=columns)
        angle_name: list[str] = []  # Name of the angles (based on type)
        angle: list[str] = []  # To save angles length
        kangle: list[str] = []  # To save angles strength
        for item in angles:
            l_line: list[str]  # breacking the line
            l_line = item.strip().split(' ')
            angle_name.append(l_line[1])
            angle.append(l_line[2])
            try:
                kangle.append(l_line[3])
            except IndexError:
                kangle.append('0.0')
        df['angle_name'] = angle_name
        df['angle'] = angle
        df['kangle'] = kangle
        return df

    def get_dihedrals(self,
                      dihedrals: list[str]  # List of dihedrals lines
                      ) -> pd.DataFrame:  # Dihderals interaction, OPLS format
        """return data for OPLS interactions"""
        columns: list[str]  # Columns for the df
        columns = ['dihedral_name', 'k1', 'k2', 'k3', 'k4']
        df = pd.DataFrame(columns=columns)
        dihedral_name: list[str] = []  # Name of the dihedrals (based on type)
        k1: list[str] = []  # To save dihedrals k1
        k2: list[str] = []  # To save dihedrals k2
        k3: list[str] = []  # To save dihedrals k3
        k4: list[str] = []  # To save dihedrals k4
        for item in dihedrals:
            l_line: list[str]  # breacking the line
            l_line = item.strip().split(' ')
            dihedral_name.append(l_line[1])
            k1.append(l_line[2])
            k2.append(l_line[3])
            k3.append(l_line[4])
            k4.append(l_line[4])
        df['dihedral_name'] = dihedral_name
        df['k1'] = k1
        df['k2'] = k2
        df['k3'] = k3
        df['k4'] = k4
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
