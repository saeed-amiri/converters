import sys
import typing
import collections
import numpy as np
import pandas as pd
import read_cvff as cvtyp
import read_lmp_data as rdlmp
from colors_text import TextColor as bcolors


class Doc:
    """convert LAMMPS data file into Material Studio:
    LAMMPS.data to .mdf and .car.
    I am doing these to use the available ForceFields, which can be
    read and used by msi2lmp tools.
    Specifically, I need this for surfactant force fields.
    I could not find a description for these two file formats.
    I used the description from:
    matsci.org/t/manual-creation-of-car-and-mdf-files-for-msi2lmp-and-emc/27494
    and available examples.
    Input:
        LAMMPS DATA file
    Output:
        Material Studio MDF and CAR files

    In .car file:
    "Upper-case atom names (“C1”, “C2”, …“C6”, and “O1”, “O2”, …) are
    unique identifiers for each atom.
    Lower-case atom names (“cp,” “c5”, “oh,” “o”,…)  are atom TYPES,
    which are used for looking up force-field parameters.
    Note that there are multiple types of carbon, oxygen, and hydrogen
    atoms corresponding to these elements in different bonding states
    and local environments. To make the file, one will have to carefu-
    lly read the description of each atom type and decide which type
    to use for all of the atoms in the molecule and what their partial
    charges should be. This is NOT trivial. (Note: You can figure out
    the partial charges by looking at the bond_increments section of
    the FRC file and sum up all of the contributions due to all of the
    bonded neighbors of each atom.)"
    """


class Car:
    """make car file by reading lammps data file"""
    def __init__(self,
                 fname: str  # Name of the input file
                 ) -> None:
        self.lmp: rdlmp.ReadData  # Whole LAMMPS data
        cvff = cvtyp.Cvff()  # cvff atom types
        self.lmp = rdlmp.ReadData(fname)
        self.to_car(cvff, fname)

    def to_car(self,
               cvff: cvtyp.Cvff,  # cvff atom types
               fname: str  # name of the input file
               ) -> None:
        """call all the methods to convert data"""
        self.df: pd.DataFrame = self.mk_df(self.lmp.Atoms_df, cvff)  # Car df
        self.write_car(fname)

    def mk_df(self,
              atoms: pd.DataFrame,  # Atoms_df of the full atom LAMMPS
              cvff: cvtyp.Cvff  # Atoms types, names, mass
              ) -> pd.DataFrame:  # Car DataFrame
        """make df in the form of the car file format"""
        columns: list[str]  # name of the df columns
        columns = ['atom_name', 'x', 'y', 'z', 'mol_name', 'mol', 'type',
                   'element', 'charge']
        df: pd.DataFrame  # Car DataFrame
        df = pd.DataFrame(columns=columns)
        df['x'] = atoms['x']  # Coordinates
        df['y'] = atoms['y']  # Coordinates
        df['z'] = atoms['z']  # Coordinates
        df['mol'] = atoms['mol']  # Index of the mol
        df['mol_name'] = ['UNK1' for _ in df.index]  # Name of the mol
        df['type'] = atoms['name']  # Type of the atom
        elements: list[str] = self.get_element(atoms, cvff)  # Symbol of atoms
        df['element'] = elements
        df['charge'] = atoms['charge']
        df['atom_name'] = self.mk_atom_name(elements)
        return df

    def get_element(self,
                    atoms: pd.DataFrame,  # Atoms_df of the full atom LAMMPS
                    cvff: cvtyp.Cvff  # Atoms types, names, mass
                    ) -> list[str]:  # Elements names
        """get the element of each type from cvff file"""
        elements: list[str]  # Elements name
        names: list[str] = list(atoms['name'])  # Names of the elements
        elements = [
            cvff.df.loc[cvff.df['Type'] == item]['Element'].item() for
            item in names]
        return elements

    def mk_atom_name(self,
                     elements: list[str]  # Element symbols of each atom
                     ) -> list[str]:  # list of the atoms names
        """retrun name for each atom based on thier element
        It first get unique elements and thier number of repetition
        Then rename them based on the elements
        """
        element_set: list[str] = self.seen_set(elements)  # Unique elements
        l_occurence: list[int] = []  # number of each occurence of each atom
        l_element: list[str] = []  # Name of the each element, this is safer
        main_ind: list[int] = []  # index of each atom in the element list
        atom_names = ['None' for _ in elements]  # final names of the atoms
        l_unq: list[str]  # name of the duplicated atoms
        for uniq in element_set:
            l_unq = []
            for k, item in enumerate(elements):
                if item == uniq:
                    l_unq.append(item)
                    main_ind.append(k)
            for ind, elem in enumerate(l_unq):
                l_occurence.append(ind+1)
                l_element.append(elem)
        for i, el in enumerate(main_ind):
            atom_names[el] = f'{l_element[i]}{l_occurence[i]}'
        return atom_names

    def seen_set(self,
                 lst: list[typing.Any]  # to drop duplicate with keping order
                 ) -> list[typing.Any]:
        """remove duplicated item with keeping order of them in the
        main list"""
        seen: set[str] = set()
        seen_add = seen.add
        return [x for x in lst if not (x in seen or seen_add(x))]

    def write_car(self,
                  fname: str  # Name of the input file
                  ) -> None:
        """write .car file with name of the input file"""
        outname: str = f'{fname.split(".")[0]}.car'  # Name of the output file
        with open(outname, 'w') as f:
            df: pd.DataFrame = self.df.copy()  # To  keep main df intact
            f.write(f'!BIOSYSM XXX X\n')
            f.write(f'PBC=OFF\n')
            f.write(f'Material Studio Generated CAR File\n')
            f.write(f'!written by: {self.__class__.__name__}\n')
            df = df.astype({'x': float, 'y': float, 'z': float})
            df.to_csv(f, sep='\t', header=None, index=False,
                      float_format='%f')
            f.write('end\n')
            f.write('end\n')
            f.write(f'\n')


class Mdf:
    """write the MDF file
    The mdf files has a following columns:
    @column 0 name -> This line Will NOT write into file
    @column 1 element
    @column 2 atom_type
    @column 3 charge_group
    @column 4 isotope
    @column 5 formal_charge
    @column 6 charge
    @column 7 switching_atom
    @column 8 oop_flag
    @column 9 chirality_flag
    @column 10 occupancy
    @column 11 xray_temp_factor
    @column 12 connections

    @molecules XXX: name
    """
    def __init__(self,
                 fname: str,  # Name of the input file
                 car_df: pd.DataFrame,  # from Car class
                 lmp: rdlmp.ReadData  # main data read by Car class
                 ) -> None:
        self.to_mdf(fname, car_df, lmp)

    def to_mdf(self,
               fname: str,  # Name of the input file
               car_df: pd.DataFrame,  # from Car class
               lmp: rdlmp.ReadData  # main data read by Car class
               ) -> None:
        """call all the methods and write the file"""
        self.mk_df(car_df, lmp)

    def mk_df(self,
              car_df: pd.DataFrame,  # from Car class
              lmp: rdlmp.ReadData  # main data read by Car class
              ) -> pd.DataFrame:  # To write into file
        """make DataFrame in the MDF format"""


if __name__ == '__main__':
    car = Car(sys.argv[1])
    mdf = Mdf(sys.argv[1], car.df, car.lmp)
