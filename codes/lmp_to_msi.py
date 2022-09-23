import pandas as pd
import sys
import read_lmp_data as rdlmp
import read_cvff as cvtyp
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
        lmp: rdlmp.ReadData  # Whole LAMMPS data
        cvff = cvtyp.Cvff()  # cvff atom types
        lmp = rdlmp.ReadData(fname)
        self.to_car(lmp, cvff)

    def to_car(self,
               lmp: rdlmp.ReadData,  # Whole LAMMPS data
               cvff # cvff atom types
               ) -> None:
        """call all the methods to convert data"""
        print(cvff.df)

    def mk_df(self,
              atoms: pd.DataFrame):
        """"""

if __name__ == '__main__':
    msi = Mdf(sys.argv[1])
