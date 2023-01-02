import pandas as pd
import read_lmp_data as relmp
from colors_text import TextColor as bcolors


class Doc:
    """read the LAMMPS data file and get force field parameters and write
        `itp` formats for GROMACS conversion"""


class ItpStyleInfos:
    """ informations in the itp file:
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
    dihedrals."""


class Itp:
    """get data from main"""
    def __init__(self,
                 lmp: relmp.ReadData  # LAMMPS data file
                 ) -> None:
        self.mk_itp(lmp)

    def mk_itp(self,
               lmp: relmp.ReadData  # LAMMPS data file
               ) -> None:
        """call functions"""
        self.mk_atoms(lmp)
        self.mk_bonds(lmp)
        self.mk_angles(lmp)
        self.mk_dihedrals(lmp)
        self.mk_pairs(lmp)

    def mk_atoms(self,
                 lmp: relmp.ReadData  # LAMMPS data file
                 ) -> pd.DataFrame:
        pass

    def mk_bonds(self,
                 lmp: relmp.ReadData  # LAMMPS data file
                 ) -> pd.DataFrame:
        pass

    def mk_angles(self,
                  lmp: relmp.ReadData  # LAMMPS data file
                  ) -> pd.DataFrame:
        pass

    def mk_dihedrals(self,
                     lmp: relmp.ReadData  # LAMMPS data file
                     ) -> pd.DataFrame:
        pass

    def mk_pairs(self,
                 lmp: relmp.ReadData  # LAMMPS data file
                 ) -> pd.DataFrame:
        pass
