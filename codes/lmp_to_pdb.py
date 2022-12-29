import pandas as pd


class Doc:
    """convert the Atoms section into PDB file
    Input:
        LAMMPS data file from lmp_itp_pdb.py:
            Atoms_df, Masses_df
    Output:
        pd.DataFrame for pdb file

    Main notes:
        ? The names of atoms should be based on the force field:
        > can be set in the input files in the Masses section after a `#.`
        ? The names of the Residues (molecules) should be set somehow
        since, in the lammps, they just separated by numbers:
        > The second name after the atom name in the Masses section
        should be the name of the residue to which the atom belongs
        and should be three or four characters long.
        - Limitations for atoms' names also should be considered.
        ? The element symbols should also be considered:
        The easiest way is to set it in the Masses section after the
        residues' names and in CAPTIAL letters.
        ? The ATOMS or HATOM:
        > It also set after the symbols in Masses section.
        ? TER: This shows that the last atom in the chain prevents the
        display of a connection to the next chain. It specified by the
        index of the aotm:
        > Bit tricky! Do it in the code by selecting the final atom as
        TER if the it is ATOM in the Masses section. For my input should
        work since the last atom is the final atom.
        ? Duplicate Atom Names: One possible editing mistake is the
        failure to uniquely name all atoms within a given residue.
        > Do it in the code by adjusting duplicate atom names by adding
        an index.

        - The following recordes will be left empty chars:
            17	    Alternate location indicator		character
            22	    Chain identifier		            character
            27	    Code for insertions of residues		character
            55-60	Occupancy	                right	real (6.2)
            61-66	Temperature factor	        right	real (6.2)
            73-76	Segment identifier¶	        left	character

    """


class PdbStyleInfo:
    """
        The PDB file consider is in standard PDB format.

        convert LAMMPS data file to a standard PDB file format based on:
    [https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html]
        A PDB file has a format as:
        ATOM 	atomic coordinate record containing the X,Y,Z
            orthogonal Å coordinates for atoms in standard residues]
            (amino acids and nucleic acids).

        HATATM 	atomic coordinate record containing the X,Y,Z orthogonal Å
        coordinates for atoms in nonstandard residues. Nonstandard residues
        include inhibitors, cofactors, ions, and solvent. The only functional
        difference from ATOM records is that HETATM residues are by default
        not connected to other residues. Note that water residues should be
        in HETATM records.

        Protein Data Bank Format for ATOM:
          Coordinate Section
            Record Type	Columns	Data 	Justification	Data Type
            ATOM 	1-4	“ATOM”	                    	character
            7-11#	Atom serial number      	right	integer
            13-16	Atom name	                left*	character
            17	    Alternate location indicator		character
            18-20§	Residue name	            right	character
            22	    Chain identifier		            character
            23-26	Residue sequence number	    right	integer
            27	    Code for insertions of residues		character
            31-38	X orthogonal Å coordinate	right	real (8.3)
            39-46	Y orthogonal Å coordinate	right	real (8.3)
            47-54	Z orthogonal Å coordinate	right	real (8.3)
            55-60	Occupancy	                right	real (6.2)
            61-66	Temperature factor	        right	real (6.2)
            73-76	Segment identifier¶	        left	character
            77-78	Element symbol              right	character
            79-80	Charge		                        character

        Protein Data Bank Format for HATATM:
          Coordinate Section
            Record Type	Columns	Data 	Justification	Data Type
            1-6	“HETATM”		character
            7-80	same as ATOM records

        #Chimera allows (nonstandard) use of columns 6-11 for the integer
            atom serial number in ATOM records, and in TER records, only the
            “TER” is required.
        *Atom names start with element symbols right-justified in columns
            13-14 as permitted by the length of the name. For example, the
            symbol FE for iron appears in columns 13-14, whereas the symbol
            C for carbon appears in column 14 (see Misaligned Atom Names).
            If an atom name has four characters, however, it must start in
            column 13 even if the element symbol is a single character
            (for example, see Hydrogen Atoms).
        §Chimera allows (nonstandard) use of four-character residue names
            occupying an additional column to the right.
        ¶Segment identifier is obsolete, but still used by some programs.
            Chimera assigns it as the atom attribute pdbSegment to allow
            command-line specification.
        The format of ecah section is (fortran style):
        Format (A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)
        """

class Pdb:
    """convert to PDB with considers all the above concerns"""
    def __init__(self,
                Masses_df: pd.DataFrame,  # df contains info about atoms
                Atoms_df: pd.DataFrame  # df contains atoms' coordinats
                ) -> None:
        """call the main functions"""
        self.get_atoms_info(Masses_df)
    
    def get_atoms_info(self,
                       Masses_df: pd.DataFrame,  # df contains info about atoms
                       ) -> pd.DataFrame:
        """get data in the Masses_df and check them"""
        print(Masses_df)