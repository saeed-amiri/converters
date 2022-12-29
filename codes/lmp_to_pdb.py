import pandas as pd


class Pdb:
    """convert the Atoms section into PDB file
    Input:
        LAMMPS data file from lmp_itp_pdb.py:
            Atoms_df, Mass_df
    Output:
        pd.DataFrame for pdb file
    """
