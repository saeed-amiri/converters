import typing
import pandas as pd
import read_param as rdprm
import clean_data as cllmp


class Doc:
    """update the DataFrame from read_param by the type of each set
    from clean_data
    input:
        DataFrames from clean_data
        DataFrames from read_param
    Output:
        File in JSON with all the interactions that can be read by
        combination codes
    """

# Rteturn unique atoms
def seen_set(
             lst: list[typing.Any]  # to drop duplicate with keping order
             ) -> list[typing.Any]:
    """remove duplicated item with keeping order of them in the
    main list"""
    seen: set[str] = set()
    seen_add = seen.add
    return [x for x in lst if not (x in seen or seen_add(x))]


class GetType:
    """set type for bonds, angles, dihedrals"""
    def __init__(self,
                 data: cllmp.CleanData  # Cleaned data to write
                 ) -> None:
        param: rdprm.ReadParam  # All the data from param file
        param = rdprm.ReadParam(fname='param')
        self.set_types(param, data)

    def set_types(self,
                  param: rdprm.ReadParam,  # ForceField from param
                  data: cllmp.CleanData  # Cleaned data to write
                  ) -> None:
        """call all the methods to set types and write them"""
        self.set_atoms(param.atoms, data.Atoms_df)
    
    def set_atoms(self,
                  atoms_prm: pd.DataFrame,  # LJ parameters for atoms
                  df: pd.DataFrame  # Atoms information
                  ) -> pd.DataFrame:  # Updated atoms parameters
        """update parameters of atoms with thier types"""
        atoms_name: list[str]  # Names of the atoms != Names in bonds, etc
        atoms_name = list(df['name'])
        atoms_name = seen_set(atoms_name)
        atoms_type: dict[str, int] = {}  # Atoms name with thier types
        for item in atoms_name:
            typ: int = df.loc[df['name'] == item]['typ']  # Type of the atom
            atoms_type[item] = seen_set(typ)[0]
        i_type: list[int] = []  # Type of each atom
        for item in atoms_prm['atom_name']:
            i_type.append(atoms_type[item])
        atoms_prm['type'] = i_type
        return atoms_prm



