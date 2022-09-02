import sys
import json_to_lmp as jlmp

class Doc:
    """convert files in other format into LAMMPS data file
    Input:
        a data file in JSON
    Output:
        a data file in LAMMPS
    """

fname: str = sys.argv[1]  # File name to convert
data = jlmp.ConvertJson(fname)
