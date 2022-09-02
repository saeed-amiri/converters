import sys
import json_to_lmp as jlmp
import write_lmp as lmpwrt


class Doc:
    """convert files in other format into LAMMPS data file
    Input:
        a data file in JSON
    Output:
        a data file in LAMMPS
    """


fname: str = sys.argv[1]  # File name to convert
data = jlmp.ConvertJson(fname)
output_fname = 'test.data'
wrt = lmpwrt.WriteLmp(data, output_fname)
wrt.write_lmp()
