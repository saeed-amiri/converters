import sys
import json_to_lmp as jlmp
import write_lmp as lmpwrt
import itp_pdb_lmp as ipl


class Doc:
    """convert files in other format into LAMMPS data file
    Input:
        a data file in JSON
        or
        a data file in itp (it needs two file with same name with .pdb)
    Output:
        a data file in LAMMPS
    """


fname: str = sys.argv[1]  # File name to convert
f_extansion: str  # Extension of the file to call proper function
f_extansion = fname.split('.')[1]
if f_extansion == 'json':
    data = jlmp.ConvertJson(fname)
if f_extansion == 'itp':
    data = ipl.ItpPdb(fname.split('.')[0])
output_fname = 'test.data'
wrt = lmpwrt.WriteLmp(data, output_fname)
wrt.write_lmp()
