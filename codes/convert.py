import sys
import json_to_lmp as jlmp
import write_lmp as lmpwrt
import itp_pdb_lmp as ipl
import clean_data as cllmp
import write_param as rdprm


class Doc:
    """convert files in other format into LAMMPS data file
    Input:
        a data file in JSON
        or
        a data file in itp (it needs two file with same name with .pdb)
        or
        a data file in data (LAMMPS full atom) to clean up
    Output:
        a data file in LAMMPS
    """


fname: str = sys.argv[1]  # File name to convert
f_extansion: str  # Extension of the file to call proper function
f_extansion = fname.split('.')[1]
output_fname = f'{fname.split(".")[0]}_lmp.data'
if f_extansion == 'json':
    data = jlmp.ConvertJson(fname)
elif f_extansion == 'itp':
    data = ipl.ItpPdb(fname.split('.')[0])
elif f_extansion == 'data':
    data = cllmp.CleanData(fname)
    output_fname = f'CTAB_lmp.data'
    prm = rdprm.WriteParam(data, output_fname)
wrt = lmpwrt.WriteLmp(data, output_fname)
wrt.write_lmp()
