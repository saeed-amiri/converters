import os
import pandas as pd


class Doc:
    """read atoms_types_cvff.data file and return it as DataFrame object
    Input:
        atoms_types_cvff.data (muyst be at same directory as this script)
    Output:
        DataFrame as an object
    """


class Cvff:
    """get atom names, element type, and Mass from the file"""
    def __init__(self) -> None:
        """read the file"""
        fname: str  # Name the file
        dir: str  # directory of the file
        f: str  # Name of the file
        dir = '/scratch/saeed/MyScripts/converters/codes'
        f = 'atom_types_cvff.data'
        fname = os.path.join(dir, f)
        self.df: pd.DataFrame = self.get_cvff(fname)  # All the atoms_type

    def get_cvff(self,
                 fname: str  # Name the file
                 ) -> pd.DataFrame:  # The main DataFrame of infos
        """read the file here"""
        line: str  # Each line of the file
        l_line: list[str]  # list of each char of each line
        l_cvff: list[list[str]] = []  # All the processed lines
        with open(fname, 'r') as f:
            while True:
                line = f.readline()
                if line.strip().startswith('!'):
                    pass
                else:
                    l_line = self.process_line(line.strip())
                    if l_line:
                        l_cvff.append(l_line)
                if not line:
                    break
        df: pd.DataFrame  # All the infos
        columns: list[str]  # Columns of the data file
        columns = ['Ver', 'Ref', 'Type', 'Mass', 'Element',
                   'Connections', 'Comment']
        df = pd.DataFrame(l_cvff, columns=columns)
        return df

    def process_line(self,
                     line: str  # read line of the file
                     ) -> list[str]:  # seperated line of the input
        """break down the line into picess"""
        l_line: list[str]  # To retirn the main data
        # First break down the line
        l_line = line.split(' ')
        l_line = [item.strip() for item in l_line if item]
        # Rejoin them so ever char seperated with one space
        line = ' '.join(l_line)
        l_line = line.split(' ', 6)  # item.index > 6 are comments
        l_line = [item.strip() for item in l_line if item]
        return l_line


if __name__ == '__main__':
    cvff = Cvff()
