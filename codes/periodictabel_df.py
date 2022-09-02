from pprint import pprint
import typing
import json
import pandas as pd


class Doc:
    """convert the JSON file of periodic table to DataFrame
    Input:
        ./periodic_table.json'
    Output:
        pandas DataFrame
    """


class PeriodicTable:
    """read the input file and convert to DataFrame"""
    def __init__(self) -> None:
        fname: str  # Name of the file, IT IS STATIC!
        fname = '/scratch/saeed/MyScripts/converters/codes/periodic_table.json'
        data: dict[str, typing.Any] = self.read_json(fname)
        self.peridic_table = self.dict_to_df(data)

    def dict_to_df(self,
                   data: dict[str, typing.Any]
                   ) -> pd.DataFrame:
        """convert dictionary to dataframe"""
        df = pd.DataFrame.from_dict(data)
        return df

    def read_json(self,
                  fname: str  # Name of the input file
                  ) -> dict[str, typing.Any]:
        """read json file as a dictionary"""
        with open(fname, 'r') as f:
            data = json.load(f)
        return data['elements'][0]


if __name__ == '__main__':
    table = PeriodicTable()
