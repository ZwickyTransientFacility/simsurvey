#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" Read and parse the simulated data. """

import pandas as pd

__all__ =["read_mockdata"]


def read_mockdata(mockdata,sep=" ",comment="#",**kwargs):
    """ read file formated as
    # name1 name2 etc.
    data1 data2 etc.
    data1 data2 etc.
    data1 data2 etc.


    Parameters:
    ----------
    sep : string, default ' ' (tab-stop)
        Delimiter to use. Regular expressions are accepted.
    comment: string, default "#"
        any line starting with this charater (single character)
        will be 
    
    **kwargs goes to pandas.read_table()
        
    
    Returns
    -------
    Pandas.DataFrame
    """
    return pd.read_table(mockdata, sep=sep, comment=comment,
                         names=open(mockdata).read().splitlines()[0].split()[1:],
                         **kwargs)
