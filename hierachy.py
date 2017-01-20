# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 16:45:49 2016

@author: Dell
"""

import seaborn as sns; sns.set()
flights = sns.load_dataset("flights")
flights = flights.pivot("month", "year", "passengers")
g = sns.clustermap(flights)