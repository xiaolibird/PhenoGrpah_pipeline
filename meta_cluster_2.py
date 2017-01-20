# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 22:14:10 2017

@author: Dell
"""
import numpy as np
import pandas as pd
import fcs_reader as fcsrd
import glob
import os

from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_rand_score
from sklearn.metrics import precision_recall_fscore_support


#from sklearn.cluster import KMeans, MiniBatchKMeans



quantile_995 = np.load("normalizer.npy")
n_healthy = 5
n_patients = 16
n_conditions = 18;
import phenograph as pg
from time import clock

#       clustering all data for one sample
sample_names = []
for i in range(n_patients):
#for i in range(1):
    fname = glob.glob('F:\\cytowork\\experiment_44185_files\\*.fcs')[n_conditions*(i+n_healthy):n_conditions*(i+n_healthy+1)]
#    create folder 
    sample_name = fname[0].split('\\')

    sample_name = sample_name[-1];
    sample_name = sample_name.split('_')
    sample_name = sample_name[0]
    sample_names.append(sample_name)        
    

flag = False    
for sample_name in sample_names:    
    if flag:
        centroids = np.load(sample_name+"\\centroids.npy")
        print(centroids)
        flag = True
    else:
        centroids = np.vstack((centroids,np.load(sample_name+"\\centroids.npy")))
        print(centroids)