# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 14:48:32 2017

@author: Dell
"""
import numpy as np
#import pandas as pd
import fcs_reader as fcsrd
import glob
import os

from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_rand_score
from sklearn.metrics import precision_recall_fscore_support


#from sklearn.cluster import KMeans, MiniBatchKMeans



quantile_995 = np.load("normalizer.npy")
n_patients = 5 + 16;
n_conditions = 18;
import phenograph as pg
from time import clock
for i in range(1):
#for i in :
    fname = glob.glob('F:\\cytowork\\experiment_44185_files\\*.fcs')[n_conditions*i:n_conditions*(i+1)]
#    create folder 
    sample_name = fname[0].split('\\')
    sample_name = sample_name[-1]
    sample_name = sample_name.split('_')
    sample_name = sample_name[0]
    os.mkdir(sample_name)
    for ii in range(len(fname)):
        if ii == 0:           
            meta, data_numpy = fcsrd.parse_fcs(fname[ii], meta_data_only=False, output_format='ndarray', reformat_meta=True)
#            surface_marker_data = data_numpy[:,surface_marker.index - 1]
            surface_marker_data = data_numpy           
        else:
            meta, data_numpy = fcsrd.parse_fcs(fname[ii], meta_data_only=False, output_format='ndarray', reformat_meta=True)      
            surface_marker_data = np.vstack((surface_marker_data,data_numpy)) 
            
#       choose channel    
        meta_data = meta['_channels_']   
        channel_names = (meta_data['$PnS'].values)
        selected_index = [] 
        for name in channel_names:
            if "CD" in name:
                selected_index.insert(0,name)                     
        selected_index.insert(0,"HLA-DR")
        selected_index = np.array(selected_index)
        surface_marker = meta_data[meta_data['$PnS'].isin(selected_index)]  
        surface_marker_data_selected = surface_marker_data[:,surface_marker.index - 1]
        
#       transform and normalize
        suface_marker_data_transformed = np.arcsinh(surface_marker_data_selected/5)
        suface_marker_data_normalized = np.divide(suface_marker_data_transformed,np.tile(np.transpose(quantile_995),(np.shape(suface_marker_data_transformed)[0],1)))
        
        
#       slice ground truth provided by author
        ground_truth = np.array(["PhenoGraph"])
        ground_truth = meta_data[meta_data['$PnS'].isin(ground_truth)]  
        ground_truth = ground_truth.index - 1
        ground_truth = surface_marker_data[:,ground_truth]
        ground_truth = np.reshape(ground_truth, (np.product(ground_truth.shape),))
        
#       clustering all data for one sample
        n_randstart = 10
        c = np.zeros((np.shape(suface_marker_data_normalized)[0],n_randstart))
        q = np.zeros((1,n_randstart))
        t = np.zeros((1,n_randstart))
        nmi = np.zeros((1,n_randstart))
        ari = np.zeros((1,n_randstart))
        fm = np.zeros((1,n_randstart))
        
        for j in range(n_randstart):
            
            start=clock()
            communities, graph, Q = pg.cluster(suface_marker_data_normalized, k=50, directed=False, prune=True, min_cluster_size=20, jaccard=True, primary_metric='euclidean', n_jobs=-1, q_tol=1e-3)
#            kmeans = KMeans(n_clusters=2, random_state=0).fit(suface_marker_data_normalized)          
            stop=clock()
            
#            Q = 1
#            communities = kmeans.labels_
#            c[:,j] = communities   #labels
            c[:,j] = communities
            q[:,j] = Q             #modularity
            t[:,j] = stop - start  #running time
            
#            calculate other validation parameters
            nmi[:,j] = normalized_mutual_info_score(ground_truth,communities)
            ari[:,j] = adjusted_rand_score(ground_truth,communities)
           
            temp = precision_recall_fscore_support(ground_truth,communities,average='weighted')
            fm[:,j] = temp[2]
            
        np.save(sample_name + "\\communities.npy",c)        
        np.save(sample_name + "\\modularity.npy",q)  
        np.save(sample_name + "\\runningtime.npy",t)
        np.save(sample_name + "\\truth.npy",ground_truth)
        np.save(sample_name + "\\nmi.npy",nmi)
        np.save(sample_name + "\\ari.npy",ari)
        np.save(sample_name + "\\fm.npy",fm)
        
        
        