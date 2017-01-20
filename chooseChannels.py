# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 20:58:22 2016

@author: Dell
"""
#==============================================================================
# choose surface markers
#==============================================================================
import numpy as np
#import pandas as pd
import fcs_reader as fcsrd
import glob
import os


#fname = glob.glob('F:\\cytowork\\experiment_44185_files\\*.fcs')[0:18*5]
#for i in range(len(fname)):
#    if i == 0:
#        
#        meta, data_numpy = fcsrd.parse_fcs(fname[i], meta_data_only=False, output_format='ndarray', reformat_meta=True)
#    
#        meta_data = meta['_channels_']   
#        channel_names = (meta_data['$PnS'].values)
#        selected_index = [] 
#        for name in channel_names:
#            if "CD" in name:
#                selected_index.insert(0,name)        
#        selected_index.insert(0,"HLA-DR")
#        selected_index = np.array(selected_index)
#        
#        surface_marker = meta_data[meta_data['$PnS'].isin(selected_index)]        
#        surface_marker_data = data_numpy[:,surface_marker.index - 1]
#        
#        
#    else:
#        meta, data_numpy = fcsrd.parse_fcs(fname[i], meta_data_only=False, output_format='ndarray', reformat_meta=True)      
#        surface_marker_data = np.vstack((surface_marker_data,data_numpy[:,surface_marker.index - 1]))        
#        
#        
#        
#
#        
#suface_marker_data_transformed = np.arcsinh(surface_marker_data/5)
#quantile_995 = np.percentile(suface_marker_data_transformed,99.5,axis=0)
#suface_marker_data_normalized = np.divide(suface_marker_data_transformed,np.tile(np.transpose(quantile_995),(np.shape(suface_marker_data_transformed)[0],1)))
#
#np.save("normalizer.npy",quantile_995)

#suface_marker_data_normalized[suface_marker_data_normalized > 1] = 1
#subsampled = np.random.randint(0,len(suface_marker_data_transformed),200000)
#suface_marker_data_transformed = suface_marker_data_transformed[np.unique(subsampled),:]
#==============================================================================
# do phenograph clustering
#==============================================================================
quantile_995 = np.load("normalizer.npy")
n_patients = 5 + 16;
n_conditions = 18;
import phenograph as pg
from time import clock
#for i in range(n_patients):
for i in [0,]:
    fname = glob.glob('F:\\cytowork\\experiment_44185_files\\*.fcs')[n_conditions*i:n_conditions*(i+1)]
    sample_name = fname[0].split('\\')
    sample_name = sample_name[-1];
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
        
#        transform and normalize
        suface_marker_data_transformed = np.arcsinh(surface_marker_data_selected/5)
        suface_marker_data_normalized = np.divide(suface_marker_data_transformed,np.tile(np.transpose(quantile_995),(np.shape(suface_marker_data_transformed)[0],1)))
#       slice ground truth provided by author
        ground_truth = np.array(["PhenoGraph"])
        ground_truth = meta_data[meta_data['$PnS'].isin(ground_truth)]  
        ground_truth = ground_truth.index - 1
        ground_truth = surface_marker_data[:,ground_truth]
        #==============================================================================
        # test 1 all data points 205496
        #==============================================================================
        c = np.zeros((np.shape(suface_marker_data_normalized)[0],30))
        q = np.zeros((1,30))
        t = np.zeros((1,30))
        for j in range(30):
            start=clock()
            communities, graph, Q = pg.cluster(suface_marker_data_normalized, k=50, directed=False, prune=True, min_cluster_size=20, jaccard=True, primary_metric='euclidean', n_jobs=-1, q_tol=1e-3)
            stop=clock()
            c[:,j] = communities
            q[:,j] = Q
            t[:,j] = stop - start
            

        
        from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_rand_score,precision_recall_fscore_support
        
            normalized_mutual_info_score()
            adjusted_rand_score()
            precision_recall_fscore_support()
        
        
        np.save(sample_name + "\\communities.npy",c)        
        np.save(sample_name + "\\modularity.npy",q)  
        np.save(sample_name + "\\runningtime.npy",t)
        np.save(sample_name + "\\truth.npy",ground_truth)
#finish=clock()
#print (finish-start)/1000000
#==============================================================================
# test 2 90% of all data points 205496*90
#==============================================================================
#c = np.zeros((np.shape(suface_marker_data_normalized)[0],120))
#q = np.zeros((1,120))
#
#count = 0
#for k in (15,30,45,60):
#    for i in range(30):
#        start=clock()
#        communities, graph, Q = pg.cluster(suface_marker_data_normalized,k=k)
#        stop=clock()
#        c[:,count] = communities
#        q[:,count] = Q
#        t[:,count] = stop - start
#        count = count + 1
#        
#np.save("communities_90.npy",c)        
#np.save("modularity_90.npy",q)
#np.save("runningtime_90.npy",t) 
#==============================================================================
# do t-SNE
#==============================================================================
#import matplotlib.pyplot as plt
#from matplotlib import offsetbox
#from sklearn import manifold
#
#
#print("Computing t-SNE embedding")
#tsne = manifold.TSNE(n_components=2, init='pca', random_state=0)
##t0 = time()
#X_tsne = tsne.fit_transform(suface_marker_data_transformed)
#
#plt.scatter(X_tsne[:,0], X_tsne[:,1], c=communities,s=10,alpha=0.5)
#plt.show()
                     
#==============================================================================
# plot
#==============================================================================


#==============================================================================
# SARA
#==============================================================================
