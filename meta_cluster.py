# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 14:48:32 2017

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


for i in range(n_healthy + n_patients):
#for i in range(1):
    fname = glob.glob('F:\\cytowork\\experiment_44185_files\\*.fcs')[n_conditions*(i):n_conditions*(i+1)]
#    create folder 
    sample_name = fname[0].split('\\')

    sample_name = sample_name[-1];
    sample_name = sample_name.split('_')
    sample_name = sample_name[0]
    print(sample_name)
#    os.mkdir(sample_name)
    sel = range(n_conditions)
    for ii in sel:
        if ii == sel[0]:           
            meta, data_numpy = fcsrd.parse_fcs(fname[ii], meta_data_only=False, output_format='ndarray', reformat_meta=True)
#            surface_marker_data = data_numpy[:,surface_marker.index - 1]
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
            surface_marker_data = data_numpy[:,surface_marker.index - 1]
            surface_marker_data = np.arcsinh(surface_marker_data/5)
            surface_marker_data = np.divide(surface_marker_data,np.tile(np.transpose(quantile_995),(np.shape(surface_marker_data)[0],1)))
                     
        else:
            meta, data_numpy = fcsrd.parse_fcs(fname[ii], meta_data_only=False, output_format='ndarray', reformat_meta=True)      
            data_numpy = np.arcsinh(data_numpy[:,surface_marker.index - 1]/5)
            data_numpy = np.divide(data_numpy,np.tile(np.transpose(quantile_995),(np.shape(data_numpy)[0],1)))
            surface_marker_data = np.vstack((surface_marker_data,data_numpy)) 
            
    
#    print(surface_marker_data.shape[0])     
    
    

    
#       transform and normalize
#    suface_marker_data_transformed = np.arcsinh(surface_marker_data_selected/5)

#    suface_marker_data_normalized = np.divide(suface_marker_data_transformed,np.tile(np.transpose(quantile_995),(np.shape(suface_marker_data_transformed)[0],1)))

#    for x in range(np.shape(suface_marker_data_transformed)[0]):
#        
#        suface_marker_data_transformed[x,:] = suface_marker_data_transformed[x,:]/quantile_995
#    
#    suface_marker_data_normalized = suface_marker_data_transformed
    for ii in sel:
        if ii == sel[0]:           
            meta, data_numpy = fcsrd.parse_fcs(fname[ii], meta_data_only=False, output_format='ndarray', reformat_meta=True)    

            cluster_label_channel = data_numpy
                     
        else:
            meta, data_numpy = fcsrd.parse_fcs(fname[ii], meta_data_only=False, output_format='ndarray', reformat_meta=True)      
            cluster_label_channel = np.vstack((cluster_label_channel,data_numpy))    
    
    
#       slice ground truth provided by author
    ground_truth = np.array(["PhenoGraph"])
    ground_truth = meta_data[meta_data['$PnS'].isin(ground_truth)]  
    ground_truth = ground_truth.index - 1
    ground_truth = cluster_label_channel[:,ground_truth]
    ground_truth = np.reshape(ground_truth, (np.product(ground_truth.shape),))
        

    print("sample {} has {} cells and {} clusters in total".format(sample_name,str(surface_marker_data.shape[0]),str(sum(np.unique(ground_truth)>=0))))
        
        
        
#   get centroid cluster
    flag = False 
#    print('sample #' + sample_name + ' has ' + str(len(np.unique(ground_truth))) + ' clusters')
    for l in np.unique(ground_truth):
        if l >= 0:
            temp_c = np.median(surface_marker_data[ground_truth == l ,:],0)
            temp_a = np.sum(ground_truth==l)
            
#            print('cluster #' + str(l) + ' has ' + str(temp_a) + ' cells')
            if temp_a > 20:
                if flag == False: 
                    centroids = temp_c
                    abundance = [temp_a]
                    
                    flag = True
                else:
                    centroids = np.vstack((centroids,temp_c))
                    abundance.append(temp_a)

                
                
    del ground_truth
    del surface_marker_data           
    
    
#    save every sample centroids and cell abundance
    if os.path.exists(sample_name) == False:
        os.mkdir(sample_name)
    np.save(sample_name +"\\centroids.npy",centroids)
    np.save(sample_name +"\\abundance.npy",abundance)
    
print("Finished centroids calculating and saving results")    
#       clustering all data for one sample
sample_names = []
for i in range(n_patients):
#for i in range(1):
    fname = glob.glob('F:\\cytowork\\experiment_44185_files\\*.fcs')[n_conditions*(i+n_healthy):n_conditions*(i+1+n_healthy)]
#    create folder 
    sample_name = fname[0].split('\\')

    sample_name = sample_name[-1];
    sample_name = sample_name.split('_')
    sample_name = sample_name[0]
    sample_names.append(sample_name)        
    
# read all centroids
flag = False    
for sample_name in sample_names:    
    if flag == False:
        centroids = np.load(sample_name+"\\centroids.npy")
        
#        print(centroids.shape[0])
        flag = True
    else:
        centroids = np.vstack((centroids,np.load(sample_name+"\\centroids.npy")))
#        print(centroids.shape[0])
        
        

n_randstart = 30
c = np.zeros((np.shape(centroids)[0],n_randstart))
q = np.zeros((1,n_randstart))
t = np.zeros((1,n_randstart))
#        nmi = np.zeros((1,n_randstart))
#        ari = np.zeros((1,n_randstart))
#        fm = np.zeros((1,n_randstart))
        
for j in range(n_randstart):
    print("=======================================================================")
    start=clock()
    communities, graph, Q = pg.cluster(centroids, k=15, directed=False, prune=True, min_cluster_size=2, jaccard=True, primary_metric='euclidean', n_jobs=-1, q_tol=1e-4)
#            kmeans = KMeans(n_clusters=2, random_state=0).fit(suface_marker_data_normalized)          
    stop=clock()
    
#            communities = kmeans.labels_
#            c[:,j] = communities   #labels
    c[:,j] = communities
    q[:,j] = Q             #modularity
    t[:,j] = stop - start  #running time
            
#            calculate other validation parameters
#            nmi[:,j] = normalized_mutual_info_score(ground_truth,communities)
#            ari[:,j] = adjusted_rand_score(ground_truth,communities)
#           
#            temp = precision_recall_fscore_support(ground_truth,communities,average='weighted')
#            fm[:,j] = temp[2]
print("Finished meta Clustering and saving results")          
np.save(sample_name + "\\meta_communities.npy",c)        
np.save(sample_name + "\\meta_modularity.npy",q)  
np.save(sample_name + "\\meta_runningtime.npy",t)
#        np.save(sample_name + "\\meta_truth.npy",ground_truth)
#        np.save(sample_name + "\\meta_nmi.npy",nmi)
#        np.save(sample_name + "\\meta_ari.npy",ari)
#        np.save(sample_name + "\\meta_fm.npy",fm)
#        
#        
#print(np.unique(communities) )

communities = c[:,np.abs(q.ravel()-max(q.ravel()))<10e-9]
communities = communities[:,0]
communities = communities.ravel()

    
mc_centroids = np.zeros((int(max(communities)),centroids.shape[1]))
for label in np.unique(communities):
    if label >= 0 :
        mc_centroids[label-1,:] = np.median(centroids[communities.ravel()==label,:],0)


import seaborn as sns
ax = sns.heatmap(mc_centroids.T,vmin=0.05,vmax=0.8,linewidths=0.005, linecolor="black",xticklabels=['MC{0}'.format(i+1) for i in range(int(max(communities))) ],yticklabels=selected_index,cmap='Spectral_r')      