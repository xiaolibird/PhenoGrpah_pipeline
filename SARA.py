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

#from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_rand_score
#from sklearn.metrics import precision_recall_fscore_support


#from sklearn.cluster import KMeans, MiniBatchKMeans


from scipy import stats

def sara(basal,sti,niter=100):
    if len(basal) == 0 or len(sti) == 0:
        return 0
    D, p_val = stats.ks_2samp(basal,sti)
#    if p_val < 0.05:
#        print(p_val)
#    pool = np.hstack((basal,sti)) 
    pool = np.concatenate((basal,sti))
    xedges = np.linspace(min(pool),max(pool),100)
    
    hist_b, _ = np.histogram(basal,xedges)
    hist_s, _ = np.histogram(sti,xedges)
    
    hist_b = np.cumsum(hist_b)/len(basal)
    hist_s = np.cumsum(hist_s)/len(sti)
    
    emd = np.sum(np.abs(hist_s - hist_b))
    
    return np.sign(np.median(sti)-np.median(basal)) * emd * (1-p_val)
#quantile_995 = np.load("normalizer.npy")
n_healthy = 5
n_patients = 16
n_conditions = 18
import phenograph as pg
from time import clock

row_name = []
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
    
#   confirm control and compare
    seq = np.array(range(len(fname)))
    basal = np.array([2,3])
#    sti = seq.diff(basal)
    
#   get colum information from 1st basal fcs
    meta, data_numpy = fcsrd.parse_fcs(fname[basal[0]], meta_data_only=False, output_format='ndarray', reformat_meta=True)
    meta_data = meta['_channels_']   
    channel_names = (meta_data['$PnS'].values)
    selected_index = [] 
    for name in channel_names:
        if "p" in name and name != "PhenoGraph":
            
            selected_index.insert(0,name)                     

    selected_index = np.array(selected_index[1:])
    phospho_marker = meta_data[meta_data['$PnS'].isin(selected_index)]
    
    phospho_marker['$PnS']
    phospho_marker.index - 1

    basal_phospho_marker_data = np.arcsinh(data_numpy[:,phospho_marker.index - 1]/5)
#   slice ground truth provided by author
    basal_label = np.array(["PhenoGraph"])
    basal_label = meta_data[meta_data['$PnS'].isin(basal_label)]  
    basal_label = basal_label.index - 1
    basal_label_data = data_numpy[:,basal_label]
    
#    ground_truth = np.reshape(ground_truth, (np.product(ground_truth.shape),))
#   stack the rest of basal fcs       
    for b in basal[1:]:
        meta, data_numpy = fcsrd.parse_fcs(fname[b], meta_data_only=False, output_format='ndarray', reformat_meta=True)
        basal_phospho_marker_data = np.vstack((basal_phospho_marker_data,np.arcsinh(data_numpy[:,phospho_marker.index - 1]/5)))
        basal_label_data = np.vstack((basal_label_data,data_numpy[:,basal_label]))
#   clustering label
    basal_label_data = np.reshape(basal_label_data, (np.product(basal_label_data.shape),))
        
    sara_score = []
    n_sti = 0
    n_label = np.sum(np.unique(basal_label_data) >= 0)
    n_marker = len(phospho_marker.index)
    col_name = []
    
    for ii in seq: #16
        if sum( basal == ii ) == 0: 
#           exclude basal itself
            sti_name = fname[ii].split('\\')
            sti_name = sti_name[-1];
            sti_name = sti_name.split('_')
            sti_name = sti_name[2]
            if sti_name == 'Basal1' or sti_name == 'Basal2':
                sti_name = 'BEZ-235'

            n_sti = n_sti + 1
            meta, data_numpy = fcsrd.parse_fcs(fname[ii], meta_data_only=False, output_format='ndarray', reformat_meta=True)  
            sti_phospho_marker_data = np.arcsinh(data_numpy[:,phospho_marker.index - 1]/5)
            sti_label_data = data_numpy[:,basal_label]
            sti_label_data = np.reshape(sti_label_data, (np.product(sti_label_data.shape),))
            
#            sort label
            
            for marker in range(len(phospho_marker.index)): #15
                col_name.append(sti_name + "->" + channel_names[phospho_marker.index - 1][marker])
                
                for l in np.unique(basal_label_data): #33
                    if l >= 0 and sum(basal_label_data == l) >20:
                        
                        row_name.append(sample_name + "#C" + str(l))
                        temp_sara = sara(basal_phospho_marker_data[basal_label_data == l,marker],sti_phospho_marker_data[sti_label_data == l,marker])                       
                        print("sampe {} in cluster #{}  marker {} under perturbation {} has SARA {}".format(sample_name, str(int(l)),channel_names[phospho_marker.index - 1][marker],sti_name,temp_sara))  
                        sara_score.append(temp_sara)
                 
                 
#                    print("this label finished=======================================================")
#                print("this marker finished=======================================================")
#            print("this stimulation finished=======================================================")
#    sara_score = np.reshape(np.array(sara_score),(len(np.unique(sti_label_data)),len(basal_label),len))            
    sara_score = np.array(sara_score)
    sara_socre_r = np.zeros(len(sara_score))
    sara_socre_r = np.reshape(sara_score,[n_label,n_sti*n_marker])
    sara_score = np.reshape(sara_score,[n_sti, n_marker, n_label])
    
    for dd in range(n_label):
        sara_socre_r[dd,:] = (sara_score[:,:,dd].ravel())
                   
#    save every sample SARA scores
    if os.path.exists(sample_name) == False:
        os.mkdir(sample_name)
        
#    row_name = np.unique(row_name)
#    col_name
    np.save(sample_name +"\\sara_scores.npy",sara_socre_r)
#    np.save(sample_name +"\\abundance.npy",abundance)

for i in range(n_healthy + n_patients):
#for i in range(1):
    fname = glob.glob('F:\\cytowork\\experiment_44185_files\\*.fcs')[n_conditions*(i):n_conditions*(i+1)]
#    create folder 
    sample_name = fname[0].split('\\')

    sample_name = sample_name[-1];
    sample_name = sample_name.split('_')
    sample_name = sample_name[0]
    print(sample_name)
    
    gate = []
    
    if i == 0:
        sara_score = np.load(sample_name +"\\sara_scores.npy")
        centroids = np.load(sample_name +"\\centroids.npy")
        gate.append(sara_score.shape[0])
    else:
        temp = np.load(sample_name +"\\sara_scores.npy")
        gate.append(temp.shape[0])
        sara_score = np.vstack((sara_score, np.load(sample_name +"\\sara_scores.npy")))
        centroids = np.vstack((centroids, np.load(sample_name +"\\centroids.npy")))
        
    centroids = np.array(centroids)
    shape = sara_score.shape
    sara_score = stats.zscore(sara_score.flatten())
    sara_score = np.reshape(sara_score,shape)


import seaborn as sns
ax = sns.heatmap(centroids,vmin=0.1,vmax=1,yticklabels=False,xticklabels=False,cmap='Spectral_r')
