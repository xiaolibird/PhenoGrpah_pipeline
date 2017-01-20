# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 15:40:36 2016

@author: Dell
"""

import numpy as np
#import pandas as pd
import fcs_reader as fcsrd
import glob
from itertools import product
from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_rand_score, adjusted_mutual_info_score, silhouette_score, homogeneity_completeness_v_measure
import itertools

def loadfcs():
    fname = glob.glob('*.fcs')[0]
    meta, data_numpy = fcsrd.parse_fcs(fname, meta_data_only=False, output_format='ndarray', reformat_meta=True)   
    data = meta['_channels_']
    a = (data['$PnS'].values)
    index = []
    
    for name in a:
        if "CD" in name:
            index.insert(0,name)
            
    index.insert(0,"abTCR")
    index.insert(0,"MHCII")
    index.insert(0,"B220")
    index.insert(0,"Areg")
    index.insert(0,"Ly6G-C")
    index.insert(0,"Ly6G")
    index.insert(0,"Ly6C")
    index.insert(0,"BST2")
    index.insert(0,"Ter119")
    index.insert(0,"F4-80")
    index.insert(0,"gdTCR")
    index.insert(0,"SiglecF")
    
    index = np.array(index)
    
    surface_marker = data[data['$PnS'].isin(index)]
    
    surface_marker_data = data_numpy[:,surface_marker.index - 1]
    
    suface_marker_data_transformed = np.arcsinh(5*surface_marker_data)
    
    quantile_995 = np.percentile(suface_marker_data_transformed,99.5,axis=0)
    
    suface_marker_data_normalized = np.divide(suface_marker_data_transformed,np.tile(np.transpose(quantile_995),(np.shape(suface_marker_data_transformed)[0],1)))
    suface_marker_data_normalized[suface_marker_data_normalized > 1] = 1
    return suface_marker_data_normalized


def cal_pairwise(cluster_channel):
    ari = []
    ami = []
    nmi = []
    
    rows, cols = np.shape(cluster_channel)
    
    rare_indexes = []
    for c in range(cols):
        mem = []
        for n in np.unique(cluster_channel[:,c]):
            n_cell = sum(cluster_channel[:,c] == n)
            if n_cell/len(cluster_channel[:,c]) < 0.01:
                mem.append(n)
        new_col = [False] * len(cluster_channel[:,c])
        for r in mem:        
            new_col += cluster_channel[:,c] == r
        rare_indexes.append(new_col)
        
    combs = tuple(itertools.combinations(range(cols),2))
    
    for c in combs:
        i,j = c
        
        ii = rare_indexes[i] + rare_indexes[j]
        print(i,j,'\n')
        t = adjusted_rand_score(cluster_channel[:,i][ii], cluster_channel[:,j][ii])
        ari.append(t)
        t = adjusted_mutual_info_score(cluster_channel[:,i][ii], cluster_channel[:,j][ii])
        ami.append(t)
        t = normalized_mutual_info_score(cluster_channel[:,i][ii], cluster_channel[:,j][ii])
        nmi.append(t)
    
    return ari,ami,nmi


def run():
#    data = loadfcs()
    labels = np.load("communities.npy")
    labels = labels.astype(np.int)
    ari,ami,nmi = cal_pairwise(labels[:,0:30])
    
    np.save("ari_k_15",ari)
    np.save("ami_k_15",ami)
    np.save("nmi_k_15",nmi)
    
    ari,ami,nmi = cal_pairwise(labels[:,30:60])
    
    np.save("ari_k_30",ari)
    np.save("ami_k_30",ami)
    np.save("nmi_k_30",nmi)
    
    ari,ami,nmi = cal_pairwise(labels[:,60:90])
    
    np.save("ari_k_45",ari)
    np.save("ami_k_45",ami)
    np.save("nmi_k_45",nmi)
    
    ari,ami,nmi = cal_pairwise(labels[:,90:120])
    
    np.save("ari_k_60",ari)
    np.save("ami_k_60",ami)
    np.save("nmi_k_60",nmi)
#sht = []
#entrp1 = []
#entrp2 = []
#entrp3 = []




if __name__ == '__main__':
    run()
#def run():
#    data = loadfcs()
#    labels = load()
#    print(data)
