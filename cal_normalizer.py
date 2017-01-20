# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 18:28:51 2017

@author: Dell
"""
import numpy as np
#import pandas as pd
import fcs_reader as fcsrd
import glob


fname = glob.glob('F:\\cytowork\\experiment_44185_files\\*.fcs')[0:18*5]
for i in range(len(fname)):
    if i == 0:
        
        meta, data_numpy = fcsrd.parse_fcs(fname[i], meta_data_only=False, output_format='ndarray', reformat_meta=True)
    
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
        
        
    else:
        meta, data_numpy = fcsrd.parse_fcs(fname[i], meta_data_only=False, output_format='ndarray', reformat_meta=True)      
        surface_marker_data = np.vstack((surface_marker_data,data_numpy[:,surface_marker.index - 1]))        
        
        
        

        
suface_marker_data_transformed = np.arcsinh(surface_marker_data/5)
quantile_995 = np.percentile(suface_marker_data_transformed,99.5,axis=0)
suface_marker_data_normalized = np.divide(suface_marker_data_transformed,np.tile(np.transpose(quantile_995),(np.shape(suface_marker_data_transformed)[0],1)))

np.save("normalizer.npy",quantile_995)