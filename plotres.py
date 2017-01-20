# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 19:01:19 2016

@author: Dell
"""

import numpy as np

ami_k_15 = np.load("ami_k_15.npy")
ami_k_30 = np.load("ami_k_30.npy")
ami_k_45 = np.load("ami_k_45.npy")
ami_k_60 = np.load("ami_k_60.npy")

nmi_k_15 = np.load("nmi_k_15.npy")
nmi_k_30 = np.load("nmi_k_30.npy")
nmi_k_45 = np.load("nmi_k_45.npy")
nmi_k_60 = np.load("nmi_k_60.npy")

ari_k_15 = np.load("ari_k_15.npy")
ari_k_30 = np.load("ari_k_30.npy")
ari_k_45 = np.load("ari_k_45.npy")
ari_k_60 = np.load("ari_k_60.npy")

b = np.concatenate((ami_k_15,ami_k_30,ami_k_45,ami_k_60,nmi_k_15,nmi_k_30,nmi_k_45,nmi_k_60,ari_k_15,ari_k_30,ari_k_45,ari_k_60),axis=0)
#b = np.transpose(b)
b.shape=(12,435)
b = np.transpose(b)
import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
sns.set(style="whitegrid")
df = pd.DataFrame(b,columns=["AMI k=15","AMI k=30","AMI k=45","AMI k=60","NMI k=15","NMI k=30","NMI k=45","NMI k=60","ARI k=15","ARI k=30","ARI k=45","ARI k=60"])
# Load the example dataset of brain network correlations
#df = sns.load_dataset("brain_networks", header=[0, 1, 2], index_col=0)

# Pull out a specific subset of networks
#used_networks = [1, 3, 4, 5, 6, 7, 8, 11, 12, 13, 16, 17]
#used_columns = (df.columns.get_level_values("network")
#                          .astype(int)
#                          .isin(used_networks))
#df = df.loc[:, used_columns]

# Compute the correlation matrix and average over networks
#corr_df = df.corr().groupby(level="network").mean()
#corr_df.index = corr_df.index.astype(int)
#corr_df = corr_df.sort_index().T

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(11, 6))

# Draw a violinplot with a narrower bandwidth than the default
#sns.violinplot(data=df, palette="Spectral_r", bw=.2, cut=1, linewidth=1)
sns.boxplot(data=df, palette="Spectral_r",  linewidth=1.2)

# Finalize the figure
#ax.set(ylim=(-.7, 1.05))
#sns.despine(left=True, bottom=True)