# % of weighted habitat area connected by inter patch movement vs. PLD

import os
import pandas as pd
import geopandas as gp
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from math import log10, floor
import statsmodels


###############################################

dirs = [
    'seagrass_20200228_SS201701',
    'seagrass_20200309_SS201705',
    'seagrass_20200309_SS201708',
    'seagrass_20200310_SS201101',
    'seagrass_20200310_SS201105',
    'seagrass_20200310_SS201108',
    'seagrass_20200327_SS201401',
    'seagrass_20200327_SS201405',
    'seagrass_20200327_SS201408',
    ]

root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_runs_cluster\seagrass'
plds = ['01', '03', '07', '21', '60']
overall_indices = r'conefor\conefor_connectivity_pld{}\overall_indices.txt'

eca_intra = 20.11273
# calculated by squaring the area of each patch, adding, then taking the square root

df_indices = pd.DataFrame()
for pld in plds:
    for dir in dirs:
        indices = os.path.join(root, dir, overall_indices.format(pld))
        df_ind = pd.read_csv(indices, sep='\t', names=['conefor_index', 'value'])
        df_ind['pld'] = int(pld)
        df_ind = df_ind.pivot(index='pld', columns='conefor_index', values='value').reset_index()
        df_ind['datepld'] = dir[20:] + pld
        df_indices = df_indices.append(df_ind)

df_indlog = df_indices
df_indlog['pldlog'] = np.log(df_indices.pld) # log pld
df_indlog['season'] = df_indlog.datepld.str[4:6].astype(int) # color by season
df_indlog['season'] = df_indlog.season.replace([1,5,8],['winter', 'spring', 'summer'])

# create new column and calculate % increase
df_indlog['ec_perc_inc'] = ((df_indlog['EC(PC)'] - eca_intra) /  df_indlog['EC(PC)']) *100

sns.set()
sns.set_style('white')
sns.set_context('paper')
colors = ['#377eb8', '#4daf4a', '#ff7f00']
sns.set_palette(colors)

f = sns.lmplot(x="pldlog", y="ec_perc_inc", data=df_indlog, fit_reg=False, hue='season', x_jitter=0.02, legend='full')
f._legend.remove()
sns.lineplot(x='pldlog', y='ec_perc_inc', data=df_indlog, ci=95, err_style='band', color='grey')
# band is a 95% confidence interval for the means
f.set(xlabel='ln PD (days)', ylabel='% of weighted habitat area connected \n by inter patch movement')
sns.despine(top=False, right=False)
f.savefig(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\publications_figures\chap1\fig04_ECAnumLOG_perc.svg')

# interpretation:
# possible max EC value ~39-45. This is when every patch either gives all of its
# particles out equally, or it gives all of its particles to the largest patch, 
# and the largest patch gives all of its particles to 2nd largest patch. Both
# scenarios are obviously unrealistic, especially since some patches dont have
# enough particles to connect to every other patch anyways.

# Notes on y-axis label:
# For someone that knows the PC
# metric it might take more explaining since it simplifies it a bit:
# Since the maximum of ECA in my directed graph is double the ECAintra,
# then the percent of ECA attributable to interconnectivity is equivalent to the
# percent of habitat area (weighted) that is connected by dispersal (interpatch)

# the key here is that it is not the COUNT of habitat connected. It is the
# equivalent AREA connected. This is another great reason for why we should add
# area as a attribute. Saying that x% of meadows are connected may not matter if
# that percentage only includes small meadows and all the large meadows are
# isolated.

########################