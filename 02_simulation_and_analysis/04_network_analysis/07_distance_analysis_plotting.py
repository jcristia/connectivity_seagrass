# plot connection probability vs. distance

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp

# euclidean distances for all combinations of nodes
distances = r"C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\distance_analysis\distance_analysis_mapping\euc_lines_ALL.csv"
# overall averaged established connections
conns = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_runs_cluster\seagrass\output_figs_SALISHSEA_ALL\connectivity_average_ALL.shp'

# There were a few meadows that did not have their ocean distances calculated in the previous script, but they did not throw an error so I did not know it at the time.
# Identify the uIDs of the meadows that did not have their distances calculated. There should be 5.
# It looks like they weren't calculated because they were land locked by the rasterization, but still in water.
# I have enough points to establish the general relationship.
# There are 5 in total, randomly distributed in space.

# calc frequency and look for singletons
dist_df = pd.read_csv(distances)
dist_df_freq = dist_df.origin_id.value_counts()

# uIDs not considered:
uIDs_not = [720, 192, 569, 762, 867]
# NOTE: 720 was not calculated at all, so it won't have a self connection

# drop these values and self connections from the dataframe
dist_df = dist_df[dist_df.DestID != dist_df.origin_id]

# access my overall averaged connection dataset
# drop self connections
# drop any lines that are to/from the ones I couldn't calculate
conns_df = gpd.read_file(conns)
conns_df = conns_df.drop(['geometry', 'dateSEA', 'dateYR', 'date'], 1)
conns_df = conns_df[conns_df.from_id != conns_df.to_id]
conns_df = conns_df[~conns_df.from_id.isin(uIDs_not)]
conns_df = conns_df[~conns_df.to_id.isin(uIDs_not)]
conns_df = pd.DataFrame(conns_df) # geopandas to pandas

# use time_int to categorize by PLD
conns_df.loc[conns_df.timeintavg < 2880, 'pld'] = 60
conns_df.loc[conns_df.timeintavg < 1008, 'pld'] = 21
conns_df.loc[conns_df.timeintavg < 336, 'pld'] = 7
conns_df.loc[conns_df.timeintavg < 144, 'pld'] = 3
conns_df.loc[conns_df.timeintavg < 48, 'pld'] = 1

# pandas merge, keep all records from distance dataframe
df_merge = dist_df.merge(conns_df, how='left', left_on=['origin_id', 'DestID'], right_on=['from_id', 'to_id'])

# for any distance combinations that do not have a connection strength, fill as 0
df_merge.probavgm = df_merge.probavgm.fillna(0)
# convert to km
df_merge['distkm'] = (df_merge.Shape_Leng)/1000.0

# create dfs with and without zeros and also a big one with duplicates of connections so that I can compare with and without zero values in one plot with seaborn
df_merge['withzeros'] = 'yes'
df_nozero = df_merge[df_merge.probavgm > 0]
df_nozero.withzeros = 'no'
df_concat = pd.concat([df_merge, df_nozero]).reset_index(drop=True)




############
# PLOTTING
############


df_nozero['probperc'] = df_nozero.probavgm * 100
df_nozero['probperclog'] = np.log10(df_nozero.probperc)
def func(x, a, b):
    return a * x**(b)
popt, pcov = curve_fit(func, df_nozero.distkm, df_nozero.probperclog)
# "Use non-linear least squares to fit a function, f, to data."
print(popt) # to see a,b

px = np.linspace(0,200, 300)
py = a * px**(b)
nom = unp.nominal_values(py)
std = unp.std_devs(py)

sns.set()
sns.set_style('white')
sns.set_context('paper')
f = sns.lmplot(
    x="distkm", 
    y="probperclog", 
    data=df_nozero, 
    hue='pld',
    hue_order=[60,21,7,3,1], 
    scatter=True, 
    fit_reg=False, 
    scatter_kws={"s": 1, 'alpha':1},
    legend=True,
    legend_out=False,
    ) # plot points
plt.plot(px, nom, 'dimgray') # plot fitted curve

# The access to the legend is different with lmplot.
# I need to order my legend in a certain way, so I need to do this manually
# unfortunately. 
pal = sns.color_palette() # I'm just using the default
pal_hex = pal.as_hex()[:5]
pal_hex.reverse()
handles = []
labels = ['1', '3', '7', '21', '60']
import matplotlib.lines as mlines
for h, l in zip(pal_hex, labels):
    blue_line = mlines.Line2D([], [], color=h, linestyle='None', marker='o', markersize=2, label=l)
    handles.append(blue_line)
plt.legend(title='PD (days)', frameon=False, handles=handles)

f.set(xlim=(0,200))
f.set(xlabel='Distance (km)', ylabel=r'log$_{10}$ Connection probability (%)')
f.savefig(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\publications_figures\chap1\fig08_conn_v_dist_log.svg')
###########

# get r squared
residuals = df_nozero.probperclog- func(df_nozero.distkm, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((df_nozero.probperclog-np.mean(df_nozero.probperclog))**2)
r_squared = 1 - (ss_res / ss_tot)
print(r_squared)

# EQUATION
# y = -0.52(x^0.38)
# r2 = 0.44

