# compare runs with increasing particle count

# release from 12 meadows
# 84 releases
# base amount of particles: 629 per release
# 6 runs with increasing particle count: 1x, 2x, 4x, 8x, 16x, 32x

# comparison:
# the increase in total number of connections by particle count
# proportional increase in number of outgoing connections per patch
# strength of new connections at each time step

import os
import pandas as pd
import geopandas as gp
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

fig_save = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\particle_num'
root_folder = r'D:\Hakai\script_runs\seagrass\seagrass_20200205'
file = r'seagrass_20200205_{}\outputs\shp\connectivity_sg1_pld45.shp'
factor = [1, 2, 4, 8, 16, 32]
conn = []
for f in range(1,len(factor)+1):
    conn.append(os.path.join(root_folder, file.format(f)))
   

########
# increase in total number of connections compared to factor
########

fcount = []
for shp in conn:
    gdf = gp.read_file(shp)
    l = len(gdf)
    fcount.append(l)
sb_sc = sns.regplot(x=factor, y=fcount, order=2, ci=None)
sb_sc.set(xlabel = "Particle factor", ylabel = "Connection count")
plt.show()
fig = sb_sc.get_figure()
fig.savefig(os.path.join(fig_save, 'connections_total'))


########
# average increase in number of outgoing connections per patch
########

for shp, f in zip(conn, factor):
    gdf = gp.read_file(shp)
    gdf = gdf.groupby('from_id').agg(
        fromid_count = ('from_id','count')
        ).reset_index()
    gdf = gdf.astype({'from_id':int})
    gdf = gdf.rename(columns={'fromid_count':'count_' + str(f)})
    if f==1:
        gdf_all = gdf.astype({'from_id':int})
    else:
        gdf_all = pd.merge(gdf_all, gdf, 'outer', on='from_id')

# Plots:
# average % increase of connections per step for all meadows (e.g. count_2 - count_1 / count_1), so just 4 data points total
cols = gdf_all.columns[1:]
gdf_diff = pd.DataFrame()
for i in range(len(cols)-1):
    gdf_diff["step{}".format(i)] = (gdf_all[cols[i+1]] - gdf_all[cols[i]]) / gdf_all[cols[i]]

gdf_t = gdf_diff.transpose().reset_index()
gdf_diff_mean = gdf_diff.mean().to_frame().reset_index()
gdf_t['mean'] = gdf_diff_mean[0]
# NOTE: this isn't the best way to display the mean on the graph, but it will do for now
gdf_m = gdf_t.melt('index', var_name='cols', value_name='vals')
sb_ln = sns.pointplot('index', 'vals', hue='cols', data=gdf_m, palette=['grey','grey','grey','grey','grey','grey','grey','grey','grey','grey','grey','grey', 'red'])
sb_ln.set(xlabel = "step x 2", ylabel = "prop increase of connections per meadow")
sb_ln.legend().remove()
plt.show()
fig = sb_ln.get_figure()
fig.savefig(os.path.join(fig_save, 'connections_prop_increase'))



########
# strength of new connections at each time step
########

for shp, f in zip(conn, factor):
    gdf = gp.read_file(shp)
    gdf = gdf[['from_id', 'to_id', 'prob']]
    gdf['factor'] = f
    if f==1:
        gdf_all = gdf
    else:
        gdf_all = gdf_all.append(gdf)

gdf_all = gdf_all.astype({'from_id':int, 'to_id':int}) # there still a mix of datatypes in the columns for some reason. This was super important to do or else the code below didn't recognize duplicates.
gdf_all = gdf_all.set_index(['factor', 'from_id', 'to_id'])
gdf_us = gdf_all.unstack(level='factor')

gdf_std = gdf_all.groupby(['from_id', 'to_id']).agg(
    freq=('prob','count'),
    prob_std=('prob','std')
    )
gdf_std

# from gdf_us, get connections that did not exist in the previous factor
df = pd.DataFrame()
cols = gdf_us.columns
for i in range(1,len(cols[1:])+1):
    df['step_' + str(i)] = np.where(gdf_us[cols[i-1]].isna(), gdf_us[cols[i]], np.NaN)
df.mean()
np.nanmean(df['step_5'])
df.median() # this shows what I am looking for: a clear drop between 2 and 4x and then things stay low
np.nanmedian(df['step_5'])  # I don't think the nanmean makes a difference
#hist = df.hist()
#plt.show()

df = df.reset_index()
df_m = df.melt('index', var_name='cols', value_name='vals')
fig, ax = plt.subplots()
ax.set(yscale="log")
sb_bx = sns.boxplot('cols', 'vals', ax=ax, data=df_m)
sb_bx.set(xlabel = "simulation", ylabel = "strength of new connections")
plt.show()
fig = sb_bx.get_figure()
fig.savefig(os.path.join(fig_save, 'connections_new_strength'))