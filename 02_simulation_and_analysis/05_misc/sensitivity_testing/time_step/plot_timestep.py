

##################################################

# distribution of January data:
import numpy as np
import netCDF4 as nc
file_salish_cluster = r'D:\Hakai\models\salishsea\salishseacast_20160101_20160301\forcing\SalishSea_1h_20160101_20160301_opendrift.nc'
d = nc.Dataset(file_salish_cluster)
u = d.variables["u"]
v = d.variables["v"]
current = np.sqrt(u[:]**2 + v[:]**2)
np.median(current[:])
np.quantile(current, 0.75)
np.quantile(current, 0.99)
'''
99% of values are below 1.14m/s.
With a 440m grid size, that means that anything less than ~383 seconds time step could still be accurate.
'''

# runtime:
# 5s: 111m
# 15s: 58m
# 30s: 34m
# 60s: 18m
# 90s: 14m
# 180s: 12m
# 300s: 11m



##################################################

import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import seaborn as sns


# Compare particle tracks over 3 days with different time steps.

# Comparison of first particle
fig, ax = plt.subplots(figsize=(12, 8))
ax.set_title(f'Seagrass Jan 01 2016')
for ts, color in zip([5, 15, 30, 60, 90, 180, 300], ['c', 'r', 'k', 'g', 'm', 'b', 'y']):
    ts_str = str(ts)
    ds = xr.open_dataset(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\time_step\working_simulation\outputs\seagrass12particles_' + ts_str + 's.nc')
    ax.plot(ds.lon[0, :], ds.lat[0, :], 'o-', color=color, fillstyle='none', label=f'opendrift {ts}s',)
ax.legend()
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.show()
fig.savefig(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\time_step\OpenDrift_interpolation_1particle.pdf', bbox_inches='tight')

# Comparison of all 12 particles
fig, axs = plt.subplots(4,3)
axs = axs.reshape(12)
fig.subplots_adjust(hspace=1)
fig.suptitle('Seagrass Jan 01 2016')
for ax, drifter in zip(axs, range(12)):
    ax.set_title(f'Drifter {drifter}')
    for ts, color in zip([5, 15, 30, 60, 90, 180, 300], ['c', 'r', 'k', 'g', 'm', 'b', 'y']):
        ts_str = str(ts)
        ds = xr.open_dataset(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\time_step\working_simulation\outputs\seagrass12particles_' + ts_str + 's.nc')
        ax.plot(
            ds.lon[drifter, :], ds.lat[drifter, :], label=f'{ts}s')
        ax.plot(ds.lon[drifter, 0].item(), ds.lat[drifter, 0].item(), marker='*', ms=10)
ax.legend()
plt.show()
fig.savefig(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\time_step\OpenDrift_interpolation_12particles.pdf', bbox_inches='tight')


#############################################

# compare particle difference to a baseline of 5s (assuming it is the most accurate)

dts = ['5', '15', '30', '60', '90', '180', '300']
fig, ax = plt.subplots(figsize=(10, 10))
ds = xr.open_dataset(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\time_step\working_simulation\outputs\seagrass12particles_5s.nc')
time, lon, lat = [ds[key] for key in ['time', 'lon', 'lat']]

for dt in dts:
    ds = xr.open_dataset(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\time_step\working_simulation\outputs\seagrass12particles_' + dt + 's.nc')
    d = np.sqrt((ds.lon[:, 1:] - lon[:, 1:])**2 + (ds.lat[:, 1:] - lat[:, 1:])**2)
    #d = d / np.sqrt((ds.lon[:, 1:] - ds.lon[:, 0])**2 + (ds.lat[:, 1:] - ds.lat[:, 0])**2)
    ax.plot(time[1:], d.mean(axis=0), label=f'{dt} s')
ax.legend()
ax.set_title(f'Average separation distance for all particles when compared to 5s time step')
plt.show()
fig.savefig(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\time_step\OpenDrift_interpolation_12particles_distancediff.pdf', bbox_inches='tight')


#############################################