import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import seaborn as sns



##################################################

# Compare particle tracks over 3 days with different window resolutions.


# Comparison of first particle
fig, ax = plt.subplots(figsize=(12, 8))
ax.set_title(f'Seagrass Jan 01 2016')
for res, color in zip([0.01, 0.005, 0.004, 0.001, 0.0005], ['c', 'r', 'k', 'g', 'm']):
    res_str = str(res).split('.')[1]
    ds = xr.open_dataset(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\working_simulation\outputs\seagrass12particles_' + res_str + '.nc')
    ax.plot(ds.lon[0, :], ds.lat[0, :], 'o-', color=color, fillstyle='none', label=f'opendrift dx={res}deg',)
#ax.set_xlim([ds.lon[0].min().item() - 0.01, ds.lon[0].max().item() + 0.01])
#ax.set_ylim([ds.lat[0].min().item() - 0.01, ds.lat[0].max().item() + 0.01])
ax.legend()
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.show()
fig.savefig(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\OpenDrift_interpolation_1particle.pdf', bbox_inches='tight')

# Comparison of all 12 particles
fig, axs = plt.subplots(4,3)
axs = axs.reshape(12)
fig.subplots_adjust(hspace=1)
fig.suptitle('Seagrass Jan 01 2016')
for ax, drifter in zip(axs, range(12)):
    ax.set_title(f'Drifter {drifter}')
    for res, color in zip([0.01, 0.005, 0.004, 0.001, 0.0005], ['c', 'r', 'k', 'g', 'm']):
        res_str = str(res).split('.')[1]
        ds = xr.open_dataset(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\working_simulation\outputs\seagrass12particles_' + res_str + '.nc')
        ax.plot(
            ds.lon[drifter, :], ds.lat[drifter, :], label=f'dx={res}deg')
        ax.plot(ds.lon[drifter, 0].item(), ds.lat[drifter, 0].item(), marker='*', ms=10)
ax.legend()
plt.show()
fig.savefig(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\OpenDrift_interpolation_12particles.pdf', bbox_inches='tight')



#############################################

# compare particle difference to a baseline of 0.005 (assuming it is the most accurate since it is closest to the real resolution)


dts = ['01', '005', '004', '001', '0005']
fig, ax = plt.subplots(figsize=(10, 10))
ds = xr.open_dataset(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\working_simulation\outputs\seagrass12particles_005.nc')
time, lon, lat = [ds[key] for key in ['time', 'lon', 'lat']]

for dt in dts:
    ds = xr.open_dataset(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\working_simulation\outputs\seagrass12particles_' + dt + '.nc')
    d = np.sqrt((ds.lon[:, 1:] - lon[:, 1:])**2 + (ds.lat[:, 1:] - lat[:, 1:])**2)
    #d = d / np.sqrt((ds.lon[:, 1:] - ds.lon[:, 0])**2 + (ds.lat[:, 1:] - ds.lat[:, 0])**2)
    ax.plot(time[1:], d.mean(axis=0), label=f'{dt} deg')
ax.legend()
ax.set_title(f'Average separation distance for all particles when compared to 0.005 window resolution')
plt.show()
fig.savefig(r'C:\Users\jcristia\Documents\GIS\MSc_Projects\Hakai\scripts_dev_scratch\sensitivity_testing\window_res\OpenDrift_interpolation_12particles_distancediff.pdf', bbox_inches='tight')


#############################################