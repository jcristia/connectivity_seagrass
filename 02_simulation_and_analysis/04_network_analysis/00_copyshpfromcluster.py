# create text to copy back shapefiles from cluster so that I can just paste it all in at once

import os
import pandas as pd

path_local = r'Hakai\script_runs\seagrass\{}\seagrass_{}\outputs'
sg_folder_count = 9
path_remote = 'jcristia@zoology.ubc.ca:flex/runs/{}/seagrass_{}/outputs/shp'
runs = [
    'seagrass_20200327_SS201401',
    'seagrass_20200327_SS201405',
    'seagrass_20200327_SS201408'
    ]

# create paths
paths_l = []
paths_r = []
for run in runs:
    for s in range(1, sg_folder_count+1):
        l = path_local.format(run, s)
        paths_l.append(l)
        r = path_remote.format(run, s)
        paths_r.append(r)

# create commands
commands = ['cd d:']
for pl,pr in zip(paths_l, paths_r):
    command = '''scp -r {} {}'''.format(pr, pl.replace('\\', '/'))
    commands.append(command)
df_commands = pd.DataFrame(commands)
    
# output dataframe to txt, without index or headers
df_commands.to_csv(r'commands.txt', header=None, index=None)

# text file is saved where this script is
# open, copy text and run in bash shell