# Determine PLD intervals for biophysical model

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv('plds.csv')

# drop species that do not have PLD data
df = df.dropna(subset=['pld_days'])

# for any species with a PLD of zero, we will consider it capable of rafting on eelgrass. We are using a rafting duration of 21 days based on studies of how long eelgrass stays buoyant.
# select rows where pld_days is 0 and just the pld_days column, change to 21
df.loc[(df.pld_days == 0), 'pld_days'] = 21

# add in row for a pld of 1 day. I want a low value for any species that will have chance dispersal for a few hours.
df = df.append({'common': '1daygeneral', 'pld_days': 1}, ignore_index=True)


# plot ln(pld) against pld
df['pld_days_ln'] = np.log(df['pld_days'])

dims = (8.5, 5.5)
sns.set()
sns.set_style('white')
sns.set_context('paper')
f, ax = plt.subplots(figsize=dims) 
sxy = sns.scatterplot(x='pld_days', y='pld_days_ln', data=df, ax=ax)
sxy.set(xlabel='PD (days)', ylabel='ln PD (days)')
ax.xaxis.set_ticks(np.arange(0, np.around(max(df['pld_days']), -1) + 10, 10)) # -1 rounds it to nearest 10
ax.set_xlim(0)
ax.set_ylim(0)

# add in points for ln(x), where x is an integer
yint_mid = [0, 1, 2, 3, 4]
xint_ln_mid = np.exp(yint_mid)
labels = ['1 day', '~3 days', '~7 days', '~3 weeks', '~ 2 months']
plt.scatter(xint_ln_mid, yint_mid, color='red', zorder=10, clip_on=False) # zorder and clip make sure point on axis is placed on top
for i, txt in enumerate(labels):
    ax.annotate(txt, (xint_ln_mid[i], yint_mid[i]),
                textcoords="offset points", # how to position text
                xytext=(17,9), # distance from point to text
                ha='center',
                weight='bold') # horizontal alignment

# create mid point lines for bins
yint = [0.5, 1.5, 2.5, 3.5]
xint_ln = np.exp(yint)
plt.vlines(xint_ln, 0, yint, linestyle='dashed', linewidth=1, color='grey')
plt.hlines(yint, 0, xint_ln, linestyle='dashed', linewidth=1, color='grey')

#leave for reference to know how to fill in between lines:
#ax.fill_between(xint_ln, 0, yint, where=xint_ln<xint_ln[4], color='red')

plt.show()

f.savefig('pds.svg', bbox_inches='tight')

# interpretation:
#There is a lot of variation in my PLD values. Some short PLDs, and some that are very large. Since I am not running a simulation for each individual, and since there is not a lot of info on exact PLDs, then I will want to create bins of PLDs and match species to those. However, to create bins of equal size without having to do too many to cover the full range, I will need to create bins on a log scale. Therefore, I will use a bin size of e^x. So EXP(0,1,2,3,4). This equates to PLDs of (1, 2.7, 7.4, 20.1, 54.6). I will round these to easily referenced times (1, 3, 7, 21, 60).
#This is convenient because my PLD data is somewhat clustered around these values, and now I have a very easily understood bin that is supported mathematically and biologically.
#I will max at 60 days partly for computational efficiency reasons, but also because in the Salish Sea almost everything will have stranded or made it to the model boundary by this time.

