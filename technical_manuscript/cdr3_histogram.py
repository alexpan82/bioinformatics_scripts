# Takes a mixcr clones.txt file and 'heavy' or 'light' as input
# Creates a histogram counting the number of (heavy or light chain) clones vs read coverage

import matplotlib.pyplot as plt
import pandas as pd
from sys import argv

script, mixcr_input, chain = argv


df = pd.read_csv(mixcr_input, sep="\t")

# Filter by chain
if argv[2] == 'light':
	chain_df = df.loc[df.allVHitsWithScore.str.contains('IGK') | df.allVHitsWithScore.str.contains('IGL')].copy()

elif argv[2] == 'heavy':
	chain_df = df.loc[df.allVHitsWithScore.str.contains('IGH')].copy()

else:
	print("Please input heavy or light as 2nd argument")
	quit()

del df
df = chain_df
del chain_df			

plt.rcParams.update({'font.size': 12})
hist = df['cloneCount'].hist(bins=int(list(df['cloneCount'])[0]/2))
plt.title(mixcr_input.split('.')[0] + ' Number %s clones: %s' % (chain, df.shape[0]))
#plt.grid(False)
plt.ylim([0, 30])
fig = hist.get_figure()
fig.savefig(mixcr_input.split('.')[0] + '.%s.hist.png' % chain)
