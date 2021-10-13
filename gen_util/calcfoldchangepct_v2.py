# Cacluates the log2 foldchange between different columns of a counts table
# Ex: python calcfoldchangepct.py counts.tsv

from sys import argv
from collections import defaultdict
import pandas as pd
import numpy as np

#script, counts_file, sections = argv
script, counts_file = argv

if counts_file.endswith('csv'):
	sep = ","
else:
	sep = "\t"
	
# Read counts table into pd.dataframe
counts_df = pd.read_csv(counts_file, sep=sep, index_col=0)

# Log normalize counts
counts_df = counts_df.transform(lambda x: np.log2(x + 1))

# Compare columns and append fold changes to new dataframe
foldchange_df = pd.DataFrame()
num_cols = len(counts_df.columns)

for n in range(0, num_cols-1, 2):
	col_name = counts_df.columns[n+1] + "_vs_" + counts_df.columns[n]
	foldchange_df[col_name] = counts_df[counts_df.columns[n+1]] - counts_df[counts_df.columns[n]]

foldchange_df.index = counts_df.index

# Print
foldchange_df.to_csv('l2fc_' + counts_file, sep="\t")
