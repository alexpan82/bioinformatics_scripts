# This script takes coverage output from multiple instances of RNASeQC (*Coverage*.txt)
# Merges them and plots with confidance interval
# Must specifiy with to normalize height or not
# x axis already assumed to be normalized

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#sns.set(style="darkgrid")

parser = argparse.ArgumentParser()
parser.add_argument("--files","-f", type = str, nargs='+', dest = "infiles", help = "takes coverage output from multiple instances of RNASeQC (*Coverage*.txt)")
parser.add_argument("-n", "--normalize", type = str, dest = "normalize", help = "Normalize y-axis before combining")
args = parser.parse_args()

# Append all dataframes into 1 df
li = []

for filename in args.infiles:
	df = pd.read_csv(filename, index_col=None, header=0)
	df['Position'] = np.arange(len(df))
	df['Norm'] = df[df.columns[0]] / df[df.columns[0]].sum()
	df['Count'] = df[df.columns[0]]
	# Dropping 1st column bcuz I want to change its name
	df = df.drop(df.columns[[0]], axis=1) 
	li.append(df)

frame = pd.concat(li, axis=0, ignore_index=True)

#sns.relplot(x="Position", y="Count", kind="line", ci='sd', data=frame)
sns.relplot(x="Position", y="Count", kind="line", ci=None, data=frame)
plt.show()
