import matplotlib.pyplot as plt
import numpy as np
from sys import argv

script, rnaseqc_file, output = argv

hist = []

with open(rnaseqc_file, 'r') as f:
	for line in f:
		line = line.strip()
		try:
			hist.append(float(line))
		except:
			pass

y = np.array(hist)
x = np.array(range(0, len(hist)))

plt.plot(x, y)
plt.title(output)
plt.tight_layout()
plt.savefig(output)
