# Plots the dropout by High/Med/Low expression relative to a reference
# Usage: python countByExpLevel.py A.txt Ref.txt
# Both files need to have gene names as the 1st column and a column named 'TPM'

from sys import argv
from collections import defaultdict

script, inputfile, reference = argv


