from sys import argv
import pandas as pd
from datetime import date
script= argv[0]
BL3f = argv[1]
BL5f = argv[2]
Y53f = argv[3]
Y55f = argv[4]
Patient = argv[5]

def intersections(lst1, lst2):
    return list(set(lst1) & set(lst2))

BL3 = pd.read_csv(BL3f, sep = '\t', index_col=0)

Y53 = pd.read_csv(Y53f, sep = '\t', index_col=0)

BL5 = pd.read_csv(BL5f, sep = '\t', index_col=0)

Y55 = pd.read_csv(Y55f, sep = '\t', index_col=0)
left = Y53f.split('_')[1]
right = BL3f.split('_')[1]

mayr = pd.read_csv('/fs/scratch/PAS0472/osu9900/mrd_lcrna_0825020/output/B_CELLS/comparisons_spliced_without_UTRs_lufei/41586_2018_465_MOESM3_ESM.csv', index_col=3)

print("Table for " + Patient +"_spliced_withoutUTRs")


today = date.today()
print('Date generated: ' + str(today))
print("\n")
print(right+' 3\' trunction: ' + str(len(set(BL3.index.tolist()))))
print(left+' 3\' trunction: ' + str(len(set(Y53.index.tolist()))))
print(right+' 5\' trunction: ' + str(len(set(BL5.index.tolist()))))
print(left+' 5\' trunction: ' + str(len(set(Y55.index.tolist()))))

print("\n")
print(right+' 3\' truncation overlaps with mayr list')
print('length:' + str(len(set(BL3.index.tolist()).intersection(set(mayr.index.tolist())))))
print(','.join(set(BL3.index.tolist()).intersection(set(mayr.index.tolist()))))

print("\n")
print(left+' 3\' truncation overlaps with mayr list')
print('length:' + str(len(set(Y53.index.tolist()).intersection(set(mayr.index.tolist())))))
print(','.join(set(Y53.index.tolist()).intersection(set(mayr.index.tolist()))))

print("\n")
print(right+' 5\' truncation overlaps with mayr list')
print('length:' + str(len(set(BL5.index.tolist()).intersection(set(mayr.index.tolist())))))
print(','.join(set(BL5.index.tolist()).intersection(set(mayr.index.tolist()))))

print("\n")
print(left+' 5\' truncation overlaps with mayr list')
print('length:' + str(len(set(Y55.index.tolist()).intersection(set(mayr.index.tolist())))))

print(','.join(set(Y55.index.tolist()).intersection(set(mayr.index.tolist()))))