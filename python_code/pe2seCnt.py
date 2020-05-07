import os
import sys
from collections import defaultdict

pe_cnt = sys.argv[1]

signatures = dict()
for line in open(pe_cnt):
    if line.startswith("#"):
        continue
    columns = line.split()
    columns = [col.strip() for col in columns]
    left_sig = columns[0].split("^")[0]
    right_sig = columns[0].split("^")[1]
    uniq_cnt = int(columns[2])
    if left_sig in signatures:
        signatures[left_sig] += uniq_cnt
    else:
        signatures[left_sig] = uniq_cnt
    if right_sig in signatures:
        signatures[right_sig] += uniq_cnt
    else:
        signatures[right_sig] = uniq_cnt

for s in signatures:
    print s+"\t"+str(signatures[s])+"\t"+str(signatures[s])


