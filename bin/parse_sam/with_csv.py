__author__ = 'langmead'

import csv
import numpy

mapqs = []
zt_zs = []

with open('input.sam', 'rb') as fh:
    reader = csv.reader(fh, delimiter='\t', quotechar=None)
    while True:
        try:
            row = reader.next()
        except StopIteration:
            break
        if row[0][0] == '@':
            continue
        zt_zs.append(row[-1][5:])
        mapqs.append(row[4])

print(zt_zs[:10])
print(mapqs[:10])
