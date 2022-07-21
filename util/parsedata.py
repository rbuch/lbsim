#! /usr/bin/env python3

import matplotlib.pyplot as plt

import sys

lbMap = {}

currentLB = None

for filename in sys.argv[1:]:
    with open(filename, 'r') as f:
        data = f.readlines()
        for line in data:
            if 'time' in line:
                currentLB = line.split(':')[0].split()[-1]
            elif 'Ratio' in line:
                ratio = float(line.strip().split('=')[1][:-1])
                assert(currentLB is not None)
                if currentLB not in lbMap:
                    lbMap[currentLB] = [ratio]
                else:
                    lbMap[currentLB].append(ratio)

print(lbMap)

plt.style.use('seaborn-colorblind')
fig, ax = plt.subplots()
ax.set_title('LB Performance')
ax.set_ylabel('Max/Avg Ratio')

labels = [key for key in lbMap][1:]
data = [lbMap[key] for key in lbMap][1:]

ax.boxplot(data, labels=labels)
plt.show()
