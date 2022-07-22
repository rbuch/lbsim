#! /usr/bin/env python3

import matplotlib.pyplot as plt
import os
from pathlib import Path
import sys

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--inputdir", help="Root processor directory)", required=True)
parser.add_argument("-o", "--output", help="Graph output filename")

args = parser.parse_args()


objToResults = {}

currentLB = None

rootdir = args.inputdir

if not os.path.isdir(rootdir):
    print("Error: given argument {} is not a directory".format(rootdir))
    sys.exit(1)

p = Path(rootdir)
numPes = int(p.stem)

objtiers = sorted([x for x in p.iterdir()])


for objtier in objtiers:
    numObjs = int(objtier.stem)
    lbMap = {}
    for filename in objtier.iterdir():
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
    objToResults[numObjs] = lbMap

print(objToResults)

# trigger core fonts for PDF backend
plt.rcParams["pdf.use14corefonts"] = True

plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica"
})

plt.style.use('seaborn-colorblind')
#fig, axs = plt.subplots(len(objToResults), sharey=True)
fig, axs = plt.subplots(len(objToResults))


for numObjs, i in zip(objToResults, range(len(objToResults))):
    lbMap = objToResults[numObjs]
    ax = axs[i]
    ax.set_title('LB Performance ({} PEs, {} Objs)'.format(numPes, numObjs))
    ax.set_ylabel('Max/Avg Ratio')
    labels = [key for key in lbMap][1:]
    labels = list(map(lambda x: '$k$d '+x[3:]+"-norm" if 'rkd' in x else x, labels))
    labels = list(map(lambda x: x.replace('Inf', '$\\infty$'), labels))
    data = [lbMap[key] for key in lbMap][1:]

    ax.boxplot(data, labels=labels, whis=(0,100))

plt.show()

#fig.tight_layout()
fig.savefig(args.output, format='pdf', bbox_inches='tight')
