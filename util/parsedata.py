#! /usr/bin/env python3

import matplotlib.pyplot as plt
import os
from pathlib import Path
import sys

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--inputdir", help="Root processor directory)", required=True)
parser.add_argument("-o", "--output", help="Graph output filename stem")

args = parser.parse_args()


MAX = 'Maximum'
TOTAL = 'Total'

objToResults = {MAX: {}, TOTAL: {}}

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
    lbMap = {MAX: {}, TOTAL: {}}
    for filename in objtier.iterdir():
        with open(filename, 'r') as f:
            data = f.readlines()
            for line in data:
                if 'time' in line:
                    currentLB = line.split(':')[0].split()[-1]
                elif 'Ratio' in line:
                    ratioString = line.strip().split('(')[1][:-1]
                    totalRatio, maxRatio = [float(x.split('=')[1]) for x in
                                       ratioString.split(',')]
                    assert(currentLB is not None)
                    if currentLB not in lbMap[MAX]:
                        lbMap[MAX][currentLB] = [maxRatio]
                        lbMap[TOTAL][currentLB] = [totalRatio]
                    else:
                        lbMap[MAX][currentLB].append(maxRatio)
                        lbMap[TOTAL][currentLB].append(totalRatio)
    objToResults[MAX][numObjs] = lbMap[MAX]
    objToResults[TOTAL][numObjs] = lbMap[TOTAL]

#print(objToResults)

# trigger core fonts for PDF backend
plt.rcParams["pdf.use14corefonts"] = True

plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica"
})

plt.style.use('seaborn-colorblind')
#fig, axs = plt.subplots(len(objToResults), sharey=True)

for kind in objToResults:
    dataset = objToResults[kind]

    fig, axs = plt.subplots(len(dataset))

    for numObjs, i in zip(dataset, range(len(dataset))):
        lbMap = dataset[numObjs]
        ax = axs[i]
        ax.set_title('LB Performance ({}, {} PEs, {} Objs)'.format(kind, numPes, numObjs))
        ax.set_ylabel('Max/Avg Ratio')
        labels = [key for key in lbMap][1:]
        labels = list(map(lambda x: '$k$d '+x[3:]+"-norm" if 'rkd' in x else x, labels))
        labels = list(map(lambda x: x.replace('Inf', '$\\infty$'), labels))
        data = [lbMap[key] for key in lbMap][1:]

        ax.boxplot(data, labels=labels, whis=(0,100))

    #plt.show()

    #fig.tight_layout()
    fig.set_size_inches(9, len(dataset) * 3.6)
    fig.savefig("{}-{}.pdf".format(args.output, kind), format='pdf', bbox_inches='tight')
