#! /usr/bin/env python

from __future__ import print_function
import sys
import xml.etree.ElementTree as ET
import collections

usage = "Usage: {0} <DemultiplexingStats.xml>".format(sys.argv[0])

if len(sys.argv) != 2:
    print(usage, file=sys.stderr)
    sys.exit(2)

#print("xml file is:", sys.argv[1])

parser = ET.XMLParser(encoding="utf-8")
tree = ET.parse(sys.argv[1], parser=parser)
root = tree.getroot()

def chunkIt(seq, num):
    avg = len(seq)/float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

# Get Samples
sample_list = []
for s in root.iter('Sample'):
    n = s.get('name')
    if not n.startswith("unknown"):
        if not n.startswith("all"):
            sample_list.append(n)

# Get barcode counts
barcode_list = collections.defaultdict(list)
count_list = collections.defaultdict(list)

for s in sample_list:
    counts = []
    for r in root.findall("./Flowcell/Project/*[@name='"+s+"']/Barcode"):
        if not r.get('name').startswith("all"):
            barcode_list[s] += {r.get('name')}
    for c in root.findall("./Flowcell/Project/*[@name='"+s+"']/Barcode/Lane/BarcodeCount"):
        counts.append(c.text)
	countb = sum(map(int, chunkIt(counts, 2)[0]))
	count_list[s] = countb

scounts=sum(count_list.values())

# Get the number of reads
countall = []
for a in root.findall("./Flowcell/Project/*[@name='all']/Barcode/Lane/BarcodeCount"):
    countall.append(a.text)

count_list['UNDETERMINED'] = sum(map(int, chunkIt(countall, 2)[0]))- scounts

# Return table
sample_list.append('UNDETERMINED')
barcode_list['UNDETERMINED'] = 'NONE'
for s in sample_list:
    print(s, "\t", "".join(barcode_list[s]), "\t", count_list[s])

