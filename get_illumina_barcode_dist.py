#! /usr/bin/env python

from __future__ import print_function
import sys
import os
from subprocess import Popen, PIPE
import csv

usage = "Usage: {0} <demux-stats.html>".format(sys.argv[0])

if len(sys.argv) != 2:
    print(usage, file=sys.stderr)
    sys.exit(2)

parse_demux = os.path.join( os.path.dirname(os.path.realpath(__file__)), 'parse_demux_stats.py')

python = '/g/solexa/bin/software/Python-2.7.3/python'
l = os.path.split(os.path.split(os.path.dirname(sys.argv[1]))[0])[1]
lane = l.strip("Unaligned_lane")
folder = os.path.split(os.path.split(os.path.dirname(sys.argv[1]))[0])[0]
barcode_fn = os.path.join(folder, "log_pipeline","lane" + lane + "_barcodes.txt" )

python = '/g/solexa/bin/software/Python-2.7.3/python'
l = os.path.split(os.path.split(os.path.dirname(sys.argv[1]))[0])[1]
lane = l.strip("Unaligned_lane")
folder = os.path.split(os.path.split(os.path.dirname(sys.argv[1]))[0])[0]
barcode_fn = os.path.join(folder, "log_pipeline","lane" + lane + "_barcodes.txt" )
#print(barcode_fn)
illum_bc = {}
with open(barcode_fn) as f_in:
    reader = csv.reader(f_in, delimiter='\t')
    for row in reader:
        #if row[0].startswith('illumina') :
        illum_bc[row[1]] = row[0]

prog_out = Popen([python, parse_demux, sys.argv[1]], stdout=PIPE)
reader = csv.reader(prog_out.stdout, delimiter='\t')
for row in reader:
    if row[0] == lane:
        if row[2] == 'NoIndex':
            bc_id = "NONE"
            bc_seq = "NONE"
        elif row[2] == 'Undetermined':
            bc_id = "UNDETERMINED"
            bc_seq = "NONE"
        else:
            bc_id = illum_bc[row[2]] if row[2] in illum_bc else "UNKNOWN"
            bc_seq = row[2]
        print(bc_id, bc_seq, row[3], sep='\t')
