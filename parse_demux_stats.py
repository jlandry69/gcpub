#! /usr/bin/env python

from __future__ import print_function
import sys
import re
from bs4 import BeautifulSoup

usage = "Usage: {0} <demux-stats.html>".format(sys.argv[0])

if len(sys.argv) != 2:
    print(usage, file=sys.stderr)
    sys.exit(2)

with open(sys.argv[1]) as f_in:
    soup = BeautifulSoup(f_in.read())

tables = soup.findAll("table")
headers = [h.string for h in tables[0].findAll('th')]

rows = tables[1].findAll('tr')
for r in rows:
    cols = [c.string for c in r.findAll('td')]
    d = dict(zip(headers, cols))
    print(d['Lane'], d['Sample ID'], d['Index'], d['# Reads'].replace(',', ''), sep='\t')
