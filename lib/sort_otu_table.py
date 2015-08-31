#!/usr/bin/env python

import sys
import csv
import re
import os
import select

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]

def main(x):
    r = csv.DictReader(x, delimiter='\t')
    new_order = sorted(r.fieldnames, key=natural_sort_key)
    otu = new_order.index('OTUId')
    new_order.insert(0, new_order.pop(otu))
    w = csv.DictWriter(sys.stdout, new_order, delimiter='\t')
    w.writeheader()
    for a in r:
        w.writerow(a)
    
if not select.select([sys.stdin,],[],[],0.0)[0]:
    if len(sys.argv) < 2:
        print "no input specified"
        os._exit(1)
    else:
        input = open(sys.argv[1], 'rb')
        main(input)
        input.close()
else:
    input = sys.stdin.readlines()
    main(input)
