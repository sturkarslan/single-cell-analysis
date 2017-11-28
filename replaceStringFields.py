#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 10:18:11 2017

@author: sturkars
"""

variantsfile = "/Volumes/omics4tb/sturkarslan/dvh-mutation-verifications/dvh-UA3-152-03and09_allsamples_mutation_counts_verified.txt"
outfile = "/Volumes/omics4tb/sturkarslan/dvh-mutation-verifications/dvh-UA3-152-03and09_allsamples_mutation_counts_verified_curated2.txt"

f = open(variantsfile, 'r')
h = open(outfile, 'w')

for line in f:
    #print(line)
    fields = line.split("\t")
    fields[-1] = fields[-1].replace('\n','')
    fields = ["0" if "|0|0|OK" in x else x for x in fields]
    fields = ["1" if "|1|1|OK" in x else x for x in fields]
    fields = ["3" if "|0|3|CHECK" in x else x for x in fields]
    fields = ["3" if "|1|3|CHECK" in x else x for x in fields]
    fields = ["3" if "3|OK" in x else x for x in fields]
    fields = ["3" if "3|CHECK" in x else x for x in fields]
    line2write = "\t".join(fields)
    h.write(line2write)
    h.write("\n")
h.close()    
