#!/bin/bash

#$ -N scite-03and091cells
#$ -o /proj/omics4tb/sturkarslan/dvh-mutation-verifications/scitelog-dvh-0309_1cell.txt
#$ -e /proj/omics4tb/sturkarslan/dvh-mutation-verifications/scitelog-dvh-0309_1cell.txt
#$ -P Bal_sturkars
#$ -pe serial 16
#$ -q baliga
#$ -S /bin/bash
#$ -m abe
#$ -M sturkarslan@systemsbiology.org

cd /users/sturkars/github/SCITE/

./scite -i /proj/omics4tb/sturkarslan/dvh-mutation-verifications/dvh-UA3-152-03and09_SC_curated_1cells_mutation_matrix.txt -names /proj/omics4tb/sturkarslan/dvh-mutation-verifications/dvh-UA3-152-03and09_SC_curated_1cells_mutation_names.txt -n 131 -m 188 -r 1 -l 900000 -fd 6.04e-5 -ad 0.21545 0.21545 -a -e 0.7
