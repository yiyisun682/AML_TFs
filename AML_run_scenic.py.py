#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## python version 3.8.0
## Construct the transcriptional regulatory networks of 
## Acute myeloid leukemia using pySCENIC

import sys,os,re
####################################
#command: python .py
####################################
scenic='~/miniconda3/bin/pyscenic'
TF_file='~/pyscenic/hg_TFs_2.txt'
annotation='~/pyscenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
feather='~/pyscenic/hg19-500bp-upstream-10species.mc9nr.feather'


os.system(scenic+' grn -o TFtarget.tab --num_workers 5 A_matrix.csv '+TF_file)
os.system('cat TFtarget.tab | tr "\t" "," > TFtarget.csv')
os.system(scenic+' ctx -o regulons.csv --num_workers 5 --annotations_fname '+annotation+' --expression_mtx_fname A_matrix.csv TFtarget.csv '+feather)
os.system(scenic+' aucell -o aucell.csv --num_workers 5 A_matrix.csv regulons.csv')
print('done!')

