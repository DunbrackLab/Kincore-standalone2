#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 11:57:13 2020
Last edited on Wed Aug 30, 2023 by RLD
@author: vivekmodi
"""

import os, subprocess

# delete all HMM files except identified group
def delete_files(pdbfilename, model_id, chain_id, identified_group):
    model_id=str(int(model_id))
    filename=f'{pdbfilename[0:-4]}_{chain_id}.fasta'
    if os.path.isfile(filename): os.remove(filename)

    for group in ('AGC','CAMK','CK1','CMGC','NEK','RGC','STE','TKL','TYR',
                  'HASP','WNK','BUB','ULK','OTHER','PINK1','EIF2AK41',
                  'POMK','PEAK','SCYL','PXK','STK31','MAP3K123','MOS','RPS6KC1','TBCK',
                  'PKDCC','RNASEL','TP53RK','PIK3R4','PAN3', 'CDC7', 'ATM'): 
        
        filename=f'{pdbfilename[0:-4]}_{model_id}_{chain_id}_{group}.hmmer.txt'
        if os.path.isfile(filename):   os.remove(filename)
