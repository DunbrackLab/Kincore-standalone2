#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 08:43:59 2020
Last edited on Wed Aug 30, 2023 by RLD
@author: vivekmodi
"""

from Bio import SearchIO
from os.path import exists

def identify_group(pdbfilename,index,conf_df):
    model_id=int(conf_df.at[index,'Model_id'])
    chain_id=conf_df.at[index,'Chain_id']
    if chain_id==" ": chain_id=""
    
    scores=dict()
    
    for group in ('AGC','CAMK','CK1','CMGC','NEK','RGC','STE','TKL','TYR','HASP',
                  'WNK','BUB','ULK','OTHER','PINK1','EIF2AK41','POMK','MAP3K123',
                  'PEAK','SCYL','PXK','STK31','MOS','RPS6KC1','TBCK',
                  'PKDCC','RNASEL','TP53RK','PIK3R4','PAN3','CDC7','ATM'): 
        hmmerfilename=f'{pdbfilename[0:-4]}_{model_id}_{chain_id}_{group}.hmmer.txt'
        if exists(hmmerfilename):
            hmm_result=   SearchIO.read(hmmerfilename, format='hmmer3-text')
            scores[group]=0.0
            for hits in hmm_result:
                for hsp in hits:
                    scores[group]+=hsp.bitscore  # sum up bit scores

    # Find group with highest bit score
    maxScore=0.0
    maxGroup='OTHER'
    for group in scores:
        if group=="MAP3K123" and scores[group]<400: continue  # only use for MAP3K12 and MAP3K13, score ~ 465
        if scores[group]>maxScore:
            maxScore=scores[group]
            maxGroup=group
            conf_df.at[index,'Group']=group
            conf_df.at[index,'Score']=maxScore

    return conf_df
