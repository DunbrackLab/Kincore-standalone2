#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 08:29:49 2020
Last edited on Wed Aug 30, 2023 by RLD
@author: vivekmodi
"""
import subprocess


def run_hmmsearch(pwd,pdbfilename,index,conf_df):
    model_id=int(conf_df.at[index,'Model_id'])
    chain_id=conf_df.at[index,'Chain_id']
    if chain_id==" ": chain_id=""  # for MD trajectories
    
    # run all of them for "OTHER" kinases except ones that have their own hmms

    cmd=""
    groupcount=0
    for group in ('AGC_','CAMK_','CK1_','CMGC_','NEK_','RGC_','STE_','TKL_','TYR_'):
        if group in pdbfilename:
            groupname=group.replace("_","")
            groupcount += 1
            cmd += f'hmmsearch  -o {pdbfilename[0:-4]}_{model_id}_{chain_id}_{groupname}.hmmer.txt {pwd}/{groupname}.hmm {pdbfilename[0:-4]}_{chain_id}.fasta;'

    if "MAP3K12" in pdbfilename or "MAP3K13" in pdbfilename:
        groupname="MAP3K123"
        cmd += f'hmmsearch  -o {pdbfilename[0:-4]}_{model_id}_{chain_id}_{groupname}.hmmer.txt {pwd}/{groupname}.hmm {pdbfilename[0:-4]}_{chain_id}.fasta;'
        groupcount += 1
        
    if groupcount==0:
        for groupname in ('AGC','CAMK','CK1','CMGC','NEK','RGC','OTHER','STE','TKL','TYR',
                          'WNK','ULK','BUB','POMK','PINK1','EIF2AK41','HASP','PEAK','SCYL',
                          'EIF2AK41','STK31','PXK','MOS','RPS6KC1','TBCK','MAP3K123',
                          'PKDCC','RNASEL','TP53RK','PIK3R4','PAN3','CDC7','ATM'):
            cmd += (f'hmmsearch  -o {pdbfilename[0:-4]}_{model_id}_{chain_id}_{groupname}.hmmer.txt {pwd}/{groupname}.hmm {pdbfilename[0:-4]}_{chain_id}.fasta; ')

    process=subprocess.Popen(cmd,shell=True)
    process.communicate()
    process.wait()
    return
