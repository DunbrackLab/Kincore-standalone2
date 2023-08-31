#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 14:13:52 2020
Last edited on Wed Aug 30, 2023 by RLD
@author: vivekmodi
"""

import gzip
from Bio import PDB

def extract_ligands(pdbfilename,index,df,structure):
    omitligands=('1PE','1PE','2HT','2PE','5LS','5LS','5TH','5TK','5TM','7PE','ACE','ACT','ACY','AF3','ALA','ALY','ARS','AU','BA','BCT','BEN','BME','BOG','BR',\
                 'BU1','BU3','BUD','BWB','CA','CAC','CAF','CAS','CAS','CD','CIT','CL','CME','CME','CO','CO3','CS','CSD','CSO','CSO','CSS','CSS','CSX','CXM',\
                 'CXS','CXS','DIO','DMS','DOD','DTD','DTT','DTT','DTV','DVT','DVT','EDO','EMC','EOH','EPE','EPE','FLC','FMT','GBL','GG5','GLC','GLC','GOL','HC4',\
                 'HG','HSJ','IMD','IOD','IPA','IPH','K','KCX','LGY','MES','MG8','MGF','MK8','MLA','MLI','MLY','MOH','MPD','MRD','MSE','MXE','MYR','NA',\
                 'NEP','NH4','NI','NO3','OCS','OCT','OCY','P4G','P6G','PDX','PEG','PG0','PG4','PGE','PGF','PGO','PHU','PO4','PPI','PSE','PTR','PUP','PZO','S26',\
                 'SBT','SCN','SCS','SEP','SEP','SGM','SIN','SO3','SO4','SR','SRT','SVQ','T8L','TAM','TAR','TCE','TFA','TLA','TLA','TMA','TPO','TRS','UNX','VO4',\
                 'YT3','Z4K','ZN','P4C')
    modified_aa=('ACE','ALY','AME','BWB','CAF','CAS','CME','COM','CSD','CSO','CSS','CSX','CXM','CY0','CYO','KCX','LGY','MHO','MK8','MLY','MSE','NEP','NMM',\
                 'OCS','OCY','PHD','PTR','SCS','SEP','T8L','TPO','UNK')
    
    model_id=df.at[index,'Model_id']
    chain_id=df.at[index,'Chain_id']

    ligandname=''
    ligandpresent=0
    ligandlist=list();

    for model in structure:
        for chain in model:
            if int(model.id)==int(model_id) and chain.id==chain_id:
                for residue in chain:
                        if 'H_' in residue.id[0]:
                            ligand_id=residue.id[1]
                            ligandname=str(residue.id[0][2:].strip(" "))

                            if ligandname in omitligands:
                                continue
                            if ligandname in modified_aa:
                                continue


                            else:
                                ligandcount=0
                                for atom in residue:       #iterate over ligand atoms
                                    if ligandcount==0:     #if the same ligand is not counted before
                                        for residue2 in chain:   #iterate over protein residues
                                            if PDB.is_aa(residue2):
                                                for atom2 in residue2:
                                                    # if ligand atom makes contact with residue within kinase domain;
                                                    # Some PDB files have spaces around atom name so use strip()
                                                    distance=residue[atom.fullname.strip()]-residue2[atom2.fullname.strip()]
                                                    #to make sure atoms for the same ligand are not counted twice
                                                    if distance<=4 and ligandcount==0:                             
                                                        ligandpresent=1
                                                        ligandcount=1
                                                        ligandlist.append(f'{ligandname}:{ligand_id}')     #format ['ATP:300']
                                                        df.at[index,'Ligand']=','.join(ligandlist)




    if ligandpresent==0:
        ligandname='No_ligand'
        ligandlist.append(ligandname)
        df.at[index,'Ligand']=','.join(ligandlist)

    return df
