#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 08:02:41 2020
Last edited on Wed Aug 30, 2023 by RLD                                                                                                                                                                                                                                   
@author: vivekmodi
"""
from Bio import PDB
import gzip
import sys
import pandas as pd

omitligands=('1PE', '1PE', '2HT', '2PE', '5LS', '5LS', '5TH', '5TK', '5TM', '7PE',
             'ACE', 'ACT', 'ACY', 'AF3', 'ALA', 'ALY', 'ARS', 'AU',
             'BA',  'BCT', 'BEN', 'BME', 'BOG', 'BR',  'BU1', 'BU3', 'BUD', 'BWB',
             'CA',  'CAC', 'CAF', 'CAS', 'CAS', 'CD',  'CIT', 'CL',  'CME', 'CO',  'CO3',
             'CS',  'CSD', 'CSO', 'CSS', 'CSX', 'CXM', 'CXS', 'CXS',
             'DIO', 'DMS', 'DOD', 'DTD', 'DTT', 'DTV', 'DVT', 'DVT',
             'EDO', 'EMC', 'EOH', 'EPE', 'FLC', 'FMT', 'GBL', 'GG5', 'GLC', 'GLC', 'GOL',
             'HC4', 'HG',  'HSJ', 'IMD', 'IOD', 'IPA', 'IPH', 'K',   'KCX', 'LGY',
             'MES', 'MG',  'MG8', 'MGF', 'MK8', 'MLA', 'MLI', 'MLY', 'MN',  'MOH', 'MPD', 'MRD', 'MSE', 'MXE', 'MYR',
             'NA',  'NEP', 'NH4', 'NI',  'NO3', 'OCS', 'OCT', 'OCY', 'P4G', 'P6G', 'PDX', 'PEG', 'PG0',
             'PG4', 'PGE', 'PGF', 'PGO', 'PHU', 'PO4', 'PPI', 'PSE', 'PTR', 'PUP', 'PZO',
             'S26', 'SBT', 'SCN', 'SCS', 'SEP', 'SGM', 'SIN', 'SO3', 'SO4', 'SR',  'SRT','SVQ',
             'T8L', 'TAM', 'TAR', 'TCE', 'TFA', 'TLA', 'TMA', 'TPO', 'TRS', 'UNX', 'VO4',
             'YT3', 'Z4K', 'ZN',  'P4C')

def compute_distance_from_rre4(structure, model_id, chain_id, ligandname, ligandid, rre4num):
    min_rre4=999

    for model in structure:
        for chain in model:
            model_id=int(model_id)
            if int(model.id)==int(model_id) and chain.id==chain_id:
                for residue1 in chain:
                    if residue1.id[0]==('H_'+ligandname) and residue1.id[1]==int(ligandid):
                        for atom1 in residue1:
                            if atom1.element!='H':

                                for residue2 in chain:
                                    if residue2.get_id()[1]==rre4num:    #Only rre4num
                                        if PDB.is_aa(residue2):  #Only protein atoms or modified residues
                                            for atom2 in residue2:
                                                if atom2.element!='H' and atom2.fullname.strip() not in ('CA','O','N','C'):   #Side chain contact
                                                    distance=residue1[atom1.fullname.strip()]-residue2[atom2.fullname.strip()]
                                                    if min_rre4 > distance:
                                                        min_rre4 = distance
    return min_rre4

def compute_distance_from_hinge(structure,model_id,chain_id,ligandname,ligandid,hinge1):
    min_hinge=999;
    for model in structure:
        for chain in model:
            if int(model.id)==int(model_id) and chain.id==chain_id:
                for residue1 in chain:
                    if residue1.id[0]==('H_'+ligandname) and residue1.id[1]==int(ligandid):
                        for atom1 in residue1:
                            if atom1.element!='H':

                                for residue2 in chain:
                                    if residue2.get_id()[1]>=hinge1 and residue2.get_id()[1]<=hinge1+2:        #compute distance from hinge
                                        if PDB.is_aa(residue2):
                                            for atom2 in residue2:
                                                if atom2.element!='H' and (atom2.fullname.strip()=='O' or atom2.fullname.strip()=='N'):    #Main chain contact
                                                    distance=residue1[atom1.fullname.strip()]-residue2[atom2.fullname.strip()]
                                                    if min_hinge>distance:
                                                        min_hinge=distance
    return min_hinge

def compute_distance_from_pocket_residues(structure,model_id,chain_id,ligandname,ligandid,
                                          type2_resi,back_pocket1_num,back_pocket2_num,back_pocket3_num,front_pocket_num,xdfg_num,spatial):
    dfgcontact=0;dfgoutcontact=0;distance=dict()
    contact_list=list();contact_list_front=list()
    backpocket_count=dict();frontpocket_count=dict()
    backpocket_count[ligandname+':'+ligandid]=0
    frontpocket_count[ligandname+':'+ligandid]=0

    for model in structure:
        for chain in model:
            if int(model.id)==int(model_id) and chain.id==chain_id:
                for residue1 in chain:
                    if residue1.id[0]==('H_'+ligandname) and residue1.id[1]==int(ligandid):
                        for atom1 in residue1:
                            if atom1.element!='H':
                                for residue2 in chain:
                                    if not PDB.is_aa(residue2): continue
                                    res_id=residue2.get_id()[1]

                                    if (res_id in back_pocket1_num or
                                        res_id in back_pocket2_num or
                                        res_id in back_pocket3_num or
                                        int(res_id)==int(xdfg_num) or
                                        int(res_id)==int(xdfg_num+1) or
                                        int(res_id)==int(xdfg_num+2) or
                                        res_id in type2_resi):

                                        for atom2 in residue2:
                                            if atom2.element!='H':
                                                distance[res_id]=round(float((residue1[atom1.fullname.strip()]-residue2[atom2.fullname.strip()])),1)

                                                #Contact with Type2 pocket present
                                                if res_id in type2_resi and distance[res_id]<=4.5:
                                                    dfgoutcontact+=1

                                                #Contact with X-D mainchain
                                                if (distance[res_id]<=4 and
                                                    (int(res_id)==int(xdfg_num) or int(res_id)==int(xdfg_num+1))and
                                                    (atom2.fullname.strip()=='O' or atom2.fullname.strip()=='N') and
                                                    dfgcontact==0 and res_id not in contact_list):
                                                    dfgcontact=1
                                                    backpocket_count[ligandname+':'+ligandid]+=1
                                                    contact_list.append(res_id)

                                                #Contact with Phe sidechain
                                                elif (distance[res_id]<=4 and
                                                      int(res_id)==int(xdfg_num+2) and
                                                      res_id not in contact_list and
                                                      (atom2.fullname.strip()!='O' and atom2.fullname.strip()!='N' and atom2.fullname.strip()!='CA') and
                                                      (spatial!='DFGout' and spatial!='DFGinter')):
                                                    backpocket_count[ligandname+':'+ligandid]+=1
                                                    contact_list.append(res_id)

                                                #Contact with non-XDF backpocket residues
                                                elif (distance[res_id]<=4 and
                                                      int(res_id)!=int(xdfg_num) and
                                                      int(res_id)!=int(xdfg_num+1) and
                                                      int(res_id)!=int(xdfg_num+2) and
                                                      res_id not in contact_list):
                                                    backpocket_count[ligandname+':'+ligandid]+=1
                                                    contact_list.append(res_id)

                                                    #Contact with N-ter of C-helix (Frontpocket)
                                                    if res_id in front_pocket_num and res_id not in contact_list_front:
                                                        frontpocket_count[ligandname+':'+ligandid]+=1
                                                        contact_list_front.append(res_id)


    return frontpocket_count, backpocket_count, dfgoutcontact, distance

#If two chains in the same PDB have Type1 and Type1.5 labels then keep only Type1.5 for both the chains
def correct_chain_diff_in_ligand_type_labels(df):   
    for i in df.index:
        group1=df.at[i,'Group']
        model1=df.at[i,'Model_id']
        chain1=df.at[i,'Chain_id']
        ligand_label1=df.at[i,'Ligand_label']
        ligand_name1=df.at[i,'Ligand']

        for j in df.index:
            group2=df.at[j,'Group']
            model2=df.at[j,'Model_id']
            chain2=df.at[j,'Chain_id']
            ligand_label2=df.at[j,'Ligand_label']
            ligand_name2=df.at[j,'Ligand']

            if group1==group2 and model1==model2 and chain1!=chain2:    #This condition will still be true if the two chains are from different proteins, which will be wrong!
                if ',' in ligand_name1:
                    for position1,ligand_n1 in enumerate(ligand_name1.split(',')):
                        for position2,ligand_n2 in enumerate(ligand_name2.split(',')):
                            if ligand_n1==ligand_n2:
                                if 'Type1'==ligand_label1.split(',')[position1] and 'Type1.5' in ligand_label2.split(',')[position2]:
                                    df.at[i,'Ligand_label']=ligand_label2
                else:
                    if ligand_name1==ligand_name2:
                        if 'Type1'==ligand_label1 and 'Type1.5' in ligand_label2:
                            df.at[i,'Ligand_label']=ligand_label2
    return df


def classify_ligands(pdbfilename,index,df,structure):
        df.at[index,'Ligand_label']='None'    #make default label None
        if 'No_ligand' in df.at[index,'Ligand']:
            return df
        model_id=df.at[index,'Model_id']
        chain_id=str(df.at[index,'Chain_id'])
        spatial=df.at[index,'Spatial_label']
        xdfg_resi=df.at[index,'XDFG_num']
        rre4num=df.at[index,'Glu4_num']
        hinge1=df.at[index,'Hinge1_num']
        type2_resi=df.at[index,'Type2_resi_num']
        back_pocket1_num=df.at[index,'Back_pocket1_num']
        back_pocket2_num=df.at[index,'Back_pocket2_num']
        back_pocket3_num=df.at[index,'Back_pocket3_num']
        front_pocket_num=df.at[index,'Front_pocket_num']
        ligand_label=list()


        ligandlist=df.at[index,'Ligand'].split(',')
        for items in ligandlist:
            ligandname=items.split(':')[0]
            ligandid=items.split(':')[1]

            #Identify allosteric ligands
            min_rre4=compute_distance_from_rre4(structure,model_id,chain_id,ligandname,ligandid,rre4num)
            min_hinge=compute_distance_from_hinge(structure,model_id,chain_id,ligandname,ligandid,hinge1)

            #Contacts with pocket residues
            (frontpocket_count, backpocket_count, dfgoutcontact,distance)=compute_distance_from_pocket_residues(structure,
            model_id,chain_id,ligandname,ligandid,type2_resi,back_pocket1_num,back_pocket2_num,back_pocket3_num,front_pocket_num,xdfg_resi,spatial)


            if min_rre4!=999 or min_hinge!=999:
                if min_rre4>=6.5 and min_hinge>=6.5:
                    ligand_label.append('Allosteric')
                elif min_hinge>=6 and backpocket_count[ligandname+':'+ligandid]>=3:
                    ligand_label.append('Type3')
                elif backpocket_count[ligandname+':'+ligandid]>=3 and frontpocket_count[ligandname+':'+ligandid]==0 and dfgoutcontact==0:
                    ligand_label.append('Type1.5_Back')
                elif backpocket_count[ligandname+':'+ligandid]>=3 and frontpocket_count[ligandname+':'+ligandid]>=1 and dfgoutcontact==0:
                    ligand_label.append('Type1.5_Front')
                elif backpocket_count[ligandname+':'+ligandid]>=3 and dfgoutcontact>=1 and spatial=='DFGout':
                    ligand_label.append('Type2')
                else:
                    ligand_label.append('Type1')

        df.at[index,'Ligand_label']=','.join(ligand_label)


        df=correct_chain_diff_in_ligand_type_labels(df)
        return df
