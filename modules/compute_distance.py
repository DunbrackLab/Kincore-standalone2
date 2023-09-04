#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 09:58:48 2020
Last edited on Wed Aug 30, 2023 by RLD                                                                                                                                                                                                                                    
@author: vivekmodi
"""

import numpy as np
from Bio import PDB

def compute_distance(pdbfilename,index,conf_df,structure):
    conf_df.at[index,'Glu4-Phe-dis']=999.0
    conf_df.at[index,'Lys-Phe-dis']=999.0
    conf_df.at[index,'Lys-Glu-dis']=999.0
    conf_df.at[index,'LysNZ-GluOE-dis']=999.0
    conf_df.at[index,'DFG6-XHRD-dis']=999.0
    restype_atom_dict={'F':'CZ',  'R':'CZ',  'L':'CG',  'P':'CG', 'N':'CG',
                       'M':'CE',  'S':'OG',  'H':'NE2', 'V':'CB', 'A':'CB',
                       'W':'CZ3', 'Y':'OH',  'G':'CA',  'C':'SG', 'D':'CG',
                       'E':'OE1', 'Q':'OE1', 'I':'CD1', 'K':'NZ', 'T':'OG1'}
    try:
        phe_atom_type=restype_atom_dict[conf_df.at[index,'Phe_restype']]
    except:
        return conf_df

    #Distance Glu4-Phe
    if conf_df.at[index,'Glu4_restype'] != "-" and conf_df.at[index,'Phe_restype'] != "-": 
        conf_df.at[index,'Glu4-Phe-dis']=distance_atoms(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],
                                                        conf_df.at[index,'Glu4_num'],conf_df.at[index,'Phe_num'],'CA',phe_atom_type,structure)
    #Distance Spine1
    if conf_df.at[index,'HRD_restype'] != "-" and conf_df.at[index,'Phe_restype'] != "-": 
        conf_df.at[index,'Spine1-dis']=distance_sidechains(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],
                                                        conf_df.at[index,'HRD_num'],conf_df.at[index,'Phe_num'],structure)
    #Distance Spine2
    if conf_df.at[index,'Phe_restype'] != "-" and conf_df.at[index,'Glu4_restype'] != "-": 
        conf_df.at[index,'Spine2-dis']=distance_sidechains(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],
                                                        conf_df.at[index,'Phe_num'],conf_df.at[index,'Glu4_num'],structure)
    #Distance Spine3
    if conf_df.at[index,'Glu4_restype'] != "-" and conf_df.at[index,'HPN7_restype'] != "-": 
        conf_df.at[index,'Spine3-dis']=distance_sidechains(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],
                                                        conf_df.at[index,'Glu4_num'],conf_df.at[index,'HPN7_num'],structure)
    # Spine distance is maximum of spine1, spine2, spine3
    maxspine=-1.0
    if conf_df.at[index,'Spine1-dis']>maxspine: maxspine=conf_df.at[index,'Spine1-dis']
    if conf_df.at[index,'Spine2-dis']>maxspine: maxspine=conf_df.at[index,'Spine2-dis']
    if conf_df.at[index,'Spine3-dis']>maxspine: maxspine=conf_df.at[index,'Spine3-dis']
    conf_df.at[index,'Spine-dis']=maxspine

    #Distance APE9-Arg
    if conf_df.at[index,'APE9_restype'] != "-" and conf_df.at[index,'Arg_restype'] != "-": 
        conf_df.at[index,'APE9-Arg-dis']=distance_atoms(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],
                                                         conf_df.at[index,'APE9_num'],conf_df.at[index,'Arg_num'],'CA','O',structure)

    #Distance Lys-Phe
    if conf_df.at[index,'Lys_restype'] != "-" and conf_df.at[index,'Phe_restype'] != "-": 
        conf_df.at[index,'Lys-Phe-dis']=distance_atoms(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],
                                                       conf_df.at[index,'Lys_num'],conf_df.at[index,'Phe_num'],'CA',phe_atom_type,structure)     

    if conf_df.at[index,'Lys_restype'] != "-" and conf_df.at[index,'Glu_restype'] != "-": 
        #Distance CB-CB of Lys-Glu
        Katom="CB"
        if conf_df.at[index,'Lys_restype']=="G": Katom="CA"
        Eatom="CB"
        if conf_df.at[index,'Glu_restype']=="G": Eatom="CA"
        conf_df.at[index,'Lys-Glu-dis']=distance_atoms(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],
                                                       conf_df.at[index,'Lys_num'],conf_df.at[index,'Glu_num'],Katom,Eatom,structure)      

        #Distance Lys-Glu salt bridge over all hydrogen bonding atoms
        #Have to calculate 4 distances in some kinases (Arg-Glu, Arg-Asp) and sometimes 2 (Lys-Glu) and sometimes 1
        dist1=999.0
        dist2=999.0
        dist3=999.0
        dist4=999.0

        Katom1="CB"
        Katom2="CB"
        Eatom1="CB"
        Eatom2="CB"
    
        if conf_df.at[index,'Lys_restype']=="K":  Katom1="NZ";  Katom2="NZ"
        if conf_df.at[index,'Lys_restype']=="R":  Katom1="NH1"; Katom2="NH2"
        if conf_df.at[index,'Lys_restype']=="Y":  Katom1="OH";  Katom2="OH"
        if conf_df.at[index,'Lys_restype']=="H":  Katom1="ND1"; Katom2="NE2"
        if conf_df.at[index,'Lys_restype']=="Q":  Katom1="OE1"; Katom2="NE2"
        if conf_df.at[index,'Lys_restype']=="N":  Katom1="OD1"; Katom2="ND2"
        if conf_df.at[index,'Lys_restype']=="S":  Katom1="OG";  Katom2="OG"
        if conf_df.at[index,'Lys_restype']=="T":  Katom1="OG1"; Katom2="OG1"
        if conf_df.at[index,'Lys_restype']=="R":  Katom1="NH1"; Katom2="NH2"
        if conf_df.at[index,'Lys_restype']=="C":  Katom1="SG";  Katom2="SG"
        if conf_df.at[index,'Lys_restype']=="G":  Katom1="N";   Katom2="O"
        
        if conf_df.at[index,'Glu_restype']=="E":  Eatom1="OE1"; Eatom2="OE2"
        if conf_df.at[index,'Glu_restype']=="Q":  Eatom1="OE1"; Eatom2="NE2"
        if conf_df.at[index,'Glu_restype']=="D":  Eatom1="OD1"; Eatom2="OD2"
        if conf_df.at[index,'Glu_restype']=="N":  Eatom1="OD1"; Eatom2="ND2"
        if conf_df.at[index,'Glu_restype']=="R":  Eatom1="NH1"; Eatom2="NH2"
        if conf_df.at[index,'Glu_restype']=="K":  Eatom1="NZ";  Eatom2="NZ"
        if conf_df.at[index,'Glu_restype']=="H":  Eatom1="ND1"; Eatom2="NE2"
        if conf_df.at[index,'Glu_restype']=="Y":  Eatom1="OH";  Eatom2="OH"
        if conf_df.at[index,'Glu_restype']=="S":  Eatom1="OG";  Eatom2="OG"
        if conf_df.at[index,'Glu_restype']=="T":  Eatom1="OG1"; Eatom2="OG1"
        if conf_df.at[index,'Glu_restype']=="W":  Eatom1="NE1"; Eatom2="NE1"
        if conf_df.at[index,'Glu_restype']=="G":  Eatom1="N";   Eatom2="O"
    
        dist1=distance_atoms(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],\
                             conf_df.at[index,'Lys_num'],conf_df.at[index,'Glu_num'],Katom1,Eatom1,structure)      
        dist2=distance_atoms(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],\
                             conf_df.at[index,'Lys_num'],conf_df.at[index,'Glu_num'],Katom1,Eatom2,structure)      
        dist3=distance_atoms(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],\
                             conf_df.at[index,'Lys_num'],conf_df.at[index,'Glu_num'],Katom2,Eatom1,structure)      
        dist4=distance_atoms(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],\
                             conf_df.at[index,'Lys_num'],conf_df.at[index,'Glu_num'],Katom2,Eatom2,structure)      
        
        mindist=999.0
        if dist1<mindist: mindist=dist1
        if dist2<mindist: mindist=dist2
        if dist3<mindist: mindist=dist3
        if dist4<mindist: mindist=dist4
        conf_df.at[index,'LysNZ-GluOE-dis']=mindist
    
    #Distance DFG6-XHRD, added RLD 12/24/22
    #Distance is shorter of two possible hbonds
    if conf_df.at[index,'XHRD_restype'] != "-" and conf_df.at[index,'DFG6_restype'] != "-": 
        mindist=999.0
        dist1=999.0
        dist2=999.0
        dist1=distance_atoms(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],\
                             conf_df.at[index,'XHRD_num'],conf_df.at[index,'DFG6_num'],'N','O',structure)      # backbone-backbone hydrogen bond
        dist2=distance_atoms(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],\
                             conf_df.at[index,'XHRD_num'],conf_df.at[index,'DFG6_num'],'O','N',structure)      # backbone-backbone hydrogen bond
        if dist1<mindist: mindist=dist1
        if dist2<mindist: mindist=dist2
        conf_df.at[index,'DFG6-XHRD-dis']=mindist


    return conf_df


def distance_atoms(pdbfilename, model_id, chain_id, res1, res2, atm1, atm2, structure):
    atom_present=0
    res1=int(res1)
    res2=int(res2)

    for model in structure:
        if int(model.id)!=int(model_id): continue
        for chain in model:
            if chain.id!=chain_id: continue

            # Fixed insertion code problem. Only count residues with insertion codes if they offset the sequence numbering
            # For example, if residue numbering is 101 102A 103, then 102A is not an extra residue in sequence and doesn't count
            #         but  if residue numbering is 101 101A 102 103, then 101A is an extra residue in sequence and offsets the sequence
            # RLD 1/3/23

            insertion_num=0    
            lastresiduenum="none"

            for residue in chain:
                if not PDB.is_aa(residue): continue  # skip non polypeptide residues
                if residue.get_id()[2] != ' ' and residue.get_id()[1] == lastresiduenum:  # only include duplicate residue numbers
                    insertion_num+=1
                lastresiduenum=residue.get_id()[1]

                if int(residue.id[1]) == res1-insertion_num:
                    if residue.has_id(atm1):
                        residue1=residue
                        atom_present=atom_present+1

                if int(residue.id[1])==res2-insertion_num:
                    if residue.has_id(atm2):
                        residue2=residue
                        atom_present=atom_present+1


    if atom_present==2:
        distance=np.round((residue1[atm1]-residue2[atm2]),2)
        return distance
    else:
        return 999.0


def distance_sidechains(pdbfilename, model_id, chain_id, res1, res2, structure):
    atom_present=0
    res1=int(res1)
    res2=int(res2)

    residue1=0
    residue2=0
    for model in structure:
        if int(model.id)!=int(model_id): continue
        for chain in model:
            if chain.id!=chain_id: continue

            # Fixed insertion code problem. Only count residues with insertion codes if they offset the sequence numbering
            # For example, if residue numbering is 101 102A 103, then 102A is not an extra residue in sequence and doesn't count
            #         but  if residue numbering is 101 101A 102 103, then 101A is an extra residue in sequence and offsets the sequence
            # RLD 1/3/23

            insertion_num=0    
            lastresiduenum="none"

            for residue in chain:
                if not PDB.is_aa(residue): continue  # skip non polypeptide residues
                if residue.get_id()[2] != ' ' and residue.get_id()[1] == lastresiduenum:  # only include duplicate residue numbers
                    insertion_num+=1
                lastresiduenum=residue.get_id()[1]

                if int(residue.id[1]) == res1-insertion_num:
                    residue1=residue

                if int(residue.id[1]) == res2-insertion_num:
                    residue2=residue


    if residue1 != 0 and residue2 != 0:
        mindist=999.0
        minatom1name="none"
        minatom2name="none"
        for atom1 in residue1:
            atom1name=atom1.get_name()
            if atom1.element=="H": continue
            if residue1.get_resname() == "GLY" and atom1name in ("N","C","O"): continue
            if residue1.get_resname() != "GLY" and atom1name in ("N","CA","C","O"): continue

            for atom2 in residue2:
                atom2name=atom2.get_name()
                if atom2.element=="H": continue
                if residue2.get_resname() == "GLY" and atom2name in ("N","C","O"): continue
                if residue2.get_resname() != "GLY" and atom2name in ("N","CA","C","O"): continue
                distance=np.round((residue1[atom1name]-residue2[atom2name]),2)
                if distance<mindist:
                    mindist=distance
                    minatom1name=atom1name
                    minatom2name=atom2name
        return mindist
    else:
        return 999.0
