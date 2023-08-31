#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri July 31 17:35:32 2020
Last edited on Wed Aug 30, 2023 by RLD                                                                                                                                                                                                                                    @author: vivekmodi
"""

import pandas as pd
import sys, os, gzip, argparse
import multiprocessing
from Bio import PDB
from warnings import simplefilter
import multiprocessing

def active_labels(index, conf_df):
    
    # calculate ActLoop positions, ActLoopNT and ActLoopCT
    group=conf_df.at[index, 'Group']
    hbondcutoff=3.6
    ape9cutoff=6.0
    if group == "TYR": ape9cutoff=8.0

    dis_xhrd=conf_df.at[index, 'DFG6-XHRD-dis']
    if dis_xhrd<=hbondcutoff:
        conf_df.at[index, 'ActLoopNT_label']='ActLoopNT-in'
    elif dis_xhrd<900:
        conf_df.at[index, 'ActLoopNT_label']='ActLoopNT-out'
    else:
        conf_df.at[index, 'ActLoopNT_label']="None"
                
    dis_ape9=conf_df.at[index, 'APE9-Arg-dis']
    if dis_ape9 <= ape9cutoff:
        conf_df.at[index, 'ActLoopCT_label']='ActLoopCT-in'
    elif dis_ape9 < 900:
        conf_df.at[index, 'ActLoopCT_label']='ActLoopCT-out'
    else:
        conf_df.at[index, 'ActLoopCT_label']="None"

    # BUB, HASP, PKDCC, TP53RK exception, no APE motif;
    # PEAK3 has group=PEAK and score~170 and no APE motif (PEAK1 and PRAG1, >500)
    if (group=="BUB" or group=="HASP" or group=="PKDCC" or
        group=="TP53RK" or group=="PAN3" or group=="RNASEL" or
        (group=="PEAK" and conf_df.at[index, 'Score']<300)):
        conf_df.at[index, 'ActLoopCT_label']="None"
        
    # calculate salt bridge position
    dis_saltbridge=conf_df.at[index, 'LysNZ-GluOE-dis']
    KEtype=conf_df.at[index, 'Lys_restype'] + conf_df.at[index, 'Glu_restype']

    if KEtype == 'KE' or KEtype == 'RE' : KEcutoff=hbondcutoff
    if KEtype == 'KD' or KEtype == 'RD' : KEcutoff=hbondcutoff+1.5
    if KEtype == 'KN' or KEtype == 'RN':  KEcutoff=10.0  # no saltbridge if Asn

    if KEtype == 'KE' or KEtype == 'RE' or KEtype == 'KD' or KEtype == 'RD' or KEtype == 'KN' or KEtype == 'RN':
        if dis_saltbridge<=KEcutoff:
            conf_df.at[index, 'SaltBr_label']='SaltBr-in'
        elif dis_saltbridge<900:
            conf_df.at[index, 'SaltBr_label']='SaltBr-out'
        else:
            conf_df.at[index, 'SaltBr_label']="None"
    else:
        conf_df.at[index, 'SaltBr_label']="None"

    # for MOS and PINK1 where Glu4=ALA
    if conf_df.at[index,'Glu4_restype']=="A":
        conf_df.at[index,'Spine-dis']=conf_df.at[index,'Spine1-dis']
                
    if conf_df.at[index,'Spine-dis']<0:
        conf_df.at[index,'Spine_label']='None'
    elif conf_df.at[index,'Spine-dis']<=5.0:
        conf_df.at[index,'Spine_label']="Spine-in"
    elif conf_df.at[index,'Spine-dis']<900:
        conf_df.at[index,'Spine_label']="Spine-out"
    else:
        conf_df.at[index,'Spine_label']="None"
        
    # Determine activity state: active kinase must be DFGin-BLAminus-Actloop-in with Saltbridge<= 3.6 Angstroms
    ActivityState="Inactive"
    if (conf_df.at[index, 'Spatial_label']=="DFGin" and 
        conf_df.at[index, 'Dihedral_label']=="BLAminus" and
        conf_df.at[index, 'SaltBr_label'] == "SaltBr-in" and
        conf_df.at[index, 'ActLoopNT_label']=="ActLoopNT-in" and
        conf_df.at[index, 'ActLoopCT_label']=="ActLoopCT-in"): ActivityState="Active"

    if ( (conf_df.at[index, 'Spatial_label']=="DFGin" and
          conf_df.at[index, 'Dihedral_label']=="BLAminus") and
         (conf_df.at[index,'SaltBr_label']=="None" or
          conf_df.at[index,'ActLoopNT_label']=="None" or
          conf_df.at[index,'ActLoopCT_label']=="None")):       ActivityState="None"

                 
    # skip saltbridge and Chelix for WNK kinases, MAP3K12 and MAP3K13 (Glu position is Asp)
    if group == "WNK" or group=="MAP3K123":   
        ActivityState="Inactive"
        if (conf_df.at[index, 'Spatial_label']=="DFGin" and 
            conf_df.at[index, 'Dihedral_label']=="BLAminus" and
            conf_df.at[index, 'ActLoopNT_label']=="ActLoopNT-in" and
            conf_df.at[index, 'ActLoopCT_label']=="ActLoopCT-in"): ActivityState="Active"

        if ( (conf_df.at[index, 'Spatial_label']=="DFGin" and
              conf_df.at[index, 'Dihedral_label']=="BLAminus") and
             (conf_df.at[index,'ActLoopNT_label']=="None" or
              conf_df.at[index,'ActLoopCT_label']=="None" )): ActivityState="None"
                 

    # skip ActLoopCT, no APE motif for BUB, HASP, PKDCC, TP53RK, PAN3, RNASEL, PEAK kinases
    if (group == "BUB" or group=="HASP" or group=="PKDCC" or
        group=="TP53RK" or group=="PAN3" or group=="RNASEL" or   
        (group=="PEAK" and conf_df.at[index, 'Score']<300)):
        ActivityState="Inactive"
        if (conf_df.at[index, 'Spatial_label']=="DFGin" and 
            conf_df.at[index, 'Dihedral_label']=="BLAminus" and
            conf_df.at[index, 'SaltBr_label'] == "SaltBr-in" and
            conf_df.at[index, 'ActLoopNT_label']=="ActLoopNT-in"): ActivityState="Active"
                 
        if ( (conf_df.at[index, 'Spatial_label']=="DFGin" and
              conf_df.at[index, 'Dihedral_label']=="BLAminus") and
             ( conf_df.at[index, 'SaltBr_label'] == "None" or
               conf_df.at[index,'ActLoopNT_label']=="None")):   ActivityState="None"

    conf_df.at[index,'Activity_label']=ActivityState
    return conf_df
