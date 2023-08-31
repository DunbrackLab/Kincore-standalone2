#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 12:27:13 2020
Last edited on Wed Aug 30, 2023 by RLD
@author: vivekmodi
"""

def assign_default_values(index,conf_df):

    conf_df.at[index,'Activity_label']='None'
    conf_df.at[index,'ActLoopNT_label']='None'
    conf_df.at[index,'ActLoopCT_label']='None'
    conf_df.at[index,'Spine_label']='None'
    conf_df.at[index,'SaltBr_label']='None'
    conf_df.at[index,'Chelix_label']='None'
    conf_df.at[index,'Spatial_label']='None'
    conf_df.at[index,'Dihedral_label']='None'
    conf_df.at[index,'Ligand']='None'
    conf_df.at[index,'Ligand_label']='None'

    conf_df.at[index,'Sequence']='-'
    conf_df.at[index,'First_res']=0
    conf_df.at[index,'Group']='None'

    conf_df.at[index,'DFG6-XHRD-dis']=999
    conf_df.at[index,'Spine1-dis']=999
    conf_df.at[index,'Spine2-dis']=999
    conf_df.at[index,'Spine3-dis']=999
    conf_df.at[index,'Spine-dis']=999
    conf_df.at[index,'Glu4-Phe-dis']=999
    conf_df.at[index,'APE9-Arg-dis']=999
    conf_df.at[index,'Lys-Phe-dis']=999
    conf_df.at[index,'Lys-Glu-dis']=999
    conf_df.at[index,'LysNZ-GluOE-dis']=999

    conf_df.at[index,'Lys_num']=0
    conf_df.at[index,'real_Lys_num']=0
    conf_df.at[index,'Lys_restype']='-'

    conf_df.at[index,'Glu_num']=0
    conf_df.at[index,'real_Glu_num']=0
    conf_df.at[index,'Glu_restype']='-'

    conf_df.at[index,'Glu4_num']=0
    conf_df.at[index,'real_Glu4_num']=0
    conf_df.at[index,'Glu4_restype']='-'

    conf_df.at[index,'APE_num']=0
    conf_df.at[index,'real_APE_num']=0
    conf_df.at[index,'APE_restype']='-'

    conf_df.at[index,'APE9_num']=0
    conf_df.at[index,'real_APE9_num']=0
    conf_df.at[index,'APE9_restype']='-'

    conf_df.at[index,'DFG6_num']=0
    conf_df.at[index,'real_DFG6_num']=0
    conf_df.at[index,'DFG6_restype']='-'

    conf_df.at[index,'HPN7_num']=0   # for spine calculation
    conf_df.at[index,'real_HPN7_num']=0
    conf_df.at[index,'HPN7_restype']='-'

    conf_df.at[index,'HRD_num']=0
    conf_df.at[index,'real_HRD_num']=0
    conf_df.at[index,'HRD_restype']='-'

    conf_df.at[index,'Arg_num']=0
    conf_df.at[index,'real_Arg_num']=0
    conf_df.at[index,'Arg_restype']='-'

    conf_df.at[index,'XHRD_num']=0
    conf_df.at[index,'real_XHRD_num']=0
    conf_df.at[index,'XHRD_restype']='-'

    conf_df.at[index,'XDFG_num']=0
    conf_df.at[index,'real_XDFG_num']=0
    conf_df.at[index,'XDFG_restype']='-'

    conf_df.at[index,'Asp_num']=0
    conf_df.at[index,'real_Asp_num']=0
    conf_df.at[index,'Asp_restype']='-'

    conf_df.at[index,'Phe_num']=0
    conf_df.at[index,'real_Phe_num']=0
    conf_df.at[index,'Phe_restype']='-'

    conf_df.at[index,'Hinge1_num']=0

    conf_df.at[index,'XDFG_Phi']=999
    conf_df.at[index,'XDFG_Psi']=999
    conf_df.at[index,'Asp_Phi']=999
    conf_df.at[index,'Asp_Psi']=999
    conf_df.at[index,'Asp_Chi1']=999
    conf_df.at[index,'Asp_Chi2']=999
    conf_df.at[index,'Phe_Phi']=999
    conf_df.at[index,'Phe_Psi']=999
    conf_df.at[index,'Phe_Chi1']=999
    conf_df.at[index,'Phe_Chi2']=999
    conf_df.at[index,'Gly_Phi']=999
    conf_df.at[index,'Gly_Psi']=999

    return conf_df
