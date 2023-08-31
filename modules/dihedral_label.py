#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Last edited on Wed Aug 30, 2023 by RLD

import math, numpy as np

def dihedral_label(index,conf_df,cutoff):
    dfginter={'BABtrans':(-80.20,128.22,-117.47,23.76,-85.16,133.21,181.42)}
    dfgout={'BBAminus':(-138.56,-176.12,-144.35,103.66,-82.59,-9.03,290.59)}
    
    dfgin_minus={'BLAminus':(-128.64,178.67,61.15,81.21,-96.89,20.53,289.12),
                 'ABAminus':(-111.82,-7.64,-141.55,148.01,-127.79,23.32,296.17),
                 'BLBminus':(-134.79,175.48,60.44,65.35,-79.44,145.34,287.56)}
    
    dfgin_plus={'BLAplus':(-119.24,167.71,58.94,34.08,-89.42,-8.54,55.63),
                'BLBplus':(-125.28,172.53,59.98,32.92,-85.51,145.28,49.01)}
    
    dfgin_trans={'BLBtrans':(-106.16,157.24,69.37,21.33,-61.73,134.56,215.23)}
    

    xdfg_phi=float(conf_df.at[index,'XDFG_Phi'])
    xdfg_psi=float(conf_df.at[index,'XDFG_Psi'])
    dfg_asp_phi=float(conf_df.at[index,'Asp_Phi'])
    dfg_asp_psi=float(conf_df.at[index,'Asp_Psi'])
    dfg_phe_phi=float(conf_df.at[index,'Phe_Phi'])
    dfg_phe_psi=float(conf_df.at[index,'Phe_Psi'])
    dfg_phe_chi1=float(conf_df.at[index,'Phe_Chi1'])
    
    conf_df.at[index,'Dihedral_label']='None'    #Default label is None
    conf_df.at[index,'Dihedral_dis']=999.0
    conf_df.at[index,'Dihedral_dis_NoChi1']=999.0
    
    if conf_df.at[index,'Spatial_label']=='None':
        conf_df.at[index,'Dihedral_label']='None'
        return conf_df
    
    if xdfg_phi==999 or xdfg_psi==999 or dfg_asp_phi==999 or dfg_asp_psi==999 or dfg_phe_phi==999 or  dfg_phe_psi==999 or dfg_phe_chi1==999:
        conf_df.at[index,'Dihedral_label']='None'
        return conf_df
    
    if conf_df.at[index,'Spatial_label']=='DFGin':            
        if dfg_phe_chi1>240 and dfg_phe_chi1<=360:
            conf_df=cosine_dis_without_chi1(conf_df,index,dfgin_minus,cutoff)
        if dfg_phe_chi1>0 and dfg_phe_chi1<=120:
            conf_df=cosine_dis_without_chi1(conf_df,index,dfgin_plus,cutoff)
        if dfg_phe_chi1>120 and dfg_phe_chi1<=240:
            conf_df=cosine_dis_without_chi1(conf_df,index,dfgin_trans,cutoff)
            
    if conf_df.at[index,'Spatial_label']=='DFGinter':            
        if dfg_phe_chi1>120 and dfg_phe_chi1<=240:
            conf_df=cosine_dis_without_chi1(conf_df,index,dfginter,cutoff)
        
    if conf_df.at[index,'Spatial_label']=='DFGout':     
        if dfg_phe_chi1>240 and dfg_phe_chi1<=360:
            conf_df=cosine_dis_without_chi1(conf_df,index,dfgout,cutoff)

    return conf_df
    
def cosine_dis_without_chi1(df,i,spatial,cutoff):
    x_dfg_phi=float(df.at[i,'XDFG_Phi'])
    x_dfg_psi=float(df.at[i,'XDFG_Psi'])
    dfg_asp_phi=float(df.at[i,'Asp_Phi']);
    dfg_asp_psi=float(df.at[i,'Asp_Psi'])
    dfg_phe_phi=float(df.at[i,'Phe_Phi'])
    dfg_phe_psi=float(df.at[i,'Phe_Psi'])
    min_spatial=999
    
    for clusters in spatial:    
        cosine_dis=(2/6)*((1-math.cos(math.radians(x_dfg_phi-float(spatial[clusters][0]))))+(1-math.cos(math.radians(x_dfg_psi-float(spatial[clusters][1]))))+
                          (1-math.cos(math.radians(dfg_asp_phi-float(spatial[clusters][2]))))+(1-math.cos(math.radians(dfg_asp_psi-float(spatial[clusters][3]))))+
                          (1-math.cos(math.radians(dfg_phe_phi-float(spatial[clusters][4]))))+(1-math.cos(math.radians(dfg_phe_psi-float(spatial[clusters][5])))))
        
        if cosine_dis<=min_spatial:
            df.at[i,'Dihedral_dis_NoChi1']=np.round(cosine_dis,2)
            min_spatial=cosine_dis
        
            if cosine_dis<=cutoff:
                df.at[i,'Dihedral_label']=clusters      #Only Dihedral column name is used for final labeling without chi1
        
    return df
