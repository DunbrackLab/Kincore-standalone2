#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 09:14:31 2020
Last edited on Wed Aug 30, 2023 by RLD
@author: vivekmodi
"""

from Bio import SearchIO
from Bio import PDB


# get author residue number for each important residue
def get_real_residue_number(pdbfilename, model_id, chain_id, res1, structure):
    res1=int(res1)

    for model in structure:
        if int(model.id)!=int(model_id): continue
        for chain in model:
            if chain.id!=chain_id: continue
            
            insertion_num=0    
            lastresiduenum="none"

            for residue in chain:
                if not PDB.is_aa(residue): continue  # skip non polypeptide residues
                if residue.get_id()[2] != ' ' and residue.get_id()[1] == lastresiduenum:  # only include duplicate residue numbers
                    insertion_num+=1
                lastresiduenum=residue.get_id()[1]

                if int(residue.id[1]) == res1-insertion_num:
                    return int(residue.id[1])
    return 0   # if residue not found       


# Identify specific residues for structure classification
def identify_residues(pdbfilename,index,conf_df, s2cdict):
    model_id=int(conf_df.at[index,'Model_id'])
    chain_id=conf_df.at[index,'Chain_id']
    group=conf_df.at[index,'Group']

    # first lines from CMGC.hmm with low gap penalties
    # from alignment of consensus sequences (manually edited)
    # beta strand B3 lys (K)
    #  lys=   {'AGC':30, 'BUB':30, 'CAMK':30, 'CK1':30, 'CMGC':30, 'EIF2AK41':35, 'HASP':28, 'MAP3K123':35, 'MOS':36,
    lys=   {'AGC':30, 'BUB':30, 'CAMK':30, 'CK1':30, 'CMGC':30, 'EIF2AK41':35, 'HASP':28, 'MAP3K123':35, 'MOS':36,
            'NEK':30, 'OTHER':30, 'PAN3':50, 'PEAK':41, 'PIK3R4':35, 'PINK1':72, 'PKDCC':36, 'POMK':35,
            'PXK':42, 'RGC':29, 'RNASEL':40, 'RPS6KC1':33, 'SCYL':39, 'STE':30, 'STK31':57, 'TBCK':12, 'TKL':28,
            'TP53RK':29, 'TYR':34, 'ULK':31, 'WNK':31, 'CDC7':41, 'ATM':182}

    # Chelix glu (E of RRE)
    #  glu=   {'AGC':49, 'BUB':39, 'CAMK':48, 'CK1':44, 'CMGC':47, 'EIF2AK41':70, 'HASP':52, 'MAP3K123':44, 'MOS':55, 
    glu=   {'AGC':49, 'BUB':39, 'CAMK':48, 'CK1':44, 'CMGC':45, 'EIF2AK41':70, 'HASP':52, 'MAP3K123':44, 'MOS':55, 
            'NEK':48, 'OTHER':47, 'PAN3':65, 'PEAK':57, 'PIK3R4':55, 'PINK1':93, 'PKDCC':69, 'POMK':49, 
            'PXK':60, 'RGC':45, 'RNASEL':52, 'RPS6KC1':47, 'SCYL':52, 'STE':47, 'STK31':72, 'TBCK':41, 'TKL':46, 
            'TP53RK':53, 'TYR':51, 'ULK':48, 'WNK':49, 'CDC7':55, 'ATM':190}

    # hrd motif (H of HRD)
    # hrd=   {'AGC':122, 'BUB':124, 'CAMK':121, 'CK1':118, 'CMGC':122, 'EIF2AK41':147, 'HASP':164, 'MAP3K123':117, 'MOS':146, 
    hrd=   {'AGC':122, 'BUB':124, 'CAMK':121, 'CK1':118, 'CMGC':128, 'EIF2AK41':147, 'HASP':164, 'MAP3K123':117, 'MOS':146, 
            'NEK':125, 'OTHER':123, 'PAN3':168, 'PEAK':185, 'PIK3R4':126, 'PINK1':214, 'PKDCC':146, 'POMK':129, 
            'PXK':149, 'RGC':120, 'RNASEL':131, 'RPS6KC1':597, 'SCYL':121, 'STE':121, 'STK31':151, 'TBCK':114, 'TKL':122, 
            'TP53RK':129, 'TYR':126, 'ULK':121, 'WNK':128, 'CDC7':126, 'ATM':337}

    # dfg motif (D of DFG)
    # asp=    {'AGC':142, 'BUB':155, 'CAMK':143, 'CK1':141, 'CMGC':142, 'EIF2AK41':167, 'HASP':204, 'MAP3K123':137, 'MOS':166, 
    asp=    {'AGC':142, 'BUB':155, 'CAMK':143, 'CK1':141, 'CMGC':150, 'EIF2AK41':167, 'HASP':204, 'MAP3K123':137, 'MOS':166, 
             'NEK':145, 'OTHER':144, 'PAN3':188, 'PEAK':230, 'PIK3R4':146, 'PINK1':238, 'PKDCC':165, 'POMK':154, 
             'PXK':168, 'RGC':140, 'RNASEL':151, 'RPS6KC1':617, 'SCYL':141, 'STE':141, 'STK31':171, 'TBCK':134, 'TKL':142, 
             'TP53RK':154, 'TYR':146, 'ULK':143, 'WNK':149, 'CDC7':147, 'ATM':358}
    
    # beginning of beta strand B4 (H of HPN)
    # hpn=   {'AGC':58,'BUB':55,'CAMK':57,'CK1':54,'CMGC':56,'EIF2AK41':79,'HASP':70, 'MAP3K123':53, 'MOS':63,
    hpn=   {'AGC':58,'BUB':55,'CAMK':57,'CK1':54,'CMGC':59,'EIF2AK41':79,'HASP':70, 'MAP3K123':53, 'MOS':63,
            'NEK':57,'OTHER':56,'PINK1':125,'PAN3':74,'PEAK':67,'PIK3R4':63,'PKDCC':78, 'POMK':58,
            'PXK':69,'RGC':54,'RNASEL':62, 'RPS6KC1':53, 'SCYL':61,'STE':56,'STK31':86,'TBCK':50,'TKL':55,
            'TP53RK':62,'TYR':60,'ULK':57,'WNK':58, 'CDC7':65, 'ATM':217}

    # middle of beta strand B5 (Y of YLV; Y207 of "YLI" in AURKA)
    # ylv=   {'AGC':75, 'BUB':72, 'CAMK':74, 'CK1':71, 'CMGC':75, 'EIF2AK41':100, 'HASP':117, 'MAP3K123':70, 'MOS':85,
    ylv=   {'AGC':75, 'BUB':72, 'CAMK':74, 'CK1':71, 'CMGC':80, 'EIF2AK41':100, 'HASP':117, 'MAP3K123':70, 'MOS':85,
            'NEK':74, 'OTHER':74, 'PINK1':169, 'PAN3':95, 'PEAK':131,'PIK3R4':80, 'PKDCC':100, 'POMK':73,
            'PXK':86, 'RGC':71, 'RNASEL':79, 'RPS6KC1':70, 'SCYL':78,'STE':73, 'STK31':105, 'TBCK':68, 'TKL':72,
            'TP53RK':79, 'TYR':77, 'ULK':74, 'WNK':79, 'CDC7':82, 'ATM':231}

    # end of E-helix (Y of YLH)
    #  ylh=   {'AGC':114,'BUB':116,'CAMK':113,'CK1':110,'CMGC':114,'EIF2AK41':139,'HASP':155,'MAP3K123':109, 'MOS':138,
    ylh=   {'AGC':114,'BUB':116,'CAMK':113,'CK1':110,'CMGC':120,'EIF2AK41':139,'HASP':155,'MAP3K123':109, 'MOS':138,
            'NEK':117,'OTHER':115,'PINK1':206,'PAN3':160,'PEAK':177,'PIK3R4':118,'PKDCC':135,'POMK':118,
            'PXK':141,'RGC':111,'RNASEL':122, 'RPS6KC1':589, 'SCYL':112,'STE':113,'STK31':143,'TBCK':106,'TKL':112,
            'TP53RK':121,'TYR':118,'ULK':113,'WNK':118, 'CDC7':118, 'ATM':333}

    # end of activation loop (E of APE)
    #  ape={'AGC':170, 'BUB':184, 'CAMK':169, 'CK1':175, 'CMGC':168, 'EIF2AK41':192, 'HASP':238, 'MAP3K123':163, 'MOS':196, 
    ape={'AGC':170, 'BUB':184, 'CAMK':169, 'CK1':175, 'CMGC':175, 'EIF2AK41':192, 'HASP':238, 'MAP3K123':163, 'MOS':196,
         'NEK':172, 'OTHER':173, 'PAN3':205, 'PEAK':256, 'PIK3R4':180, 'PINK1':271, 'PKDCC':186, 'POMK':184, 
         'PXK':183, 'RGC':168, 'RNASEL':170, 'RPS6KC1':641, 'SCYL':171, 'STE':168, 'STK31':198, 'TBCK':160, 'TKL':173, 
         'TP53RK':176, 'TYR':175, 'ULK':169, 'WNK':174, 'CDC7':334, 'ATM':373}
    
    # from Kincore paper
    # ATP binding region—hinge  211–213
    # Back pocket—C-helix       166–193
    # partial regions of B4     196-204
    # partial region of B5      205–207
    # DFG                       273–275
    # Type 2-only pocket        184, 188, 247, 254
    # back_pocket1 = range(lys+4,hpn+4), Align col 106-185, AURKA 166-194 (not 193 as in paper)
    # type2_resi (Glu+3,Glu+7,ylh,hrd) 
    # back_pocket2: Align col 187-196  = backpocket1C+2, backpocket2N+9
    # back_pocket3 = range(ylv-2, ylv+1, Align col 420-423 
    # frontpocket = lys+4, Glu-2,  Align col 106-144  
    # xdf = XDF,XDF+3, Align col 1337-1340 
    ape9={}
    glu4={}
    xdfg={}
    phe={}
    gly={}
    dfg6={}
    xhrd={}
    arg={}
    hinge1={}
    back_pocket1={}
    back_pocket2={}
    back_pocket3={}
    type2_resi={}
    front_pocket={}
    xdf_residues={}
    for family in glu:
        ape9[family] =         ape[family]-8
        glu4[family] =         glu[family]+4
        arg[family] =          hrd[family]+1
        xhrd[family] =         hrd[family]-1
        xdfg[family] =         asp[family]-1
        phe[family] =          asp[family]+1
        gly[family] =          asp[family]+2
        dfg6[family] =         asp[family]+5
        hinge1[family]=         ylv[family]+4
        type2_resi[family] =   (glu[family]+3, glu[family]+7, ylh[family], hrd[family])
        back_pocket1[family] = range(lys[family]+4, hpn[family]+4)
        back_pocket2[family] = range(hpn[family]+6, hpn[family]+15)
        back_pocket3[family] = range(ylv[family]-2, ylv[family]+1)
        front_pocket[family] = range(lys[family]+4, glu[family]-2)
        xdf_residues[family] = range(asp[family]-1, asp[family]+2)

    type2_resi_num=set()
    back_pocket1_num=set()
    back_pocket2_num=set()
    back_pocket3_num=set()
    front_pocket_num=set()
    xdf_residues_num=set()

    #this is the way list is assigned to pandas, need to do it only for the first index where the column is defined
    conf_df.at[index,'Type2_resi_num']   = [type2_resi_num]   
    conf_df.at[index,'Back_pocket1_num'] = [back_pocket1_num]
    conf_df.at[index,'Back_pocket2_num'] = [back_pocket2_num]
    conf_df.at[index,'Back_pocket3_num'] = [back_pocket3_num]
    conf_df.at[index,'Front_pocket_num'] = [front_pocket_num]
    conf_df.at[index,'XDF_residues_num'] = [xdf_residues_num]

    first_res=conf_df.at[index,'First_res']
    if s2cdict is not None: first_res=0
    hmmerfilename=f'{pdbfilename[0:-4]}_{model_id}_{chain_id}_{group}.hmmer.txt'
    hmm_result=SearchIO.read(hmmerfilename, format='hmmer3-text')

    for hits in hmm_result:     #extract hit from alignment in HMM output file
        for hsps in hits:
            if hsps.bitscore<5: continue   # skip minor alignments which are wrong

            col_num=0
            hmm_index=hsps.query_start
            hit_index=hsps.hit_start+first_res-1
            alignment_length=hsps.aln.get_alignment_length()

            for hmm_res in hsps.aln[0]:
                col_num=col_num+1
                if hmm_res!='.':
                    hmm_index=hmm_index+1
                if hsps.aln[1][col_num-1]!='-':
                    hit_index=hit_index+1

                    s2cnum=hit_index
                    if s2cdict is not None: s2cnum=s2cdict[hit_index+1]

                if hmm_res ==".": continue
                
                if hmm_index==lys[group]:
                    conf_df.at[index,'Lys_num']=s2cnum
                    conf_df.at[index,'Lys_restype']=hsps.aln[1][col_num-1].upper()

                if hmm_index==glu[group]:
                    conf_df.at[index,'Glu_num']=s2cnum
                    conf_df.at[index,'Glu_restype']=hsps.aln[1][col_num-1].upper()
                    if conf_df.at[index,'Glu_restype'] == "-":
                        if hsps.aln[1][col_num-2].upper() != "-":
                            conf_df.at[index,'Glu_num']=s2cnum-1
                            conf_df.at[index,'Glu_restype']=hsps.aln[1][col_num-2].upper()
                            print("# ",pdbfilename,"Glu shifted left one residue", conf_df.at[index,'Glu_num'], conf_df.at[index,'Glu_restype'])
                        elif hsps.aln[1][col_num].upper() != "-":
                            conf_df.at[index,'Glu_num']=s2cnum+1
                            conf_df.at[index,'Glu_restype']=hsps.aln[1][col_num].upper()
                            print("# ",pdbfilename,"Glu shifted right one residue", conf_df.at[index,'Glu_num'], conf_df.at[index,'Glu_restype'])

                if hmm_index==glu4[group]:
                    conf_df.at[index,'Glu4_num']=s2cnum
                    conf_df.at[index,'Glu4_restype']=hsps.aln[1][col_num-1].upper()
                    if conf_df.at[index,'Glu4_restype'] == "-":
                        if hsps.aln[1][col_num].upper() != "-":
                            conf_df.at[index,'Glu4_num']=s2cnum+1
                            conf_df.at[index,'Glu4_restype']=hsps.aln[1][col_num].upper()
                            print("# ",pdbfilename,"Glu4 shifted right one residue", conf_df.at[index,'Glu4_num'], conf_df.at[index,'Glu4_restype'])
                        elif hsps.aln[1][col_num-2].upper() != "-":
                            conf_df.at[index,'Glu4_num']=s2cnum-1
                            conf_df.at[index,'Glu4_restype']=hsps.aln[1][col_num-2].upper()
                            print("# ",pdbfilename,"Glu4 shifted left one residue", conf_df.at[index,'Glu4_num'], conf_df.at[index,'Glu4_restype'])

                # RLD added these  12/22; ape9 find by traversing hit sequence and counting residues
                if hmm_index==ape[group]:
                    conf_df.at[index,'APE_num']=s2cnum
                    conf_df.at[index,'APE_restype']=hsps.aln[1][col_num-1].upper()

                    rescount=0
                    for i in range(0,20):
                        if col_num-1-i < 0: continue
                        if hsps.aln[1][col_num-1-i] != "-":
                            rescount = rescount+1
                        else:
                            continue

                        if rescount==9 and conf_df.at[index,'APE9_num'] == 0: 
                            conf_df.at[index,'APE9_num']=s2cnum-rescount+1
                            conf_df.at[index,'APE9_restype']=hsps.aln[1][col_num-1-i].upper()
                            break
                        
                if hmm_index==hrd[group]:
                    conf_df.at[index,'HRD_num']=s2cnum
                    conf_df.at[index,'HRD_restype']=hsps.aln[1][col_num-1].upper()

                if hmm_index==xhrd[group]:
                    conf_df.at[index,'XHRD_num']=s2cnum
                    conf_df.at[index,'XHRD_restype']=hsps.aln[1][col_num-1].upper()

                if hmm_index==arg[group]:
                    conf_df.at[index,'Arg_num']=s2cnum
                    conf_df.at[index,'Arg_restype']=hsps.aln[1][col_num-1].upper()

                if hmm_index==xdfg[group]:
                    conf_df.at[index,'XDFG_num']=s2cnum
                    conf_df.at[index,'XDFG_restype']=hsps.aln[1][col_num-1].upper()

                if hmm_index==asp[group]:
                    conf_df.at[index,'Asp_num']=s2cnum
                    conf_df.at[index,'Asp_restype']=hsps.aln[1][col_num-1].upper()

                    rescount=0
                    for i in range(0,20):
                        if  col_num+i > alignment_length: continue
                        if hsps.aln[1][col_num-1+i] != "-": rescount = rescount+1

                        if rescount==6 and conf_df.at[index,'DFG6_num'] == 0: 
                            conf_df.at[index,'DFG6_num']=s2cnum+rescount-1
                            conf_df.at[index,'DFG6_restype']=hsps.aln[1][col_num-1+i].upper()
                            break

                if hmm_index==hpn[group]:
                    rescount=0
                    for i in range(0,20):
                        if  col_num+i > alignment_length: continue
                        if hsps.aln[1][col_num-1+i] != "-": rescount = rescount+1

                        if rescount==7 and conf_df.at[index,'HPN7_num'] == 0: 
                            conf_df.at[index,'HPN7_num']=s2cnum+rescount-1
                            conf_df.at[index,'HPN7_restype']=hsps.aln[1][col_num-1+i].upper()
                            break
                        
                if hmm_index==phe[group]:
                    conf_df.at[index,'Phe_num']=s2cnum
                    conf_df.at[index,'Phe_restype']=hsps.aln[1][col_num-1].upper()

                if hmm_index==gly[group]:
                    conf_df.at[index,'Gly_num']=s2cnum
                    conf_df.at[index,'Gly_restype']=hsps.aln[1][col_num-1].upper()

                if hmm_index==hinge1[group]:
                    conf_df.at[index,'Hinge1_num']=int(s2cnum)

                if hmm_index in type2_resi[group]:
                    type2_resi_num.add(int(s2cnum))

                if hmm_index in back_pocket1[group]:
                    back_pocket1_num.add(int(s2cnum))

                if hmm_index in back_pocket2[group]:
                    back_pocket2_num.add(int(s2cnum))

                if hmm_index in back_pocket3[group]:
                    back_pocket3_num.add(int(s2cnum))

                if hmm_index in xdf_residues[group]:
                    xdf_residues_num.add(int(s2cnum))

                if hmm_index in front_pocket[group]:
                    front_pocket_num.add(int(s2cnum))
               
    conf_df.at[index,'Type2_resi_num']  =type2_resi_num   
    conf_df.at[index,'Back_pocket1_num']=back_pocket1_num
    conf_df.at[index,'Back_pocket2_num']=back_pocket2_num
    conf_df.at[index,'Back_pocket3_num']=back_pocket3_num
    conf_df.at[index,'Front_pocket_num']=front_pocket_num
    conf_df.at[index,'XDF_residues_num']=xdf_residues_num
    return conf_df

