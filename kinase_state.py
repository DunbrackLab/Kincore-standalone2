#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri July 31 17:35:32 2020
Last edited on Wed August 30, 2023 by RLD
@author: vivekmodi
"""

import pandas as pd
import sys, os, gzip, argparse
import multiprocessing
from Bio import PDB
from warnings import simplefilter
import multiprocessing

from modules.active_labels          import active_labels
from modules.assign_default_values  import assign_default_values
from modules.chelix                 import chelix_conformation
from modules.classify_ligands       import classify_ligands
from modules.compute_dihedrals      import compute_dihedrals
from modules.compute_distance       import compute_distance
from modules.delete_files           import delete_files
from modules.dihedral_label         import dihedral_label
from modules.extract_ligands        import extract_ligands
from modules.extract_sequence       import extract_seq
from modules.identify_group         import identify_group
from modules.identify_residues      import identify_residues, get_real_residue_number
from modules.identify_restypes      import identify_restypes
from modules.run_hmmsearch          import run_hmmsearch
from modules.spatial_label          import spatial_label

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

#Function to read input file list (filename.txt)
def read_inputlist_pool(listfilename, hmm_loc):
    
    fileset=set()
    INPUTFILE=open(listfilename,'r')
    maxlen=0
    for line in INPUTFILE:
        line=line.strip()
        linelist=line.split()
        name=linelist[0]  # first item on line
        fileset.add(name)
        line=line.replace('.pdb','')
        line=line.replace('.cif.gz','')
        line=line.replace('.cif','')
        if len(line)>maxlen: maxlen=len(line)
    INPUTFILE.close()
    
    with multiprocessing.Pool() as pool:
        results = pool.map(identify_state, [(filename, hmm_loc, maxlen) for filename in fileset])

    return results
    
#Function to read input file list (filename.txt)
def read_inputlist(listfilename, hmm_loc):
    
    fileset=set()
    INPUTFILE=open(listfilename,'r')
    maxlen=0
    for line in INPUTFILE:
        line=line.strip()
        linelist=line.split()
        name=linelist[0]  # first item on line
        fileset.add(name)
        line=line.replace('.pdb','')
        line=line.replace('.cif.gz','')
        line=line.replace('.cif','')
        if len(line)>maxlen: maxlen=len(line)
    INPUTFILE.close()
    
    for filenane in fileset:
        (name,string)=identify_state((filename, hmm_loc, maxlen))
        print(string)
    
#Identify the state of a kinase structure in pdbfilename
def identify_state(args):
    (pdbfilename, pwd, maxfilenamelen)=args
    conf_df=pd.DataFrame()
    chain_list=list()

    if len(pdbfilename)<3: return
    #Check if the file is compressed
    if pdbfilename[-3:].lower()=='.gz':
        handle=gzip.open(pdbfilename, 'rt')
    else:
        handle=open(pdbfilename, 'r')

    # Biopython parser
    if pdbfilename[-4:].lower()=='.cif' or pdbfilename[-7:].lower()==".cif.gz":
        parser=PDB.MMCIFParser(QUIET=True)
    elif pdbfilename[-4:].lower()=='.pdb' or pdbfilename[-7:].lower()==".pdb.gz":
        parser=PDB.PDBParser(QUIET=True)
    else:
        return pdbfilename, "NotPDB","NotPDB"
    
    structure=parser.get_structure(pdbfilename, handle)
    string=""
    errorstring=""
    index=-1
    for model in structure:
        for chain in model:

            # add chains to list; skip chains that are too short
            if len(chain.get_list())<=200:  # mininum length of kinase chain
                errorstring += f'# {pdbfilename} Model {int(model.id)}, Chain {chain.id} is too short to be a kinase (length {len(chain.get_list())}).\n'
                continue 

            index+=1
            chain_list.append(chain.id)

            # default values
            conf_df.at[index, 'Model_id']=int(float(model.id))
            conf_df.at[index, 'Chain_id']=chain.id
            conf_df=assign_default_values(index, conf_df)

            # get sequence and set up s2cdict (sequence-coordinate correspondence)
            s2cdict=None
            conf_df,s2cdict=extract_seq(pdbfilename, index, conf_df)

            # run hmmsearch to identify kinase family 
            # if no group is identified, delete HMM files and go to next chain/model
            run_hmmsearch(pwd, pdbfilename, index, conf_df)
            conf_df=identify_group(pdbfilename, index, conf_df)
            group=conf_df.at[index,'Group']
            if group=='None':
                delete_files(pdbfilename, conf_df.at[index, 'Model_id'],
                             conf_df.at[index, 'Chain_id'], conf_df.at[index, 'Group']) 
                errorstring += f'# {pdbfilename} Model {int(model.id)}, Chain {chain.id} is probably not a protein kinase.\n'
                continue

            # identify key residues from HMM alignment and s2cdict
            # if any of K, E, X, D, F are not identified: print errors but continue since some features may still be calculated
            conf_df=identify_residues(pdbfilename, index, conf_df, s2cdict)
            for item in ('Lys_restype','Glu4_restype','Glu_restype','XDFG_restype','Phe_restype','Gly_restype',
                         'Asp_restype','XHRD_restype','HRD_restype','DFG6_restype', 'HPN7_restype',
                         'APE_restype', 'APE9_restype'):
                (resname,restype)=item.split("_")
                if conf_df.at[index, item]== 'X': 
                    errorstring += (f'# {pdbfilename} '
                                    f'Group {conf_df.at[index,"Group"]} '
                                    f'Score {conf_df.at[index,"Score"]:0.1f} '
                                    f'Model {int(conf_df.at[index,"Model_id"])}, '
                                    f'Chain {conf_df.at[index,"Chain_id"]} {resname} residue is missing.\n')

                if conf_df.at[index, item]== '-': 
                    errorstring += (f'# {pdbfilename} '
                                    f'Group {conf_df.at[index,"Group"]} '
                                    f'Score {conf_df.at[index,"Score"]:0.1f} '
                                    f'Model {int(conf_df.at[index,"Model_id"])}, '
                                    f'Chain {conf_df.at[index,"Chain_id"]} {resname} residue may be misaligned.\n')
                    
            # calculate distances, dihedrals, spatial label, dihedral label, chelix conformation, active labels, ligand labels
            conf_df= compute_distance(pdbfilename, index, conf_df, structure)
            conf_df= compute_dihedrals(pdbfilename, index, conf_df, structure)
            conf_df= spatial_label(index, conf_df)
            conf_df= dihedral_label(index, conf_df, 0.45)
            conf_df= chelix_conformation(index, conf_df)
            conf_df= active_labels(index, conf_df)
            conf_df= extract_ligands(pdbfilename, index, conf_df, structure)
            conf_df= classify_ligands(pdbfilename, index, conf_df, structure)

            # get PDB residue numbers for identified residues
            model_id=conf_df.at[index,"Model_id"]
            chain_id=conf_df.at[index,"Chain_id"]

            for resname in ("XDFG","Lys","Glu","Glu4","Asp","Phe","Gly","HRD","XHRD","Arg","APE","APE9","DFG6","HPN7"):
                res1=int(conf_df.at[index,resname+"_num"])
                if s2cdict is not None:
                    conf_df.at[index, f'real_{resname}_num'] = conf_df.at[index, f'{resname}_num']
                else:
                    conf_df.at[index, "real_" + resname + "_num"]=get_real_residue_number(pdbfilename, model_id, chain_id, res1, structure)
                conf_df.at[index, resname+"_string"]=  f"{resname}.{conf_df.at[index, resname+'_restype']}{int(conf_df.at[index, 'real_'+resname+'_num'])}"

            # Determine final group from hmm group
            finalgroup=group
            if group != "None":
                if group=="CDC7": finalgroup="OTHER"
                if group=="PAN3": finalgroup="OTHER"
                if group=="TP53RK": finalgroup="OTHER"
                if group=="PIK3R4": finalgroup="OTHER"
                if group=="RNASEL": finalgroup="OTHER"
                if group=="PKDCC": finalgroup="OTHER"
                if group=="BUB": finalgroup="OTHER"
                if group=="ULK": finalgroup="OTHER"
                if group=="WNK": finalgroup="OTHER"
                if group=="PINK1": finalgroup="OTHER"
                if group=="EIF2AK41": finalgroup="OTHER"
                if group=="PEAK": finalgroup="OTHER"
                if group=="POMK": finalgroup="OTHER"
                if group=="HASP": finalgroup="OTHER"
                if group=="SCYL": finalgroup="OTHER"
                if group=="STK31": finalgroup="OTHER"
                if group=="PXK": finalgroup="OTHER"
                if group=="MOS": finalgroup="OTHER"
                if group=="TBCK": finalgroup="OTHER"
                if group=="RPS6KC1": finalgroup="OTHER"
                if group=="MAP3K123": finalgroup="TKL"
                if conf_df.at[index, 'Score']<190.0:
                    finalgroup="OTHER"  # some OTHER kinases have higher scores for CAMK, CMGC, etc than the OTHER HMM

            # process filename
            shortfilename=pdbfilename
            for extension in (".cif.gz", ".pdb.gz", ".cif", ".pdb"):
                shortfilename=shortfilename.replace(extension,"")
            shortlist=shortfilename.split("/")
            shortfilename=shortlist[-1]
            shortlist=shortfilename.split("_")
            if len(shortlist)>1:
                kinase=shortlist[0] + "_" + shortlist[1]
            else:
                kinase=shortlist[0]
            if int(maxfilenamelen)>3:
                filenamelen=maxfilenamelen
            else:
                filenamelen=len(shortfilename)

            shortstate=f'{conf_df.at[index,"Activity_label"]}_{conf_df.at[index,"Spatial_label"]}_{conf_df.at[index,"Dihedral_label"]}_{conf_df.at[index,"SaltBr_label"]}_{conf_df.at[index,"ActLoopNT_label"]}_{conf_df.at[index,"ActLoopCT_label"]}'
            shortstate=shortstate.replace("ActLoop","")
            shortstate=shortstate.replace("SaltBr","SB")
            shortstate=shortstate.replace("-","")
            
            # output string
            introstring= (f'{kinase:14} '
                          f'{shortfilename:{filenamelen}} '
                          f'{conf_df.at[index, "Model_id"]:2.0f} '
                          f'{conf_df.at[index, "Chain_id"]:<4} '
                          f'{shortstate:<44}  '
                          )

            string =string +  (
                f'{introstring}'
                f'{"Labels":9} '
                f'family {finalgroup:<5} '
                f'hmm {group:<8} '
                f'score {conf_df.at[index,"Score"]:6.1f}   ' 
                f'{conf_df.at[index, "Activity_label"]:<8} '
                f'{conf_df.at[index, "Spatial_label"]:<8} '
                f'{conf_df.at[index, "Dihedral_label"]:<8} '
                f'{conf_df.at[index, "Chelix_label"]:<11} '
                f'{conf_df.at[index, "SaltBr_label"]:<11} '
                f'{conf_df.at[index, "ActLoopNT_label"]:<14} '
                f'{conf_df.at[index, "ActLoopCT_label"]:<14} '
                f'{conf_df.at[index, "Spine_label"]:<10}\n'

                f'{introstring}'
                f'{"Residues":9} '
                f'{conf_df.at[index, "Lys_string"]:<9} '
                f'{conf_df.at[index, "Glu_string"]:<9} '
                f'{conf_df.at[index, "Glu4_string"]:<9} '
                f'{conf_df.at[index, "HPN7_string"]:<10} '
                f'{conf_df.at[index, "XHRD_string"]:<10} '
                f'{conf_df.at[index, "HRD_string"]:<9} '
                f'{conf_df.at[index, "Arg_string"]:<9} '
                f'{conf_df.at[index, "XDFG_string"]:<10} '
                f'{conf_df.at[index, "Asp_string"]:<9} '
                f'{conf_df.at[index, "Phe_string"]:<9} '
                f'{conf_df.at[index, "Gly_string"]:<9} '
                f'{conf_df.at[index, "DFG6_string"]:<9} '
                f'{conf_df.at[index, "APE9_string"]:<9} '
                f'{conf_df.at[index, "APE_string"]:<9}\n'

                f'{introstring}'
                f'{"Distances":9} '
                f'Glu4_Phe {conf_df.at[index, "Glu4-Phe-dis"]:>6.2f} '
                f'Lys_Phe {conf_df.at[index, "Lys-Phe-dis"]:>6.2f} '
                f'Lys_Glu {conf_df.at[index, "Lys-Glu-dis"]:>6.2f} '
                f'SaltBr {conf_df.at[index, "LysNZ-GluOE-dis"]:>6.2f} '
                f'DFG6_XHRD {conf_df.at[index, "DFG6-XHRD-dis"]:>6.2f} '
                f'APE9_Arg {conf_df.at[index, "APE9-Arg-dis"]:>6.2f} '
                f'Spine  {conf_df.at[index, "Spine1-dis"]:>6.2f} '
                f'{conf_df.at[index, "Spine2-dis"]:>6.2f} '
                f'{conf_df.at[index, "Spine3-dis"]:>6.2f} '
                f'{conf_df.at[index, "Spine-dis"]:>6.2f}\n'

                f'{introstring}'
                f'{"Dihedrals":9} '
                f'X {conf_df.at[index, "XDFG_Phi"]:>7.2f} '
                f'{conf_df.at[index, "XDFG_Psi"]:>7.2f} '
                f'D {conf_df.at[index, "Asp_Phi"]:>7.2f} '
                f'{conf_df.at[index, "Asp_Psi"]:>7.2f} '
                f'{conf_df.at[index, "Asp_Chi1"]:>7.2f} '
                f'{conf_df.at[index, "Asp_Chi2"]:>7.2f} '
                f'F {conf_df.at[index, "Phe_Phi"]:>7.2f} '
                f'{conf_df.at[index, "Phe_Psi"]:>7.2f} '
                f'{conf_df.at[index, "Phe_Chi1"]:>7.2f} '
                f'{conf_df.at[index, "Phe_Chi2"]:>7.2f} '
                f'G {conf_df.at[index, "Gly_Phi"]:>7.2f} '
                f'{conf_df.at[index, "Gly_Psi"]:>7.2f}\n'

                f'{introstring}'
                f'{"Ligands":9} '
                f'{conf_df.at[index, "Ligand"]}    '
                f'{conf_df.at[index, "Ligand_label"]}'
                f'\n'
            )                   

            # delete hmmer files
            delete_files(pdbfilename, conf_df.at[index, 'Model_id'], conf_df.at[index, 'Chain_id'], conf_df.at[index, 'Group'])
    return pdbfilename,string, errorstring
                        

# main
if __name__ == '__main__':
    #__file__ is where kinasestate.py is; HMMs are a subdirectory of that

    hmm_loc=os.path.dirname(os.path.realpath(__file__))+'/HMMs'    
    parser=argparse.ArgumentParser()
    parser.add_argument('PDB', help='PDB file: filename.cif, filename.pdb, filename.cif.gz, filename.pdb.gz, or filename.txt. Text file is for a list of PDB filenames.')
    
    args=parser.parse_args()
    
    pdbfilename=args.PDB
    lenpdbfilename=len(pdbfilename)

    # file of pdbfilenames
    if '.txt' in pdbfilename:
        results=read_inputlist_pool(pdbfilename, hmm_loc)
        ERROR=open('error.log','w')
        for item in sorted(results):
            (pdbfilename, string, errorstring)=item
            if len(string)==0: continue
            if string.endswith('\n'):
                print(string,end='')
            else:
                print(string)
            print("")
            if "#" in errorstring: ERROR.write(errorstring)

    elif len(pdbfilename)>3:
        (filename,result,errorstring)=identify_state((pdbfilename, hmm_loc, lenpdbfilename))
        if result.endswith('\n'):
            print(result,end='')
        else:
            print(result)
        if "#" in errorstring:
            ERROR=open('error.log','w')
            ERROR.write(errorstring)
