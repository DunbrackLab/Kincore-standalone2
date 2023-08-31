#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Last edited on Wed Aug 30, 2023 by RLD

from Bio.PDB import MMCIFParser, PDBParser, MMCIF2Dict
from Bio.SeqUtils import seq1, seq3
import os, sys, gzip, re
from Bio import PDB

# get file handle for gzipped or not gzipped files
def open_file(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

# get sequence from PDB format files from SEQRES lines
def get_sequence_from_seqres(filename, chain_id):
    sequence = []
    with open_file(filename) as file:
        for line in file:
            if line.startswith("SEQRES") and line.split()[2] == chain_id:
                seqres = line.split()[4:]
                for aa in seqres:
                    try:
                        sequence.append(seq1(aa))
                    except KeyError:
                        sequence.append("X")
    return "".join(sequence) if sequence else None

# get sequence from coordinates in case file does not have sequence records
def get_sequence_from_coordinates(structure, chain_id, model_id):
    sequence = []
    chain = structure[model_id][chain_id]
    prev_residue_num = None
    rescount=0
    first_res=0
    for residue in chain:
        if PDB.is_aa(residue):
            curr_residue_num = residue.get_id()[1]
            if int(curr_residue_num)>50000: continue  # for PDBrenum files
            if prev_residue_num is not None:
                missing_residues = residue.get_id()[1] - prev_residue_num - 1
                sequence.extend(["X"] * missing_residues)
            sequence.append(seq1(residue.get_resname()))
            prev_residue_num = residue.get_id()[1]
            rescount += 1
            if rescount==1: first_res = residue.get_id()[1]
    return ''.join(sequence), first_res

# main routine to get sequence from file of different types
def get_sequence_from_file(filename, chain_id, model_id):

    # get filename extension (.cif or .pdb)
    filenamelist=filename.split(".")
    if filenamelist[-1]=="gz":
        gzipped=True
        file_extension=filenamelist[-2]
    else:
        gzipped=False
        file_extension=filenamelist[-1]
        
    # Parse the structure file
    if file_extension == "cif":
        with open_file(filename) as file:
            mmcif_dict = MMCIF2Dict.MMCIF2Dict(file)
        parser = MMCIFParser()
        parser.QUIET = True
        with open_file(filename) as file:
            structure = parser.get_structure("structure_id", file)
    elif file_extension == "pdb":
        parser = PDBParser(QUIET=True)
        with open_file(filename) as file:
            structure = parser.get_structure("structure_id", file)
    else:
        raise ValueError("Invalid file format. Please use .cif or .pdb file.")

    # Extract the protein sequence
    protein_sequence=None
    first_res=None
    if file_extension == "cif" and "_entity_poly.pdbx_seq_one_letter_code" in mmcif_dict:
        chain_list=mmcif_dict["_entity_poly.pdbx_strand_id"]
        protein_sequence_list=mmcif_dict["_entity_poly.pdbx_seq_one_letter_code"]
        for i in range(len(chain_list)):
            chains=chain_list[i].split(",")
            if chain_id in chains:
                protein_sequence=protein_sequence_list[i]
        protein_sequence = ''.join(protein_sequence.splitlines())
        protein_sequence = re.sub('\([A-Z0-9][A-Z0-9][A-Z0-9]\)',"X",protein_sequence) # modified residues replaced with X
    elif file_extension == "pdb":
        protein_sequence = get_sequence_from_seqres(filename, chain_id)

    protein_sequence_from_coor,first_res = get_sequence_from_coordinates(structure, chain_id, model_id)

    # get s2c information (sequence to residue number correspondence)
    s2cdict=None
    if file_extension == "cif" and "_pdbx_poly_seq_scheme.pdb_seq_num" in mmcif_dict:
        s2cdict={}
        chainidlist=mmcif_dict["_pdbx_poly_seq_scheme.pdb_strand_id"]
        seqidlist=mmcif_dict["_pdbx_poly_seq_scheme.seq_id"]
        pdbseqnumlist=mmcif_dict["_pdbx_poly_seq_scheme.pdb_seq_num"]
        for i in range(len(chainidlist)):
            if chainidlist[i] != chain_id: continue
            seqid=int(seqidlist[i])
            pdbseqnum=int(pdbseqnumlist[i])
            s2cdict[seqid]=pdbseqnum

    frag_index=0
    if protein_sequence is None:
        protein_sequence = protein_sequence_from_coor
    else:
        frag=protein_sequence_from_coor[0:10]
        frag=frag.replace("X",".")
        frag_search=re.search(frag,protein_sequence)
        if frag_search is not None:
            frag_index=frag_search.start()
        else:
            frag=protein_sequence_from_coor[0:5]
            frag=frag.replace("X",".")
            frag_search=re.search(frag,protein_sequence)
            if frag_search is not None:
                frag_index=frag_search.start()
    if frag_index >=0: first_res = first_res - frag_index
    else: print("error", filename, frag, protein_sequence_from_coor, protein_sequence)
    return protein_sequence, first_res, s2cdict

# extract sequence from PDB file in pdb or cif format
def extract_seq(pdbfilename,index,conf_df):
    model_id=conf_df.at[index,'Model_id']
    chain_id=conf_df.at[index,'Chain_id']
    if chain_id==" ": chain_id=""
    protein_sequence, first_res, s2cdict = get_sequence_from_file(pdbfilename, chain_id, model_id)

    conf_df.at[index,'Sequence']=protein_sequence
    conf_df.at[index,'First_res']=first_res
    fhandle_outputseq=open(f"{pdbfilename[0:-4]}_{chain_id}.fasta",'w')
    fhandle_outputseq.write(f">{pdbfilename[0:-4]}_{chain_id}\n")
    fhandle_outputseq.write(protein_sequence + "\n")
    fhandle_outputseq.close()
    return conf_df, s2cdict
