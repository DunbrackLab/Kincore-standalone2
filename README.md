# Kincore-standalone
# Use the following command to create a virtual enviroment and install the necessary packages.

$git clone https://github.com/DunbrackLab/Kincore-standalone2

$cd Kincore-standalone2

# (or "cd Kincore-standalone2-main" if you downloaded the zip file from github)

$conda create --name 'kincore-standalone2' python=3.8 pandas numpy biopython hmmer --channel conda-forge --channel bioconda

$conda activate kincore-standalone2

$python3 kinase_state.py -h

# To run Kincore-standalone2 for a single protein structure, use these commands:
$python3 kinase_state.py 1ol5.cif
$python3 kinase_state.py 1ol5.pdb
$python3 kinase_state.py 1ol5.cif.gz
$python3 kinase_state.py 1ol5.pdb.gz

# To run Kincore-standalone2 on a set of structures, place the filenames (or full paths plus the filenames) in a file with the extension ".txt" and give this command:
$python3 kinase_state.py list.txt > output.txt

# The output of Kincore-standalone2 consists of 5 lines for each kinase structure:
$python3 kinase_state.py 1ol5.cif
1ol5  0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Labels    family CAMK  hmm CAMK     score  271.0   Active   DFGin    BLAminus Chelix-in   SaltBr-in   ActLoopNT-in   ActLoopCT-in   Spine-in  
1ol5  0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Residues  Lys.K162  Glu.E181  Glu4.Q185 HPN7.L196  XHRD.I253  HRD.H254  Arg.R255  XDFG.A273  Asp.D274  Phe.F275  Gly.G276  DFG6.V279 APE9.G291 APE.E299 
1ol5  0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Distances Glu4_Phe   5.58 Lys_Phe  15.10 Lys_Glu   9.45 SaltBr   2.91 DFG6_XHRD   2.93 APE9_Arg   3.68 Spine    3.59   3.51   3.91   3.91
1ol5  0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Dihedrals X -139.20 -172.08 D   52.22   80.16 -163.54   -3.74 F  -90.83   24.66  282.61   72.55 G  -49.87  -44.18
1ol5  0 A    Active_DFGin_BLAminus_SBin_NTin_CTin          Ligands   ADP:1388,MG:1389,MG:1390,MG:1394    Type1,Allosteric,Allosteric,Allosteric

# Each line starts with the filename (with the extension and pathname removed), the model number, and chain. In the example, these are "1ol5  0 A". This is followed by a long string that contains the state of the kinase in abbreviated format, consisting of the Activity_label, the Spatial_label, the Dihedral_label, the SaltBr_label, the ActLoopNT_label, and the ActLoopCT_label, e.g. Active_DFGin_BLAminus_SBin_NTin_CTin.

# The first line contains the main kinase family (AGC, CAMK, CK1, CMGC, NEK, OTHER, RGC, STE, TKL, TYR) as well as the HMM model and HMM score used to identify the family. There are extra HMMs for unusual kinases, such as BUB, PEAK, and TP53RK that are members of the OTHER family. After the score follows the conformational labels including the Spine_label at the end of the line.
Labels    family CAMK  hmm CAMK     score  271.0   Active   DFGin    BLAminus Chelix-in   SaltBr-in   ActLoopNT-in   ActLoopCT-in   Spine-in

# The second line provides the residue identifications used to make the conformational assignments. Each is named followed by a dot and then the residue type in one letter code and residue number. The residues in order are: the (1) Lys and (2) Glu of the N-terminal domain saltbridge; (3) the Glu4 residue (4 residues after the saltbridge Glu) used to calculate the spatial label; (4) the HPN7 residue (7 residues after the conserved HPN motif in the loop following the C-helix, HPNxxxX) which is the 4th spine residue;  the (5) XHRD, (6) His of HRD, and (7) Arg of HRD; the XDFG motif residues consisting of (8) X, (9) Asp, (10) Phe, (11) Gly, and (12) DFG6 residue, which is used to determine the ActLoopNT state along with the XHRD residue; the (13) APE9 and APE residues. The APE9 residue is used to determine the ActLoopCT state along with the HRD-Arg residue.
Residues  Lys.K162  Glu.E181  Glu4.Q185 HPN7.L196  XHRD.I253  HRD.H254  Arg.R255  XDFG.A273  Asp.D274  Phe.F275  Gly.G276  DFG6.V279 APE9.G291 APE.E299 

# The third line contains the distances measured to determine the conformational state of the kinase. These distances are: (1) the Glu4-CA/Phe-CZ and (2) the Lys-CA/Phe-CZ distances used to determine the DFGin, DFGinter, or DFGout label; (3) the Lys-CB/Glu-CB distance that determines the Chelix-in or Chelix-out label; (4) the Lys-NZ/Glu-OE distance that determines whether the saltbridge exists (SaltBr-in) or not (SaltBr-out); (5) the DFG6-XHRD hydrogen bond distance used to determine the ActLoopNT status; (6) the APE9-HRDArg distance useed to determine the ActLooPCT status; and (7) the three spine distances followed by the maximum spine distance. The three spine distances are between the side chain pairs: HRD-His/DFG-Phe, DFG-Phe/Glu4, and Glu4/HPN7.
Distances Glu4_Phe   5.58 Lys_Phe  15.10 Lys_Glu   9.45 SaltBr   2.91 DFG6_XHRD   2.93 APE9_Arg   3.68 Spine    3.59   3.51   3.91   3.91

# The fourth line contains the backbone dihedral angles of the XDFG motif and the side-chain dihedrals of the Asp and Phe residues.
Dihedrals X -139.20 -172.08 D   52.22   80.16 -163.54   -3.74 F  -90.83   24.66  282.61   72.55 G  -49.87  -44.18

# The fifth line contains the ligands and their types. The ligands and ligand numbers are given first separated by commas, followed by the ligand types, separated by commas. 
Ligands   ADP:1388,MG:1389,MG:1390,MG:1394    Type1,Allosteric,Allosteric,Allosteric
