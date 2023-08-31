# Kincore-standalone
# Use the following command to create a virtual enviroment and install the necessary packages.

$git clone https://github.com/DunbrackLab/Kincore-standalone2

$cd Kincore-standalone2

$conda create --name 'kincore-standalone2' python=3.8 pandas numpy biopython hmmer --channel conda-forge --channel bioconda

$conda activate kincore-standalone2

$python3 kinase_state.py -h
