Run the files in this directory using these steps:

1 Use this command to make a list of files

  ls -1 *.cif > test.txt

2 Type this to set up Kincore-standalone2

  conda activate Kincore-standalone2

3 Type this command to run the test structures

  python3 path_to/kinase_state.py test.txt > test2.out

4 Compare the results with the file in this folder, test.out

  diff test.out test2.out

