# Steps to Reproduce Benchmark Study

## Clone this repository
## Install required libraries for algorithms 
### Install Clustal Omega 
1. Download the Source code .tar.gz : http://www.clustal.org/omega/
2. follow the steps in this install file: http://www.clustal.org/omega/INSTALL
### Install KAlign
1. follow the installation steps here: https://github.com/TimoLassmann/kalign
### Install MAFFT 
1. follow the installation steps here: https://mafft.cbrc.jp/alignment/software/
### Install MUSCLE
1. follow the installation steps here: https://github.com/rcedgar/muscle?tab=readme-ov-file


## Modify the paths in the algorithm 
in the run_msa function, ensure that you modify the paths for each algorithms to your respective path based on where you have installed tha algorithm and its required libraries

## Run the benchmark algorithm
cd into the benchmarkStudyCS4775 folder and run the following terminal command: python3 benchmark.py

## Run algorithm on short sequences and long sequences
Modify the paths of the reference, output, and input files in the main function. 
Include the following pathe for long sequences:
1. msf_ref = os.path.join(working_directory, f"reference_files/test{i}_ref.msf")
2. input_sequences = os.path.join(working_directory, f"input_files/test{i}_input.fasta")
3. reference_alignment = os.path.join(working_directory, f"reference_files/test{i}_ref.fasta")
4. output_alignment_{algorithm} = os.path.join(working_directory, f"output_files/{algorithm}/test{i}_output.fasta")

Include the following paths for short sequences:
1. msf_ref = os.path.join(working_directory, f"short_reference_files/test{i}_ref.msf")
2. input_sequences = os.path.join(working_directory, f"short_input_files/test{i}_input.fasta")
3. reference_alignment = os.path.join(working_directory, f"short_reference_files/test{i}_ref.fasta")
4. output_alignment_{algorithm} = os.path.join(working_directory, f"short_output_files/{algorithm}/test{i}_output.fasta")



## Check results 
Results of the benchmark study will be written to the results.txt file

