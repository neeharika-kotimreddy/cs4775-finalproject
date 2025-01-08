# Benchmarking Study: Evaluation of Multiple Sequence Alignment Algorithms

This project evaluates four prominent Multiple Sequence Alignment (MSA) algorithms: MUSCLE, MAFFT, Clustal Omega, and Kalign. The benchmark provides insights into their performance based on speed, accuracy, and Sum-of-Pairs (SP) scores, aiding researchers in selecting the most appropriate algorithm for various biological sequence alignment tasks. This is my final project for the class CS 4775: Computational Genetics and Genomics, which I took the fall 2023 semester.

## Description

Multiple Sequence Alignment (MSA) plays a critical role in computational biology, facilitating evolutionary analysis and protein structure prediction. This project benchmarks four popular MSA algorithms by analyzing their performance on sequence datasets of varying lengths, sourced from BAliBASE Reference Sets 1 and 9. The study examines the computational speed, alignment accuracy, and SP scores of each algorithm and compares their efficiency on short versus long sequences. The goal is to guide researchers in choosing the optimal alignment algorithm for their specific data requirements.

## Getting Started

### Dependencies

Before using the project, ensure the following are installed:
- **Python 3.x** with the following libraries:
  - `os`
  - `Biopython`
  - `time`
- External tools:
  - [Clustal Omega](http://www.clustal.org/omega/)
  - [KAlign](https://github.com/TimoLassmann/kalign)
  - [MAFFT](https://mafft.cbrc.jp/alignment/software/)
  - [MUSCLE](https://github.com/rcedgar/muscle?tab=readme-ov-file)


### Installing

1. Clone the repository.
2. Install required libraries and external tools listed above.
3. Modify the `run_msa` function in `benchmark.py` to reflect the correct paths to your installed tools.

### Executing the program
1. Navigate to the project directory and run 
     ```bash
     python3 benchmark.py
2. To analyze both short and long sequences, modify the paths in the main function of benchmark.py as follows:
  Long Sequences:
   ```bash
      1. msf_ref = os.path.join(working_directory, f"reference_files/test{i}_ref.msf")
      2. input_sequences = os.path.join(working_directory, f"input_files/test{i}_input.fasta")
      3. reference_alignment = os.path.join(working_directory, f"reference_files/test{i}_ref.fasta")
      4. output_alignment_{algorithm} = os.path.join(working_directory, f"output_files/{algorithm}/test{i}_output.fasta")
  
  Short Sequences:
  ```bash
      1. msf_ref = os.path.join(working_directory, f"short_reference_files/test{i}_ref.msf")
      2. input_sequences = os.path.join(working_directory, f"short_input_files/test{i}_input.fasta")
      3. reference_alignment = os.path.join(working_directory, f"short_reference_files/test{i}_ref.fasta")
      4. output_alignment_{algorithm} = os.path.join(working_directory, f"short_output_files/{algorithm}/test{i}_output.fasta")
      
      
      
    ## Check results 
Results of the benchmark study will be written to the results.txt file

