import os
import sys
import time
import subprocess
import psutil
from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline, TCoffeeCommandline
from Bio import AlignIO
import psutil
import os
import matplotlib.pyplot as plt
import tracemalloc

"""run python3 benchmark.py"""
def calculate_sp_score(alignment):
    sp_score = 0
    for i in range(len(alignment[0])):  # Iterate over each column
        column = alignment[:, i]
        for j in range(len(column)):
            for k in range(j + 1, len(column)):
                if column[j] == column[k]:
                    sp_score += 1  # Increase score for matching pairs
    return sp_score




def run_msa(algorithm, input_file, output_file):
    """Runs the specified algorithm on the input file and returns the execution time"""
    start_time = time.time()
    tracemalloc.start()

    if algorithm == "kalign":
        kalign_path = "/Users/likitag/Downloads/kalign/build/src/kalign"

        kalign_command = [kalign_path, '-i', input_file, '-o', output_file]

        try:
            subprocess.run(kalign_command, check=True)
            print("Kalign alignment completed successfully.")
            # Calculate the SP score for your alignment

        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running Kalign: {e}")

    if algorithm == "clustalo":
        clustal_path = "/Users/likitag/Downloads/clustal-omega-1.2.4/src/clustalo"
        cmd = [
            clustal_path, 
            '-i', input_file, 
            '-o', output_file, 
            '--force' 
        ]


        try:
            subprocess.run(cmd, check=True)

        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running {algorithm}: {e}")
            return None


    if algorithm == "muscle":
        muscle_path = "/Users/likitag/Downloads/muscle-5.1.0/src/Darwin/muscle"

        command = [
            muscle_path,
            "-align", input_file,
            "-output", output_file
        ]

        try:
            process = psutil.Process(os.getpid())
            mem_before = process.memory_info().rss / 1024 / 1024
            subprocess.run(command, check=True)
            mem_after = process.memory_info().rss / 1024 / 1024
            print(f"Memory used: {mem_after - mem_before} MB")
            print("Muscle alignment completed successfully.")

        except subprocess.CalledProcessError as e:
            return f"Error running Muscle: {e}"

    if algorithm == "mafft":
        command = "/Users/likitag/Downloads/mafft/core/mafft " + input_file + " > " + output_file
        os.system(command)

    _ , peakMemory = tracemalloc.get_traced_memory()
    


    execution_time = (time.time() - start_time) * 1000
    tracemalloc.stop()
    
    alignment = AlignIO.read(output_file, "fasta")
    sp_score = calculate_sp_score(alignment)


    return execution_time, peakMemory, sp_score

def accuracy_comparison(output_alignment, reference_alignment): 
    """Given the output alignment and the benchamrk reference alignment, will calculate an accuracy score."""
    total_accuracy = 0

    #convert each fasta file to a dictionary mapping 
    ref_dic = fasta_to_dict(reference_alignment)
    out_dic = fasta_to_dict(output_alignment)

    #calculate accuracy for each sequence, and add to total 
    for seq in out_dic: 
        seq1 = out_dic[seq]
        seq2 = ref_dic[seq]
        
        matches = sum(c1 == c2 for c1, c2 in zip(seq1, seq2))
        total = max(len(seq1), len(seq2))

        if total > 0: 
            total_accuracy+=(matches/total)

    #compute average accuracy across all sequences (should be between 0 and 1)
    return total_accuracy / len(out_dic)



def fasta_to_dict(file_path):
    """Convert a FASTA file to a dictionary where keys are sequence IDs and values are sequences."""
    sequences = {}
    current_id = None

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                current_id = line[1:].split()[0] 
                sequences[current_id] = ""
            elif current_id is not None:
                sequences[current_id] += line.strip()

    return sequences

def msf_to_fasta(input_file, output_file):
    """Convert a MSF file to a FASTA file, since the Balibase database has reference files in msf format, but we want to convert it to fasta format """
    
    alignment = AlignIO.read(input_file, "msf")
    AlignIO.write(alignment, output_file, "fasta")

    return

def main():
    working_directory = os.getcwd()

    muscle_times = []
    mafft_times = []
    clustalo_times = []
    kalign_times = []
    

    muscle_accuracies = []
    mafft_accuracies = []
    clustalo_accuracies = []
    kalign_accuracies = []

    muscle_memories = []
    mafft_memories = []
    clustalo_memories = []
    kalign_memories = []

    muscle_sp_scores = []
    mafft_sp_scores = []
    clustalo_sp_scores = []
    kalign_sp_scores = []

    for i in range(1, 10):
        msf_ref = os.path.join(working_directory, f"reference_files/test{i}_ref.msf")
        input_sequences = os.path.join(working_directory, f"input_files/test{i}_input.fasta")
        reference_alignment = os.path.join(working_directory, f"reference_files/test{i}_ref.fasta")
        msf_to_fasta(msf_ref, reference_alignment)

        #kalign alignment 
        output_alignment_kalign = os.path.join(working_directory, f"output_files/kalign/test{i}_output.fasta")
        time_kalign, memory_kalign, sp_score_kalign = run_msa("kalign", input_sequences, output_alignment_kalign)

        # print(f"KAlign execution time for test{i}: {time_kalign}")

        accuracy_score_kalign = accuracy_comparison(output_alignment_kalign, reference_alignment)

        # print(f"KAlign Accuracy for test{i}: {accuracy_score_kalign}")

        #clustal omega alignment 
        output_alignment_clustalo = os.path.join(working_directory, f"output_files/clustalo/test{i}_output.fasta")
        time_clustalo, memory_clustalo, sp_score_clustalo = run_msa("clustalo", input_sequences, output_alignment_clustalo)
        # print(f"Clustal Omega execution time for test{i}: {time_clustalo}")

        accuracy_score_clustalo = accuracy_comparison(output_alignment_clustalo, reference_alignment)
        # print(f"Clustal Omega Accuracy for test{i}: {accuracy_score_clustalo}")

        # Muscle Alignment
        output_alignment_muscle = os.path.join(working_directory, f"output_files/muscle/test{i}_output.fasta")
        time_muscle, memory_muscle, sp_score_muscle = run_msa("muscle", input_sequences, output_alignment_muscle)
        #print(f"MUSCLE execution time for test{i}: {time_muscle}")

        accuracy_score_muscle = accuracy_comparison(output_alignment_muscle, reference_alignment)
        #print(f"MUSCLE Accuracy for test{i}: {accuracy_score_muscle}")

        # MAFFT Alignment
        output_alignment_mafft = os.path.join(working_directory, f"output_files/mafft/test{i}_output.fasta")
        time_mafft, memory_mafft, sp_score_mafft = run_msa("mafft", input_sequences, output_alignment_mafft)
        #print(f"MAFFT execution time for test{i}: {time_mafft}")

        accuracy_score_mafft = accuracy_comparison(output_alignment_mafft, reference_alignment)
        #print(f"MAFFT Accuracy for test{i}: {accuracy_score_mafft}")

        muscle_times.append(time_muscle)
        mafft_times.append(time_mafft)
        clustalo_times.append(time_clustalo)
        kalign_times.append(time_kalign)

        muscle_accuracies.append(accuracy_score_muscle)
        mafft_accuracies.append(accuracy_score_mafft)
        clustalo_accuracies.append(accuracy_score_clustalo)
        kalign_accuracies.append(accuracy_score_kalign)

        muscle_memories.append(memory_muscle)
        mafft_memories.append(memory_mafft)
        clustalo_memories.append(memory_clustalo)
        kalign_memories.append(memory_kalign)

        muscle_sp_scores.append(sp_score_muscle)
        mafft_sp_scores.append(sp_score_mafft)
        clustalo_sp_scores.append(sp_score_clustalo)
        kalign_sp_scores.append(sp_score_kalign)

    print("MUSCLE times:")
    print(muscle_times)
    muscle_avg_time = sum(muscle_times)/len(muscle_times)
    print("MUSCLE accuracies:")
    print(muscle_accuracies)
    muscle_avg_accuracy = sum(muscle_accuracies)/len(muscle_accuracies)
    print("MUSCLE memories:")
    print(muscle_memories)
    muscle_avg_memory = sum(muscle_memories)/len(muscle_memories)
    print("MUSCLE SP scores:")
    print(muscle_sp_scores)
    muscle_avg_sp_score = sum(muscle_sp_scores)/len(muscle_sp_scores)

    print("MAFFT times:")
    print(mafft_times)
    mafft_avg_time = sum(mafft_times)/len(mafft_times)
    print("MAFFT accuracies:")
    print(mafft_accuracies)
    mafft_avg_accuracy = sum(mafft_accuracies)/len(mafft_accuracies)
    print("MAFFT memories:")
    print(mafft_memories)
    mafft_avg_memory = sum(mafft_memories)/len(mafft_memories)
    print("MAFFT SP scores:")
    print(mafft_sp_scores)
    mafft_avg_sp_score = sum(mafft_sp_scores)/len(mafft_sp_scores)

    print("Clustal Omega times:")
    print(clustalo_times)
    clustalo_avg_time = sum(clustalo_times)/len(clustalo_times)
    print("Clustal Omega accuracies:")
    print(clustalo_accuracies)
    clustalo_avg_accuracy = sum(clustalo_accuracies)/len(clustalo_accuracies)
    print("Clustal Omega memories:")
    print(clustalo_memories)
    clustalo_avg_memory = sum(clustalo_memories)/len(clustalo_memories)
    print("Clustal Omega SP scores:")
    print(clustalo_sp_scores)
    clustalo_avg_sp_score = sum(clustalo_sp_scores)/len(clustalo_sp_scores)

    print("kalign times:")
    print(kalign_times)
    kalign_avg_time = sum(kalign_times)/len(kalign_times)
    print("KAlign accuracies:")
    print(kalign_accuracies)
    kalign_avg_accuracy = sum(kalign_accuracies)/len(kalign_accuracies)
    print("KAlign memories:")
    print(kalign_memories)
    kalign_avg_memory = sum(kalign_memories)/len(kalign_memories)
    print("KAlign SP scores:")
    print(kalign_sp_scores)
    kalign_avg_sp_score = sum(kalign_sp_scores)/len(kalign_sp_scores)


    results_file = os.path.join(working_directory, "results.txt")

    # Open the file and write the data
    with open(results_file, 'w') as file:
        file.write("muscle avg time: " + str(muscle_avg_time) + "\n")
        file.write("muscle avg accuracy: " + str(muscle_avg_accuracy) + "\n")
        file.write("muscle avg memory: " + str(muscle_avg_memory) + "\n")
        file.write("muscle avg sp score: " + str(muscle_avg_sp_score) + "\n")

        file.write("mafft avg time: " + str(mafft_avg_time) + "\n")
        file.write("mafft avg accuracy: " + str(mafft_avg_accuracy) + "\n")
        file.write("mafft avg memory: " + str(mafft_avg_memory) + "\n")
        file.write("mafft avg sp score: " + str(mafft_avg_sp_score) + "\n")

        file.write("clustalo avg time: " + str(clustalo_avg_time) + "\n")
        file.write("clustalo avg accuracy: " + str(clustalo_avg_accuracy) + "\n")
        file.write("clustalo avg memory: " + str(clustalo_avg_memory) + "\n")
        file.write("clustalo avg sp score: " + str(clustalo_avg_sp_score) + "\n")

        file.write("kalign avg time: " + str(kalign_avg_time) + "\n")
        file.write("kalign avg accuracy: " + str(kalign_avg_accuracy) + "\n")
        file.write("kalign avg memory: " + str(kalign_avg_memory) + "\n")
        file.write("kalign avg sp score: " + str(kalign_avg_sp_score) + "\n")


    print("Results written to:", results_file)







if __name__ == '__main__':
    main()


