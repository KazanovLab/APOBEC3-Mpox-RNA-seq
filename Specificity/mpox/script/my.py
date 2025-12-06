import pandas as pd
from scipy.stats import pearsonr
import numpy as np
import math
import os
import sys
import argparse


parser = argparse.ArgumentParser(prog="'script_for_MK.py'", description="Python script for specificity matrix construction and comparison.")
parser.add_argument('-i', '--input', help="Enter the whole path to input TXT file with list of positions", type=str)
parser.add_argument('-g', '--genome', help="Enter the whole path to genome FASTA file.", type=str)
parser.add_argument('-m', '--matrix', help='Enter the whole path to literature matrix TXT file.', type=str)
#parser.add_argument('-lf', '--len_frame', help="Choose the length of the frame: 11 (-5 to +5), 7 (-3 to 3), 5 (-2 to 2), 3 (-1 to 1).", type=str)
#parser.add_argument('-o', '--output', help="Enter the whole path to output TXT file with user matrix.", type=str)

args = parser.parse_args()
input_file = args.input
genome_file = args.genome
matrix_file = args.matrix
#frame = args.len_frame
#output_file = args.output

# Get list of positions #
df = pd.read_csv(input_file, dtype={'pos':int})
list_of_pos = df['pos'].tolist()

# Get sequence of genome #
with open(genome_file, 'r') as file:
    data = file.readlines()[1:]
    genome_seq = ''.join([s.strip() for s in data])

# Define frame #
frame = ['-5', '-4', '-3', '-2', '-1', '0', '+1', '+2', '+3', '+4', '+5']
'''
if frame == '11':
    frame = ['-5', '-4', '-3', '-2', '-1', '0', '+1', '+2', '+3', '+4', '+5']
elif frame == '7':
    frame = ['-3', '-2', '-1', '0', '+1', '+2', '+3']
elif frame == '5':
    frame = ['-2', '-1', '0', '+1', '+2']
elif frame == '3':
    frame = ['-1', '0', '1']
else:
    print('Enter the length of the frame correctly!')
    sys.exit()
'''
### STEP 1 - USER MATRIX CONSTRUCTION ###
# Function to create user matrix #
def create_matrix(genome, list_pos, frame):
    def get_reverse_complement(seq):
        new_seq = ''
        for n in seq:
            if n == 'A':
                new_seq += 'T'
            elif n == 'C':
                new_seq += 'G'
            elif n == 'G':
                new_seq += 'C'
            else:
                new_seq += 'A'
        return new_seq[::-1]
    
    
    # Matrix construction #
    data = pd.DataFrame(index=['A', 'C', 'G', 'T'], columns=frame, data=0)
    len_frame = len(frame)
    for p in list_pos:
        subseq = genome[(p-(math.floor(len_frame / 2)+1)):(p + math.floor(len_frame / 2))]
        mutation_pos = genome[p - 1]
        if mutation_pos in ['A', 'T']:
            print(f"Position {p} is {mutation_pos} nucleotide in genome! It will not be taken into account when constructing the matrix.")
            continue
        elif mutation_pos == 'G':
            subseq = get_reverse_complement(subseq)
        
        for s, i in zip(subseq, frame):
            data.loc[s, i] += 1

    # Normalize matrix #
    print(data)
    for i in frame:
        data[i] = round(data[i] / data[i].sum(), 3)
    
    return data

user_matrix = create_matrix(genome_seq, list_of_pos, frame)
print(user_matrix)

### STEP 2 - MATRIX COMPARISON ###
literature_matrix = pd.read_csv(matrix_file, sep='\t', index_col=0)
print(literature_matrix)
def compare_two_matrix_by_Pearson(m1, m2):
    m1 = m1[m2.columns]
    vm1 = m1.values.flatten(order='F') # to vector by column to column 
    vm2 = m2.values.flatten(order='F') # to vector by column to column
    print(vm1)
    print(vm2)
    comparison_result = pearsonr(vm1, vm2)
    
    return comparison_result

#print(f"Answer:\n{compare_two_matrix_by_Pearson(user_matrix, literature_matrix)}")

#flist = ["S1A.M11_A3A_Gordenin.txt","S1D.M11_A3A_GordeninTaylor.txt","S1B.M11_A3B_Gordenin.txt","S1E.M11_A3B_GordeninTaylor.txt","hA3B_found1.processed.txt","hA3F_found1.processed.txt","hA3F_found3.processed.txt","hA3F_found4.processed.txt","hA3G_found1.processed.txt","hA3G_found4.processed.txt"]

#for i in range(len(flist)):
#    for j in range(len(flist)):
#        if i > j:
#            continue
#        f1 = flist[i]
#        f2 = flist[j]
#        path1 = "/Users/mar/BIO/PROJECTS/MPOX/Specificity/MK/literature_matrices/"+f1
#        path2 = "/Users/mar/BIO/PROJECTS/MPOX/Specificity/MK/literature_matrices/"+f2
#        m1 = pd.read_csv(path1, sep='\t', index_col=0)
#        m2 = pd.read_csv(path2, sep='\t', index_col=0)
#        print(f1)
#        print(f2)
#        print(f"Answer:\n{compare_two_matrix_by_Pearson(m1, m2)}")
    
flist = ["S1A.M11_A3A_Gordenin.txt","S1D.M11_A3A_GordeninTaylor.txt","S1B.M11_A3B_Gordenin.txt","S1E.M11_A3B_GordeninTaylor.txt","hA3B_found1.processed.txt","hA3F_found1.processed.txt","hA3G_found1.processed.txt","hA3G_found4.processed.txt"]

for i in range(len(flist)):
    f1 = flist[i]
    path1 = "/Users/mar/BIO/PROJECTS/MPOX/Specificity/MK/literature_matrices/"+f1
    m1 = pd.read_csv(path1, sep='\t', index_col=0)
    print(f1)
    print(f"Answer:\n{compare_two_matrix_by_Pearson(user_matrix, m1)}")
    








