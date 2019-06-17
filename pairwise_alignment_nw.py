#!/usr/bin/env python
# Description = Computes Needleman-Wunsch Pairwise Alignment Algorithm

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def create_substution_matrix(residue_list, match_score=1,
                             mismatch_score=-1):
    '''
    This function creates a substituion matrix for residues:

    Arguments:
    - residue_list: A list of amino acid or dna/rna residues.
    - match_score: An integer number indicating match score. (default score = 1)
    - mismatch_score: An integer number indicating mismatch score. (default score = -1)
    '''
    scoring_matrix = pd.DataFrame(index=residue_list, columns=residue_list)
    scoring_matrix = scoring_matrix.fillna(0)
    for residue_col in residue_list:
        for residue_row in residue_list:
            if residue_col == residue_row:
                scoring_matrix.loc[residue_col, residue_row] = match_score
            else:
                scoring_matrix.loc[residue_col, residue_row] = mismatch_score
    print(scoring_matrix)
    print('\n')
    return scoring_matrix

def create_dynamic_prog_matrix(seq1, seq2, scoring_matrix, gap_penalty=-1):
    '''
    This function creates a dynamic matrix for two sequences:

    Arguments:
    - seq1: First amino acid sequence
    - seq2: Second amino acid sequence
    - scoring_matrix: Scoring matrix for amino acid
    - gap_penalty: An integer number indicating gap penalty/score. (default score = -1)
    '''
    index_list = [0]+list(seq1)
    column_list = [0]+list(seq2)
    dp_matrix = pd.DataFrame(index=index_list, columns=column_list)
    dp_matrix = dp_matrix.fillna(0)
    for i, residue_seq1 in enumerate(list(seq1)):
        dp_matrix.iloc[i+1, 0] = (i+1)*-1
        for j, residue_seq2 in enumerate(list(seq2)):
            dp_matrix.iloc[0, j+1] = (j+1)*-1
            dp_matrix.loc[residue_seq1,
                          residue_seq2] = scoring_matrix.loc[
                              residue_seq1, residue_seq2]
    scored_dp_matrix = _calculate_alignment_scores(dp_matrix, gap_penalty)
    return scored_dp_matrix

def _calculate_alignment_scores(dp_matrix, gap_penalty):
    for i, rows in enumerate(dp_matrix.index.values):
        if not rows == 0:
            for j, cols in enumerate(dp_matrix.columns.values):
                if  not cols == 0:
                    current_score = dp_matrix.iloc[i, j]
                    left_score = dp_matrix.iloc[i, j-1] + gap_penalty
                    up_score = dp_matrix.iloc[i-1, j] + gap_penalty
                    diag_score = dp_matrix.iloc[i-1, j-1] + current_score
                    high_score = max([left_score, up_score, diag_score])
                    dp_matrix.iloc[i, j] = high_score
    print(dp_matrix)
    print('\n')
    return(dp_matrix)

def trace_best_alignment(scored_dp_matrix, match_score=1,
                         mismatch_score=-1, gap_penalty=-1):
    '''
    This function traces back the best alignment.
    Diagonal arrow is a match or mismatch. Horizontal arrows introduce
    gap ("-") in the row and vertical arrows introduce gaps in the column.

    Arguments:
    - scored_dp_matrix: scored matrix for two sequences.
    - match_score: An integer number indicating match score. (default score = 1)
    - mismatch_score: An integer number indicating mismatch score. (default score = -1)
    - gap_penalty: An integer number indicating gap penalty/score. (default score = -1)
    '''
    i = len(scored_dp_matrix.index.values)-1
    j = len(scored_dp_matrix.columns.values)-1
    row_residue_list = []
    col_residue_list = []
    print("Trackback type:")
    while i > 0 and j > 0:
        current_score = scored_dp_matrix.iloc[i, j]
        left_score = scored_dp_matrix.iloc[i, j-1]
        up_score = scored_dp_matrix.iloc[i-1, j]
        diag_score = scored_dp_matrix.iloc[i-1, j-1]
        row_val = scored_dp_matrix.index.values[i]
        col_val = scored_dp_matrix.columns.values[j]
        trackback_type = ""  
        if i > 1 and j > 1 and (current_score == diag_score + match_score and row_val == col_val):
            trackback_type = "diagonal_match"
            row_val = scored_dp_matrix.index.values[i]
            col_val = scored_dp_matrix.columns.values[j]
            i -= 1
            j -= 1
        elif i > 1 and j > 1 and (current_score == diag_score + mismatch_score and row_val != col_val):
            trackback_type = "diagonal_mismatch"
            row_val = scored_dp_matrix.index.values[i]
            col_val = scored_dp_matrix.columns.values[j]
            i -= 1
            j -= 1
        elif i > 0 and (current_score == up_score + gap_penalty):
            trackback_type = "up"
            row_val = scored_dp_matrix.index.values[i]
            col_val = '-'
            i -= 1
        elif j > 0 and (current_score == left_score + gap_penalty):
            trackback_type = "left"
            col_val = scored_dp_matrix.columns.values[j]
            row_val = '-'
            j -= 1
        elif i == 1 and j == 1:
            row_val = scored_dp_matrix.index.values[i]
            col_val = scored_dp_matrix.columns.values[j]
            i -= 1
            j -= 1
        else:
            row_val = scored_dp_matrix.index.values[i]
            col_val = scored_dp_matrix.columns.values[j]
            i -= 1
            j -= 1
        print(trackback_type)
        row_residue_list.append(row_val)
        col_residue_list.append(col_val)
    col_seq = ''.join(map(str, col_residue_list[::-1]))
    row_seq = ''.join(map(str, row_residue_list[::-1]))
    return col_seq, row_seq


## Edited the previous function trace_best_alignment 
## by removing match and mismatch score

def trace_best_alignment_with_blosum(scored_dp_matrix, gap_penalty=-4):
    '''
    This function traces back the best alignment.
    Diagonal arrow is a match or mismatch. Horizontal arrows introduce
    gap ("-") in the row and vertical arrows introduce gaps in the column.

    Arguments:
    - scored_dp_matrix:   scored matrix for two sequences.
    - gap_penalty:        An integer number indicating gap penalty/score. 
                          (default gap penalty in NCBI BLOSUM62 = -4)
    '''
    i = len(scored_dp_matrix.index.values)-1
    j = len(scored_dp_matrix.columns.values)-1
    row_residue_list = []
    col_residue_list = []
    match_positions = []
    print("\nTrackback type:")
    while i > 0 and j > 0:
        current_score = scored_dp_matrix.iloc[i, j]
        left_score = scored_dp_matrix.iloc[i, j-1]
        up_score = scored_dp_matrix.iloc[i-1, j]
        diag_score = scored_dp_matrix.iloc[i-1, j-1]
        row_val = scored_dp_matrix.index.values[i]
        col_val = scored_dp_matrix.columns.values[j]
        trackback_type = ""  
        if i > 1 and j > 1 and current_score == diag_score:
            trackback_type = "diagonal_match"
            row_val = scored_dp_matrix.index.values[i]
            col_val = scored_dp_matrix.columns.values[j]
            i -= 1
            j -= 1
            match_positions.append(row_val)
        elif i > 0 and (current_score == up_score + gap_penalty):
            trackback_type = "up"
            row_val = scored_dp_matrix.index.values[i]
            col_val = '-'
            i -= 1
        elif j > 0 and (current_score == left_score + gap_penalty):
            trackback_type = "left"
            col_val = scored_dp_matrix.columns.values[j]
            row_val = '-'
            j -= 1
        else:
            trackback_type = "diagonal_match"
            row_val = scored_dp_matrix.index.values[i]
            col_val = scored_dp_matrix.columns.values[j]
            i -= 1
            j -= 1
            match_positions.append(row_val)
        print(trackback_type)
        row_residue_list.append(row_val)
        col_residue_list.append(col_val)
    print("\nTotal aligned positions: {}".format(len(match_positions)))
    col_seq = ''.join(map(str, col_residue_list[::-1]))
    row_seq = ''.join(map(str, row_residue_list[::-1]))
    return col_seq, row_seq

# Format BLOSUM subsitution file to pandas matrix

def read_blosum_file_to_matrix(blosum_file):
    '''
    Creates a matrix from NCBI BLOSUM file.

    Arguments:
    - blosum_file:     provide local file BLOSUM substitution matrix.
                       (Download the current version: 
                        https://www.ncbi.nlm.nih.gov/Class/BLAST/BLOSUM62.txt)
    '''
                
    header_list= ['A', 'R', 'N', 'D', 'C', 'Q',
                      'E', 'G', 'H', 'I', 'L', 'K', 
                      'M', 'F', 'P', 'S', 'T', 'W', 
                      'Y', 'V', 'B', 'Z', 'X', '*'] # * is gap
    blosum = pd.read_csv(blosum_file, skiprows=6, delim_whitespace=True)
    blosum = blosum.replace('NaN', 0)
    blosum.columns = header_list
    blosum.index = header_list
    return blosum

###################################

# Test pairwise alignment for Protein Sequences

human_p53 = 'TFSDLWKLLPENNV'
mouse_p53 = 'SQETFSGLWKLLPP'
scoring_matrix = create_substution_matrix(amino_acid_list)
scored_dp_matrix = create_dynamic_prog_matrix(human_p53, mouse_p53,
                                                scoring_matrix)
aligned_seq1, aligned_seq2 = trace_best_alignment(scored_dp_matrix)
print('Optimal gloabl alignment of the given sequences is:\n{}\n{}'.format(
    aligned_seq1, aligned_seq2))

# Test pairwise alignment with BLOSUM 62 Metrices
blosum = read_blosum_file_to_matrix('blosum62.txt')
scored_dp_matrix_blosum = create_dynamic_prog_matrix(
                            human_p53, mouse_p53, 
                            blosum, gap_penalty=-4)
print(scored_dp_matrix_blosum)
aligned_seq1, aligned_seq2 = trace_best_alignment_with_blosum(scored_dp_matrix)
print('Optimal gloabl alignment of the given sequences is:\n{}\n{}'.format(
    aligned_seq1, aligned_seq2))
create_heatmap_from_matrix(scored_dp_matrix_blosum)

# Test pairwise alignment for DNA Sequences

dna_list = ('A', 'C', 'T', 'G', 'U')
dna1 = 'GATTACA'
dna2 = 'GCATGCU'
# scoring_matrix = create_substution_matrix(dna_list)
# scored_dp_matrix = create_dynamic_prog_matrix(dna1, dna2,

amino_acid_list = ('A', 'C', 'D', 'E', 'F',
                   'G', 'H', 'I', 'K', 'L',
                   'M', 'N', 'P', 'Q', 'R',
                   'S', 'T', 'V', 'W', 'Y')
