#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Routines for measuring propensity of a sequence to slip into the -1 reading frame and minimize slip propensity over a defined subsequence
@author: Charles Kuntz :: cpkuntz@iu.edu
"""

import sys
import os
import re
import textwrap
import string
from itertools import groupby
import numpy as np
import pandas as pd


def domain_processor(domain_sequence_raw):
    domain_sequence = re.sub("[^a-zA-Z]", "", domain_sequence_raw)
    domain_codons = textwrap.wrap(domain_sequence, 3)
    domain_length = len(domain_sequence)
    return domain_sequence, domain_codons, domain_length


def parse_df(excel_workbook):
    pd_df = pd.read_excel(excel_workbook, sheet_name=0, header=0, index_col=0)
    wt_seq_codons = list(pd_df.columns.values)
    wt_seq_str = ''.join(wt_seq_codons)
    mut_codons = list(pd_df.index.values)
    sequences_codons = []
    wt_codon_idx = 0
    while wt_codon_idx < len(wt_seq_codons):
        for mut_codon in mut_codons:
            if wt_codon_idx == 0:
                codon1 = mut_codon
                codon2 = wt_seq_codons[1]
                codon3 = wt_seq_codons[2]
            if wt_codon_idx == 1:
                codon1 = wt_seq_codons[0]
                codon2 = mut_codon
                codon3 = wt_seq_codons[2]
            if wt_codon_idx == 2:
                codon1 = wt_seq_codons[0]
                codon2 = wt_seq_codons[1]
                codon3 = mut_codon
            seq_codons = [codon1, codon2, codon3]
            sequences_codons.append(seq_codons)
        wt_codon_idx += 1

    sequences_str = []
    for sequence_codons in sequences_codons:
        sequence_str = ''.join(sequence_codons)
        sequences_str.append(sequence_str)

    df_column = flatten_df_func(pd_df)

    flat_array = df_column.values
    flat_2dlist = flat_array.tolist()
    flat_list = []
    for items in flat_2dlist:
        item = items[0]
        flat_list.append(item)
    
    fs_efficiencies = flat_list

    return sequences_str, fs_efficiencies

def measure_slipsite_window(str_9mer, scoring_sys):
    window = textwrap.wrap(str_9mer, 3)
    codons, codons_fs, pairs, pairs_fs = frameshift_slipsite(window)
    fs_score = measure_slipsite(codons, codons_fs, pairs, pairs_fs, scoring_sys)
    return fs_score


def measure_single_heptamer(str_7mer, scoring_sys):
    # convert 7mer into 9mer
    str_9mer = 'GC'+str_7mer
    fs_score = measure_slipsite_window(str_9mer, scoring_sys)
    return fs_score


def measure_multiple_sequences(list_9mers, scoring_sys):
    scores = []
    for sequence in list_9mers:
        fs_score = measure_slipsite_window(sequence, scoring_sys)
        scores.append(fs_score)
    return scores

def flatten_df_func(pd_df):
    np_array = pd_df.values
    np_array_flat = np_array.flatten('F')
    np_array_long = np.transpose(np_array_flat)
    df_array_long = pd.DataFrame(np_array_long)
    return df_array_long

def genetic_code_func():
    genetic_code  = """
        TTT F      CTT L      ATT I      GTT V
        TTC F      CTC L      ATC I      GTC V
        TTA L      CTA L      ATA I      GTA V
        TTG L      CTG L      ATG M      GTG V
        TCT S      CCT P      ACT T      GCT A
        TCC S      CCC P      ACC T      GCC A
        TCA S      CCA P      ACA T      GCA A
        TCG S      CCG P      ACG T      GCG A
        TAT Y      CAT H      AAT N      GAT D
        TAC Y      CAC H      AAC N      GAC D
        TAA *      CAA Q      AAA K      GAA E
        TAG *      CAG Q      AAG K      GAG E
        TGT C      CGT R      AGT S      GGT G
        TGC C      CGC R      AGC S      GGC G
        TGA *      CGA R      AGA R      GGA G
        TGG W      CGG R      AGG R      GGG G
        """
    codon_finder = re.compile(r'[ATCG]{3}')
    amino_acid_finder = re.compile(r'\ \w{1}[\ |\n]|\*')
    codon_list = codon_finder.findall(genetic_code)
    amino_acid_list = [x.strip() for x in amino_acid_finder.findall(genetic_code)]
    genetic_code_dict = {}
    i = 0
    while i < len(codon_list):
        genetic_code_dict[codon_list[i]] = amino_acid_list[i]
        i += 1
    return genetic_code_dict


def reverse_complement_func(forwardSequence):
    forwardAlphabet = 'AGCTagct'
    revcompAlphabet = 'TCGAtcga'
    reverseSequence = forwardSequence[::-1]
    reverseComplementSequence = reverseSequence.translate({ord(x):
        y for (x, y) in zip(forwardAlphabet, revcompAlphabet)})
    return reverseComplementSequence


def ribosome(domain_codons, genetic_code_dict):
    domain_aminoacids = []
    for codon in domain_codons:
        amino_acid = genetic_code_dict[codon]
        domain_aminoacids.append(amino_acid)
    return domain_aminoacids


def get_codon_table():
    codon_table = {
    'A': ('GCT', 'GCC', 'GCA', 'GCG'),
    'C': ('TGT', 'TGC'),
    'D': ('GAT', 'GAC'),
    'E': ('GAA', 'GAG'),
    'F': ('TTT', 'TTC'),
    'G': ('GGT', 'GGC', 'GGA', 'GGG'),
    'I': ('ATT', 'ATC', 'ATA'),
    'H': ('CAT', 'CAC'),
    'K': ('AAA', 'AAG'),
    'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
    'M': ('ATG',),
    'N': ('AAT', 'AAC'),
    'P': ('CCT', 'CCC', 'CCA', 'CCG'),
    'Q': ('CAA', 'CAG'),
    'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
    'T': ('ACT', 'ACC', 'ACA', 'ACG'),
    'V': ('GTT', 'GTC', 'GTA', 'GTG'),
    'W': ('TGG',),
    'Y': ('TAT', 'TAC'),
    '*': ('TAA', 'TAG', 'TGA'),
    }
    return codon_table


def get_codon_freq_table():
    codon_freq_table = {
    'A': (['GCT', 0.26], ['GCC', 0.40], ['GCA', 0.23], ['GCG', 0.11]),
    'C': (['TGT', 0.45], ['TGC', 0.55]),
    'D': (['GAT', 0.46], ['GAC', 0.54]),
    'E': (['GAA', 0.42], ['GAG', 0.58]),
    'F': (['TTT', 0.45], ['TTC', 0.55]),
    'G': (['GGT', 0.16], ['GGC', 0.34], ['GGA', 0.25], ['GGG', 0.25]),
    'I': (['ATT', 0.36], ['ATC', 0.48], ['ATA', 0.16]),
    'H': (['CAT', 0.41], ['CAC', 0.59]),
    'K': (['AAA', 0.42], ['AAG', 0.58]),
    'L': (['TTA', 0.07], ['TTG', 0.13], ['CTT', 0.13], ['CTC', 0.20], ['CTA', 0.07], ['CTG', 0.41]),
    'M': (['ATG', 1.0],),
    'N': (['AAT', 0.46], ['AAC', 0.54]),
    'P': (['CCT', 0.28], ['CCC', 0.33], ['CCA', 0.27], ['CCG', 0.11]),
    'Q': (['CAA', 0.25], ['CAG', 0.75]),
    'R': (['CGT', 0.08], ['CGC', 0.19], ['CGA', 0.11], ['CGG', 0.21], ['AGA', 0.20], ['AGG', 0.20]),
    'S': (['TCT', 0.18], ['TCC', 0.22], ['TCA', 0.15], ['TCG', 0.06], ['AGT', 0.15], ['AGC', 0.24]),
    'T': (['ACT', 0.24], ['ACC', 0.36], ['ACA', 0.28], ['ACG', 0.12]),
    'V': (['GTT', 0.18], ['GTC', 0.24], ['GTA', 0.11], ['GTG', 0.47]),
    'W': (['TGG', 1.0],),
    'Y': (['TAT', 0.43], ['TAC', 0.57]),
    '*': (['TAA', 0.28], ['TAG', 0.20], ['TGA', 0.52]),
    }
    return codon_freq_table


def calc_avg_freq(domain_codons):
    codon_freq_table = get_codon_freq_table()
    codon_freqs_dict = {}
    for residue, codons in codon_freq_table.items():
        for codon, freq in codons:
            codon_freqs_dict[codon] = freq
    domain_freqs = []
    for position in domain_codons:
        domain_freqs.append(codon_freqs_dict[position])
    len_domain = len(domain_codons)
    domain_freqs_sum = sum(domain_freqs)
    domain_freqs_avg = domain_freqs_sum/len_domain
    return domain_freqs_avg


def frameshift_slipsite(window):
    p_site_codon = window[1]
    p_site_anticodon = reverse_complement_func(p_site_codon)
    p_site_rev = p_site_anticodon[::-1]
    a_site_codon = window[2]
    a_site_anticodon = reverse_complement_func(a_site_codon)
    a_site_rev = a_site_anticodon[::-1]
    window_bases = ''.join(window)
    frameshifted_window = textwrap.wrap(window_bases[2:-1], 3)
    p_site_fs = frameshifted_window[0]
    a_site_fs = frameshifted_window[1]
    pairs_fs = [''.join([p_site_rev[0], p_site_fs[0]]),
                ''.join([p_site_rev[1], p_site_fs[1]]),
                ''.join([p_site_rev[2], p_site_fs[2]]),
                ''.join([a_site_rev[0], a_site_fs[0]]),
                ''.join([a_site_rev[1], a_site_fs[1]]),
                ''.join([a_site_rev[2], a_site_fs[2]])]
    pairs =    [''.join([p_site_rev[0], p_site_codon[0]]),
                ''.join([p_site_rev[1], p_site_codon[1]]),
                ''.join([p_site_rev[2], p_site_codon[2]]),
                ''.join([a_site_rev[0], a_site_codon[0]]),
                ''.join([a_site_rev[1], a_site_codon[1]]),
                ''.join([a_site_rev[2], a_site_codon[2]])]
    codons = [[p_site_rev, p_site_codon], [a_site_rev, a_site_codon]]
    codons_fs = [[p_site_rev, p_site_fs], [a_site_rev, a_site_fs]]
    # p_site_triplet = [p_site_rev, p_site_fs]
    # a_site_triplet = [a_site_rev, a_site_fs]
    return codons, codons_fs, pairs, pairs_fs

def measure_slipsite(codons, codons_fs, pairs, pairs_fs, scoring_sys):
    if scoring_sys == 1:
        fs_score = measure_slipsite_simple(codons, codons_fs, pairs, pairs_fs)
    if scoring_sys == 2:
        fs_score = measure_slipsite_turner(codons, codons_fs, pairs, pairs_fs)
    return fs_score

def measure_slipsite_simple(codons, codons_fs, pairs, pairs_fs):
    simplest_scores = {'AA':3, 'GG':3, 'GA':3, 'AG':3,
                       'TT':2, 'CC':2, 'CT':2, 'TC':2,
                       'AT':-2, 'TA':-2, 'GC':-3, 'CG':-3,
                       'GT':-1, 'TG':-1, 'AC':1, 'CA':1}
    energies_zero = [simplest_scores[pair] for pair in pairs]
    energies_fs = [simplest_scores[pair] for pair in pairs_fs]
    fs_score = sum(energies_fs) - sum(energies_zero)
    return fs_score


def measure_slipsite_turner(codons, codons_fs, pairs, pairs_fs):
    import turner_energies as turner
    """
    p_site_codons = ['ABC', 'DEF']
    a_site_codons = ['UVW', 'XYZ']
    pairs = [('AB', 'DE'),
             ('BC', 'EF'),
             ('UV', 'XY'),
             ('VW', 'YZ')]
    """
    pairs_rna = [pair.replace('T', 'U') for pair in pairs]
    pairs_fs_rna = [pair.replace('T', 'U') for pair in pairs_fs]

    p_pairs = pairs_rna[0:3]
    p_pairs_fs = pairs_fs_rna[0:3]
    a_pairs = pairs_rna[3:]
    a_pairs_fs = pairs_fs_rna[3:]

    p_site_rev_raw = codons_fs[0][0]
    p_site_fs_raw = codons_fs[0][1]
    p_site_codon_raw = codons[0][1]
    a_site_rev_raw = codons_fs[1][0]
    a_site_fs_raw = codons_fs[1][1]
    a_site_codon_raw = codons[1][1]

    p_site_rev = p_site_rev_raw.replace('T', 'U')
    p_site_fs = p_site_fs_raw.replace('T', 'U')
    p_site_codon = p_site_codon_raw.replace('T', 'U')
    a_site_rev = a_site_rev_raw.replace('T', 'U')
    a_site_fs = a_site_fs_raw.replace('T', 'U')
    a_site_codon = a_site_codon_raw.replace('T', 'U')

    p_nn_pairs_fs = [(p_site_rev[0:2], p_site_fs[0:2]), # Watson-Crick positions P site
                     (p_site_rev[1:], p_site_fs[1:])]   # Wobble position P site
                   
    a_nn_pairs_fs = [(a_site_rev[0:2], a_site_fs[0:2]), # Watson-Crick positions A site
                     (a_site_rev[1:], a_site_fs[1:])]   # Wobble position A site

    p_nn_pairs =    [(p_site_rev[0:2], p_site_codon[0:2]), # Watson-Crick positions P site
                     (p_site_rev[1:], p_site_codon[1:])]   # Wobble position P site

    a_nn_pairs =    [(a_site_rev[0:2], a_site_codon[0:2]), # Watson-Crick positions A site
                     (a_site_rev[1:], a_site_codon[1:])]   # Wobble position A site

    energy_p = turner.calc_energy(p_pairs, p_nn_pairs)
    energy_a = turner.calc_energy(a_pairs, a_nn_pairs)
    energy_p_fs = turner.calc_energy(p_pairs_fs, p_nn_pairs_fs)
    energy_a_fs = turner.calc_energy(a_pairs_fs, a_nn_pairs_fs)

    energies_zero = energy_p + energy_a
    energies_fs = energy_p_fs + energy_a_fs
    fs_score = energies_fs - energies_zero
    return fs_score

def measure_slipsite_scores(domain_codons, scoring_sys):
    idx = 0
    scores = []
    windows = []
    while idx < len(domain_codons) - 3:
        window = domain_codons[idx:idx+3]
        windows.append(window)
        codons, codons_fs, pairs, pairs_fs = frameshift_slipsite(window)
        fs_score = measure_slipsite(codons, codons_fs, pairs, pairs_fs, scoring_sys)
        scores.append(fs_score)
        idx += 1
    return scores, windows


def get_GC(sequence):
    num_GC = 0
    seqlen = len(sequence)
    for base in sequence:
        if base == 'G' or base == 'C':
            num_GC += 1
    percent_GC = (num_GC/seqlen)*100
    return percent_GC


def maximize_friction_window(window_wt_codons, scoring_sys):
    codons_wt, codons_fs_wt, pairs_wt, pairs_fs_wt = frameshift_slipsite(window_wt_codons)
    window_friction_score_wt = measure_slipsite(codons_wt, codons_fs_wt, pairs_wt, pairs_fs_wt, scoring_sys)
    codon_dict = genetic_code_func()
    amino_acid_dict = get_codon_table()
    window_wt_amino_acids = ribosome(window_wt_codons, codon_dict)
    possible_codon_sequences = backtranslate_subsequence(window_wt_amino_acids, amino_acid_dict)
    window_friction_score_mut = window_friction_score_wt
    max_friction_window = window_wt_codons
    for candidate_window in possible_codon_sequences:
        codons, codons_fs, pairs, pairs_fs = frameshift_slipsite(candidate_window)
        candidate_friction_score = measure_slipsite(codons, codons_fs, pairs, pairs_fs, scoring_sys)
        if candidate_friction_score > window_friction_score_mut:
            window_friction_score_mut = candidate_friction_score
            max_friction_window = candidate_window
    delta_friction = window_friction_score_mut - window_friction_score_wt
    return max_friction_window, delta_friction


def backtranslate_subsequence(protein_seq, codon_table):
    sequences = [codon for codon in codon_table[protein_seq[0]]]
    for amino_acid in protein_seq[1:]:
        to_extend = sequences
        sequences = []
        for codon in codon_table[amino_acid]:
            for sequence in to_extend:
                sequence += codon
                sequences.append(sequence)
    seq_codons_list = []
    for sequence in sequences:
        seq_codon = textwrap.wrap(sequence, 3)
        seq_codons_list.append(seq_codon)
    return seq_codons_list


def get_candidate_windows(domain_codons, scoring_sys):
    wt_scores, wt_windows = measure_slipsite_scores(domain_codons, scoring_sys)
    overlapping_delta_friction = []
    overlapping_maxed_windows = []
    for idx in range(len(domain_codons)-3):
        window = domain_codons[idx:idx+3]
        max_window, delta_friction = maximize_friction_window(window, scoring_sys)
        overlapping_delta_friction.append(delta_friction)
        overlapping_maxed_windows.append(max_window)
    return overlapping_delta_friction, overlapping_maxed_windows, wt_scores, wt_windows


def maximize_friction(delta_friction_list, maxed_windows, wt_scores, wt_windows, domain_codons, cutoff_stdevs):
    friction_optimized_codons = []
    avg_wt_scores = np.mean(wt_scores)
    std_wt_scores = np.std(wt_scores)
    if cutoff_stdevs == 'SetManualCutoff':
        cutoff = float(input("You have opted to enter a manual floating point cutoff. Please enter now: "))
    else:
        cutoff = avg_wt_scores - (cutoff_stdevs * std_wt_scores)
    i = 0
    while i < len(wt_scores) - 5:
        if wt_scores[i] > cutoff:
            friction_optimized_codons.append(wt_windows[i][0])
            i += 1
        else:
            # 1. Set master candidate equal to candidate
            # 2. Check if idx+1 and idx+2 wt score are below  master_candidate score.
            # 3. If yes, one of these is the new master candidate. Index increment goes to appropriate position
            # 4. If no, we add codons from master_candidate to optimized sequence
            
            
            # Set index counter and master window
            #idx = 1
            candidate_window = maxed_windows[i]
            candidate_delta_friction = delta_friction_list[i]
            master_window = candidate_window
            master_delta_friction = candidate_delta_friction

            #Keep track of which window is the master with these variables
            master_window_idx0 = 1
            master_window_idx1 = 0
            master_window_idx2 = 0

            # Add first codon of master window to optimized sequence
            friction_optimized_codons.append(candidate_window[0])

            # Set variables for idx+1 and idx+2 windows
            plus1_candidate = maxed_windows[i+1]
            plus1_delta = delta_friction_list[i+1]
            plus2_candidate = maxed_windows[i+2]
            plus2_delta = delta_friction_list[i+2]
            
            # If plus 1 candidate is both less than cutoff and greater friction score than current master, it is new master
            if all([wt_scores[i+1] < cutoff, plus1_delta > master_delta_friction]):
                #idx += 1
                master_window = plus1_candidate
                master_delta_friction = plus1_delta
                master_window_idx0 = 0
                master_window_idx1 = 1
                master_window_idx2 = 0

            # If plus 2 candidate is both less than cutoff and greater friction score than current master, it is the new master
            if all([wt_scores[i+2] < cutoff, plus2_delta > master_delta_friction]):
                #idx += 2
                master_window = plus2_candidate
                master_delta_friction = plus2_delta
                master_window_idx0 = 0
                master_window_idx1 = 0
                master_window_idx2 = 1
            
            if master_window_idx0 == 1:
                friction_optimized_codons.append(candidate_window[1])
                friction_optimized_codons.append(candidate_window[2])
                i += 3

            if master_window_idx1 == 1:
                friction_optimized_codons.append(plus1_candidate[0])
                friction_optimized_codons.append(plus1_candidate[1])
                friction_optimized_codons.append(plus1_candidate[2])
                i += 4
            
            if master_window_idx2 == 1:
                # First position is second position of candidate_window, and first position of plus1_candidate. Which is the better codon to use?
                friction_optimized_codons.append(candidate_window[1])
                friction_optimized_codons.append(plus2_candidate[0])
                friction_optimized_codons.append(plus2_candidate[1])
                friction_optimized_codons.append(plus2_candidate[2])
                i += 5

    idx = len(wt_scores) - 5
    while idx < len(domain_codons):
        friction_optimized_codons.append(domain_codons[idx])
        idx += 1
    
    return friction_optimized_codons

def load_fasta_data(path_fasta):
    with open(path_fasta, "r") as fasta_handle:
        faiter = (x[1] for x in groupby(fasta_handle, lambda line: line[0] == ">"))
        for header in faiter:
            header_str = header.__next__()[0:].strip()
            seq = "".join(s.strip() for s in faiter.__next__())
            yield (header_str, seq)

def parse_fasta(path_fasta, ensembl_flag, short_transcript_filter):
    fasta_generator = load_fasta_data(path_fasta)
    fasta_list = []
    for header, sequence in fasta_generator:
        header_list = header.split()
        if ensembl_flag == True:
            ensembl_gene = header_list[3][5:]
            ensembl_transcript = header_list[0][1:]
            try:
                gene_name = header_list[6][12:]
            except IndexError:
                gene_name = ensembl_gene
            locus_coord_str = header_list[2][11:]
        else:
            gene_name = header_list[0][1:]
        sequence_str = ''.join(sequence)
        if ensembl_flag == True:
            output_line = [ensembl_gene, ensembl_transcript, gene_name, locus_coord_str, len(sequence_str), sequence_str]
        else:
            output_line = [gene_name, sequence_str]
        len_seq = len(sequence_str)
        # Basic checks on quality of transcript:
        pass_conditions = [(len_seq % 3 == 0),
                           (sequence_str[0:3] == 'ATG'),
                           (sequence_str[-3:] in ('TAG', 'TAA', 'TGA')),
                           ('N' not in sequence_str)]
        if all(pass_conditions):
            fasta_list.append(output_line)

    if ensembl_flag == True:
        stop_filtered_list = stop_codon_filter(fasta_list)

        if short_transcript_filter == True:
            output_fasta_list = fasta_transcript_filter(stop_filtered_list)
        else:
            output_fasta_list = stop_filtered_list

    else:
        output_fasta_list = fasta_list

    if ensembl_flag == True:
        ensembl_output_list = []
        for entry in output_fasta_list:
            [ensembl_gene, ensembl_transcript, gene_name, locus_coord_str, len_seq, sequence_str] = entry
            output_entry = [gene_name, sequence_str]
            ensembl_output_list.append(output_entry)
        final_output_fasta_list = ensembl_output_list
    else:
        final_output_fasta_list = output_fasta_list

    return final_output_fasta_list

def stop_codon_filter(cds_fasta_list):
    stop_filtered_list = []
    genetic_code_dict = genetic_code_func()
    for transcript_data in cds_fasta_list:
        sequence_raw = transcript_data[5]
        transcript_seq = re.sub("[^a-zA-Z]", "", sequence_raw)
        seq_codons = textwrap.wrap(transcript_seq, 3)
        seq_protein = ribosome(seq_codons, genetic_code_dict)
        seq_protein_str = ''.join(seq_protein[0:-1])
        if '*' not in seq_protein_str:
            stop_filtered_list.append(transcript_data)
        else:
            continue
    return stop_filtered_list


def fasta_transcript_filter(fasta_list):
    unique_genes = []
    filtered_list = []
    for entry in fasta_list:
        current_gene = entry[2]
        if current_gene not in unique_genes:
            unique_genes.append(current_gene)
    binned_list = []
    for idx in range(len(unique_genes)):
        binned_list.append([])
    for idx in range(len(unique_genes)):
        gene_name = unique_genes[idx]
        for entry in fasta_list:
            if entry[2] == gene_name:
                binned_list[idx].append(entry)
            else:
                continue
    for gene_bin in binned_list:
        len_longest_transcript = 0
        idx_longest_transcript = 0
        for idx in range(len(gene_bin)):
            len_transcript = gene_bin[idx][4]
            if len_transcript > len_longest_transcript:
                len_longest_transcript = len_transcript
                idx_longest_transcript = idx
        filtered_list.append(gene_bin[idx_longest_transcript])
    return filtered_list

def canonical_slipsite_counter(sequence_heptamers):
    slipsite_list = import_text_file('canonical_slipsites.txt')
    slipsite_count = 0
    for heptamer in sequence_heptamers:
        if heptamer in slipsite_list:
            slipsite_count += 1
    return slipsite_count


def sequence_friction_scan(seqname, input_sequence, scoring_sys):
    processed_sequence, codons_raw, seqlen = domain_processor(input_sequence)
    input_scores, input_windows = measure_slipsite_scores(codons_raw, scoring_sys)
    input_windows_bases = [''.join(window) for window in input_windows]
    input_heptamers = [window_bases[2:] for window_bases in input_windows_bases]
    canonical_slipsite_count = canonical_slipsite_counter(input_heptamers)
    len_seq = int(len(input_sequence)/3)
    percent_GC_wt = get_GC(input_sequence)
    avg_friction_wt = np.mean(input_scores)
    domain_freqs_avg = calc_avg_freq(codons_raw)
    output_list = [seqname+'|Length: '+str(len_seq)+'|XXX YYY Z slipsites: '+str(canonical_slipsite_count)+'|Avg fric: {:.3f}'.format(avg_friction_wt)+'|%GC: {:.1f}'.format(percent_GC_wt)+'|Avg codon freq: {:.3f}'.format(domain_freqs_avg), 'Friction scores', input_heptamers, input_scores]
    return output_list


def sequence_friction_scan_meta(seqname, input_sequence, scoring_sys, cutoff):
    processed_sequence, codons_raw, seqlen = domain_processor(input_sequence)
    input_scores, input_windows = measure_slipsite_scores(codons_raw, scoring_sys)
    input_windows_bases = [''.join(window) for window in input_windows]
    input_heptamers = [window_bases[2:] for window_bases in input_windows_bases]
    canonical_slipsite_count = canonical_slipsite_counter(input_heptamers)
    len_seq = int(len(input_sequence)/3)
    percent_GC_wt = get_GC(input_sequence)
    avg_friction_wt = np.mean(input_scores)
    num_slipsites = sum(x <= cutoff for x in input_scores)
    slipsites_per_hundred = 100 * (num_slipsites/len_seq)
    domain_freqs_avg = calc_avg_freq(codons_raw)
    scores_str_list = [str(x) for x in input_scores]
    scores_str = ','.join(scores_str_list)
    output_list = [str(seqname), str(len_seq), str(canonical_slipsite_count), str(avg_friction_wt), str(percent_GC_wt), str(domain_freqs_avg), str(num_slipsites), str(slipsites_per_hundred), scores_str]
    final_output_str = ','.join(output_list)
    return final_output_str

def print_scan_score_csv(output_list, jobname):
    score_csv = open(jobname+'_scores.csv', 'w')
    [output_tag, scorecolumn_tag, input_heptamers, input_scores] = output_list
    print(output_tag, file = score_csv)
    positions = [index+1 for index in range(len(input_scores))]
    print(','.join(['heptamer', 'position', 'friction_score']), file = score_csv)
    idx = 0
    while idx < len(positions):
        csv_line = [input_heptamers[idx], str(positions[idx]), str(input_scores[idx])]
        csv_line_str = ','.join(csv_line)
        print(csv_line_str, file = score_csv)
        idx += 1
    score_csv.close()
    
def print_scan_score_csv_meta(output_list, jobname):
    score_csv = open(jobname+'_metadata.csv', 'w')
    print(','.join(['gene name', 'length', 'XXXYYYZ slipsites', 'avg friction_score', 'percent GC', 'avg codon frequency', 'num slipsites', 'slipsites per hundred codons', 'sequence scores']), file = score_csv)
    for entry in output_list:
        print(entry, file = score_csv)
    score_csv.close()

def scan_multiple_sequences_in_fasta(path_fasta, scoring_sys, ensembl_flag, short_transcript_filter):
    fasta_list = parse_fasta(path_fasta, ensembl_flag, short_transcript_filter)
    for entry in fasta_list:
        gene_name = entry[0]
        sequence = entry[1]
        output_list = sequence_friction_scan(gene_name, sequence, scoring_sys)
        print_scan_score_csv(output_list, gene_name)


def scan_multiple_sequences_in_fasta_meta(path_fasta, scoring_sys, ensembl_flag, jobname, cutoff, short_transcript_filter):
    fasta_list = parse_fasta(path_fasta, ensembl_flag, short_transcript_filter)
    multiseq_data = []
    for entry in fasta_list:
        gene_name = entry[0]
        sequence = entry[1]
        output_list = sequence_friction_scan_meta(gene_name, sequence, scoring_sys, cutoff)
        multiseq_data.append(output_list)
    print_scan_score_csv_meta(multiseq_data, jobname)


def optimize_sequence(input_sequence, num_opts, start_in, end_in, cutoff_stdevs, scoring_sys):

    processed_sequence, codons_raw, seqlen = domain_processor(input_sequence)

    if start_in == 'StartDefault':
        start = 1
    else:
        start = start_in
    
    if end_in == 'EndDefault':
        end = len(codons_raw)
    else:
        end = end_in

    input_scores, input_windows = measure_slipsite_scores(codons_raw, scoring_sys)

    domain_codons = codons_raw[start:end]
    not_domain_fiveprime = codons_raw[0:start]
    not_domain_threeprime = codons_raw[end:]

    delta_friction_list_r1, maxed_windows_r1, wt_scores_r1, wt_windows_r1 = get_candidate_windows(domain_codons, scoring_sys)

    friction_optimized_codons_r1 = maximize_friction(delta_friction_list_r1, maxed_windows_r1, wt_scores_r1, wt_windows_r1, domain_codons, cutoff_stdevs)

    if num_opts == 1:
        friction_optimized_domain_codons_out = friction_optimized_codons_r1

    if num_opts == 2:

        delta_friction_list_r2, maxed_windows_r2, wt_scores_r2, wt_windows_r2 = get_candidate_windows(friction_optimized_codons_r1, scoring_sys)
        friction_optimized_codons_r2 = maximize_friction(delta_friction_list_r2, maxed_windows_r2, wt_scores_r2, wt_windows_r2, friction_optimized_codons_r1, cutoff_stdevs)

        friction_optimized_domain_codons_out = friction_optimized_codons_r2

    if num_opts == 3:

        delta_friction_list_r2, maxed_windows_r2, wt_scores_r2, wt_windows_r2 = get_candidate_windows(friction_optimized_codons_r1, scoring_sys)
        friction_optimized_codons_r2 = maximize_friction(delta_friction_list_r2, maxed_windows_r2, wt_scores_r2, wt_windows_r2, friction_optimized_codons_r1, cutoff_stdevs)

        delta_friction_list_r3, maxed_windows_r3, wt_scores_r3, wt_windows_r3 = get_candidate_windows(friction_optimized_codons_r2, scoring_sys)
        friction_optimized_codons_r3 = maximize_friction(delta_friction_list_r3, maxed_windows_r3, wt_scores_r3, wt_windows_r3, friction_optimized_codons_r2, cutoff_stdevs)

        friction_optimized_domain_codons_out = friction_optimized_codons_r3

    opt_codons_complete = []

    for codon in not_domain_fiveprime:
        opt_codons_complete.append(codon)

    for codon in friction_optimized_domain_codons_out:
        opt_codons_complete.append(codon)

    for codon in not_domain_threeprime:
        opt_codons_complete.append(codon)

    output_scores, output_windows = measure_slipsite_scores(opt_codons_complete, scoring_sys)

    output_seq = ''.join(opt_codons_complete)

    input_windows_bases = [''.join(window) for window in input_windows]
    input_heptamers = [window_bases[2:] for window_bases in input_windows_bases]
    output_windows_bases = [''.join(window) for window in output_windows]
    output_heptamers = [window_bases[2:] for window_bases in output_windows_bases]

    return input_scores, input_heptamers, output_scores, output_heptamers, output_seq


def prepare_output(input_seq, input_scores, output_seq, output_scores, seqname):

    len_seq = int(len(input_seq)/3)
    percent_GC_wt = get_GC(input_seq)
    percent_GC_opt = get_GC(output_seq)
    avg_friction_wt = np.mean(input_scores)
    avg_friction_opt = np.mean(output_scores)

    output_wt = [seqname+'_WT', 'Length: '+str(len_seq), 'Avg: {:.3f}'.format(avg_friction_wt), '%GC: {:.1f}'.format(percent_GC_wt), input_seq]
    output_opt = [seqname+'_opt', 'Length: '+str(len_seq), 'Avg: {:.3f}'.format(avg_friction_opt), '%GC: {:.1f}'.format(percent_GC_opt), output_seq]

    output_list = [output_wt, output_opt]

    return output_list


def print_fasta(filtered_list, outname):
     myfasta = open(outname+'.fasta', 'w')
     for record in filtered_list:
         print('>'+record[0]+'|'+record[1]+'|'+record[2]+'|'+record[3], file = myfasta)
         idx = 0
         columns = 60
         while idx < len(record[-1]):
             if len(record[-1][idx:]) > columns:
                 print(record[-1][idx:idx+columns], file = myfasta)
             else:
                 print(record[-1][idx:], file = myfasta)
             idx += columns
     myfasta.close()

def print_score_csv(input_scores, input_heptamers, output_scores, output_heptamers, jobname):
    score_csv = open(jobname+'_scores.csv', 'w')
    positions = [index+1 for index in range(len(input_scores))]
    print('seq_position,wt_heptamer,wt_score,opt_heptamer,opt_score', file = score_csv)
    idx = 0
    while idx < len(positions):
        csv_line = [str(positions[idx]), input_heptamers[idx], str(input_scores[idx]), output_heptamers[idx], str(output_scores[idx])]
        csv_line_str = ','.join(csv_line)
        print(csv_line_str, file = score_csv)
        idx += 1
    score_csv.close()

def print_dms_score_csv(list_seqs, fs_dms, scores, jobname):
    score_csv = open(jobname+'.csv', 'w')
    print('sequence,fs_efficiency,friction_score', file = score_csv)
    idx = 0
    while idx < len(list_seqs):
        csv_line = [str(list_seqs[idx]), str(fs_dms[idx]), str(scores[idx])]
        csv_line_str = ','.join(csv_line)
        print(csv_line_str, file = score_csv)
        idx += 1
    score_csv.close()

def import_fasta_seq(path_to_file):
    file_obj = open(path_to_file)
    data_raw = []
    for line in file_obj:
        data_raw.append(line.split())
    file_obj.close()
    idx = 1
    wt_sequence = ''
    while idx < len(data_raw):
        wt_seq_line = data_raw[idx][0].upper()
        wt_sequence += wt_seq_line
        idx += 1
    return wt_sequence

def import_text_file(path_to_file):
    file_obj = open(path_to_file)
    data_raw = []
    for line in file_obj:
        data_raw.append(line.split())
    file_obj.close()
    seqs_str_list = []
    for line in data_raw:
        item = line[0]
        seqs_str_list.append(item)
    return seqs_str_list


def convert_7N_to_9N(seq_7N_list):
    seq_9N_list = []
    for seq_7N in seq_7N_list:
        seq_9N = 'GC'+seq_7N
        seq_9N_list.append(seq_9N)
    return seq_9N_list

def print_7N_score_csv(list_seqs, scores, jobname):
    score_csv = open(jobname+'.csv', 'w')
    print('sequence,friction_score', file = score_csv)
    idx = 0
    while idx < len(list_seqs):
        csv_line = [str(list_seqs[idx]), str(scores[idx])]
        csv_line_str = ','.join(csv_line)
        print(csv_line_str, file = score_csv)
        idx += 1
    score_csv.close()

def interactive_configuration():
    while True:
        print("Which analysis do you want to perform?")
        print("[1] Measure friction of a sequence and get friction-minimized sequence")
        print("[2] Upload table of frameshifting ratios measured by DMS and calculate friction of each sequence")
        print("[3] Upload text file containing list of heptamers and get friction scores for each sequence")
        print("[4] Perform a series of friction scans on a fasta file containing gene sequences")
        print("[5] Perform a series of friction scans and only save metadata for each sequence")
        analysis_choice = int(input("Please enter [1], [2], [3], [4], or [5]: "))
        if analysis_choice in (1, 2, 3, 4, 5):
            break
        else:
            print("There has been an error. Please re-enter.")
    if analysis_choice == 1:
        path_fasta = str(input("Please enter a path to the FASTA file containing the sequence to be optimized: "))
        path_text = 'null'
        seqname = str(input("Please enter the name of the sequence being optimized: "))
        try:
            start = int(input("Please enter the number of the first codon in the range\n of sequence you want optimized, counting from 1\nOr hit [Enter] for beginning: ")) - 1
        except ValueError:
            start = 'StartDefault'
        try:
            end = int(input("Please enter the number of the last codon in the range\n of the sequence you want optimized, counting from 1\nOr hit [Enter] for end: ")) - 1
        except ValueError:
            end = 'EndDefault'
        while True:
            num_iterations = int(input("Please enter the number of iterations (1-3) you want to run. More iterations will optimize more codons: "))
            if num_iterations in (1,2,3):
                break
            else:
                print("Number of iterations must be 1-3. Please re-enter.")
        try:
            cutoff = float(input("Please enter a cutoff in standard deviations below the mean friction score to use as a criterion for optimization of a codon (2 works well): "))
        except ValueError:
            cutoff = 'SetManualCutoff'
        jobname = str(input("Please enter a name for this job: "))
        path_excel = 'null'
    if analysis_choice == 2:
        path_fasta = 'null'
        seqname = 'null'
        start = 0
        end = 10
        num_iterations = 0
        cutoff = 1.0
        path_text = 'null'
        path_excel = str(input("Please enter a path to the Excel sheet containing frameshifting efficiencies from DMS: "))
        jobname = str(input("Please enter a name for this job: "))
    if analysis_choice == 3:
        path_fasta = 'null'
        seqname = 'null'
        start = 0
        end = 10
        num_iterations = 0
        cutoff = 1.0
        path_excel = 'null'
        path_text = str(input("Please enter a path to a text file containing a list of heptamers: "))
        jobname = str(input("Please enter a name for this job: "))
    if analysis_choice == 4:
        path_fasta = str(input("Please enter a path to the FASTA file containing the sequences to be scanned: "))
        seqname = 'null'
        start = 0
        end = 10
        num_iterations = 0
        cutoff = 1.0
        path_excel = 'null'
        path_text = 'null'
        jobname = 'null'
    if analysis_choice == 5:
        path_fasta = str(input("Please enter a path to the FASTA file containing the sequences to be scanned: "))
        seqname = 'null'
        start = 0
        end = 10
        num_iterations = 0
        cutoff = float(input("Please enter a cutoff friction score to use as a criterion for defining and counting a slip site (2.14 works well): "))
        path_excel = 'null'
        path_text = 'null'
        jobname = str(input("Please enter a name for this job: "))
    while True:
        print("Which scoring system do you want to use")
        print("[1] Use number of hydrogen-bonds formed and broken (simple scoring)")
        print("[2] Use Turner nearest neighbor parameters")
        scoring_sys = int(input("Please enter [1] or [2]: "))
        if scoring_sys in (1, 2):
            break
        else:
            print("There has been an error. Please re-enter.")
    return analysis_choice, path_fasta, seqname, start, end, num_iterations, cutoff, jobname, path_excel, path_text, scoring_sys


if __name__ == '__main__':

    analysis_choice, path_fasta, seqname, start, end, num_opts, cutoff, jobname, path_excel, path_text, scoring_sys = interactive_configuration()

    if analysis_choice == 1:
        wt_seq = import_fasta_seq(path_fasta)
        input_scores, input_heptamers, output_scores, output_heptamers, output_seq = optimize_sequence(wt_seq, num_opts, start, end, cutoff, scoring_sys)
        output_list = prepare_output(wt_seq, input_scores, output_seq, output_scores, seqname)
        print_score_csv(input_scores, input_heptamers, output_scores, output_heptamers, jobname)
        print_fasta(output_list, jobname)

    if analysis_choice == 2:
        sequences_str, fs_efficiencies = parse_df(path_excel)
        fric_scores = measure_multiple_sequences(sequences_str, scoring_sys)
        print_dms_score_csv(sequences_str, fs_efficiencies, fric_scores, jobname)

    if analysis_choice == 3:
        seqs_7N_list = import_text_file(path_text)
        seqs_9N_list = convert_7N_to_9N(seqs_7N_list)
        fric_scores = measure_multiple_sequences(seqs_9N_list, scoring_sys)
        print_7N_score_csv(seqs_7N_list, fric_scores, jobname)


    if analysis_choice == 4:
        while True:
            print("Does the input fasta have Ensembl-format headers?")
            print("[1] Yes")
            print("[2] No, the headers just contain the gene names right after the \'>\'")
            ensembl_flag_inp = int(input("Please enter [1] or [2]: "))
            if ensembl_flag_inp in (1, 2):
                break
            else:
                print("There has been an error. Please re-enter.")
        if ensembl_flag_inp == 1:
            ensembl_flag = True
        else:
            ensembl_flag = False

        while True:
            print("Do you want to filter transcripts shorter than the longest transcript for a given gene in the fasta?")
            print("[1] Yes")
            print("[2] No")
            filter_flag_inp = int(input("Please enter [1] or [2]: "))
            if filter_flag_inp in (1, 2):
                break
            else:
                print("There has been an error. Please re-enter.")
        if filter_flag_inp == 1:
            short_transcript_filter = True
        else:
            short_transcript_filter = False

        scan_multiple_sequences_in_fasta(path_fasta, scoring_sys, ensembl_flag, short_transcript_filter)

    if analysis_choice == 5:
        while True:
            print("Does the input fasta have Ensembl-format headers?")
            print("[1] Yes")
            print("[2] No, the headers just contain the gene names right after the \'>\'")
            ensembl_flag_inp = int(input("Please enter [1] or [2]: "))
            if ensembl_flag_inp in (1, 2):
                break
            else:
                print("There has been an error. Please re-enter.")
        if ensembl_flag_inp == 1:
            ensembl_flag = True
        else:
            ensembl_flag = False

        while True:
            print("Do you want to filter transcripts shorter than the longest transcript for a given gene in the fasta?")
            print("[1] Yes")
            print("[2] No")
            filter_flag_inp = int(input("Please enter [1] or [2]: "))
            if filter_flag_inp in (1, 2):
                break
            else:
                print("There has been an error. Please re-enter.")
        if filter_flag_inp == 1:
            short_transcript_filter = True
        else:
            short_transcript_filter = False

        scan_multiple_sequences_in_fasta_meta(path_fasta, scoring_sys, ensembl_flag, jobname, cutoff, short_transcript_filter)
