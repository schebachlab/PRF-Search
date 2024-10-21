#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to search transcriptome databases (Ensembl CDS) for Harrington motifs
"""
import RNA
import re
import textwrap
import csv
import numpy as np
import pandas as pd
import random
import sys
import pickle
import gc
from collections import Counter
import multiprocessing as mp
from scipy import stats
from scipy.stats import mannwhitneyu, binomtest, ks_2samp, kstwobign, zscore, chisquare
from scipy.stats import norm
from itertools import groupby
import matplotlib.pyplot as plt
import frameshift_routines_v10 as frameshift
import time
import configparser
from Bio import Align
from Bio.Align import substitution_matrices
import difflib


def parse_blastp_output(blastp_results_path):
    """
    Parses the blastp output (TSV) and returns a dictionary where keys are Ensembl transcript ID and values are
    a list of matches, each with Uniprot Accession ID and similarity score.
    """
    blastp_matches = {}
    with open(blastp_results_path, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            ensembl_transcript_id = row[0]
            uniprot_accession_id = row[1].split('|')[1]  # Extract the Uniprot ID from the full string
            similarity_score = float(row[2])

            # Add the match to the dictionary
            if ensembl_transcript_id not in blastp_matches:
                blastp_matches[ensembl_transcript_id] = []
            blastp_matches[ensembl_transcript_id].append((uniprot_accession_id, similarity_score))
    
    return blastp_matches


def bin_transcripts_by_gene_id(filtered_list):
    """
    Bins Ensembl transcripts by their shared Ensembl Gene ID.
    
    Parameters:
    filtered_list (list): List of parsed and filtered Ensembl transcripts.
    
    Returns:
    dict: Dictionary where keys are Ensembl Gene IDs and values are lists of transcripts associated with that Gene ID.
    """
    binned_dict = {}
    
    # Iterate through the filtered_list and bin by ensembl_gene_id (which is the first item in each entry)
    for entry in filtered_list:
        ensembl_gene_id = entry[0]  # Assuming the Ensembl Gene ID is at index 0
        
        if ensembl_gene_id not in binned_dict:
            binned_dict[ensembl_gene_id] = []
        
        # Append the transcript entry to the list of transcripts for this gene
        binned_dict[ensembl_gene_id].append(entry)
    
    return binned_dict


def match_ensembl_to_uniprot(blastp_results, ensembl_gene_id_dict):
    """
    Matches Ensembl transcripts to Uniprot sequences using the parsed blastp output.
    
    Arguments:
    - filtered_list: List of Ensembl transcripts, labeled as canonical or alternative.
    - blastp_results: Dictionary of blastp results, mapping Ensembl transcript IDs to Uniprot matches.
    - ensembl_gene_id_dict: Dictionary mapping Ensembl Gene IDs to their corresponding transcripts.
    
    Returns:
    - ensembl_to_uniprot: A dictionary mapping Ensembl transcript IDs to Uniprot Accession IDs and similarity scores.
    """
    ensembl_to_uniprot = {}

    for ensembl_gene_id, transcripts in ensembl_gene_id_dict.items():
        # Keep track of the best Uniprot match (highest similarity score)
        best_match = None
        best_similarity_score = -1
        best_uniprot_id = None
        
        # First pass: identify the best Uniprot match for the canonical transcript
        for transcript in transcripts:
            transcript_id = transcript[1]  # Ensembl transcript ID

            if transcript_id in blastp_results:
                for uniprot_id, similarity_score in blastp_results[transcript_id]:
                    if similarity_score > best_similarity_score:
                        best_similarity_score = similarity_score
                        best_uniprot_id = uniprot_id
                        best_match = transcript  # Store the transcript with the best match
        
        # If a match is found, update all transcripts in this gene with the Uniprot ID of the canonical transcript
        if best_uniprot_id:
            for transcript in transcripts:
                transcript_id = transcript[1]  # Ensembl transcript ID
                # Check if transcript has an individual match
                transcript_matches = blastp_results.get(transcript_id)
                if transcript_matches:
                    for uniprot_id, similarity_score in transcript_matches:
                        if uniprot_id == best_uniprot_id:
                            ensembl_to_uniprot[transcript_id] = [uniprot_id, similarity_score]
                            break
                    else:
                        ensembl_to_uniprot[transcript_id] = [best_uniprot_id, "No match to this alternative transcript"]
                else:
                    # No match found for this transcript
                    ensembl_to_uniprot[transcript_id] = [best_uniprot_id, "No match to this alternative transcript"]
        else:
            # No match found for any transcript in this gene
            for transcript in transcripts:
                transcript_id = transcript[1]
                ensembl_to_uniprot[transcript_id] = ["No Uniprot Accession ID found", "No similarity score with any Uniprot sequence found"]

    return ensembl_to_uniprot



def process_transcripts_with_uniprot(filtered_list, blastp_results_path, tmd_data, outname):
    """
    Processes Ensembl transcripts to associate them with Uniprot sequences based on blastp results.

    Parameters:
    - filtered_list (list): The list of parsed and filtered Ensembl transcripts.
    - blastp_results_path (str): Path to the TSV file with blastp results.

    Returns:
    dict: A dictionary mapping Ensembl transcript IDs to Uniprot metadata, including Uniprot accession and similarity score.
    """
    # Step 1: Parse the blastp results
    blastp_results = parse_blastp_output(blastp_results_path)

    # Step 2a: Filter transcripts that lack a TMD
    # Skip this for now - will filter by KW-0472 later after harrington motif list has been made
    #tmd_filtered_list = filter_transcripts_with_tmd(filtered_list, tmd_data)

    # Step 2b: Bin transcripts by Ensembl Gene ID
    ensembl_gene_id_dict = bin_transcripts_by_gene_id(filtered_list)

    # Step 3: Match Ensembl transcripts to Uniprot sequences
    ensembl_uniprot_dict = match_ensembl_to_uniprot(blastp_results, ensembl_gene_id_dict)

    output_file = outname+'_ensembl_uniprot_dict_dataframe.xlsx'

    # Step 4: Write the results to an Excel file
    output_ensembl_to_uniprot_excel(ensembl_uniprot_dict, filtered_list, output_file)

    return ensembl_uniprot_dict


def output_ensembl_to_uniprot_excel(ensembl_uniprot_dict, filtered_list, output_file):
    """
    Outputs the mapping between Ensembl transcripts and Uniprot sequences to an Excel file.

    Parameters:
    - ensembl_uniprot_dict (dict): A dictionary mapping Ensembl transcript IDs to Uniprot accession IDs and similarity scores.
    - filtered_list (list): The list of parsed and filtered Ensembl transcripts.
    - output_file (str): Path to save the output Excel file.
    """
    data = []
    
    for entry in filtered_list:
        ensembl_gene_id = entry[0]  # Assuming Ensembl Gene ID is at index 0
        ensembl_transcript_id = entry[1]  # Assuming Ensembl Transcript ID is at index 1
        canonical_status = entry[-1]  # Assuming canonical status is at the last index

        # Retrieve the Uniprot data
        uniprot_data = ensembl_uniprot_dict.get(ensembl_transcript_id, ["No Uniprot Accession ID found", "No similarity score found"])
        uniprot_accession_id = uniprot_data[0]
        similarity_score = uniprot_data[1]

        # Append the row to the data list
        data.append([ensembl_transcript_id, ensembl_gene_id, canonical_status, uniprot_accession_id, similarity_score])

    # Convert to a DataFrame
    df = pd.DataFrame(data, columns=['Ensembl Transcript ID', 'Ensembl Gene ID', 'Canonical Status', 'Uniprot Accession ID', 'Similarity Score'])

    # Output the DataFrame to Excel
    df.to_excel(output_file, index=False)
    print(f"Ensembl to Uniprot mapping saved to {output_file}")


def append_uniprot_accession_to_motif(harrington_motifs, ensembl_uniprot_dict):
    """
    Appends Uniprot accession and similarity score to the Harrington motif data.

    Parameters:
    harrington_motifs (list): List of Harrington motif entries, where each entry has the Ensembl transcript ID as the first element.
    ensembl_uniprot_dict (dict): Dictionary mapping Ensembl transcript IDs to a list containing Uniprot accession and similarity score.

    Returns:
    list: Updated Harrington motifs list with Uniprot accessions appended to each entry.
    """
    updated_harrington_motifs = []

    for motif_entry in harrington_motifs:
        ensembl_transcript_id = motif_entry[0]  # Assuming the first column is Ensembl transcript ID
        uniprot_info = ensembl_uniprot_dict.get(ensembl_transcript_id)

        if uniprot_info:
            uniprot_id = uniprot_info[0] if len(uniprot_info) > 0 else 'No Uniprot sequence identified'
            similarity_score = uniprot_info[1] if len(uniprot_info) > 1 else 'N/A'
        else:
            uniprot_id = 'No Uniprot sequence identified'
            similarity_score = 'N/A'

        # Append Uniprot accession and similarity score to the motif entry
        updated_entry = motif_entry + [uniprot_id, similarity_score]
        updated_harrington_motifs.append(updated_entry)

    return updated_harrington_motifs


def genetic_code_func():
    genetic_code = """
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


def get_GC(sequence):
    num_GC = 0
    seqlen = len(sequence)
    for base in sequence:
        if base == 'G' or base == 'C':
            num_GC += 1
    percent_GC = (num_GC / seqlen) * 100
    return percent_GC


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
    domain_freqs_avg = domain_freqs_sum / len_domain
    return domain_freqs_avg


def load_fasta_data(path_fasta):
    with open(path_fasta, "r") as file_handle:
        faiter = (x[1] for x in groupby(file_handle, lambda line: line[0] == ">"))
        for header in faiter:
            header_str = header.__next__()[0:].strip()
            seq = ''.join(s.strip() for s in faiter.__next__())
            yield (header_str, seq)

def parse_fasta(path_fasta):
    fasta_generator = load_fasta_data(path_fasta)
    fasta_list = []
    for header, sequence in fasta_generator:
        header_list = header.split()
        ensembl_gene = header_list[3][5:]
        ensembl_transcript = header_list[0][1:]
        try:
            gene_name = header_list[6][12:]
        except IndexError:
            gene_name = ensembl_gene
        locus_coord_str = header_list[2][11:]
        sequence_str = ''.join(sequence)
        output_line = [ensembl_gene, ensembl_transcript, gene_name, locus_coord_str, len(sequence_str), sequence_str]
        len_seq = len(sequence_str)
        # Basic checks on quality of transcript:
        pass_conditions = [(len_seq % 3 == 0),
                           (sequence_str[0:3] == 'ATG'),
                           (sequence_str[-3:] in ('TAG', 'TAA', 'TGA')),
                           ('N' not in sequence_str)]
        if all(pass_conditions):
            fasta_list.append(output_line)
    return fasta_list

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

def fasta_transcript_filter_canonical_only(fasta_list):
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


def fasta_transcript_filter_first_version(fasta_list):
    unique_genes = []
    filtered_list_longest = list([])
    filtered_list_transcript_id_only = []
    filtered_list_shorter = list([])
    for entry in fasta_list:
        current_gene = entry[0]
        if current_gene not in unique_genes:
            unique_genes.append(current_gene)
    binned_list = []
    for idx in range(len(unique_genes)):
        binned_list.append([])
    for idx in range(len(unique_genes)):
        gene_name = unique_genes[idx]
        for entry in fasta_list:
            if entry[0] == gene_name:
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
        longest_gene_bin = list(gene_bin[idx_longest_transcript])
        filtered_list_transcript_id_only.append(longest_gene_bin[1])
        filtered_list_longest.append(longest_gene_bin)
    for gene_bin in binned_list:
        for idx in range(len(gene_bin)):
            shorter_gene_bin = list(gene_bin[idx])
            if gene_bin[idx][1] not in filtered_list_transcript_id_only:
                filtered_list_shorter.append(shorter_gene_bin)

    canonical_count = len(filtered_list_longest)
    alt_spliced_count = len(filtered_list_shorter)
    filtered_list = []
    for entry in filtered_list_longest:
        entry.append('Canonical')
        filtered_list.append(entry)
    
    for entry in filtered_list_shorter:
        entry.append('Alternative')
        filtered_list.append(entry)

    #return filtered_list, canonical_count, alt_spliced_count, filtered_list_longest, filtered_list_shorter
    return filtered_list, canonical_count, alt_spliced_count


def fasta_transcript_filter(fasta_list):
    unique_genes = []
    filtered_list = []

    canonical_count = 0
    alt_spliced_count = 0

    # Identify unique genes
    for entry in fasta_list:
        current_gene = entry[0]
        if current_gene not in unique_genes:
            unique_genes.append(current_gene)

    # Create bins for transcripts belonging to the same gene
    binned_list = [[] for _ in range(len(unique_genes))]
    for idx, gene_name in enumerate(unique_genes):
        for entry in fasta_list:
            if entry[0] == gene_name:
                binned_list[idx].append(entry)

    # Label transcripts as canonical or non-canonical
    for gene_bin in binned_list:
        # Find the longest transcript in each bin
        len_longest_transcript = max(gene_bin, key=lambda x: x[4])[4]
        
        for entry in gene_bin:
            # Label as canonical or non-canonical
            if entry[4] == len_longest_transcript:
                entry.append('Canonical')
                canonical_count += 1
            else:
                entry.append('Alternative')
                alt_spliced_count += 1
            
            filtered_list.append(entry)

    return filtered_list, canonical_count, alt_spliced_count



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
    file_obj.close()
    return seqs_str_list

def import_raw(path_to_file):
    file_obj = open(path_to_file)
    data_raw = []
    for line in file_obj:
        data_raw.append(line.split())
    file_obj.close()
    return data_raw


# Existing function to read friction data
def parse_friction_csv(path_friction_csv):
    friction_data = {}
    with open(path_friction_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            heptamer = row['sequence']
            friction_score = float(row['friction_score'])
            friction_data[heptamer] = friction_score
    return friction_data

# New function to predict secondary structure
def predict_secondary_structure(sequence):
    rna_sequence = sequence.replace('T', 'U')
    structure, mfe = RNA.fold(rna_sequence)
    return structure, mfe


def detect_pseudoknot(structure):
    return any(char in structure for char in ['[', ']', '{', '}', '<', '>'])


def detect_stem_loop(structure, min_stem_length=3, max_loop_size=8):
    stack = []
    stems = []

    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                start = stack.pop()
                stems.append((start, i))

    # Filter stems by length and check for loops
    for stem in stems:
        stem_length = stem[1] - stem[0] + 1
        if stem_length >= 2 * min_stem_length:  # Each stem has at least min_stem_length pairs
            loop_size = min(structure[stem[0]+1:stem[1]].count('.'), max_loop_size)
            if loop_size > 0:
                return True  # Found a valid stem-loop
    return False


def filter_transcripts_with_tmd(stop_filtered_list, tmd_data):
    """
    Filters stop_filtered_list to only include transcripts that have TMD data in tmd_data.

    Parameters:
    stop_filtered_list (list): List of Ensembl transcripts that passed the stop codon filter.
    tmd_data (dict): Dictionary containing TMD information, with Ensembl transcript IDs as keys.

    Returns:
    list: Filtered list of Ensembl transcripts that have corresponding TMD data.
    """
    filtered_list_with_tmd = []

    for transcript in stop_filtered_list:
        ensembl_transcript_id = transcript[1]  # Assuming the second item in each entry is the Ensembl transcript ID
        
        # Only include the transcript if it has corresponding TMD data in tmd_data
        if ensembl_transcript_id in tmd_data:
            filtered_list_with_tmd.append(transcript)
    
    return filtered_list_with_tmd


def get_canonical_count(filtered_list, tmd_data, harrington_motifs, outname):

    # Total number of transcripts
    total_canonical_count = 0
    total_alt_spliced_count = 0
    total_count = len(filtered_list)

    # Total number of transripts with TMDs
    canonical_count_with_tmds = 0
    alt_spliced_count_with_tmds = 0
    total_count_with_tmds = 0

    # Number of TMDs in topcons data
    unique_tmd_entries = []
    transcripts_with_tmds = []
    for transcript_id, topcons_info in tmd_data.items():
        if transcript_id not in transcripts_with_tmds:
            transcripts_with_tmds.append(transcript_id)

        for idx, tmd_info in enumerate(topcons_info):
            tmd_seq = tmd_info['adjusted_tmd_seq']  # Extract tmd_type
            tm_start = tmd_info['adjusted_tmd_start_position']
            tm_end = tmd_info['adjusted_tmd_end_position']
            tmd_entry = [transcript_id, tmd_seq, tm_start, tm_end]
            unique_tmd_entries.append(tmd_entry)


    num_unique_tmds = len(unique_tmd_entries)
    num_transcripts_with_tmds = len(transcripts_with_tmds)

    # Total number of unique genes in CDS database
    unique_gene_ids = []

    # Total number of unique genes in CDS database with TMDs
    unique_gene_ids_with_tmds = []

    # Total number of unique HGNC symbols in CDS database
    unique_hgnc_symbols = []

    # Total number of unique HGNC symbols in CDS database with TMDs
    unique_hgnc_symbols_with_tmds = []

    unique_transcripts_with_tmds = []

    for entry in filtered_list:
        gene_id = entry[0]
        hgnc_symbol = entry[2]
        transcript_id = entry[1]
        canon_state = entry[-1]
        if gene_id not in unique_gene_ids:
            unique_gene_ids.append(gene_id)

        if hgnc_symbol not in unique_hgnc_symbols:
            unique_hgnc_symbols.append(hgnc_symbol)

        if canon_state == 'Canonical':
            total_canonical_count += 1
        else:
            total_alt_spliced_count += 1

        if transcript_id in tmd_data:
            unique_transcripts_with_tmds.append(transcript_id)
            total_count_with_tmds += 1

            if gene_id not in unique_gene_ids_with_tmds:
                unique_gene_ids_with_tmds.append(gene_id)
            
            if hgnc_symbol not in unique_hgnc_symbols_with_tmds:
                unique_hgnc_symbols_with_tmds.append(hgnc_symbol)

            if canon_state == 'Canonical':
                canonical_count_with_tmds += 1
            else:
                alt_spliced_count_with_tmds += 1
    
    tmd_entries_in_filtered_data = []
    for tmd_entry in unique_tmd_entries:
        transcript_id = tmd_entry[0]
        if transcript_id in unique_transcripts_with_tmds:
            tmd_entries_in_filtered_data.append(tmd_entry)

    num_tmd_entries_in_filtered_data = len(tmd_entries_in_filtered_data)

    # Total number of unique genes
    num_unique_genes = len(unique_gene_ids)

    # Total number of unique genes with TMDs
    num_unique_genes_with_tmds = len(unique_gene_ids_with_tmds)

    # Total number of unique HGNC symbols
    num_unique_hgnc_symbols = len(unique_hgnc_symbols)

    # Total number of unique HGNC symbols with TMDs
    num_unique_hgnc_symbols_with_tmds = len(unique_hgnc_symbols_with_tmds)

    # Number of Harrington motifs identified
    num_harrington_motifs = len(harrington_motifs)

    # Number of Harrington motifs in alternative isoforms
    alt_motifs = 0

    # Number of Harrington motifs in canonical isoforms
    canonical_motifs = 0

    for entry in harrington_motifs:
        canon_status = entry[-2]
        if canon_status == 'Alternative':
            alt_motifs += 1
        if canon_status == 'Canonical':
            canonical_motifs += 1

    logfile = open(outname+'_logfile.txt', 'w')
    
    print(f"Total number of transcripts in Ensembl CDS database: {total_count}", file=logfile)
    print(f"Number of canonical transcripts: {total_canonical_count}", file=logfile)
    print(f"Number of alternate isoform transcripts: {total_alt_spliced_count}", file=logfile)
    print(f"Total number of transcripts in Ensembl with TMDs (filtered): {total_count_with_tmds}", file=logfile)
    print(f"Number of canonical transcripts with TMDs: {canonical_count_with_tmds}", file=logfile)
    print(f"Number of alternate isoform transcripts with TMDs: {alt_spliced_count_with_tmds}", file=logfile)
    print(f"Number of TMDs predicted by TOPCONS2 in Ensembl CDS database: {num_unique_tmds}", file=logfile)
    print(f"Number of TMDs predicted by TOPCONS2 in filtered Ensembl data: {num_tmd_entries_in_filtered_data}", file=logfile)
    print(f"Separate count of transcripts in TOPCONS data: {num_transcripts_with_tmds}", file=logfile)
    print(f"Number of unique Ensembl gene IDs encountered in Ensembl CDS database: {num_unique_genes}", file=logfile)
    print(f"Number of unique Ensembl gene IDs with TMDs encountered in Ensembl CDS database: {num_unique_genes_with_tmds}", file=logfile)
    print(f"Number of unique HGNC symbols encountered in Ensembl CDS database: {num_unique_hgnc_symbols}", file=logfile)
    print(f"Number of unique HGNC symbols with TMDs encountered in Ensembl CDS database: {num_unique_hgnc_symbols_with_tmds}", file=logfile)
    print(f"Number of Harrington motifs identified: {num_harrington_motifs}", file=logfile)
    print(f"Number of Harrington motifs identified in canonical transcripts: {canonical_motifs}", file=logfile)
    print(f"Number of Harrington motifs identified in alternative isoform transcripts: {alt_motifs}", file=logfile)

    logfile.close()

    output_list = [total_canonical_count,
                   total_alt_spliced_count,
                   total_count,
                   canonical_count_with_tmds,
                   alt_spliced_count_with_tmds,
                   total_count_with_tmds,
                   num_unique_tmds,
                   num_tmd_entries_in_filtered_data,
                   num_transcripts_with_tmds,
                   num_unique_genes,
                   num_unique_genes_with_tmds,
                   num_unique_hgnc_symbols,
                   num_unique_hgnc_symbols_with_tmds]

    return output_list



def parse_tmd_csv(path_tmd_csv):
    tmd_data = {}
    with open(path_tmd_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            transcript_id = row['seq_name']
            tmd_data.setdefault(transcript_id, []).append({
                'adjusted_tmd_seq': row['adjusted_tmd_seq'],
                'len_adjusted_tmd': int(row['len_adjusted_tmd']),
                'adjusted_tmd_start_position': int(row['adjusted_tmd_start_position']),
                'adjusted_tmd_end_position': int(row['adjusted_tmd_end_position']),
                'adjusted_tmd_delta_G': float(row['adjusted_tmd_delta_G']),
                'adjusted_tmd_upstream_loop_length': int(row['adjusted_tmd_upstream_loop_length']),
                'adjusted_tmd_downstream_loop_length': int(row['adjusted_tmd_downstream_loop_length']),
                'tmd_type': row['tmd_type']  # Add tmd_type to the data structure
            })
    return tmd_data




def motif_search(filtered_list, tmd_data, slipsite_list, friction_data, gap_boundary_near, gap_boundary_far, scoring_sys):
    genetic_code_dict = genetic_code_func()
    codon_table = get_codon_table()
    codon_freq_table = get_codon_freq_table()

    set_24_canonical_slipsites = {'AAAAAAA',
                                  'AAAAAAT',
                                  'AAAAAAC',
                                  'AAATTTA',
                                  'AAATTTT',
                                  'AAATTTC',
                                  'TTTAAAA',
                                  'TTTAAAT',
                                  'TTTAAAC',
                                  'TTTTTTA',
                                  'TTTTTTT',
                                  'TTTTTTC',
                                  'CCCAAAA',
                                  'CCCAAAT',
                                  'CCCAAAC',
                                  'CCCTTTA',
                                  'CCCTTTT',
                                  'CCCTTTC',
                                  'GGGAAAA',
                                  'GGGAAAT',
                                  'GGGAAAC',
                                  'GGGTTTA',
                                  'GGGTTTT',
                                  'GGGTTTC'}

    transcript_codons_dict = {}
    transcript_protein_dict = {}
    transcript_header_dict = {}

    passing_transcripts = []
    motifs_within_gap = 0  # Count motifs within the user-defined gap boundaries
    num_trials_all = 0 # Number of attempts to search for slippery sequence 45 codons downstream from TMDs for use in binomial test
    num_trials_canonical = 0
    num_trials_alternative = 0

    # Populate dictionaries for transcript codons, proteins, and headers
    for transcript_data in filtered_list:
        ensembl_transcript_name = transcript_data[1]
        passing_transcripts.append(ensembl_transcript_name)
        ensembl_gene_name = transcript_data[0]
        common_gene_name = transcript_data[2]
        sequence_raw = transcript_data[5]
        transcript_seq = re.sub("[^a-zA-Z]", "", sequence_raw)
        seq_codons = textwrap.wrap(transcript_seq, 3)
        seq_protein = ribosome(seq_codons, genetic_code_dict)
        locus_coord_str = transcript_data[3].split(':')
        canonical_status = transcript_data[6]  # Extract the canonical or non-canonical status
        try:
            locus_coord_start, locus_coord_end = int(locus_coord_str[2]), int(locus_coord_str[3])
        except ValueError:
            locus_coord_start, locus_coord_end = int(locus_coord_str[3]), int(locus_coord_str[4])
        header_data = [ensembl_gene_name, common_gene_name, locus_coord_start, locus_coord_end, transcript_seq, canonical_status]
        
        for codon in seq_codons:
            transcript_codons_dict.setdefault(ensembl_transcript_name, []).append(codon)
        for residue in seq_protein:
            transcript_protein_dict.setdefault(ensembl_transcript_name, []).append(residue)
        for datum in header_data:
            transcript_header_dict.setdefault(ensembl_transcript_name, []).append(datum)

    harrington_motif_indexes = []

    # Main loop to process each transcript and its TMDs
    for transcript_id in passing_transcripts:
        if transcript_id not in tmd_data:
            continue
        tmd_info_list = tmd_data[transcript_id]  # List of TMDs for the transcript
        [ensembl_gene_name, common_gene_name, locus_coord_start, locus_coord_end, transcript_seq, canonical_status] = transcript_header_dict[transcript_id]
        transcript_codons = transcript_codons_dict[transcript_id]
        transcript_residues = transcript_protein_dict[transcript_id]
        transcript_residues_str = ''.join(transcript_residues)

        for idx, tmd_info in enumerate(tmd_info_list):
            tmd_type = tmd_info['tmd_type']  # Extract tmd_type
            tm_start = tmd_info['adjusted_tmd_start_position']
            tm_end = tmd_info['adjusted_tmd_end_position']
            tm_dg = tmd_info['adjusted_tmd_delta_G']
            tm_protein = transcript_residues[tm_start-1:tm_end]
            tm_domain_str = ''.join(tm_protein)
            tm_codons = transcript_codons[tm_start-1:tm_end]
            tm_codons_str = ''.join(tm_codons)
            candidates_in_window = []


            # Check to see if there is heptamer sequence 45 codons downstream from end of TMD
            if len(transcript_codons) >= tm_end + 45 + 2:
                num_trials_all += 1 # This provides num_trials for doing the binomial test for successful slip site trial at 45 codons distance from TMD
                if canonical_status == 'Canonical':
                    num_trials_canonical += 1
                if canonical_status == 'Alternative':
                    num_trials_alternative += 1

            # Set up the search window
            i = max(0, tm_end + gap_boundary_near)  # Ensure i is within bounds
            j = min(len(transcript_codons), tm_end + gap_boundary_far)  # Ensure j is within bounds

            # Check if there is another TMD starting within the search window
            if idx < len(tmd_info_list) - 1:
                next_tm_end = tmd_info_list[idx + 1]['adjusted_tmd_end_position'] + gap_boundary_near
                if i < next_tm_end <= j:
                    j = next_tm_end  # Adjust the window end to the next TMD

            ss_search_codons = transcript_codons[i:j]
            ss_search_idx = 0
            while ss_search_idx < len(ss_search_codons) - 2:
                candidate_codons = ss_search_codons[ss_search_idx:ss_search_idx+3]
                candidate_slipsite_long = ''.join(candidate_codons)
                candidate_slipsite = candidate_slipsite_long[2:]
                
                if candidate_slipsite in slipsite_list:
                    location_ideality_abs = abs((tm_end + 45) - (ss_search_idx + i + 1))
                    location_ideality = (tm_end + 45) - (ss_search_idx + i + 1)
                    ss_coordinate_nuc = (ss_search_idx + 1 + i) * 3 + 2
                    ss_coordinate_codon = i + ss_search_idx + 1
                    gap_size = i + ss_search_idx + 1 - tm_end
                    chr_coordinate = locus_coord_start + ss_coordinate_nuc
                    gap_codons = transcript_codons[tm_end:ss_coordinate_codon-1]
                    gap_nuc_raw = ''.join(transcript_codons[tm_end:ss_coordinate_codon])
                    gap_nuc = gap_nuc_raw[:-1] if gap_nuc_raw else ''  # Ensure it's not empty
                    gap_gc_content = get_GC(gap_nuc)
                    gap_avg_codon_freq = calc_avg_freq(gap_codons)

                    # Extract downstream sequence for frameshifted peptide
                    downstream_seq_codons = transcript_codons[ss_coordinate_codon-1:]
                    downstream_seq_bases = ''.join(downstream_seq_codons)
                    
                    if len(downstream_seq_bases) >= 3:  # Ensure there's enough sequence to process
                        downstream_codons_fs = textwrap.wrap(downstream_seq_bases[2:-1], 3) # Exclude the last nucleotide
                        downstream_peptide_fs_full = ribosome(downstream_codons_fs, genetic_code_dict)
                        fs_peptide_list = []
                        for residue in downstream_peptide_fs_full:
                            if residue != '*':
                                fs_peptide_list.append(residue)
                            else:
                                fs_peptide_list.append(residue)
                                break
                        fs_peptide = ''.join(fs_peptide_list)
                        terminated = 'yes' if fs_peptide and fs_peptide[-1] == '*' else 'no'
                    else:
                        fs_peptide = ''
                        terminated = 'no'

                    # Calculate friction score
                    if candidate_slipsite in friction_data:
                        friction_score = friction_data[candidate_slipsite]
                    else:
                        friction_score = frameshift.measure_slipsite_window(candidate_slipsite_long, scoring_sys)

                    if candidate_slipsite in set_24_canonical_slipsites:
                        slipsite_strict = 'Strict'
                    else:
                        slipsite_strict = 'Loose'

                    # Predict mRNA secondary structure in a defined window
                    downstream_window_start = (ss_coordinate_codon * 3) + 5
                    downstream_window_end = min(len(transcript_codons), downstream_window_start + 60)
                    downstream_window = ''.join(transcript_codons[downstream_window_start:downstream_window_end])
                    structure, mfe = predict_secondary_structure(downstream_window)
                    stem_loop_found = detect_stem_loop(structure)
                    pseudoknot_found = detect_pseudoknot(structure)

                    candidate_entry = [transcript_id,
                                       ensembl_gene_name,
                                       common_gene_name,
                                       tm_start,
                                       tm_end,
                                       tm_dg,
                                       tm_domain_str,
                                       tm_codons_str,
                                       ss_coordinate_nuc,
                                       ss_coordinate_codon,
                                       chr_coordinate,
                                       candidate_slipsite,
                                       gap_size,
                                       location_ideality,
                                       gap_nuc,
                                       gap_gc_content,
                                       gap_avg_codon_freq,
                                       fs_peptide,
                                       terminated,
                                       transcript_seq,
                                       transcript_residues_str,
                                       friction_score,
                                       slipsite_strict,
                                       "Yes" if stem_loop_found else "No",  # Stem-loop presence
                                       "Yes" if pseudoknot_found else "No",  # Pseudoknot presence
                                       structure,
                                       canonical_status,
                                       tmd_type]  # Add tmd_type to the output

                    candidates_in_window.append(candidate_entry)
                
                    # Use gap_boundary_near and gap_boundary_far for counting motifs within the user-defined window
                    if gap_boundary_near <= gap_size <= gap_boundary_far:
                        motifs_within_gap += 1
                
                ss_search_idx += 1

            # Evaluate and select the best candidate if multiple are found
            if len(candidates_in_window) == 0:
                slipsite_str = 'None'
            elif len(candidates_in_window) == 1:
                slipsite_str = candidates_in_window[0][0]
                winning_idx = 0
            else:
                slipsite_str = candidates_in_window[0][0]
                ideality = candidates_in_window[0][11]  # location_ideality index in candidate_entry
                winning_idx = 0
                idx = 1
                while idx < len(candidates_in_window):
                    if candidates_in_window[idx][11] < ideality:
                        slipsite_str = candidates_in_window[idx][0]
                        ideality = candidates_in_window[idx][11]
                        winning_idx = idx
                    idx += 1

            if slipsite_str != 'None':
                found_motif = candidates_in_window[winning_idx]
                harrington_motif_indexes.append(found_motif)

    return harrington_motif_indexes, motifs_within_gap, transcript_codons_dict, num_trials_all, num_trials_canonical, num_trials_alternative





def plot_slipsite_tmd_distances_with_categories(distance_data, output_file):
    """
    Plots the distances between slippery sequences and their nearest upstream TMDs,
    categorized by sequence type and slip site type (strict vs loose).
    
    Parameters:
        distance_data (dict): Categorized distance data returned from analyze_slipsite_to_tmd_distances.
        output_file (str): Path to save the histogram plot.
    """
    plt.figure(figsize=(10, 6))

    # Plot histograms for each category (strict and loose, all, canonical, alternative)
    plt.hist(distance_data['strict_all'], bins=50, color='blue', alpha=0.7, label='Strict - All')
    plt.hist(distance_data['loose_all'], bins=50, color='skyblue', alpha=0.7, label='Loose - All')
    plt.hist(distance_data['strict_canonical'], bins=50, color='green', alpha=0.5, label='Strict - Canonical')
    plt.hist(distance_data['loose_canonical'], bins=50, color='lightgreen', alpha=0.5, label='Loose - Canonical')
    plt.hist(distance_data['strict_alternative'], bins=50, color='orange', alpha=0.5, label='Strict - Alternative')
    plt.hist(distance_data['loose_alternative'], bins=50, color='yellow', alpha=0.5, label='Loose - Alternative')

    plt.xlabel('Distance from nearest upstream TMD (codons)', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title('Distribution of Distances between Slippery Sequences and Nearest Upstream TMD', fontsize=16)
    plt.legend()
    plt.grid(axis='y', alpha=0.75)
    plt.savefig(output_file)
    plt.close()
    print(f"Histogram saved to {output_file}")


def incremental_histogram(background_files, bin_edges):
    """
    Builds a histogram incrementally from background data stored in multiple files.
    
    Parameters:
        background_files (list): List of file paths containing the background data.
        bin_edges (array-like): Edges of the histogram bins.
    
    Returns:
        hist (numpy array): Counts for each bin.
    """
    hist = np.zeros(len(bin_edges) - 1, dtype=np.float64)
    
    for file in background_files:
        with open(file, 'rb') as f:
            data = pickle.load(f)
            counts, _ = np.histogram(data, bins=bin_edges)
            hist += counts
            del data
            gc.collect()
    
    return hist

def mann_whitney_u_test_incremental(actual_distances, background_files, lower_range=-10, upper_range=10):
    """
    Performs the Mann-Whitney U test on the subset of distances between lower_range and upper_range codons,
    processing background data incrementally without sampling.

    Parameters:
        actual_distances (list): List of distances from actual data.
        background_files (list): List of file paths containing the background data.
        lower_range (int): Lower bound of the range (inclusive).
        upper_range (int): Upper bound of the range (inclusive).

    Returns:
        U statistic and p-value from the Mann-Whitney U test.
    """
    import gc
    from scipy.stats import mannwhitneyu

    # Extract the distances in the desired range for actual data
    actual_over_range = [d for d in actual_distances if lower_range <= d <= upper_range]

    if len(actual_over_range) == 0:
        print("No actual distances in the specified range for Mann-Whitney U Test.")
        return None, None

    # Initialize a list to store background distances within the range
    background_over_range = []

    # Process background files incrementally
    for file in background_files:
        with open(file, 'rb') as f:
            data = pickle.load(f)
            # Filter data to desired range
            filtered_data = [d for d in data if lower_range <= d <= upper_range]
            background_over_range.extend(filtered_data)
            del data
            del filtered_data
            gc.collect()

    if len(background_over_range) == 0:
        print("No background distances in the specified range for Mann-Whitney U Test.")
        return None, None

    # Perform Mann-Whitney U test
    U_statistic, p_value = mannwhitneyu(actual_over_range, background_over_range, alternative='two-sided')

    return U_statistic, p_value

def mann_whitney_u_test_incremental_with_sampling(actual_distances, background_files, lower_range=-10, upper_range=10, max_background_samples=100000):
    """
    Performs the Mann-Whitney U test on the subset of distances between lower_range and upper_range codons,
    processing background data incrementally.

    Parameters:
        actual_distances (list): List of distances from actual data.
        background_files (list): List of file paths containing the background data.
        lower_range (int): Lower bound of the range (inclusive).
        upper_range (int): Upper bound of the range (inclusive).
        max_background_samples (int): Maximum number of background samples to use.

    Returns:
        U statistic and p-value from the Mann-Whitney U test.
    """
    import gc
    import random
    from scipy.stats import mannwhitneyu

    # Extract the distances in the desired range for actual data
    actual_over_range = [d for d in actual_distances if lower_range <= d <= upper_range]

    if len(actual_over_range) == 0:
        print("No actual distances in the specified range for Mann-Whitney U Test.")
        return None, None

    # Initialize background_over_range as empty list
    background_over_range = []

    # Process background files incrementally
    for file in background_files:
        with open(file, 'rb') as f:
            data = pickle.load(f)
            # Filter data to desired range
            filtered_data = [d for d in data if lower_range <= d <= upper_range]

            # Calculate remaining samples needed
            remaining_samples = max_background_samples - len(background_over_range)
            if remaining_samples <= 0:
                break

            if len(filtered_data) <= remaining_samples:
                background_over_range.extend(filtered_data)
            else:
                # Randomly sample to limit memory usage
                background_over_range.extend(random.sample(filtered_data, remaining_samples))

            del data
            del filtered_data
            gc.collect()

    if len(background_over_range) == 0:
        print("No background distances in the specified range for Mann-Whitney U Test.")
        return None, None

    # Perform Mann-Whitney U test
    U_statistic, p_value = mannwhitneyu(actual_over_range, background_over_range, alternative='two-sided')

    return U_statistic, p_value


def mann_whitney_u_test_v42(actual_distances, background_distances, lower_range=-10, upper_range=10):
    """
    Performs the Mann-Whitney U test on the subset of distances between 40-50 codons.
    
    Parameters:
        actual_distances (list): List of distances from actual data.
        background_distances (list): List of distances from background distribution.
    
    Returns:
        U statistic and p-value from the Mann-Whitney U test.
    """
    # Extract the 40-50 codon range for both actual and background distances
    actual_over_range = [d for d in actual_distances if lower_range <= d <= upper_range]
    background_over_range = [d for d in background_distances if 40 <= d <= 50]
    
    if len(actual_over_range) > 0 and len(background_over_range) > 0:
        # Perform Mann-Whitney U test
        U_statistic, p_value = mannwhitneyu(actual_over_range, background_over_range, alternative='two-sided')
        return U_statistic, p_value
    else:
        print("Insufficient data in the 40-50 codon range for Mann-Whitney U Test.")
        return None, None

def binomial_test_incremental(actual_distances, background_files):
    """
    Performs a binomial test to assess enrichment of slippery sequences exactly 45 codons upstream of TMDs,
    processing background data incrementally.

    Parameters:
        actual_distances (list): List of distances from actual data.
        background_files (list): List of file paths containing the background data.

    Returns:
        p_value: P-value from the binomial test.
    """
    from scipy.stats import binomtest
    import gc

    # Count successes in actual data
    k_observed = actual_distances.count(0)
    n_total = len(actual_distances)

    if n_total == 0:
        print("No actual distances available for the binomial test.")
        return None

    # Count successes and total in background data incrementally
    k_expected = 0
    n_expected = 0

    for file in background_files:
        with open(file, 'rb') as f:
            data = pickle.load(f)
            k_expected += data.count(0)
            n_expected += len(data)
            del data
            gc.collect()

    if n_expected == 0:
        print("No background distances available for the binomial test.")
        return None

    p_expected = k_expected / n_expected

    if p_expected == 0:
        print("Expected probability for 45 codons is zero. Cannot perform binomial test.")
        return None

    # Perform a one-tailed binomial test
    p_value = binomtest(k_observed, n_total, p_expected, alternative='greater').pvalue

    return p_value



def ks_test_incremental_old(actual_distances, background_files, bin_edges):
    """
    Performs the Kolmogorov-Smirnov test using incremental ECDFs.

    Parameters:
        actual_distances (list): List of distances from actual data.
        background_files (list): List of file paths containing the background data.
        bin_edges (array-like): Edges of the histogram bins.

    Returns:
        d_statistic (float): KS test statistic.
    """
    import numpy as np

    # Build background histogram incrementally
    hist_background = incremental_histogram(background_files, bin_edges)
    ecdf_background = np.cumsum(hist_background) / np.sum(hist_background)

    # Build actual data histogram
    hist_actual, _ = np.histogram(actual_distances, bins=bin_edges)
    ecdf_actual = np.cumsum(hist_actual) / np.sum(hist_actual)

    # Compute the maximum difference between the ECDFs
    d_statistic = np.max(np.abs(ecdf_actual - ecdf_background))

    return d_statistic


def ks_test_incremental(actual_distances, background_files, bin_edges):
    """
    Performs the Kolmogorov-Smirnov test using incremental ECDFs.

    Parameters:
        actual_distances (list): List of distances from actual data.
        background_files (list): List of file paths containing the background data.
        bin_edges (array-like): Edges of the histogram bins.

    Returns:
        d_statistic (float): KS test statistic.
        p_value (float): KS test p-value.
    """
    import numpy as np
    from scipy.stats import ks_2samp

    # Build background histogram incrementally
    hist_background = incremental_histogram(background_files, bin_edges)
    cdf_background = np.cumsum(hist_background) / np.sum(hist_background)

    # Build actual data histogram
    hist_actual, _ = np.histogram(actual_distances, bins=bin_edges)
    cdf_actual = np.cumsum(hist_actual) / np.sum(hist_actual)

    # Compute the maximum difference between the CDFs
    d_statistic = np.max(np.abs(cdf_actual - cdf_background))

    # Estimate the effective sample sizes
    n1 = len(actual_distances)
    n2 = np.sum(hist_background)

    # Compute the p-value using the Kolmogorov distribution
    # For large sample sizes, the p-value approaches zero
    # So we can approximate the p-value or report it as approximately zero
    en = np.sqrt(n1 * n2 / (n1 + n2))
    try:
        from scipy.stats import kstwobign
        p_value = kstwobign.sf(d_statistic * en)
    except ImportError:
        # If kstwobign is not available, use an approximation
        p_value = np.exp(-2 * (d_statistic * en) ** 2)

    return d_statistic, p_value



def compare_distributions(distances_dict, background_file_prefix, output_file_prefix):
    """
    Compares the actual and background distributions using statistical tests and plots for each category.

    Parameters:
        distances_dict (dict): Dictionary with six real slip site categories.
        background_file_prefix (str): Prefix for the background distance files.
        output_file_prefix (str): Prefix for saving the comparison plots.
    """
    import os
    import pickle
    import numpy as np
    import matplotlib.pyplot as plt
    import gc
    from scipy.stats import mannwhitneyu  # Import Mann-Whitney U test

    distances_logfile = open(output_file_prefix+'_distance_distribution_logfile.txt', 'w')

    category_mapping = {
        'strict_all': 'All',
        'loose_all': 'All',
        'strict_canonical': 'Canonical',
        'loose_canonical': 'Canonical',
        'strict_alternative': 'Alternative',
        'loose_alternative': 'Alternative'
    }

    for category, actual_distances in distances_dict.items():
        background_category = category_mapping[category]
        background_files = [f for f in os.listdir() if f.startswith(f"{background_file_prefix}_{background_category}_batch_")]

        if not actual_distances or not background_files:
            print(f"Skipping category '{category}' due to empty data.", file=distances_logfile)
            continue

        print(f"Processing category: {category}", file=distances_logfile)
        print(f"Actual Distances: {len(actual_distances)}", file=distances_logfile)

        # Define bin edges for histograms
        min_distance = min(actual_distances)
        max_distance = max(actual_distances)
        num_bins = 100  # Adjust as needed
        bin_edges = np.linspace(min_distance - 1, max_distance + 1, num_bins + 1)

        # Perform Kolmogorov-Smirnov Test using ks_test_incremental
        d_statistic, p_value_ks = ks_test_incremental(actual_distances, background_files, bin_edges)
        print(f"KS Test Statistic: {d_statistic:.4f}", file=distances_logfile)
        print(f"KS Test p-value: {p_value_ks:.4e}", file=distances_logfile)

        # Build histograms for plotting
        hist_background = incremental_histogram(background_files, bin_edges)
        cdf_background = np.cumsum(hist_background) / np.sum(hist_background)

        hist_actual, _ = np.histogram(actual_distances, bins=bin_edges)
        cdf_actual = np.cumsum(hist_actual) / np.sum(hist_actual)

        # Perform Binomial Test
        p_value_binom = binomial_test_incremental(actual_distances, background_files)
        if p_value_binom is not None:
            print(f"Binomial test p-value: {p_value_binom:.4e}", file=distances_logfile)
        else:
            print("Binomial test could not be performed.", file=distances_logfile)

        # Perform Mann-Whitney U Test
        U_statistic, p_value_mw = mann_whitney_u_test_incremental(
            actual_distances,
            background_files,
            lower_range=-10,  # Adjust the range as needed
            upper_range=10,
            #max_background_samples=100000  # Adjust the sample size as needed
        )
        if p_value_mw is not None:
            print(f"Mann-Whitney U Test Statistic: {U_statistic:.4f}", file=distances_logfile)
            print(f"Mann-Whitney U test p-value: {p_value_mw:.4e}", file=distances_logfile)
        else:
            print("Mann-Whitney U test could not be performed.", file=distances_logfile)

        # Plotting
        hist_actual_normalized = hist_actual / np.sum(hist_actual)
        hist_background_normalized = hist_background / np.sum(hist_background)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        plt.figure(figsize=(12, 8))
        plt.bar(bin_centers, hist_actual_normalized, width=np.diff(bin_edges), alpha=0.7,
                label=f'{category} Slippery Sequences', color='skyblue', edgecolor='black', align='center')
        plt.bar(bin_centers, hist_background_normalized, width=np.diff(bin_edges), alpha=0.7,
                label='Random Heptamer Positions', color='salmon', edgecolor='black', align='center')
        plt.xlabel('Distance from Nearest Upstream TMD (codons)', fontsize=14)
        plt.ylabel('Normalized Frequency', fontsize=14)
        plt.title(f'Comparison of Distance Distributions: {category} vs. Background', fontsize=16)
        plt.legend(fontsize=12)
        plt.grid(axis='y', alpha=0.75)
        plt.savefig(f"{output_file_prefix}_{category}_comparison.png")
        plt.close()
        print(f"Comparison plot for {category} saved.", file=distances_logfile)

        # Clean up
        del hist_actual, hist_background, cdf_actual, cdf_background
        gc.collect()

    distances_logfile.close()




# Serial version of the background generation function - no multiprocessing
def generate_background_distribution_serial(filtered_list, tmd_data, real_slipsite_list, num_iterations=100):
    """
    Generates a background distribution of distances by sampling random heptamer positions.
    
    Parameters:
        filtered_list (list): List of filtered transcripts with sequences.
        tmd_data (dict): Dictionary containing TMD positions for each transcript.
        real_slipsite_list (list): The actual list of slippery sequences.
        num_iterations (int): Number of random heptamer sets to generate.
        
    Returns:
        list: List of distances between random heptamers and nearest upstream TMDs.
    """
    background_distances = []
    possible_nucleotides = ['A', 'C', 'G', 'T']
    num_slipsite_samples = len(real_slipsite_list)
    
    for x in range(num_iterations):
        iteration_number = x + 1

        # Generate a random set of heptamers with the same size as real_slipsite_list
        random_heptamers = [''.join(random.choices(possible_nucleotides, k=7)) for _ in range(num_slipsite_samples)]
        
        # Analyze distances for this random set
        distance_data = analyze_slipsite_to_tmd_distances(filtered_list, tmd_data, random_heptamers, outname="background_analysis")
        
        # Accumulate distances from this iteration
        background_distances.extend(distance_data['All'])  # You could store distances for Canonical/Alternative separately if needed.
    
        sys.stdout.write(f"\rCompleted iteration {iteration_number} of {num_iterations}")
        sys.stdout.flush()

    return background_distances


def analyze_slipsite_to_tmd_distances(filtered_list, tmd_data, slipsite_list, outname):
    """
    Analyzes the distances between slippery sequences and their nearest upstream TMDs, now categorized by strict vs loose slip sites.
    
    Parameters:
        filtered_list (list): List of labeled transcripts with sequences and canonical/non-canonical labels.
        tmd_data (dict): Dictionary containing TMD positions for each transcript.
        slipsite_list (list): List of known slippery sequences (loose set).
        outname (str): Base output name for files.
    
    Returns:
        dict: A dictionary containing categorized distances for strict/loose slip sites (all, canonical, alternative).
    """
    # List of 24 canonical slip sites (strict set)
    strict_slip_sites = ['AAAAAAA', 'AAAAAAT', 'AAAAAAC', 'AAATTTA', 'AAATTTT', 'AAATTTC',
                         'TTTAAAA', 'TTTAAAT', 'TTTAAAC', 'TTTTTTA', 'TTTTTTT', 'TTTTTTC',
                         'CCCAAAA', 'CCCAAAT', 'CCCAAAC', 'CCCTTTA', 'CCCTTTT', 'CCCTTTC',
                         'GGGAAAA', 'GGGAAAT', 'GGGAAAC', 'GGGTTTA', 'GGGTTTT', 'GGGTTTC']

    # Initialize distance data structure for six categories
    distance_data = {
        "strict_all": [],
        "strict_canonical": [],
        "strict_alternative": [],
        "loose_all": [],
        "loose_canonical": [],
        "loose_alternative": []
    }

    for transcript_data in filtered_list:
        ensembl_transcript_id = transcript_data[1]  # Assuming transcript ID is at index 1
        sequence_raw = transcript_data[5]  # Assuming sequence is at index 5
        canonical_status = transcript_data[-1]  # Label is now at the last index
        
        # Clean sequence and convert to codons
        sequence = ''.join(filter(str.isalpha, sequence_raw)).upper()
        codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]
        
        # Skip if transcript has no TMD data
        if ensembl_transcript_id not in tmd_data:
            continue
        
        # Get TMD positions for this transcript
        tmd_positions = []
        for tmd in tmd_data[ensembl_transcript_id]:
            tm_start = tmd['adjusted_tmd_start_position'] - 1
            tm_end = tmd['adjusted_tmd_end_position'] - 1
            tmd_positions.append((tm_start, tm_end))
        
        # Iterate over codon positions to find slippery sequences
        for i in range(len(codons) - 2):
            heptamer_long = ''.join(codons[i:i+3])
            heptamer = heptamer_long[2:]  # Extract the heptamer
            
            if heptamer in slipsite_list or heptamer in strict_slip_sites:
                slipsite_position = i  # Position of the slippery sequence start codon
                
                # Find all upstream TMDs
                upstream_tmds = [tm for tm in tmd_positions if tm[1] < slipsite_position]
                
                if not upstream_tmds:
                    continue  # No upstream TMDs; skip this slippery sequence
                
                # Calculate distances to all upstream TMDs and select the minimum (nearest)
                distances_to_tmds = [slipsite_position - tm[1] for tm in upstream_tmds]
                distance_from_harrington_ideality = [tmd_ss_distance - 45 for tmd_ss_distance in distances_to_tmds]
                distance_from_motif_ideality_upstream_only = [gap_ideality for gap_ideality in distance_from_harrington_ideality if gap_ideality > -45] # Could use -25 to allow a full TMD between slip site and next upstream TMD, or could use -45 
                abs_distance_from_harrington_ideality = [abs(distance) for distance in distance_from_motif_ideality_upstream_only]
                min_distance_idx = abs_distance_from_harrington_ideality.index(min(abs_distance_from_harrington_ideality))
                min_distance = distance_from_motif_ideality_upstream_only[min_distance_idx]
                
                # Categorize by strict vs loose slip sites
                if heptamer in strict_slip_sites:
                    distance_data["strict_all"].append(min_distance)
                    if canonical_status == 'Canonical':
                        distance_data["strict_canonical"].append(min_distance)
                    else:
                        distance_data["strict_alternative"].append(min_distance)

                # Always add to the loose+strict category
                distance_data["loose_all"].append(min_distance)
                if canonical_status == 'Canonical':
                    distance_data["loose_canonical"].append(min_distance)
                else:
                    distance_data["loose_alternative"].append(min_distance)

    # Save distances to CSV with categories
    raw_data_file = f"{outname}_v55_slipsite_tmd_distances_strict_and_loose_slipsite_distribution.csv"
    with open(raw_data_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Distance_from_TMD', 'Category'])
        
        # Write distances for all categories to the CSV
        for distance in distance_data['strict_all']:
            csv_writer.writerow([distance, 'strict_all'])
        for distance in distance_data['strict_canonical']:
            csv_writer.writerow([distance, 'strict_canonical'])
        for distance in distance_data['strict_alternative']:
            csv_writer.writerow([distance, 'strict_alternative'])
        for distance in distance_data['loose_all']:
            csv_writer.writerow([distance, 'loose_all'])
        for distance in distance_data['loose_canonical']:
            csv_writer.writerow([distance, 'loose_canonical'])
        for distance in distance_data['loose_alternative']:
            csv_writer.writerow([distance, 'loose_alternative'])

    # Return categorized distances
    return distance_data


def analyze_slipsite_to_tmd_distances_loose_only(filtered_list, tmd_data, slipsite_list, outname):
    """
    Analyzes the distances between slippery sequences and their nearest upstream TMDs.
    
    Parameters:
        filtered_list (list): List of labeled transcripts with sequences and canonical/non-canonical labels.
        tmd_data (dict): Dictionary containing TMD positions for each transcript.
        slipsite_list (list): List of known slippery sequences.
        outname (str): Base output name for files.
    
    Returns:
        dict: A dictionary containing categorized distances.
    """

    distances_all = []
    distances_canonical = []
    distances_non_canonical = []

    for transcript_data in filtered_list:
        ensembl_transcript_id = transcript_data[1]  # Assuming transcript ID is at index 1
        sequence_raw = transcript_data[5]  # Assuming sequence is at index 5
        canonical_status = transcript_data[-1]  # Label is now at the last index
        
        # Clean sequence and convert to codons
        sequence = ''.join(filter(str.isalpha, sequence_raw)).upper()
        codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]
        
        # Skip if transcript has no TMD data
        if ensembl_transcript_id not in tmd_data:
            continue
        
        # Get TMD positions for this transcript
        tmd_positions = []
        for tmd in tmd_data[ensembl_transcript_id]:
            tm_start = tmd['adjusted_tmd_start_position'] - 1
            tm_end = tmd['adjusted_tmd_end_position'] - 1
            tmd_positions.append((tm_start, tm_end))
        
        # Iterate over codon positions to find slippery sequences
        for i in range(len(codons) - 2):
            heptamer_long = ''.join(codons[i:i+3])
            heptamer = heptamer_long[2:]  # Extract the heptamer
            
            if heptamer in slipsite_list:
                slipsite_position = i  # Position of the slippery sequence start codon
                
                # Find all upstream TMDs
                upstream_tmds = [tm for tm in tmd_positions if tm[1] < slipsite_position]
                
                if not upstream_tmds:
                    continue  # No upstream TMDs; skip this slippery sequence
                
                # Calculate distances to all upstream TMDs and select the minimum (nearest)
                distances_to_tmds = [slipsite_position - tm[1] for tm in upstream_tmds]
                distance_from_harrington_ideality = [tmd_ss_distance - 45 for tmd_ss_distance in distances_to_tmds]
                distance_from_motif_ideality_upstream_only = [gap_ideality for gap_ideality in distance_from_harrington_ideality if gap_ideality > -45] # Could use -25 to allow a full TMD between slip site and next upstream TMD, or could use -45 
                abs_distance_from_harrington_ideality = [abs(distance) for distance in distance_from_motif_ideality_upstream_only]
                min_distance_idx = abs_distance_from_harrington_ideality.index(min(abs_distance_from_harrington_ideality))
                min_distance = distance_from_motif_ideality_upstream_only[min_distance_idx]

                distances_all.append(min_distance)
                
                # Append to the appropriate category based on canonical status
                if canonical_status == 'Canonical':
                    distances_canonical.append(min_distance)
                else:
                    distances_non_canonical.append(min_distance)
    
    # Save distances to CSV with categories
    raw_data_file = f"{outname}_heptamer_tmd_distances.csv"
    with open(raw_data_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Distance_from_TMD', 'Category'])
        for distance in distances_all:
            csv_writer.writerow([distance, 'All'])
        for distance in distances_canonical:
            csv_writer.writerow([distance, 'Canonical'])
        for distance in distances_non_canonical:
            csv_writer.writerow([distance, 'Alternative'])

    # Return categorized distances
    return {
        'All': distances_all,
        'Canonical': distances_canonical,
        'Alternative': distances_non_canonical
    }




# Print motif results to CSV
def print_motif_csv(output_array, outname):
    with open(outname + '.csv', 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        
        # Updated headers to match the structure of candidate_entry and include Uniprot information
        csv_writer.writerow([
            'transcript_id', 'gene_id', 'gene_name', 'tm_start', 'tm_end', 'tm_dg', 'tm_domain', 'tm_codons',
            'ss_coord_nuc', 'ss_coord_codon', 'chr_coordinate', 'slipsite_seq', 'gap_size', 'location_ideality',
            'gap_nuc', 'gap_%gc', 'gap_avg_codon_freq', 'fs_peptide', 'terminated', 'transcript_seq_nuc',
            'transcript_seq_res', 'friction_score', 'slipsite_canonical?', 'stem_loop_found', 'pseudoknot_found', 'secondary_structure',
            'canonical_status', 'tmd_type', 'uniprot_id', 'similarity_score'
        ])
        
        # Write motif entries to the CSV
        for motif in output_array:
            csv_writer.writerow(motif)




def summarize_transcripts(tmd_data, canonical_count, alt_spliced_count, harrington_motifs):
    # X1: Total number of TMDs in the dataset
    X1 = sum([len(tmd_info) for tmd_info in tmd_data.values()])

    # X2: Number of canonical sequences
    X2 = canonical_count

    # X3: Number of alternatively spliced sequences
    X3 = alt_spliced_count

    # X4: Number of identified Harrington motifs
    X4 = len(harrington_motifs)

    return X1, X2, X3, X4


def generate_random_heptamers(num_slipsite_samples, possible_nucleotides):
    """Generates a random list of heptamers of a given size."""
    return [''.join(random.choices(possible_nucleotides, k=7)) for _ in range(num_slipsite_samples)]


def draw_random_heptamers_less_efficient(num_slipsite_samples, shuffled_list_heptamers):
    random_heptamers = []
    num_unique_heptamers_drawn = 0
    while num_unique_heptamers_drawn < num_slipsite_samples:
        random_heptamer_candidate = random.choice(shuffled_list_heptamers)
        if random_heptamer_candidate not in random_heptamers:
            random_heptamers.append(random_heptamer_candidate)
            num_unique_heptamers_drawn += 1
    return random_heptamers

def draw_random_heptamers(num_slipsite_samples, shuffled_list_heptamers):

    random_heptamers = set()
    
    while len(random_heptamers) < num_slipsite_samples:
        random_heptamer_candidate = random.choice(shuffled_list_heptamers)
        random_heptamers.add(random_heptamer_candidate)
    
    return list(random_heptamers)


def analyze_iteration(filtered_list, tmd_data, num_slipsite_samples, shuffled_list_of_heptamers_all, shuffled_list_of_heptamers_canonical, shuffled_list_of_heptamers_alternative, possible_nucleotides, transcript_type):
    """
    Performs a single iteration of distance analysis for a random set of heptamers based on the transcript type.
    
    Args:
        transcript_type: Indicates whether to use 'All', 'Canonical', or 'Alternative' heptamer list.
    """
    if transcript_type == 'All':
        random_heptamers = draw_random_heptamers(num_slipsite_samples, shuffled_list_of_heptamers_all)
    elif transcript_type == 'Canonical':
        random_heptamers = draw_random_heptamers(num_slipsite_samples, shuffled_list_of_heptamers_canonical)
    elif transcript_type == 'Alternative':
        random_heptamers = draw_random_heptamers(num_slipsite_samples, shuffled_list_of_heptamers_alternative)
    else:
        raise ValueError(f"Invalid transcript type: {transcript_type}")

    # Analyze distances using the chosen random heptamers
    distance_data = analyze_slipsite_to_tmd_distances_loose_only(filtered_list, tmd_data, random_heptamers, outname="background_analysis_ignore_me")
    
    return distance_data

def generate_background_distribution(filtered_list, tmd_data, real_slipsite_list, shuffled_list_of_heptamers_all, shuffled_list_of_heptamers_canonical, shuffled_list_of_heptamers_alternative, num_iterations=1024, batch_size=128, num_processors=8):
    """
    Generates a background distribution of distances by sampling random heptamer positions for All, Canonical, and Alternative transcripts.
    """
    # Initialize the background distance dictionary
    background_distances = {
        'All': [],
        'Canonical': [],
        'Alternative': []
    }

    possible_nucleotides = ['A', 'C', 'G', 'T']
    num_slipsite_samples = len(real_slipsite_list)

    # Define transcript types to process
    transcript_types = ['All', 'Canonical', 'Alternative']

    # Calculate the number of batches
    total_batches = (num_iterations + batch_size - 1) // batch_size

    for batch_num in range(total_batches):
        current_batch_size = min(batch_size, num_iterations - batch_num * batch_size)
        print(f"Processing batch {batch_num + 1} of {total_batches}, iterations: {current_batch_size}")

        # Prepare arguments for each batch
        batch_args = []
        for _ in range(current_batch_size):
            batch_args.append((filtered_list, tmd_data, num_slipsite_samples, shuffled_list_of_heptamers_all, shuffled_list_of_heptamers_canonical, shuffled_list_of_heptamers_alternative, possible_nucleotides))

        # Process each transcript type separately
        for transcript_type in transcript_types:
            pool_args = [(filtered_list, tmd_data, num_slipsite_samples, shuffled_list_of_heptamers_all, shuffled_list_of_heptamers_canonical, shuffled_list_of_heptamers_alternative, possible_nucleotides, transcript_type) for _ in range(current_batch_size)]

            with mp.Pool(processes=num_processors) as pool:
                results = pool.starmap(analyze_iteration, pool_args)

            # Aggregate results and write to disk
            batch_distances = []
            for result in results:
                batch_distances.extend(result[transcript_type])

            # Save batch distances to disk
            batch_file = f"background_distances_{transcript_type}_batch_{batch_num + 1}.pkl"
            with open(batch_file, 'wb') as f:
                pickle.dump(batch_distances, f)

            print(f"Saved batch {batch_num + 1} distances for {transcript_type} to {batch_file}")

            del batch_distances
            del results
            gc.collect()

    # After all batches are processed, load and combine distances
    for transcript_type in transcript_types:
        all_batches = []
        for batch_num in range(total_batches):
            batch_file = f"background_distances_{transcript_type}_batch_{batch_num + 1}.pkl"
            with open(batch_file, 'rb') as f:
                batch_distances = pickle.load(f)
                all_batches.extend(batch_distances)

            del batch_distances
            gc.collect()

        background_distances[transcript_type] = all_batches
        print(f"Total distances for {transcript_type}: {len(all_batches)}")

        del all_batches
        gc.collect()

    return background_distances




# This function has been superceded by a new one but we still use it
def calculate_statistics(harrington_motifs, total_codons, gap_near, gap_far, num_trials=100000):
    # Count the number of observed Harrington motifs (X4)
    observed_count = len(harrington_motifs)

    # Randomly permute the positions of slippery sequences and recalculate counts
    random_counts = []
    for _ in range(num_trials):
        random_positions = random.sample(range(total_codons), observed_count)
        random_count = sum([1 for pos in random_positions if gap_near <= pos <= gap_far])
        random_counts.append(random_count)

    # X5: Calculate the fold difference between observed and mean expected count
    expected_mean = np.mean(random_counts) if np.mean(random_counts) > 0 else 1  # Avoid division by zero
    X5 = observed_count / expected_mean

    # X6: Calculate the p-value (probability of observing a count >= observed_count by chance)
    X6 = np.sum(np.array(random_counts) >= observed_count) / num_trials

    return X5, X6


def import_uniprot_fasta(path_fasta):
    fasta_generator = load_fasta_data(path_fasta)
    raw_fasta_output = []
    for header, sequence in fasta_generator:
        sequence_str = ''.join(sequence)
        output_line = [header, sequence_str]
        raw_fasta_output.append(output_line)

    uniprot_fasta_data = []
    for entry in raw_fasta_output:
        header = entry[0]
        uniprot_seq = entry[1]
        uniprot_id = header.split('|')[1]
        gene_name = header.split('=')[-3].split()[0]
        output_line = [uniprot_id, gene_name, uniprot_seq]
        uniprot_fasta_data.append(output_line)
    return uniprot_fasta_data


def print_results(output_data, outname):
    myfile = open(outname+'.csv', 'w')
    for line in output_data:
        line_string = [str(x) for x in line]
        csv_line = ','.join(line_string)
        print(csv_line, file = myfile)
    myfile.close()


def print_ensembl_uniprot_dict(ensembl_uniprot_dict, outname):
    ensembl_uniprot_list = [['ensembl_transcript_id', 'uniprot_accession_id', 'similarity_score']]
    for key, value in ensembl_uniprot_dict.items():
        ensembl_transcript_id = key
        [uniprot_accession_id, similarity_score] = value
        output_entry = [ensembl_transcript_id, uniprot_accession_id, similarity_score]
        ensembl_uniprot_list.append(output_entry)
    output_csv_name = outname+'_ensembl_uniprot_dict_full_blastp'
    print_results(ensembl_uniprot_list, output_csv_name)



def load_uniprot_excel_to_dict(excel_path, key_column='Entry'):
    """
    Loads a Uniprot Excel file and converts it into a dictionary with the specified key column.
    
    Parameters:
    excel_path (str): Path to the Excel file containing Uniprot data.
    key_column (str): The column to use as the key for the dictionary (e.g., 'Sequence', 'Entry').
    
    Returns:
    dict: A dictionary where keys are from the key_column and values are dictionaries of the other fields, 
          including the sequence itself.
    """
    # Load the Excel file into a Pandas DataFrame
    df = pd.read_excel(excel_path)

    # Fill any NaN values with empty strings to avoid issues with missing data
    df = df.fillna('')

    # Ensure the key_column exists in the DataFrame
    if key_column not in df.columns:
        raise ValueError(f"Key column '{key_column}' does not exist in the Excel file.")
    
    # Add the sequence as a value in the dictionary (even though it's the key)
    # Create the dictionary and make sure 'Sequence' is part of the values as well.
    uniprot_dict = df.set_index(key_column).T.to_dict('dict')

    # Explicitly add 'Sequence' into the dictionary values
    #for seq in uniprot_dict.keys():
    #    uniprot_dict[seq]['Sequence'] = seq  # Add the sequence as a value

    return uniprot_dict

def split_up_harrington_motifs(harrington_motifs, location_ideality_threshold, stem_loop_required=False):
    '''
    # Indexes of harrington_motifs here
    candidate_entry = ['transcript_id',         # 0
                       'ensembl_gene_name',     # 1
                       'common_gene_name',      # 2
                       'tm_start',              # 3
                       'tm_end',                # 4
                       'tm_dg',                 # 5
                       'tm_domain_str',         # 6
                       'tm_codons_str',         # 7
                       'ss_coordinate_nuc',     # 8
                       'ss_coordinate_codon',   # 9
                       'chr_coordinate',        # 10
                       'candidate_slipsite',    # 11
                       'gap_size',              # 12
                       'location_ideality',     # 13
                       'gap_nuc',               # 14
                       'gap_gc_content',        # 15
                       'gap_avg_codon_freq',    # 16
                       'fs_peptide',            # 17
                       'terminated',            # 18
                       'transcript_seq',        # 19
                       'transcript_res_str',    # 20
                       'friction_score',        # 21
                       'slipsite_strict',       # 22
                       'stem_loop_presence',    # 23
                       'pseudokno_presence',    # 24
                       'structure',             # 25
                       'canonical_status',      # 26
                       'tmd_type'               # 27
                       'uniprot_id'             # 28
                       'similarity_score']      # 29
    '''

    # Initialize strict subsets
    strict_all_motifs = []
    strict_canonical_motifs = []
    strict_alternative_motifs = []

    # Initialize loose subsets
    loose_all_motifs = []
    loose_canonical_motifs = []
    loose_alternative_motifs = []

    for entry in harrington_motifs:
        location_ideality = abs(int(entry[13]))
        slipsite_strict = entry[22]
        canonical_status = entry[26]
        terminated_status = entry[18]
        stem_loop_presence = entry[23]

        subset_dimensions = [slipsite_strict, canonical_status]

        if all([location_ideality <= location_ideality_threshold,
                stem_loop_required == True,
                stem_loop_presence == 'yes']) or all([stem_loop_required == False,
                                                      location_ideality <= location_ideality_threshold]):

            if subset_dimensions == ['Strict', 'Canonical']:
                strict_all_motifs.append(entry)
                strict_canonical_motifs.append(entry)
                loose_all_motifs.append(entry)
                loose_canonical_motifs.append(entry)
        
            if subset_dimensions == ['Strict', 'Alternative']:
                strict_all_motifs.append(entry)
                strict_alternative_motifs.append(entry)
                loose_all_motifs.append(entry)
                loose_alternative_motifs.append(entry)
        
            if subset_dimensions == ['Loose', 'Canonical']:
                loose_all_motifs.append(entry)
                loose_canonical_motifs.append(entry)
        
            if subset_dimensions == ['Loose', 'Alternative']:
                loose_all_motifs.append(entry)
                loose_alternative_motifs.append(entry)


    motif_subsets = [strict_all_motifs, strict_canonical_motifs, strict_alternative_motifs, loose_all_motifs, loose_canonical_motifs, loose_alternative_motifs]

    return motif_subsets


def filter_non_canonical_transcripts(filtered_list):
    filtered_list_canonical_only = []
    filtered_list_alternative_only = []
    for entry in filtered_list:
        canonical_status = entry[-1]
        if canonical_status == 'Canonical':
            filtered_list_canonical_only.append(entry)
        if canonical_status == 'Alternative':
            filtered_list_alternative_only.append(entry)
    return filtered_list_canonical_only, filtered_list_alternative_only


def extract_heptamers_for_prf(sequence_str):
    """
    Extract heptamers that could be involved in -1 PRF from a nucleotide sequence in codon space.
    The heptamer starts at the third base of the first codon and spans 7 nucleotides.

    Args:
        sequence_str: Nucleotide sequence (CDS) as a string.
        
    Returns:
        List of heptamers that could be involved in -1 PRF.
    """
    # Step 1: Split the nucleotide sequence into codons
    seq_codons = textwrap.wrap(sequence_str, 3)
    
    seq_heptamers = []
    
    # Step 2: Build the heptamers by starting at the 3rd base of the first codon
    for i in range(len(seq_codons) - 2):  # Need at least 3 codons to form a heptamer
        codon1 = seq_codons[i][2:]    # Take the 3rd base of the first codon
        codon2 = seq_codons[i+1]      # Take all of the second codon
        codon3 = seq_codons[i+2]      # Take all of the third codon
        
        # Combine them to form a heptamer (7 nucleotides)
        heptamer = codon1 + codon2 + codon3
        
        seq_heptamers.append(heptamer)
    
    return seq_heptamers




def calculate_heptamer_frequencies_prf(filtered_list, tmd_data):
    """
    Calculate the frequency of heptamers that could be involved in -1 PRF in the CDS sequences.

    Args:
        filtered_list: List of CDS sequences from the transcriptome.
        
    Returns:
        heptamer_frequencies: A dictionary of heptamer frequencies.
    """
    heptamer_counts = Counter()
    total_heptamers = 0
    
    list_of_heptamers = []

    filtered_list_with_tmds = []
    for entry in filtered_list:
        transcript_id = entry[1]
        if transcript_id in tmd_data:
            filtered_list_with_tmds.append(entry)

    # Iterate through each transcript's CDS sequence
    for transcript_data in filtered_list_with_tmds:
        sequence_str = transcript_data[5]  # CDS sequence (nucleotide space)
        
        # Extract heptamers for -1 PRF
        seq_heptamers = extract_heptamers_for_prf(sequence_str)
        for heptamer in seq_heptamers:
            list_of_heptamers.append(heptamer)
        
        # Count the heptamers
        for heptamer in seq_heptamers:
            heptamer_counts[heptamer] += 1
            total_heptamers += 1
    
    set_unique_heptamers = list(set(list_of_heptamers))
    shuffled_list_of_heptamers = list_of_heptamers.copy()
    random.shuffle(shuffled_list_of_heptamers)

    # Convert counts to frequencies
    heptamer_frequencies = {heptamer: count / total_heptamers for heptamer, count in heptamer_counts.items()}
    
    return heptamer_frequencies, set_unique_heptamers, shuffled_list_of_heptamers




def calculate_p_null_slippery(heptamer_frequencies, slippery_sequences):
    """
    Calculate the probability of finding a slippery sequence by summing their frequencies.

    Args:
        heptamer_frequencies: A dictionary of heptamer frequencies.
        slippery_sequences: A set of slippery sequences to consider.
        
    Returns:
        p_null: The probability of observing a slippery sequence at random.
    """
    p_null = sum(heptamer_frequencies.get(seq, 0) for seq in slippery_sequences)
    return p_null



def run_binomial_test(motif_subset, p_null, num_trials):
    """
    Perform a binomial test on the motif subset to check for enrichment at the exact gap size (e.g., 45 codons).
    
    Args:
        motif_subset: List of motifs (each motif has gap_size in position [12])
        gap_size_target: The gap size you are testing for (default is 45 codons).
        p_null: The null probability, representing the expected probability of finding a slippery sequence
                by chance based on heptamer frequencies.
                
    Returns:
        p_value: The p-value of the binomial test.
        observed_count: The number of successes (motifs at the target gap size).
        total_trials: The total number of trials (motifs).
    """
    observed_count = sum(1 for motif in motif_subset if motif[12] == 45)  # gap_size == 45 codons
    
    # Perform the binomial test
    p_value = binomtest(observed_count, num_trials, p_null, alternative='greater')
    
    return p_value, observed_count


def analyze_slipsite_to_tmd_distances_contingency_table(filtered_list, tmd_data, slipsite_list, all_heptamers):
    """
    Analyzes the distances between slippery sequences and their nearest upstream TMDs.
    
    Parameters:
        filtered_list (list): List of labeled transcripts with sequences and canonical/non-canonical labels.
        tmd_data (dict): Dictionary containing TMD positions for each transcript.
        slipsite_list (list): List of known slippery sequences.
        outname (str): Base output name for files.
    
    Returns:
        dict: A dictionary containing categorized distances.
    """

    non_slippery_heptamers = []
    for heptamer in all_heptamers:
        if heptamer not in slipsite_list:
            non_slippery_heptamers.append(heptamer)

    counts = {
        'slippery_ideal': 0,          # A
        'non_slippery_ideal': 0,      # B
        'slippery_other': 0,          # C
        'non_slippery_other': 0       # D
    }

    for transcript_data in filtered_list:
        ensembl_transcript_id = transcript_data[1]  # Assuming transcript ID is at index 1
        sequence_raw = transcript_data[5]  # Assuming sequence is at index 5
        canonical_status = transcript_data[-1]  # Label is now at the last index

        # Clean sequence and convert to codons
        sequence = ''.join(filter(str.isalpha, sequence_raw)).upper()
        codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]

        # Skip if transcript has no TMD data
        if ensembl_transcript_id not in tmd_data:
            continue

        # Get TMD positions for this transcript
        tmd_positions = []
        for tmd in tmd_data[ensembl_transcript_id]:
            tm_start = tmd['adjusted_tmd_start_position'] - 1
            tm_end = tmd['adjusted_tmd_end_position'] - 1
            tmd_positions.append((tm_start, tm_end))

        # Iterate over codon positions to find slippery sequences
        for i in range(len(codons) - 2):
            heptamer_long = ''.join(codons[i:i+3])
            heptamer = heptamer_long[2:]  # Extract the heptamer
            
            heptamer_slippery = heptamer in slipsite_list
            heptamer_position = i  # Position of the slippery sequence start codon
                
            # Find all upstream TMDs
            upstream_tmds = [tm for tm in tmd_positions if tm[1] < heptamer_position]
                
            if not upstream_tmds:
                continue  # No upstream TMDs; skip this slippery sequence
                
            # Calculate distances to all upstream TMDs and select the minimum (nearest)
            distances_to_tmds = [heptamer_position - tm[1] for tm in upstream_tmds]
            distance_from_harrington_ideality = [tmd_ss_distance - 45 for tmd_ss_distance in distances_to_tmds]
            distance_from_motif_ideality_upstream_only = [gap_ideality for gap_ideality in distance_from_harrington_ideality if gap_ideality > -45] # Could use -25 to allow a full TMD between slip site and next upstream TMD, or could use -45 
            abs_distance_from_harrington_ideality = [abs(distance) for distance in distance_from_motif_ideality_upstream_only]
            min_distance_idx = abs_distance_from_harrington_ideality.index(min(abs_distance_from_harrington_ideality))
            min_distance = distance_from_motif_ideality_upstream_only[min_distance_idx]


            if min_distance == 0:
                if heptamer_slippery:
                    counts['slippery_ideal'] += 1  # A
                else:
                    counts['non_slippery_ideal'] += 1  # B
            else:
                if heptamer_slippery:
                    counts['slippery_other'] += 1  # C
                else:
                    counts['non_slippery_other'] += 1  # D

    return counts



def analyze_slipsite_to_tmd_distances_basic(filtered_list, tmd_data, slipsite_list):
    """
    Analyzes the distances between slippery sequences and their nearest upstream TMDs.
    
    Parameters:
        filtered_list (list): List of labeled transcripts with sequences and canonical/non-canonical labels.
        tmd_data (dict): Dictionary containing TMD positions for each transcript.
        slipsite_list (list): List of known slippery sequences.
        outname (str): Base output name for files.
    
    Returns:
        dict: A dictionary containing categorized distances.
    """

    distances = []
    num_ideal_distances = 0

    for transcript_data in filtered_list:
        ensembl_transcript_id = transcript_data[1]  # Assuming transcript ID is at index 1
        sequence_raw = transcript_data[5]  # Assuming sequence is at index 5
        canonical_status = transcript_data[-1]  # Label is now at the last index

        # Clean sequence and convert to codons
        sequence = ''.join(filter(str.isalpha, sequence_raw)).upper()
        codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]

        # Skip if transcript has no TMD data
        if ensembl_transcript_id not in tmd_data:
            continue

        # Get TMD positions for this transcript
        tmd_positions = []
        for tmd in tmd_data[ensembl_transcript_id]:
            tm_start = tmd['adjusted_tmd_start_position'] - 1
            tm_end = tmd['adjusted_tmd_end_position'] - 1
            tmd_positions.append((tm_start, tm_end))

        # Iterate over codon positions to find slippery sequences
        for i in range(len(codons) - 2):
            heptamer_long = ''.join(codons[i:i+3])
            heptamer = heptamer_long[2:]  # Extract the heptamer
            
            if heptamer in slipsite_list:
                slipsite_position = i  # Position of the slippery sequence start codon
                
                # Find all upstream TMDs
                upstream_tmds = [tm for tm in tmd_positions if tm[1] < slipsite_position]
                
                if not upstream_tmds:
                    continue  # No upstream TMDs; skip this slippery sequence
                
                # Calculate distances to all upstream TMDs and select the minimum (nearest)
                distances_to_tmds = [slipsite_position - tm[1] for tm in upstream_tmds]
                distance_from_harrington_ideality = [tmd_ss_distance - 45 for tmd_ss_distance in distances_to_tmds]
                distance_from_motif_ideality_upstream_only = [gap_ideality for gap_ideality in distance_from_harrington_ideality if gap_ideality > -45] # Could use -25 to allow a full TMD between slip site and next upstream TMD, or could use -45 
                abs_distance_from_harrington_ideality = [abs(distance) for distance in distance_from_motif_ideality_upstream_only]
                min_distance_idx = abs_distance_from_harrington_ideality.index(min(abs_distance_from_harrington_ideality))
                min_distance = distance_from_motif_ideality_upstream_only[min_distance_idx]
                min_distance_abs = abs(min_distance)
                if min_distance == 0:
                    num_ideal_distances += 1

                #distances.append(min_distance_abs)
                distances.append(min_distance)
    
    # Return categorized distances
    return distances, num_ideal_distances


def process_iteration_for_tmd_shuffling(args):
    filtered_list, tmd_data, slipsite_list = args
    shuffled_tmd_data = shuffle_tmd_positions(filtered_list, tmd_data)
    shuffled_distances, num_ideal_distances = analyze_slipsite_to_tmd_distances_basic(filtered_list, shuffled_tmd_data, slipsite_list)
    return shuffled_distances, num_ideal_distances



def run_tmd_shuffling_test_parallel(filtered_list, tmd_data, slipsite_list, num_iterations=128, num_processes=8):
    """
    Perform the TMD shuffling test in parallel and compare distances to the original TMD distances.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary of TMD positions for each transcript.
        slipsite_list: List of known slippery sequences.
        num_iterations: Number of shuffling iterations to run (default is 1000).
        num_processes: Number of parallel processes (default is 4).
    
    Returns:
        shuffled_distances: A list of distances between shuffled TMDs and slip sites (one list for each iteration).
    """
    
    # Prepare the arguments for each process
    args = [(filtered_list, tmd_data, slipsite_list) for _ in range(num_iterations)]

    # Use multiprocessing Pool to parallelize the iterations
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_iteration_for_tmd_shuffling, args) # This needs to return num_ideal_distances too - but I need to incorporate a new list containing this number for each iteration

    # Unpack the results
    shuffled_distances_list = []
    num_ideal_distances_list = []

    for shuffled_distances, num_ideal_distances in results:
        shuffled_distances_list.append(shuffled_distances)
        num_ideal_distances_list.append(num_ideal_distances)

    return shuffled_distances_list, num_ideal_distances_list




import random

def shuffle_tmd_positions(filtered_list, tmd_data):
    """
    Shuffle the positions of TMDs within each transcript while maintaining the constraints:
    - Preserve the number of TMDs.
    - Preserve the length of each TMD.
    - Preserve the order of TMDs.
    - Ensure TMDs do not overlap.
    - Maintain at least minimal_spacing codons between TMDs.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary containing TMD positions for each transcript.
        minimal_spacing: Minimum number of codons required between TMDs (default is 3).
        
    Returns:
        shuffled_tmd_data: A new dictionary of shuffled TMD positions for each transcript.
    """
    # SET MINIMAL SPACING HERE
    minimal_spacing = 3
    # ADJUST MINIMAL SPACING IF DESIRED

    shuffled_tmd_data = {}
    
    # Build a dictionary with transcript sequences (for convenience)
    transcript_sequences = {entry[1]: entry[5] for entry in filtered_list}
    
    # Iterate over each transcript
    for transcript_id, tmd_list in tmd_data.items():
        if transcript_id not in transcript_sequences:
            continue
        
        sequence = transcript_sequences[transcript_id]
        seq_length = len(sequence) // 3  # Convert nucleotide sequence length to codon length
        num_tmds = len(tmd_list)
        
        # Get the lengths of each TMD
        tmd_lengths = [
            tmd['adjusted_tmd_end_position'] - tmd['adjusted_tmd_start_position'] + 1
            for tmd in tmd_list
        ]
        total_tmd_length = sum(tmd_lengths)
        
        # Calculate the minimum total gap required between TMDs
        total_min_gap = (num_tmds - 1) * minimal_spacing
        total_gap_space = seq_length - total_tmd_length
        
        # Check if it's possible to place TMDs under constraints
        extra_space = total_gap_space - total_min_gap
        if extra_space < 0:
            # Cannot fit TMDs under constraints
            continue
        
        # Distribute extra_space randomly among the gaps
        num_gaps = num_tmds + 1  # Gaps before, between, and after TMDs
        minimal_gaps = [0] + [minimal_spacing] * (num_tmds - 1) + [0]
        
        # Randomly distribute extra_space among the gaps
        delta_gaps = [0] * num_gaps
        if extra_space > 0:
            random_points = sorted(random.sample(range(extra_space + num_gaps - 1), num_gaps - 1))
            delta_gaps = [random_points[0]] + \
                         [random_points[i] - random_points[i - 1] - 1 for i in range(1, num_gaps - 1)] + \
                         [extra_space + num_gaps - 2 - random_points[-1]]
        
        # Calculate the actual gaps
        gaps = [minimal_gaps[i] + delta_gaps[i] for i in range(num_gaps)]
        
        # Determine the start positions for each TMD
        start_positions = []
        current_position = gaps[0]
        for i in range(num_tmds):
            tmd_start = current_position
            tmd_end = tmd_start + tmd_lengths[i] - 1
            start_positions.append((tmd_start, tmd_end, tmd_list[i]['tmd_type']))
            current_position = tmd_end + 1 + gaps[i + 1]
        
        # Build the new TMD list with shuffled positions
        new_tmd_list = []
        for tmd_info, (new_start, new_end, tmd_type) in zip(tmd_list, start_positions):
            new_tmd_list.append({
                'adjusted_tmd_start_position': new_start,
                'adjusted_tmd_end_position': new_end,
                'tmd_type': tmd_type,
            })
        
        # Add the shuffled TMD positions to the new dictionary
        shuffled_tmd_data[transcript_id] = new_tmd_list
    
    return shuffled_tmd_data




def run_tmd_shuffling_test(filtered_list, tmd_data, slipsite_list, num_iterations=128):
    """
    Perform the TMD shuffling test and compare distances to the original TMD distances.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary of TMD positions for each transcript.
        slipsite_list: List of known slippery sequences.
        num_iterations: Number of shuffling iterations to run (default is 1000).
    
    Returns:
        shuffled_distances: A list of distances between shuffled TMDs and slip sites (one list for each iteration).
    """
    # Initialize storage for shuffled distances across iterations
    shuffled_distances_list = []
    ideal_distances_list = []

    for _ in range(num_iterations):
        # Shuffle the TMD positions
        shuffled_tmd_data = shuffle_tmd_positions(filtered_list, tmd_data)
        
        # Calculate distances for the shuffled TMDs
        shuffled_distances, num_ideal_distances = analyze_slipsite_to_tmd_distances_basic(filtered_list, shuffled_tmd_data, slipsite_list)
        
        # Store the distances for this iteration
        shuffled_distances_list.append(shuffled_distances) # Is it okay to append a list? This will now be a list of lists - is that okay?
        ideal_distances_list.append(num_ideal_distances)
    
    return shuffled_distances_list, ideal_distances_list



def compare_distances(real_distances, shuffled_distances_list):
    """
    Compare the real distances (absolute deviations from 45) to the shuffled distances and compute Z-scores or two-sided p-values.
    
    Args:
        real_distances: The list of distances (absolute deviations from 45 codons) between real TMDs and slip sites.
        shuffled_distances_list: A list of distance lists, one for each shuffling iteration.
    
    Returns:
        z_score: Z-score comparing the real data to the shuffled data.
        p_value: Two-sided p-value indicating whether the real distances are significantly different from shuffled.
    """
    # Compute the mean distance for the real data
    real_mean_distance = sum(real_distances) / len(real_distances)
    
    # Compute the mean distance for each shuffled iteration
    shuffled_mean_distances = [sum(distances) / len(distances) for distances in shuffled_distances_list]
    
    # Compute Z-score comparing real distances to shuffled distances
    z_score = zscore([real_mean_distance] + shuffled_mean_distances)[0]
    p_value = 2 * norm.sf(abs(z_score))
    
    # Compute two-sided p-value using permutation test
    extreme_count = sum(1 for d in shuffled_mean_distances if abs(d) >= abs(real_mean_distance))
    p_value_permutation = extreme_count / len(shuffled_mean_distances) # Always gives a value of zero; not enough precision

    return z_score, p_value


def compare_tmd_shuffle_to_observed(filtered_list, tmd_data, slipsite_list, num_iterations=100):
    """
    Compare real slipsite-to-TMD distances to shuffled TMD distances and compute Z-scores or p-values.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary containing TMD positions for each transcript.
        slipsite_list: List of known slippery sequences.
        num_iterations: Number of shuffling iterations to run (default is 1000).
    
    Returns:
        z_score: Z-score comparing real distances to shuffled distances.
        p_value: Two-sided p-value indicating whether real distances are significantly different from shuffled.
    """
    
    # Step 1: Calculate real distances between slip sites and TMDs
    print("Calculating real slipsite-to-TMD distances...")
    real_distances, num_ideal_distances_real = analyze_slipsite_to_tmd_distances_basic(filtered_list, tmd_data, slipsite_list)
    
    # Step 2: Run the TMD shuffling test to generate shuffled distances
    print(f"Running TMD shuffling test with {num_iterations} iterations...")
    shuffled_distances_list, num_ideal_distances_shuffled = run_tmd_shuffling_test_parallel(filtered_list, tmd_data, slipsite_list, num_iterations=128, num_processes=8)
    
    # Step 3: Compare real distances to shuffled distances using Z-score and p-value
    print("Comparing real distances to shuffled distances...")
    z_score, p_value = compare_distances(real_distances, shuffled_distances_list)
    
    return z_score, p_value, num_ideal_distances_real, num_ideal_distances_shuffled



## Random heptamer distance test


def compare_distances_heptamers(real_distances, sampled_distances_list):
    """
    Compare the real distances to the sampled distances and compute Z-scores or two-sided p-values.
    
    Args:
        real_distances: The list of real distances between TMDs and slippery sequences.
        sampled_distances_list: A list of distance lists, one for each sampling iteration.
    
    Returns:
        z_score: Z-score comparing real data to sampled data.
        p_value: Two-sided p-value indicating whether real distances are significantly different from sampled.
    """
    # Compute the mean distance for the real data
    real_mean_distance = sum(real_distances) / len(real_distances)
    
    # Compute the mean distance for each sampling iteration
    sampled_mean_distances = [sum(distances) / len(distances) for distances in sampled_distances_list]
    
    # Compute Z-score comparing real distances to sampled distances
    z_score = zscore([real_mean_distance] + sampled_mean_distances)[0]
    p_value = 2 * norm.sf(abs(z_score))
    
    # Compute two-sided p-value using permutation test
    extreme_count = sum(1 for d in sampled_mean_distances if abs(d) >= abs(real_mean_distance))
    # No matter how many iterations I performed this p_value was 0 so I implemented a different p-value calculation based on the Z score
    p_value_permutation = extreme_count / len(sampled_mean_distances)

    return z_score, p_value

def process_iteration_for_heptamer_sampling(args):
    filtered_list, tmd_data, slipsite_list, shuffled_list_of_heptamers = args
    num_slipsite_samples = len(slipsite_list)
    
    # Calculate heptamer frequencies and generate the list of heptamers
    #_, _, shuffled_list_of_heptamers = calculate_heptamer_frequencies_prf(filtered_list, tmd_data)
    
    # Randomly sample heptamers in proportion to their frequency
    random_heptamers = draw_random_heptamers(num_slipsite_samples, shuffled_list_of_heptamers)
    
    # Calculate distances for the randomly sampled heptamers
    sampled_heptamer_distances, num_ideal_distances = analyze_slipsite_to_tmd_distances_basic(filtered_list, tmd_data, random_heptamers)
    
    return sampled_heptamer_distances, num_ideal_distances


def run_heptamer_sampling_test_parallel(filtered_list, tmd_data, slipsite_list, shuffled_heptamer_list, num_iterations=128, num_processes=8):
    """
    Perform a heptamer sampling test in parallel and compare their distances to TMDs.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary containing TMD positions for each transcript.
        slipsite_list: List of known slippery sequences.
        num_iterations: Number of sampling iterations to run (default is 1000).
        num_processes: Number of parallel processes (default is 4).
    
    Returns:
        sampled_heptamer_distances: A list of distances between sampled heptamers and TMDs (one list for each iteration).
    """
    
    # Prepare the arguments for each process
    args = [(filtered_list, tmd_data, slipsite_list, shuffled_heptamer_list) for _ in range(num_iterations)]

    # Use multiprocessing Pool to parallelize the iterations
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_iteration_for_heptamer_sampling, args)

    # Unpack the results
    sampled_heptamer_distances_list = []
    num_ideal_distances_list = []

    for sampled_heptamer_distances, num_ideal_distances in results:
        sampled_heptamer_distances_list.append(sampled_heptamer_distances)
        num_ideal_distances_list.append(num_ideal_distances)


    return sampled_heptamer_distances_list, num_ideal_distances_list


def compare_heptamer_sampling_to_observed(filtered_list, tmd_data, slipsite_list, shuffled_heptamer_list, num_iterations=128):
    """
    Compare real slipsite-to-TMD distances to distances from random heptamer sampling.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary containing TMD positions for each transcript.
        slipsite_list: List of known slippery sequences.
        num_iterations: Number of sampling iterations to run (default is 1000).
    
    Returns:
        z_score: Z-score comparing real distances to sampled distances.
        p_value: P-value indicating whether the real distances are significantly different from sampled distances.
    """
    
    # Step 1: Calculate real distances between slip sites and TMDs
    print("Calculating real slipsite-to-TMD distances...")
    real_distances, num_ideal_distances_real = analyze_slipsite_to_tmd_distances_basic(filtered_list, tmd_data, slipsite_list)
    
    # Step 2: Run the heptamer sampling test to generate random distances
    print(f"Running heptamer sampling test with {num_iterations} iterations...")
    sampled_heptamer_distances_list, num_ideal_distances_sampled = run_heptamer_sampling_test_parallel(filtered_list, tmd_data, slipsite_list, shuffled_heptamer_list, num_iterations=128, num_processes=8)
    
    # Step 3: Compare real distances to random sampled heptamer distances using Z-score and p-value
    print("Comparing real distances to sampled heptamer distances...")
    z_score, p_value = compare_distances(real_distances, sampled_heptamer_distances_list)
    
    return z_score, p_value, num_ideal_distances_real, num_ideal_distances_sampled



def run_new_statistical_tests(harrington_motifs_with_uniprot, tmd_data, filtered_list, slipsite_list, num_trials_all, num_trials_canonical, num_trials_alternative, output_file, num_processes=8):

    import csv

    strict_slip_sites = ['AAAAAAA', 'AAAAAAT', 'AAAAAAC', 'AAATTTA', 'AAATTTT', 'AAATTTC',
                         'TTTAAAA', 'TTTAAAT', 'TTTAAAC', 'TTTTTTA', 'TTTTTTT', 'TTTTTTC',
                         'CCCAAAA', 'CCCAAAT', 'CCCAAAC', 'CCCTTTA', 'CCCTTTT', 'CCCTTTC',
                         'GGGAAAA', 'GGGAAAT', 'GGGAAAC', 'GGGTTTA', 'GGGTTTT', 'GGGTTTC']

    # Get subsets of the Harrington motif dataset
    [strict_all_motifs, 
     strict_canonical_motifs, 
     strict_alternative_motifs, 
     loose_all_motifs, 
     loose_canonical_motifs, 
     loose_alternative_motifs] = split_up_harrington_motifs(harrington_motifs_with_uniprot, 20, False)

    # Get the subset of transcripts that are canonical only and alternative only
    filtered_list_canonical_only, filtered_list_alternative_only = filter_non_canonical_transcripts(filtered_list)

    # Fetch heptamer frequencies for the transcriptome
    heptamer_frequencies_all, set_unique_heptamers_all, shuffled_list_of_heptamers_all = calculate_heptamer_frequencies_prf(filtered_list, tmd_data)
    heptamer_frequencies_canonical, set_unique_heptamers_canonical, shuffled_list_of_heptamers_canonical = calculate_heptamer_frequencies_prf(filtered_list_canonical_only, tmd_data)
    heptamer_frequencies_alternative, set_unique_heptamers_alternative, shuffled_list_of_heptamers_alternative = calculate_heptamer_frequencies_prf(filtered_list_alternative_only, tmd_data)

    # Convert heptamer sets into lists for use in chi squared test
    all_heptamers_all = list(set_unique_heptamers_all)
    all_heptamers_canonical = list(set_unique_heptamers_canonical)
    all_heptamers_alternative = list(set_unique_heptamers_alternative)
    
    # Calculate p_null for each combination of heptamer frequencies and slip site lists
    p_null_strict_all = calculate_p_null_slippery(heptamer_frequencies_all, strict_slip_sites)
    p_null_strict_canonical = calculate_p_null_slippery(heptamer_frequencies_canonical, strict_slip_sites)
    p_null_strict_alternative = calculate_p_null_slippery(heptamer_frequencies_alternative, strict_slip_sites)

    p_null_loose_all = calculate_p_null_slippery(heptamer_frequencies_all, slipsite_list)
    p_null_loose_canonical = calculate_p_null_slippery(heptamer_frequencies_canonical, slipsite_list)
    p_null_loose_alternative = calculate_p_null_slippery(heptamer_frequencies_alternative, slipsite_list)

    # Initialize a list to store results
    results = []

    # Binomial tests - NB: These binomial tests are going to have p-values of 1 since the null distribution and the test distribution are the same due to me changing how this test was done. 
    # A new and better binomial test in which the null distribution is random heptamer-TMD distances and the test distribution is slippery heptamer - TMD distances was implemented in another part of the script.
    binom_results = [
        ('strict_all', run_binomial_test(strict_all_motifs, p_null_strict_all, num_trials_all)),
        ('strict_canonical', run_binomial_test(strict_canonical_motifs, p_null_strict_canonical, num_trials_canonical)),
        ('strict_alternative', run_binomial_test(strict_alternative_motifs, p_null_strict_alternative, num_trials_alternative)),
        ('loose_all', run_binomial_test(loose_all_motifs, p_null_loose_all, num_trials_all)),
        ('loose_canonical', run_binomial_test(loose_canonical_motifs, p_null_loose_canonical, num_trials_canonical)),
        ('loose_alternative', run_binomial_test(loose_alternative_motifs, p_null_loose_alternative, num_trials_alternative))
    ]

    for category, (binom_p_value, binom_observed_count) in binom_results:
        results.append({
            'test_type': 'Binomial - do not use this version of binomial test (null distribution calculation was changed to heptamer frequencies drawn from real transcripts which will result in p-value of 1); remains for legacy reasons only',
            'category': category,
            'p_value': binom_p_value,
            'observed_count': binom_observed_count
        })


    # Initialize dictionaries to store ideal distances
    num_ideal_distances_real_tmd = {}
    num_ideal_distances_shuffled_tmd = {}

    # TMD Shuffling Test
    for category, data in [('strict_all', (filtered_list, tmd_data, strict_slip_sites)),
                           ('strict_canonical', (filtered_list_canonical_only, tmd_data, strict_slip_sites)),
                           ('strict_alternative', (filtered_list_alternative_only, tmd_data, strict_slip_sites)),
                           ('loose_all', (filtered_list, tmd_data, slipsite_list)),
                           ('loose_canonical', (filtered_list_canonical_only, tmd_data, slipsite_list)),
                           ('loose_alternative', (filtered_list_alternative_only, tmd_data, slipsite_list))]:
        z_score, p_value, num_ideal_distances_real, num_ideal_distances_shuffled = compare_tmd_shuffle_to_observed(
            *data, num_iterations=1024)
        results.append({
            'test_type': 'TMD Shuffling',
            'category': category,
            'z_score': z_score,
            'p_value': p_value
        })
        num_ideal_distances_real_tmd[category] = num_ideal_distances_real
        num_ideal_distances_shuffled_tmd[category] = num_ideal_distances_shuffled


    # Initialize dictionaries to store ideal distances
    num_ideal_distances_real_heptamer = {}
    num_ideal_distances_sampled_heptamer = {}

    # Heptamer Sampling Test
    for category, data, shuffled_heptamer_list in [
            ('strict_all', (filtered_list, tmd_data, strict_slip_sites), shuffled_list_of_heptamers_all),
            ('strict_canonical', (filtered_list_canonical_only, tmd_data, strict_slip_sites), shuffled_list_of_heptamers_canonical),
            ('strict_alternative', (filtered_list_alternative_only, tmd_data, strict_slip_sites), shuffled_list_of_heptamers_alternative),
            ('loose_all', (filtered_list, tmd_data, slipsite_list), shuffled_list_of_heptamers_all),
            ('loose_canonical', (filtered_list_canonical_only, tmd_data, slipsite_list), shuffled_list_of_heptamers_canonical),
            ('loose_alternative', (filtered_list_alternative_only, tmd_data, slipsite_list), shuffled_list_of_heptamers_alternative)
        ]:
        z_score, p_value, num_ideal_distances_real, num_ideal_distances_sampled = compare_heptamer_sampling_to_observed(
            *data, shuffled_heptamer_list, num_iterations=1024)
        results.append({
            'test_type': 'Random Heptamer Sampling',
            'category': category,
            'z_score': z_score,
            'p_value': p_value
        })
        num_ideal_distances_real_heptamer[category] = num_ideal_distances_real
        num_ideal_distances_sampled_heptamer[category] = num_ideal_distances_sampled

    # Chi-Squared Test
    chi_squared_results = []

    for category, data, slip_sites, all_heptamers in [
            ('strict_all', filtered_list, strict_slip_sites, all_heptamers_all),
            ('strict_canonical', filtered_list_canonical_only, strict_slip_sites, all_heptamers_canonical),
            ('strict_alternative', filtered_list_alternative_only, strict_slip_sites, all_heptamers_alternative),
            ('loose_all', filtered_list, slipsite_list, all_heptamers_all),
            ('loose_canonical', filtered_list_canonical_only, slipsite_list, all_heptamers_canonical),
            ('loose_alternative', filtered_list_alternative_only, slipsite_list, all_heptamers_alternative)
        ]:
        counts = analyze_slipsite_to_tmd_distances_contingency_table(
            data, tmd_data, slip_sites, all_heptamers)
        
        stats_results = perform_statistical_analysis(counts)
        
        chi_squared_results.append({
            'test_type': 'Chi-Squared',
            'category': category,
            'p_value': stats_results['p_value'],
            'chi2_stat': stats_results['chi2_stat'],
            'degrees_of_freedom': stats_results['dof'],
            'odds_ratio': stats_results['odds_ratio'],
            'ci_lower': stats_results['ci_lower'],
            'ci_upper': stats_results['ci_upper'],
            'A': counts['slippery_ideal'],       # Adding count A
            'B': counts['non_slippery_ideal'],   # Adding count B
            'C': counts['slippery_other'],       # Adding count C
            'D': counts['non_slippery_other']    # Adding count D
        })

    # Write to CSV
    fieldnames = ['test_type', 'category', 'p_value', 'observed_count', 'z_score', 'chi2_stat', 'degrees_of_freedom', 'odds_ratio', 'ci_lower', 'ci_upper', 'A', 'B', 'C', 'D']
    with open(output_file+'_chisquared_test_output_results.tsv', 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()
        for result in chi_squared_results:
            # Fill missing fields with None
            for field in fieldnames:
                if field not in result:
                    result[field] = None
            writer.writerow(result)


    print(f"Chi squared test results written to {output_file+'_chisquared_test_output_results.tsv'}")
    # Additional logging or saving of intermediate data can be added here



    categories = ['strict_all', 'strict_canonical', 'strict_alternative', 
                  'loose_all', 'loose_canonical', 'loose_alternative']

    # Determine the number of iterations
    num_iterations_tmd_shuffle = len(next(iter(num_ideal_distances_shuffled_tmd.values())))

    # Write to CSV
    with open(output_file+'_tmd_shuffling_ideal_distances.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Iteration'] + categories)
        for i in range(num_iterations_tmd_shuffle):
            row = [i + 1]
            for category in categories:
                row.append(num_ideal_distances_shuffled_tmd[category][i])
            writer.writerow(row)

    # Determine the number of iterations
    num_iterations_heptamer_sampling = len(next(iter(num_ideal_distances_sampled_heptamer.values())))

    # Write to CSV
    with open(output_file+'_heptamer_sampling_ideal_distances.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Iteration'] + categories)
        for i in range(num_iterations_heptamer_sampling):
            row = [i + 1]
            for category in categories:
                row.append(num_ideal_distances_sampled_heptamer[category][i])
            writer.writerow(row)

    with open(output_file+'_ideal_distances_real_log.txt', 'w') as logfile:
        logfile.write('Ideal distances in unshuffled data (TMD Shuffling Test):\n')
        for category, count in num_ideal_distances_real_tmd.items():
            logfile.write(f"{category}: {count}\n")
        logfile.write('\nIdeal distances in unshuffled data (Heptamer Sampling Test):\n')
        for category, count in num_ideal_distances_real_heptamer.items():
            logfile.write(f"{category}: {count}\n")


    with open(output_file+'_other_new_statistical_tests_categories.tsv', 'w', newline='') as tsvfile:
        fieldnames = ['test_type', 'category', 'p_value', 'observed_count', 'z_score']
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()
        for result in results:
            writer.writerow(result)

    print(f"Results written to {output_file}")
    print("Ideal distances data saved to CSV and log files.")


    

    return shuffled_list_of_heptamers_all, shuffled_list_of_heptamers_canonical, shuffled_list_of_heptamers_alternative



def perform_statistical_analysis(counts):
    """
    Performs chi-squared test and calculates odds ratio with confidence intervals from counts.

    Args:
        counts (dict): Dictionary with counts in the contingency table.

    Returns:
        dict: Dictionary with chi-squared statistic, p-value, degrees of freedom, odds ratio, 95% confidence interval.
    """
    import numpy as np
    import scipy.stats as stats
    from statsmodels.stats.contingency_tables import Table2x2

    A = counts['slippery_ideal']
    B = counts['non_slippery_ideal']
    C = counts['slippery_other']
    D = counts['non_slippery_other']

    # Continuity correction to avoid divide-by-zero error
    if A == 0 or B == 0 or C == 0 or D == 0:
        A += 0.5
        B += 0.5
        C += 0.5
        D += 0.5

    # Construct contingency table
    contingency_table = np.array([[A, B],
                                  [C, D]])

    # Perform chi-squared test without Yates' correction
    chi2_stat, p_value, dof, expected = stats.chi2_contingency(contingency_table, correction=False)
    
    # Calculate odds ratio and 95% confidence interval
    table2x2 = Table2x2(contingency_table)
    odds_ratio = table2x2.oddsratio
    ci_low, ci_upp = table2x2.oddsratio_confint()
    #p_value_oddsratio = table2x2.test_oddsratio_null().pvalue
    p_value_oddsratio = table2x2.oddsratio_pvalue()

    results = {
        'chi2_stat': chi2_stat,
        'p_value': p_value,
        'dof': dof,
        'odds_ratio': odds_ratio,
        'ci_lower': ci_low,
        'ci_upper': ci_upp,
        'p_value_oddsratio': p_value_oddsratio
    }

    return results



def check_duplicate_columns(excel_path):
    df = pd.read_excel(excel_path)
    column_names = df.columns
    duplicate_columns = column_names[column_names.duplicated()].unique()
    print(f"Column names: ", column_names)
    if len(duplicate_columns) > 0:
        print(f"Duplicate columns found: {duplicate_columns}")
    else:
        print("No duplicate columns found.")


def run_heptamer_sampling_test_serial(filtered_list, tmd_data, slipsite_list, num_iterations=100):
    """
    Perform a heptamer sampling test by randomly sampling 7N sequences and comparing their distances to TMDs.
    
    Args:
        filtered_list: List of labeled transcripts with sequences.
        tmd_data: Dictionary containing TMD positions for each transcript.
        slipsite_list: List of known slippery sequences.
        num_iterations: Number of sampling iterations to run (default is 1000).
    
    Returns:
        shuffled_distances: A list of distances between sampled heptamers and TMDs (one list for each iteration).
    """
    # Get the total number of slippery sequences to sample
    num_slipsite_samples = len(slipsite_list)
    
    # Calculate heptamer frequencies and generate the list of heptamers
    _, _, shuffled_list_of_heptamers = calculate_heptamer_frequencies_prf(filtered_list, tmd_data)
    
    # Initialize storage for sampled heptamer distances across iterations
    sampled_heptamer_distances_list = []

    # Initialize list of ideal distances across iterations
    num_ideal_distances_list = []

    for _ in range(num_iterations):
        # Randomly sample heptamers in proportion to their frequency
        random_heptamers = draw_random_heptamers(num_slipsite_samples, shuffled_list_of_heptamers)
        
        # Calculate distances for the randomly sampled heptamers
        sampled_heptamer_distances, num_ideal_distances = analyze_slipsite_to_tmd_distances_basic(filtered_list, tmd_data, random_heptamers)
        
        # Store the distances for this iteration
        sampled_heptamer_distances_list.append(sampled_heptamer_distances)

        # Store the number of ideal distances in this iteration
        num_ideal_distances_list.append(num_ideal_distances)
    
    return sampled_heptamer_distances_list, num_ideal_distances_list



def configparser_function(config_file="config.ini"):
    """
    Loads configuration from an ini file using configparser.
    """
    config = configparser.ConfigParser()
    config.read(config_file)
    
    # Paths
    path_fasta = config.get('Paths', 'path_fasta')
    path_uniprot = config.get('Paths', 'path_uniprot')
    blastp_results_path = config.get('Paths', 'blastp_results_path')
    path_tmd_csv = config.get('Paths', 'path_tmd_csv')
    path_slipsites = config.get('Paths', 'path_slipsites')
    path_friction_csv = config.get('Paths', 'path_friction_csv')
    
    # Parameters
    gap_near = config.getint('Parameters', 'gap_near')
    gap_far = config.getint('Parameters', 'gap_far')
    outname = config.get('Parameters', 'outname')
    scoring_sys = config.getint('Parameters', 'scoring_sys')  # Will return 1 or 2
    
    # Return all values in the same order as interactive_configuration
    return path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, gap_near, gap_far, outname, scoring_sys


# Main configuration function
def interactive_configuration():
    path_fasta = input("Please enter a path to an Ensembl-format CDS database: ")
    path_uniprot = input("Please enter a path to an Excel sheet containing Uniprot sequence data: ")
    blastp_results_path = input("Please enter a path to the output of the blastp search matching Ensembl transcripts to Uniprot sequences: ")
    path_tmd_csv = input("Please enter a path to the von Heijne optimized TMD data CSV: ")
    path_slipsites = input("Please enter a path to a list of slip sites to search for: ")
    path_friction_csv = input("Please enter a path to the friction scores CSV: ")
    gap_near = int(input("Please enter a lower bound on the distance range in the sequence upstream from the slippery sequence to be searched (in codons): "))
    gap_far = int(input("Please enter an upper bound on the distance range in the sequence upstream from the slippery sequence to be searched (in codons): "))
    outname = input("Please enter a name for the output: ")
    while True:
        print("Which scoring system do you want to use")
        print("[1] Use number of hydrogen-bonds formed and broken (simple scoring)")
        print("[2] Use Turner nearest neighbor parameters")
        scoring_sys = int(input("Please enter [1] or [2]: "))
        if scoring_sys in (1, 2):
            break
        else:
            print("There has been an error. Please re-enter.")
    return path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, gap_near, gap_far, outname, scoring_sys



def main():
    
    if len(sys.argv) < 2:
        print("Usage: python script.py <config_file.ini>")
        config_mode_flag = int(input("Do you ant to revert to interactive mode? [1] yes or [2] no"))
        if config_mode_flag == 1:
            path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, gap_near, gap_far, outname, scoring_sys = interactive_configuration()
        else:
            sys.exit(1)
    
    config_file = sys.argv[1]

    path_fasta, path_uniprot, blastp_results_path, path_tmd_csv, path_slipsites, path_friction_csv, gap_near, gap_far, outname, scoring_sys = configparser_function(config_file)

    fasta_list = parse_fasta(path_fasta)
    uniprot_dict = load_uniprot_excel_to_dict(path_uniprot)
    stop_filtered_fasta_list = stop_codon_filter(fasta_list)
    slipsite_list = import_text_file(path_slipsites)
    tmd_data = parse_tmd_csv(path_tmd_csv)
    friction_data = parse_friction_csv(path_friction_csv)

    print("Processing transcripts and matching to Uniprot sequences...")
    filtered_list, canonical_count, alt_spliced_count = fasta_transcript_filter(stop_filtered_fasta_list)
    filtered_list_canonical_only, filtered_list_alternative_only = filter_non_canonical_transcripts(filtered_list)
    ensembl_uniprot_dict = process_transcripts_with_uniprot(filtered_list, blastp_results_path, tmd_data, outname)
    print_ensembl_uniprot_dict(ensembl_uniprot_dict, outname)

    # Perform the Harrington motif search
    print("Search for Harrington motifs...")
    harrington_motifs, motifs_within_gap, transcript_codons_dict, num_trials_all, num_trials_canonical, num_trials_alternative = motif_search(filtered_list, tmd_data, slipsite_list, friction_data, gap_near, gap_far, scoring_sys)
    
    print("Associaing Harrington motifs with Uniprot sequences")
    harrington_motifs_with_uniprot = append_uniprot_accession_to_motif(harrington_motifs, ensembl_uniprot_dict)


    print("Running new stats tests including Chi Squared test")
    # Get data structures for new statistical analyses
    # Note that the heptamer frequency function uses only canonical sequences to avoid skewing the distribution of heptamers
    # with Ensembl genes with a lot of alternative isoform transcripts
    shuffled_list_of_heptamers_all, shuffled_list_of_heptamers_canonical, shuffled_list_of_heptamers_alternative = run_new_statistical_tests(harrington_motifs_with_uniprot, tmd_data, filtered_list, slipsite_list, num_trials_all, num_trials_canonical, num_trials_alternative, outname+'_new_stats_tests')

    print("Finished running bank of new tests!")

    # This prints a log file with some basic information on the shape of the data
    output_list = get_canonical_count(filtered_list, tmd_data, harrington_motifs, outname)

    # These little stats were an early effort to get some sense of the shape of the data - feel free to comment out
    X1, X2, X3, X4 = summarize_transcripts(tmd_data, canonical_count, alt_spliced_count, harrington_motifs)
    X5, X6 = calculate_statistics(harrington_motifs, sum([len(codons) for codons in transcript_codons_dict.values()]), gap_near, gap_far)

    # This prints out the list of motifs identified
    print_motif_csv(harrington_motifs_with_uniprot, outname)
    
    print(f"Total canonical transcripts: {canonical_count}")
    print(f"Total alternatively-spliced transcripts: {alt_spliced_count}")
    print(f"Number of motifs within the gap range ({gap_near}-{gap_far} codons): {motifs_within_gap}")
    print(f"X1: {X1}, X2: {X2}, X3: {X3}, X4: {X4}")
    print(f"Basic X5: {X5}-fold difference from random, Basic X6: p-value = {X6}")

    print("Obtaining TMD--slippery sequence distance distribution...")
    # Analyze slipsite to TMD distances and categorize by transcript type
    distance_data = analyze_slipsite_to_tmd_distances(filtered_list, tmd_data, slipsite_list, outname)

    # Generate and plot actual distances
    histogram_output_file = f"{outname}_slipsite_tmd_distance_histogram.png"
    plot_slipsite_tmd_distances_with_categories(distance_data, histogram_output_file)

    # Generate background distribution
    print("Obtaining TMD--random heptamer background distribution...")
    background_distances = generate_background_distribution(filtered_list, tmd_data, slipsite_list, shuffled_list_of_heptamers_all, shuffled_list_of_heptamers_canonical, shuffled_list_of_heptamers_alternative, num_iterations=1024, batch_size=128, num_processors=8)

    print("TMD--random heptamer background distribution generated")

    # Compare distributions for all categories
    background_file_prefix = "background_distances"
    output_file_prefix = f"{outname}_distance_distribution_comparison"

    compare_distributions(distance_data, background_file_prefix, output_file_prefix)

    print("Analysis complete!")
    
if __name__ == '__main__':

    main()
