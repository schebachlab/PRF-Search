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
import multiprocessing as mp
import math
from scipy import stats
from scipy.stats import ks_2samp
from scipy.stats import binom_test
from scipy.stats import mannwhitneyu
from scipy.stats import poisson
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from itertools import groupby
import matplotlib.pyplot as plt
import frameshift_routines_v10 as frameshift
import time
import configparser
from Bio import Align
from Bio.Align import substitution_matrices
from Bio import SeqIO
import difflib

from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy

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


def parse_motif_csv_notdone(path_friction_csv):
    friction_data = {}
    with open(path_friction_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            heptamer = row['sequence']
            friction_score = float(row['friction_score'])
            friction_data[heptamer] = friction_score
    return friction_data


def import_csv(path_csv):
    file_obj = open(path_csv, "r")
    data_raw = []
    for line in file_obj:
        line = line.rstrip()
        line = line.split(',')
        data_raw.append(line)
    return data_raw


def import_tsv(path_tsv):
    file_obj = open(path_tsv, "r")
    data_raw = []
    for line in file_obj:
        line = line.rstrip()
        line = line.split('\t')
        data_raw.append(line)
    return data_raw


def fetch_membrane_proteins_lists(uniprot_data_raw):
    key_word_membrane_list = []
    key_word_transmembrane_helix_list = []
    all_human_reviewed_list = []
    disease_association_list = []

    uniprot_data_list = uniprot_data_raw[1:]
    for entry in uniprot_data_list:
        uniprot_accession_id = entry[0]
        uniprot_keywords_string = entry[12]
        uniprot_keywords_list = uniprot_keywords_string.split(';')
        
        all_human_reviewed_list.append(uniprot_accession_id)

        if 'Membrane' in uniprot_keywords_list:
            key_word_membrane_list.append(uniprot_accession_id)
        if 'Transmembrane helix' in uniprot_keywords_list:
            key_word_transmembrane_helix_list.append(uniprot_accession_id)
        if 'Disease variant' in uniprot_keywords_list:
            disease_association_list.append(uniprot_accession_id)

    return key_word_membrane_list, key_word_transmembrane_helix_list, disease_association_list, all_human_reviewed_list


def parse_go_terms(uniprot_data_raw):
    # uniprot_headers = uniprot_data_raw[0]
    # 0 Entry
    # 1 Reviewed
    # 2 Entry Name
    # 3 Protein names
    # 4 Gene Names
    # 5 Organism
    # 6 Length
    # 7 Gene Ontology (biological process)
    # 8 Gene Ontology (cellular component)
    # 9 Involvement in disease
    # 10 Transmembrane
    # 11 Gene Ontology (molecular function)
    # 12 Keywords

    # Unique_go_terms
    go_term_description_dict_all = {}
    go_term_description_dict_biological_process = {}
    go_term_description_dict_cellular_component = {}
    go_term_description_dict_molecular_function = {}

    # Gene:GO term dictionaries
    gene_ontology_biological_process_dict = {}
    gene_ontology_cellular_component_dict = {}
    gene_ontology_molecular_function_dict = {}

    # Other output dictionaries
    uniprot_keywords_dict = {}
    disease_involvement_dict = {}
    protein_name_dict = {}

    # Loop to iterate through uniprot items
    uniprot_data_list = uniprot_data_raw[1:]
    for entry in uniprot_data_list:

        # Fetch basic terms from list
        uniprot_accession_id = entry[0]
        uniprot_entry_name = entry[2]
        descriptive_protein_name = entry[3]
        gene_names_string = entry[4]
        gene_names_list = gene_names_string.split(' ')

        # Define regex pattern to separate GO terms from their descriptions
        pattern = r'^(.*)\s+\[(GO:\d+)\]$'

        # Process biological process GO terms
        go_biological_process_string = entry[7]
        go_biological_process_raw_list = go_biological_process_string.split(';')
        go_biological_process_list = [x.strip() for x in go_biological_process_raw_list]
        go_ids_biological_process = []
        go_descriptions_biological_process = []
        for go_string in go_biological_process_list:
            match_object = re.search(pattern, go_string)
            if match_object:
                go_description = match_object.group(1)
                go_descriptions_biological_process.append(go_description)
                go_term_id = match_object.group(2)
                go_ids_biological_process.append(go_term_id)
                # Add GO terms and their descriptions to the dictionary for all GO terms
                if go_term_id not in go_term_description_dict_all:
                    go_term_description_dict_all[go_term_id] = []
                    go_term_description_dict_all[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_all[go_term_id]:
                        go_term_description_dict_all[go_term_id].append(go_description)
                # Add GO terms and their descriptions to the dictionary for biological process GO terms
                if go_term_id not in go_term_description_dict_biological_process:
                    go_term_description_dict_biological_process[go_term_id] = []
                    go_term_description_dict_biological_process[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_biological_process[go_term_id]:
                        go_term_description_dict_biological_process[go_term_id].append(go_description)



        # These dictionaries will return a list: index 0 original string, index 1: go term ids, index 2: descriptions of go terms
        gene_ontology_biological_process_dict[uniprot_accession_id] = [go_biological_process_list, go_ids_biological_process, go_descriptions_biological_process]


        # Process cellular component GO terms
        go_cellular_component_string = entry[8]
        go_cellular_component_raw_list = go_cellular_component_string.split(';')
        go_cellular_component_list = [x.strip() for x in go_cellular_component_raw_list]
        go_ids_cellular_component = []
        go_descriptions_cellular_component = []
        for go_string in go_cellular_component_list:
            match_object = re.search(pattern, go_string)
            if match_object:
                go_description = match_object.group(1)
                go_descriptions_cellular_component.append(go_description)
                go_term_id = match_object.group(2)
                go_ids_cellular_component.append(go_term_id)
                # Add GO terms and their descriptions to the dictionary for all terms
                if go_term_id not in go_term_description_dict_all:
                    go_term_description_dict_all[go_term_id] = []
                    go_term_description_dict_all[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_all[go_term_id]:
                        go_term_description_dict_all[go_term_id].append(go_description)
                # Add GO terms and their descriptions to the dictionary for cellular component GO terms
                if go_term_id not in go_term_description_dict_cellular_component:
                    go_term_description_dict_cellular_component[go_term_id] = []
                    go_term_description_dict_cellular_component[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_cellular_component[go_term_id]:
                        go_term_description_dict_cellular_component[go_term_id].append(go_description)
        # These dictionaries will return a list: index 0 original string, index 1: go term ids, index 2: descriptions of go terms
        gene_ontology_cellular_component_dict[uniprot_accession_id] = [go_cellular_component_list, go_ids_cellular_component, go_descriptions_cellular_component]


        # Process molecular function GO terms
        go_molecular_function_string = entry[11]
        go_molecular_function_raw_list = go_molecular_function_string.split(';')
        go_molecular_function_list = [x.strip() for x in go_molecular_function_raw_list]
        go_ids_molecular_function = []
        go_descriptions_molecular_function = []
        for go_string in go_molecular_function_list:
            match_object = re.search(pattern, go_string)
            if match_object:
                go_description = match_object.group(1)
                go_descriptions_molecular_function.append(go_description)
                go_term_id = match_object.group(2)
                go_ids_molecular_function.append(go_term_id)
                # Add GO terms and their descriptions to the dictionary for all terms
                if go_term_id not in go_term_description_dict_all:
                    go_term_description_dict_all[go_term_id] = []
                    go_term_description_dict_all[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_all[go_term_id]:
                        go_term_description_dict_all[go_term_id].append(go_description)
                # Add GO terms and their descriptions to the dictioanry for molecular function GO terms
                if go_term_id not in go_term_description_dict_molecular_function:
                    go_term_description_dict_molecular_function[go_term_id] = []
                    go_term_description_dict_molecular_function[go_term_id].append(go_description)
                else:
                    if go_description not in go_term_description_dict_molecular_function[go_term_id]:
                        go_term_description_dict_molecular_function[go_term_id].append(go_description)

        # These dictionaries will return a list: index 0 original string, index 1: go term ids, index 2: descriptions of go terms
        gene_ontology_molecular_function_dict[uniprot_accession_id] = [go_molecular_function_list, go_ids_molecular_function, go_descriptions_molecular_function]


        # Process Uniprot Keywords
        uniprot_keywords_string = entry[12]
        uniprot_keywords_list = uniprot_keywords_string.split(';')
        uniprot_keywords_dict[uniprot_accession_id] = uniprot_keywords_list

        # Process Uniprot annotations regarding disease involvement
        uniprot_disease_involvement_string = entry[9]
        uniprot_disease_involvement_list = uniprot_disease_involvement_string.split(';')

        if uniprot_disease_involvement_list[0] == '':
            num_associated_diseases = 'Number of associated diseases: 0'
            disease_involvement_output_entry = [num_associated_diseases, 'No disease involvement']
        else:
            num_associated_diseases = 'Number of associated diseases: ' + str(len(uniprot_disease_involvement_list))
            disease_involvement_output_entry = [num_associated_diseases]
            for entry in uniprot_disease_involvement_list:
                disease_involvement_output_entry.append(entry)

        disease_involvement_dict[uniprot_accession_id] = disease_involvement_output_entry

        # Process Uniprot names and protein descriptions
        protein_name_entry = [uniprot_entry_name, descriptive_protein_name, gene_names_string, gene_names_list]
        protein_name_dict[uniprot_accession_id] = protein_name_entry

    all_outputs = (go_term_description_dict_all,
                   go_term_description_dict_biological_process,
                   go_term_description_dict_cellular_component,
                   go_term_description_dict_molecular_function,
                   gene_ontology_biological_process_dict,
                   gene_ontology_cellular_component_dict,
                   gene_ontology_molecular_function_dict,
                   uniprot_keywords_dict,
                   disease_involvement_dict,
                   protein_name_dict)
    
    return all_outputs


def fetch_uniprot_id_from_motif_subset(harrington_motif_subset):

    subset_uniprot_ids_list = []
    for entry in harrington_motif_subset:
        uniprot_id = entry[28]
        subset_uniprot_ids_list.append(uniprot_id)
    
    subset_uniprot_ids_set = set(subset_uniprot_ids_list)
    subset_unique_uniprot_ids_list = list(subset_uniprot_ids_set)

    return subset_unique_uniprot_ids_list



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
        location_ideality = int(entry[13])
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


def parse_ensembl_uniprot_dict(path_ensembl_uniprot_dict):
    ensembl_uniprot_dict = {}
    with open(path_ensembl_uniprot_dict, newline = '') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            transcript_id = row['ensembl_transcript_id']
            uniprot_accession_id = row['uniprot_accession_id']
            ensembl_uniprot_dict[transcript_id] = uniprot_accession_id
    return ensembl_uniprot_dict


def map_all_predicted_tmd_transcripts_to_uniprot_ids(path_tmds='topcons_ensembl_output.csv', path_ensembl_uniprot_dict='motif_search_output_v32_full_blastp_fishers_ensembl_uniprot_dict_full_blastp.csv'):
    tmd_data = parse_tmd_csv(path_tmds)
    uniprot_ensembl_dict = parse_ensembl_uniprot_dict(path_ensembl_uniprot_dict)

    uniprot_accession_ids_with_tmds = []

    for ensembl_transcript_id, uniprot_accession_id in uniprot_ensembl_dict.items():
        if ensembl_transcript_id in tmd_data:
            if uniprot_accession_id not in uniprot_accession_ids_with_tmds:
                uniprot_accession_ids_with_tmds.append(uniprot_accession_id)
    return uniprot_accession_ids_with_tmds


def print_results(output_data, outname):
    myfile = open(outname+'.csv', 'w')
    for line in output_data:
        line_string = [str(x) for x in line]
        csv_line = ','.join(line_string)
        print(csv_line, file = myfile)
    myfile.close()


def get_go_counts(study_set, background_set, go_term_dict):
    """
    Get counts for GO term enrichment analysis for Fisher's Exact Test.
    """
    # GO term counts
    go_term_counts = {}
    
    # Iterate over background to count GO terms
    for uniprot_id in background_set:
        if uniprot_id in go_term_dict:
            for go_term in go_term_dict[uniprot_id][1]:  # Access GO term IDs from the dictionary
                if go_term not in go_term_counts:
                    go_term_counts[go_term] = {'a': 0, 'b': 0, 'c': 0, 'd': 0}
                if uniprot_id in study_set:
                    go_term_counts[go_term]['a'] += 1  # Study set with GO term
                else:
                    go_term_counts[go_term]['b'] += 1  # Background with GO term

    # Calculate c and d
    for go_term, counts in go_term_counts.items():
        counts['c'] = len(study_set) - counts['a']  # Study set without GO term
        counts['d'] = len(background_set) - len(study_set) - counts['b']  # Background without GO term

    return go_term_counts


def run_fisher_exact_test(go_term_counts, go_term_description_dict_all, go_category, gene_set):
    """
    Perform Fisher's Exact Test for each GO term and collect p-values.
    """
    p_values = []
    go_terms = []
    log_odds_ratios = []
    go_descriptions = []
    go_categories = []
    gene_set_list = []
    
    for go_term, counts in go_term_counts.items():
        # Build the contingency table for Fisher's Exact Test
        contingency_table = [[counts['a'], counts['b']],
                             [counts['c'], counts['d']]]
        
        # Perform Fisher's Exact Test
        odds_ratio, p_value = stats.fisher_exact(contingency_table, alternative='greater')

        # Compute log odds ratio and safely handle cases where odds_ratio is zero or undefined
        if odds_ratio > 0:
            log_odds_ratio = math.log(odds_ratio)
        else:
            log_odds_ratio = None

        #log_odds_ratio = log_odds_ratio_prior # Is this what I want?

        p_values.append(p_value)
        go_terms.append(go_term)
        log_odds_ratios.append(log_odds_ratio)
        go_description = go_term_description_dict_all[go_term]
        go_descriptions.append(go_description)
        go_categories.append(go_category)
        gene_set_list.append(gene_set)
    
    return go_terms, p_values, log_odds_ratios, go_descriptions, go_categories, gene_set_list


def apply_fdr_correction(p_values, alpha=0.05):
    """
    Apply FDR correction using the Benjamini-Hochberg method.
    """
    _, p_adjusted, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')
    return p_adjusted



def gene_ontology_enrichment_function(motif_subset_uniprot_list, background_uniprot_list, gene_ontology_dict, go_term_description_dict_all, go_category, gene_set):
    """
    Run a Fisher's exact test and p-value FDR correction using a subset of Uniprot accession IDs associated
    with Harrington motifs, an appropriate set of background genes (membrane proteins, all reviewed, etc),
    and a dictionary associating Uniprot IDs with GO terms
    """

    # Get the GO term counts for Fisher's exact test
    go_term_counts = get_go_counts(motif_subset_uniprot_list, background_uniprot_list, gene_ontology_dict)

    # Run Fisher's exact test for category of GO terms and subset of Harrington motifs
    go_terms, p_values, log_odds_ratios, go_descriptions, go_categories, gene_set_list = run_fisher_exact_test(go_term_counts, go_term_description_dict_all, go_category, gene_set)

    # Apply FDR correction
    p_adjusted = apply_fdr_correction(p_values)

    # Save the results as a list
    output_list = list(zip(go_terms, go_descriptions, log_odds_ratios, p_values, p_adjusted, go_categories, gene_set_list))

    return output_list



def main_go_enrichment_function(path_uniprot_tsv, path_harrington_motifs):

    # Load raw uniprot data
    uniprot_data_raw = import_tsv(path_uniprot_tsv)
    uniprot_data = uniprot_data_raw[1:]
    uniprot_go_dicts = parse_go_terms(uniprot_data_raw)

    # Unwrap GO dictionaries
    (go_term_description_dict_all, 
     go_term_description_dict_biological_process, 
     go_term_description_dict_cellular_component, 
     go_term_description_dict_molecular_function, 
     gene_ontology_biological_process_dict, 
     gene_ontology_cellular_component_dict, 
     gene_ontology_molecular_function_dict, 
     uniprot_keywords_dict, 
     disease_involvement_dict, 
     protein_name_dict) = uniprot_go_dicts
    
    
    # Load Harrington motif data
    harrington_motifs_raw = import_csv(path_harrington_motifs)
    harrington_motifs = harrington_motifs_raw[1:]

    [strict_all_motifs, 
     strict_canonical_motifs, 
     strict_alternative_motifs, 
     loose_all_motifs, 
     loose_canonical_motifs, 
     loose_alternative_motifs] = split_up_harrington_motifs(harrington_motifs, 20, False)
    
    # Get lists of Uniprot accession IDs for background
    (key_word_membrane_list, 
     key_word_transmembrane_helix_list_not_used, # Change this variable name so it is not passed to gene_ontology_enrichment_function
     disease_association_list, 
     all_human_reviewed_list) = fetch_membrane_proteins_lists(uniprot_data_raw)
    
    # Get Uniprot IDs of Harrington motif subsets
    strict_all_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(strict_all_motifs)
    strict_canonical_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(strict_canonical_motifs)
    strict_alternative_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(strict_alternative_motifs)
    loose_all_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(loose_all_motifs)
    loose_canonical_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(loose_canonical_motifs)
    loose_alternative_motifs_uniprot_ids = fetch_uniprot_id_from_motif_subset(loose_alternative_motifs)

    # Re-assign background list of Uniprot sequences
    uniprot_accession_ids_with_tmds = map_all_predicted_tmd_transcripts_to_uniprot_ids('topcons_ensembl_output.csv', 'motif_search_output_v32_full_blastp_fishers_ensembl_uniprot_dict_full_blastp.csv')
    key_word_transmembrane_helix_list = uniprot_accession_ids_with_tmds

    # Run Fisher's exact test on all six motif subsets 

    # Initialize a dictionary to hold DataFrames
    dataframes = {}

    # Define categories and GO categories
    categories = [
        ('Strict_All', strict_all_motifs_uniprot_ids),
        ('Strict_Canonical', strict_canonical_motifs_uniprot_ids),
        ('Strict_Alternative', strict_alternative_motifs_uniprot_ids),
        ('Loose_All', loose_all_motifs_uniprot_ids),
        ('Loose_Canonical', loose_canonical_motifs_uniprot_ids),
        ('Loose_Alternative', loose_alternative_motifs_uniprot_ids)
    ]

    go_categories = [
        ('Biological_Process', gene_ontology_biological_process_dict),
        ('Cellular_Component', gene_ontology_cellular_component_dict),
        ('Molecular_Function', gene_ontology_molecular_function_dict)
    ]

    # Loop over categories and GO categories to create DataFrames
    for gene_set_name, motif_uniprot_ids in categories:
        for go_category_name, go_dict in go_categories:
            output_list = gene_ontology_enrichment_function(
                motif_uniprot_ids, key_word_transmembrane_helix_list, go_dict, go_term_description_dict_all, go_category_name, gene_set_name
            )
            # Convert output_list to DataFrame
            df = pd.DataFrame(output_list, columns=['go_term', 'go_description', 'log_odds_ratio', 'p_value', 'p_adjusted', 'go_category', 'gene_set'])
            # Create a key for the dataframes dictionary
            key = f"{gene_set_name}_{go_category_name}"
            dataframes[key] = df

    # Write all DataFrames to an Excel file with separate sheets
    with pd.ExcelWriter('go_enrichment_results_no_stemloop_required.xlsx') as writer:
        for sheet_name, df in dataframes.items():
            # Ensure sheet names are less than 31 characters and unique
            sheet_name = sheet_name[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)




if __name__ == '__main__':

    main_go_enrichment_function('uniprotkb_AND_reviewed_true_AND_model_o_2024_10_02.tsv', 'motif_search_output_v31_full_blastp_fishers.csv')