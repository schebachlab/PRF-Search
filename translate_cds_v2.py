#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to search transcriptome databases (Ensembl CDS) for Harrington motifs

@author: Charles Kuntz :: cpkuntz@iu.edu
Adapted from Perl codes by Christopher Hemmerich :: chemmerich@indiana.edu
"""

import sys
import os
import re
import textwrap
import string


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


def ribosome(domain_codons, genetic_code_dict):
    domain_aminoacids = []
    for codon in domain_codons:
        amino_acid = genetic_code_dict[codon]
        domain_aminoacids.append(amino_acid)
    return domain_aminoacids


def load_fasta_data(path_fasta):
    with open(path_fasta, "r") as fasta_handle:
        header = None
        sequence = []
        for line in fasta_handle:
            line = line.strip()
            if line.startswith(">"):
                if header is None:
                    header = line
                else:
                    yield header, sequence
                    header = line
                    sequence = []
            else:
                seq_on_line = list(line)
                sequence.extend(seq_on_line)

def parse_fasta_old(path_fasta):
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
        pass_conditions = [(len_seq % 3 == 0),
                           (sequence_str[0:3] == 'ATG'),
                           (sequence_str[-3:] in ('TAG', 'TAA', 'TGA')),
                           ('N' not in sequence_str)]
        if all(pass_conditions):
            fasta_list.append(output_line)
    return fasta_list


def parse_fasta(path_fasta):
    fasta_generator = load_fasta_data(path_fasta)
    fasta_list = []
    all_pass_fasta_list = []
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

        if (len_seq % 3 == 0) and ('N' not in sequence_str):
            all_pass_fasta_list.append(output_line)

    return fasta_list, all_pass_fasta_list

def translate_fasta(cds_fasta_list):
    genetic_code_dict = genetic_code_func()
    translated_data = []
    for transcript_data in cds_fasta_list:
        ensembl_transcript_name = transcript_data[1]
        sequence_raw = transcript_data[5]
        transcript_seq = re.sub("[^a-zA-Z]", "", sequence_raw)
        seq_codons = textwrap.wrap(transcript_seq, 3)
        seq_protein = ribosome(seq_codons, genetic_code_dict)
        seq_protein_str = ''.join(seq_protein[0:-1])
        translated_entry = [ensembl_transcript_name, seq_protein_str]
        if '*' not in seq_protein_str:
            translated_data.append(translated_entry)
        else:
            continue
    return translated_data


def print_fasta(translated_data, outname):
     myfasta = open(outname+'.fasta', 'w')
     for record in translated_data:
         print('>'+record[0], file = myfasta)
         idx = 0
         columns = 60
         while idx < len(record[1]):
             if len(record[1][idx:]) > columns:
                 print(record[1][idx:idx+columns], file = myfasta)
             else:
                 print(record[1][idx:], file = myfasta)
             idx += columns
     myfasta.close()

def interactive_configuration():
    path_cds_fasta = str(input("Please enter a path to an Ensembl-format CDS database to be translated: "))

    while True:
        good_transcript_switch = int(input("Do you want to only search the transcripts with start and stop codons? [1] yes or [2] no: "))
        if good_transcript_switch not in [1, 2]:
            print("Must enter [1] or [2]. Please re-enter.")
        else:
            break


    outname = str(input("Please enter a name for the output fasta: "))
    return path_cds_fasta, good_transcript_switch, outname

if __name__ == '__main__':

    path_cds_fasta, good_transcript_switch, outname = interactive_configuration()

    fasta_list_good, fasta_list_all = parse_fasta(path_cds_fasta)

    if good_transcript_switch == 1:
        cds_fasta_list = fasta_list_good
    else:
        cds_fasta_list = fasta_list_all

    translated_data = translate_fasta(cds_fasta_list)

    print_fasta(translated_data, outname)



