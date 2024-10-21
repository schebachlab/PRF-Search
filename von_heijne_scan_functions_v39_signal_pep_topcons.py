#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to adjust TMDs predicted by DeepTMHMM using the von Heijne biological scale

@author: Charles Kuntz :: cpkuntz@purdue.edu
Adapted from Perl codes from Gunnar von Heijne's lab and works on TMD predictions from DeepTMHMM
https://dgpred.cbr.su.se/
https://dtu.biolib.com/DeepTMHMM
"""

import argparse
import os
import math
import numpy as np
import csv
import re
import sys


def read_fasta(fasta_file):
    sequences = {}
    with open(fasta_file, 'r') as f:
        current_header = None
        current_sequence = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    seq_id = current_header.split()[0]  # Use only the first part of the header
                    sequences[seq_id] = ''.join(current_sequence)
                current_header = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_header:
            seq_id = current_header.split()[0]  # Use only the first part of the header
            sequences[seq_id] = ''.join(current_sequence)
    return sequences

def parse_topcons_output(path_topcons):
    topcons_predictions = {}
    topcons_output_raw = []
    with open(path_topcons, 'r') as f:
        all_lines = f.read()
        topcons_output_raw.append(all_lines)

    topcons_data_raw = topcons_output_raw[0]
    topcons_data_raw_entries = topcons_data_raw.split('##############################################################################\n')
    entries = []
    for entry in topcons_data_raw_entries[2:-1]:
        lines = entry.strip().split('\n')
        entries.append(lines)

    topcons_data = []
    for entry in entries:
        seq_id_string = entry[1]
        # uniprot_id = seq_id_string.split('|')[1]
        seq_id = seq_id_string.split(' ')[2]

        entry_seq = entry[4]
        topology_string = entry[8]
        # entry_data = [seq_id, uniprot_id, entry_seq, topology_string]
        entry_data = [seq_id, entry_seq, topology_string]
        topcons_data.append(entry_data)


    for entry in topcons_data:
        seq_id = entry[0]

        # topology_str = entry[3]
        topology_str = entry[2]
        segments = []
        current_segment_type = None
        current_segment_start = None

        for i, char in enumerate(topology_str):
            if char != current_segment_type:
                if current_segment_type is not None:
                    segments.append((current_segment_type, current_segment_start, i))
                current_segment_type = char
                current_segment_start = i + 1

        if current_segment_type is not None:
            segments.append((current_segment_type, current_segment_start, len(topology_str)))

        topology_predictions = []
        topcons_predictions[seq_id] = []
        for segment in segments:
            feature_type, start, end = segment

            if feature_type in ('M', 'S'):
                previous_segment = topology_predictions[-1] if topology_predictions else None
                next_segment = segments[segments.index(segment) + 1] if segments.index(segment) + 1 < len(segments) else None

                if previous_segment and previous_segment[0] == 'i':
                    topology = 'i->o'
                elif previous_segment and previous_segment[0] == 'o':
                    topology = 'o->i'
                elif next_segment and next_segment[0] == 'i':
                    topology = 'o->i'
                elif next_segment and next_segment[0] == 'o':
                    topology = 'i->o'
                else:
                    topology = 'unknown'

                if feature_type == 'M':
                    output_feature = 'TMhelix'
                else:
                    output_feature = 'signal'

                topology_predictions.append((output_feature, start, end, topology))
                topcons_predictions[seq_id].append((output_feature, start, end, topology))

    return topcons_predictions

def parse_gff3(gff3_file):
    tmhmm_predictions = {}
    with open(gff3_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            seq_id = fields[0].split()[0]  # Use only the first part of the header
            feature_type = fields[1]
            start = int(fields[2])
            end = int(fields[3])

            if feature_type in ('TMhelix', 'signal'):
                if seq_id not in tmhmm_predictions:
                    tmhmm_predictions[seq_id] = []

                # Determine topology based on neighboring segments
                previous_feature = next(
                    (x for x in reversed(tmhmm_predictions[seq_id]) if x[0] in ['inside', 'outside']), None)
                next_line = next(f, None)
                if next_line:
                    next_fields = next_line.strip().split('\t')
                    next_feature = next_fields[1] if len(next_fields) > 1 else None
                else:
                    next_feature = None

                if previous_feature and previous_feature[0] == 'inside':
                    topology = 'i->o'
                elif previous_feature and previous_feature[0] == 'outside':
                    topology = 'o->i'
                elif next_feature == 'inside':
                    topology = 'o->i'
                elif next_feature == 'outside':
                    topology = 'i->o'
                else:
                    topology = 'unknown'

                tmhmm_predictions[seq_id].append((feature_type, start, end, topology))

    return tmhmm_predictions


def parse_tmhmm_output(path_tmhmm_out):
    tmhmm_predictions = {}
    with open(path_tmhmm_out, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split('\t')
            if len(fields) < 5:
                continue

            seq_id = fields[0].split()[0]  # Use only the first part of the header
            num_helices = int(fields[4].split('=')[1])
            topo_string = fields[-1].split('=')[1]

            if num_helices > 0:
                list_tuples_str = re.findall(r'(\d+)-(\d+)([io])', topo_string)
                for current_tuple in list_tuples_str:
                    tm_start = int(current_tuple[0])
                    tm_end = int(current_tuple[1])
                    topo = current_tuple[2]

                    if seq_id not in tmhmm_predictions:
                        tmhmm_predictions[seq_id] = []

                    # Determine topology based on the current segment
                    if topo == 'i':
                        topology = 'o->i'
                    elif topo == 'o':
                        topology = 'i->o'
                    else:
                        topology = 'unknown'

                    tmhmm_predictions[seq_id].append(('TMhelix', tm_start, tm_end, topology))

    return tmhmm_predictions


def get_global_variables():
    profiles_dict = {
        'A': [0.1267255, 0.0215152],
        'B': [1.5556351, 0.0133523],
        'C': [-0.0765051, 0.0994228],
        'D': [1.7939795, 0.0172643],
        'E': [1.4193720, 0.0089351],
        'F': [-0.2766953, 0.0010297],
        'G': [0.4813492, 0.0047210],
        'H': [1.1998590, 0.0080127],
        'I': [-0.4597384, 0.0181495],
        'K': [1.8485768, 0.0218446],
        'L': [-0.4282992, 0.0023804],
        'M': [-0.0774786, 0.0984413],
        'N': [1.3266132, 0.0092375],
        'P': [1.0860888, 0.0100568],
        'Q': [1.3336109, 0.0111996],
        'R': [1.6492534, 0.0512044],
        'S': [0.7023921, 0.0077661],
        'T': [0.5266550, 0.0311973],
        'U': [-0.0774786, 0.0984413],
        'V': [-0.2447218, 0.0979201],
        'W': [0.2909390, 0.0189282, -0.5479140, 0.0930222, 6.4736619],
        'X': [0.6348884, 0.0180273],
        'Y': [0.6275249, 0.0103896, -0.5744404, 0.0947821, 6.9164963],
        'Z': [1.3761092, 0.0099898]
    }
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    constants = {'c0': 0.27045, 'c1': 9.29274167549645, 'c2': -0.64513139783394, 'c3': 0.00822196628688}
    return profiles_dict, amino_acids, constants


def pos_spec_contrib(aa, i, L, profiles):
    pos = 9 * (2 * (i - 1) / (L - 1) - 1)
    if aa in ['W', 'Y']:
        return (profiles[aa][0] * math.exp(-profiles[aa][1] * pos**2) +
                profiles[aa][2] * (math.exp(-profiles[aa][3] * (pos - profiles[aa][4])**2) +
                                   math.exp(-profiles[aa][3] * (pos + profiles[aa][4])**2)))
    else:
        return profiles[aa][0] * math.exp(-profiles[aa][1] * pos**2)


def calculate_dg(seq, L, profiles, constants):
    seq_length = len(seq)
    results = []
    c0, c1, c2, c3 = constants['c0'], constants['c1'], constants['c2'], constants['c3']
    piover180 = math.pi / 180

    for k in range(seq_length - L + 1):
        DG_sum, DG_sin_sum, DG_cos_sum = 0, 0, 0
        for i in range(1, L + 1):
            aa = seq[k + i - 1]
            DG = pos_spec_contrib(aa, i, L, profiles)
            DG_sum += DG
            DG_sin_sum += DG * math.sin(100 * (i - 1) * piover180)
            DG_cos_sum += DG * math.cos(100 * (i - 1) * piover180)
        helix_DG = DG_sum + c0 * math.sqrt(DG_sin_sum**2 + DG_cos_sum**2) + c1 + c2 * L + c3 * L**2
        results.append((k + 1, helix_DG))  # Adjusted to store 1-based position
    return results


def precalculate_dg(sequences, profiles, constants, Lmin, Lmax):
    dg_values = {}
    total_sequences = len(sequences)
    processed_sequences = 0

    for seq_id, seq in sequences.items():
        seq_length = len(seq)
        dg_values[seq_id] = {L: {} for L in range(Lmin, Lmax + 1)}  # Initialize dictionary for each length L

        for L in range(Lmin, Lmax + 1):
            dg_results = calculate_dg(seq, L, profiles, constants)
            for k, helix_DG in dg_results:
                dg_values[seq_id][L][k] = helix_DG  # Store delta G value in dictionary

        processed_sequences += 1
        progress = (processed_sequences / total_sequences) * 100
        sys.stdout.write(f"\rProcessed {processed_sequences}/{total_sequences} sequences ({progress:.2f}%)")
        sys.stdout.flush()

    print("\nPrecalculation complete.")
    return dg_values


def write_dg_for_sequences(dg_values, sequences, dg_outfile, Lmin, Lmax):
    with open(dg_outfile, 'w') as f:
        for seq_id, seq in sequences.items():
            print(f"Writing delta G values for sequence: {seq_id}")
            f.write(f"# {seq_id}\n")
            f.write("# position\tresidue\t" + "\t".join([f"L{L}" for L in range(Lmin, Lmax + 1)]) + "\n")

            seq_length = len(seq)
            for k in range(seq_length):
                line = [str(k + 1), seq[k]]
                for L in range(Lmin, Lmax + 1):
                    if k + 1 <= seq_length - L + 1:
                        # Use dictionary lookup for precomputed delta G value
                        dg_value = dg_values[seq_id][L].get(k + 1, None)
                        if dg_value is not None:
                            line.append(f"{dg_value:.3f}")
                        else:
                            line.append("")
                    else:
                        line.append("")
                f.write("\t".join(line) + "\n")
    print("\nDone.")


def calculate_dg_for_segment(seq, start, end, profiles, constants):
    DG_sum, DG_sin_sum, DG_cos_sum = 0, 0, 0
    L = end - start + 1
    for j in range(1, L + 1):
        aa = seq[start + j - 2]  # Adjust indexing for 0-based Python indexing
        DG = pos_spec_contrib(aa, j, L, profiles)
        DG_sum += DG
        DG_sin_sum += DG * math.sin(100 * (j - 1) * math.pi / 180)
        DG_cos_sum += DG * math.cos(100 * (j - 1) * math.pi / 180)
    helix_DG = DG_sum + constants['c0'] * math.sqrt(DG_sin_sum**2 + DG_cos_sum**2) + constants['c1'] + constants['c2'] * L + constants['c3'] * L**2
    return helix_DG

def resolve_TMs(potential_tmds, min_loop_length, num_tmds):
    resolved_tmds = []
    used_positions = set()

    for key in sorted(potential_tmds, key=lambda x: potential_tmds[x][2]):
        seg_start, seg_end, helix_DG = potential_tmds[key]
        if all(seg_start > e + min_loop_length or seg_end < s - min_loop_length for s, e, _ in resolved_tmds):
            resolved_tmds.append((seg_start, seg_end, helix_DG))
            used_positions.update(range(seg_start, seg_end + 1))
            if len(resolved_tmds) == num_tmds:
                break

    return resolved_tmds


def optimize_tmhmm_predictions(sequences, tmhmm_predictions, profiles, constants, adjustment_range, Lmin, Lmax, min_loop_length, dg_values):
    optimized_predictions = []

    for seq_id, seq in sequences.items():
        print(f"Processing sequence: {seq_id}")
        if seq_id not in tmhmm_predictions:
            print(f"No TMHMM predictions for sequence {seq_id}")
            continue

        features = [f for f in tmhmm_predictions[seq_id] if f[0] in ('TMhelix', 'signal')]
        potential_tmds = {}

        for feature in features:
            start, end = feature[1], feature[2]
            len_tmd_candidate = end - start + 1
            print(f"Feature: start={start}, end={end}, length={len_tmd_candidate}")
            center = (start + end) // 2

            for L in range(max(Lmin, len_tmd_candidate - adjustment_range), min(Lmax, len_tmd_candidate + adjustment_range) + 1):
                if L not in dg_values[seq_id]:
                    continue

                for seg_start, helix_DG in dg_values[seq_id][L].items():
                    seg_end = seg_start + L - 1

                    if abs(center - ((seg_start + seg_end) // 2)) > adjustment_range:
                        continue

                    key = f"{L}_{seg_start}"
                    potential_tmds[key] = (seg_start, seg_end, helix_DG)

        num_tmds = len(features)
        resolved_tmds = resolve_TMs(potential_tmds, min_loop_length, num_tmds)

        if len(resolved_tmds) != num_tmds:
            print(f"Warning: Number of resolved TMDs ({len(resolved_tmds)}) does not match TMHMM predictions ({num_tmds}).")

        # Sort resolved_tmds by their start positions
        resolved_tmds.sort(key=lambda x: x[0])

        for i, (best_start, best_end, best_dg) in enumerate(resolved_tmds):
            original_feature = features[i]
            feature_type = original_feature[0]
            start, end = original_feature[1], original_feature[2]
            original_dg = calculate_dg_for_segment(seq, start, end, profiles, constants)
            adjusted_upstream_loop = best_start - 1 if i == 0 else best_start - resolved_tmds[i - 1][1] - 1
            adjusted_downstream_loop = len(seq) - best_end if i == len(resolved_tmds) - 1 else resolved_tmds[i + 1][0] - best_end - 1

            adjusted_upstream_loop = max(adjusted_upstream_loop, 0)
            adjusted_downstream_loop = max(adjusted_downstream_loop, 0)

            original_len_tmd = end - start + 1
            adjusted_len_tmd = best_end - best_start + 1

            optimized_predictions.append({
                'seq_name': seq_id,
                'tmhmm_predicted_tmd_seq': seq[start - 1:end],
                'len_tmhmm_predicted_tmd': original_len_tmd,
                'tmhmm_predicted_tmd_topology': original_feature[3],
                'tmhmm_predicted_tmd_start_position': start,
                'tmhmm_predicted_tmd_end_position': end,
                'tmhmm_predicted_tmd_delta_G': original_dg,
                'tmhmm_predicted_tmd_upstream_loop_length': original_feature[1] - 1 if i == 0 else original_feature[1] - features[i - 1][2] - 1,
                'tmhmm_predicted_tmd_downstream_loop_length': len(seq) - original_feature[2] if i == len(features) - 1 else features[i + 1][1] - original_feature[2] - 1,
                'adjusted_tmd_seq': seq[best_start - 1:best_end],
                'len_adjusted_tmd': adjusted_len_tmd,
                'adjusted_tmd_start_position': best_start,
                'adjusted_tmd_end_position': best_end,
                'adjusted_tmd_delta_G': best_dg,
                'adjusted_tmd_upstream_loop_length': adjusted_upstream_loop,
                'adjusted_tmd_downstream_loop_length': adjusted_downstream_loop,
                'tmd_type': feature_type
            })

    return optimized_predictions




def parse_mutations(mutation_file):
    def is_mutation_format(mutation):
        return len(mutation) >= 3 and mutation[0].isalpha() and mutation[-1].isalpha() and mutation[1:-1].isdigit()

    mutations_dict = {}
    with open(mutation_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t' if mutation_file.endswith('.tsv') else ',')
        headers = next(reader, None)  # Read the first row as headers
        #if not is_mutation_format(headers[1]):
        #    headers = next(reader, None)  # Skip headers if the second field is not in mutation format

        for row in reader:
            if len(row) != 2:
                continue
            seq_id, mutation = row
            if seq_id not in mutations_dict:
                mutations_dict[seq_id] = []
            mutations_dict[seq_id].append(mutation)
    return mutations_dict


def apply_mutation(seq, mutation):
    pos = int(mutation[1:-1]) - 1
    wt_residue = mutation[0]
    mut_residue = mutation[-1]
    if seq[pos] != wt_residue:
        raise ValueError(f"Wild-type residue {wt_residue} does not match sequence residue {seq[pos]} at position {pos + 1}")
    return seq[:pos] + mut_residue + seq[pos + 1:]


def calculate_dg_for_sequence(seq, profiles, constants, Lmin, Lmax):
    seq_length = len(seq)
    dg_values = {}

    for L in range(Lmin, Lmax + 1):
        if L not in dg_values:
            dg_values[L] = {}
        for k in range(seq_length - L + 1):
            DG_sum, DG_sin_sum, DG_cos_sum = 0, 0, 0
            for i in range(1, L + 1):
                aa = seq[k + i - 1]
                DG = pos_spec_contrib(aa, i, L, profiles)
                DG_sum += DG
                DG_sin_sum += DG * math.sin(100 * (i - 1) * math.pi / 180)
                DG_cos_sum += DG * math.cos(100 * (i - 1) * math.pi / 180)
            helix_DG = DG_sum + constants['c0'] * math.sqrt(DG_sin_sum**2 + DG_cos_sum**2) + constants['c1'] + constants['c2'] * L + constants['c3'] * L**2
            dg_values[L][k + 1] = helix_DG

    return dg_values



def calculate_mutation_impact(seq_id, seq, mutations, tmhmm_predictions, profiles, constants, adjustment_range, Lmin, Lmax, min_loop_length, dg_values):
    mutation_results = []

    # Use precomputed ΔG values for wild-type sequence
    wt_predictions = optimize_tmhmm_predictions({seq_id: seq}, {seq_id: tmhmm_predictions[seq_id]}, profiles, constants, adjustment_range, Lmin, Lmax, min_loop_length, dg_values)

    for mutation in mutations:
        mut_seq = apply_mutation(seq, mutation)
        # Dynamically calculate ΔG values for mutated sequence
        mut_dg_values = {seq_id: calculate_dg_for_sequence(mut_seq, profiles, constants, Lmin, Lmax)}

        mut_predictions = optimize_tmhmm_predictions({seq_id: mut_seq}, {seq_id: tmhmm_predictions[seq_id]}, profiles, constants, adjustment_range, Lmin, Lmax, min_loop_length, mut_dg_values)

        for wt_pred, mut_pred in zip(wt_predictions, mut_predictions):
            wt_dg = wt_pred['adjusted_tmd_delta_G']
            mut_dg = mut_pred['adjusted_tmd_delta_G']
            ddg = mut_dg - wt_dg
            mutation_results.append({
                'seq_name': seq_id,
                'mutation': mutation,
                'wildtype_tmd_seq': wt_pred['adjusted_tmd_seq'],
                'len_wildtype_tmd': wt_pred['len_adjusted_tmd'],
                'wildtype_tmd_start': wt_pred['adjusted_tmd_start_position'],
                'wildtype_tmd_end': wt_pred['adjusted_tmd_end_position'],
                'wildtype_tmd_topology': wt_pred['tmhmm_predicted_tmd_topology'],
                'wildtype_tmd_upstream_loop_length': wt_pred['adjusted_tmd_upstream_loop_length'],
                'wildtype_tmd_downstream_loop_length': wt_pred['adjusted_tmd_downstream_loop_length'],
                'wildtype_tmd_delta_G': wt_dg,
                'mutated_tmd_seq': mut_pred['adjusted_tmd_seq'],
                'len_mutated_tmd': mut_pred['len_adjusted_tmd'],
                'mutated_tmd_start': mut_pred['adjusted_tmd_start_position'],
                'mutated_tmd_end': mut_pred['adjusted_tmd_end_position'],
                'mutated_tmd_topology': mut_pred['tmhmm_predicted_tmd_topology'],
                'mutated_tmd_upstream_loop_length': mut_pred['adjusted_tmd_upstream_loop_length'],
                'mutated_tmd_downstream_loop_length': mut_pred['adjusted_tmd_downstream_loop_length'],
                'mutated_tmd_delta_G': mut_dg,
                'delta_delta_G': ddg
            })

    return mutation_results



def main():

    parser = argparse.ArgumentParser(description='Optimize TMD predictions using von Heijne delta G scale and analyze mutations.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file')
    parser.add_argument('prediction_file', type=str, help='Input file with TOPCONS2, classic TMHMM, or DeepTMHMM predictions (TOPCONS2 text output, GFF3 for DeepTMHMM or classic TMHMM text output)')
    parser.add_argument('format', type=str, choices=['deep_tmhmm', 'classic_tmhmm', 'topcons'], help='Format of the prediction file (deep_tmhmm, classic_tmhmm, or topcons)')
    parser.add_argument('-o', '--outfile', type=str, default='output.csv', help='Output CSV file')
    parser.add_argument('-d', '--dg_outfile', type=str, default='dg_output.txt', help='Output file for delta G values')
    parser.add_argument('-a', '--adjustment_range', type=int, default=10, help='Adjustment range for optimizing TMDs')
    parser.add_argument('-lmin', '--Lmin', type=int, default=15, help='Minimum length of the sliding window')
    parser.add_argument('-lmax', '--Lmax', type=int, default=30, help='Maximum length of the sliding window')
    parser.add_argument('-m', '--min_loop_length', type=int, default=3, help='Minimum loop length between neighboring helices')
    parser.add_argument('--mutation_file', type=str, help='Optional file with mutations to analyze')
    parser.add_argument('--skip_dg_output', action='store_true', help='Skip writing the delta G output file')
    args = parser.parse_args()

    sequences = read_fasta(args.fasta_file)

    if args.format == 'deep_tmhmm':
        tmhmm_predictions = parse_gff3(args.prediction_file)
    elif args.format == 'classic_tmhmm':
        tmhmm_predictions = parse_tmhmm_output(args.prediction_file)
    else:
        tmhmm_predictions = parse_topcons_output(args.prediction_file)

    profiles, _, constants = get_global_variables()
    dg_values = precalculate_dg(sequences, profiles, constants, args.Lmin, args.Lmax)

    optimized_predictions = optimize_tmhmm_predictions(sequences, tmhmm_predictions, profiles, constants, args.adjustment_range, args.Lmin, args.Lmax, args.min_loop_length, dg_values)

    if not args.skip_dg_output and args.dg_outfile:
        print("Writing full-sequence delta G scan output...")
        write_dg_for_sequences(dg_values, sequences, args.dg_outfile, args.Lmin, args.Lmax)
        print("Finished writing full-sequence delta G scan output.")

    with open(args.outfile, 'w', newline='') as out_f:
        writer = csv.DictWriter(out_f, fieldnames=[
            'seq_name', 'tmhmm_predicted_tmd_seq', 'len_tmhmm_predicted_tmd', 'tmhmm_predicted_tmd_topology',
            'tmhmm_predicted_tmd_start_position', 'tmhmm_predicted_tmd_end_position',
            'tmhmm_predicted_tmd_delta_G', 'tmhmm_predicted_tmd_upstream_loop_length', 'tmhmm_predicted_tmd_downstream_loop_length',
            'adjusted_tmd_seq', 'len_adjusted_tmd', 'adjusted_tmd_start_position',
            'adjusted_tmd_end_position', 'adjusted_tmd_delta_G',
            'adjusted_tmd_upstream_loop_length', 'adjusted_tmd_downstream_loop_length', 'tmd_type'
        ])
        writer.writeheader()
        print("Writing output CSV..")
        for entry in optimized_predictions:
            writer.writerow(entry)

        print("Finished writing output CSV")

    if args.mutation_file:
        mutations_dict = parse_mutations(args.mutation_file)
        mutation_results = []
        for seq_id, mutations in mutations_dict.items():
            if seq_id in sequences:
                mutation_results.extend(calculate_mutation_impact(seq_id, sequences[seq_id], mutations, tmhmm_predictions, profiles, constants, args.adjustment_range, args.Lmin, args.Lmax, args.min_loop_length, dg_values))

        mutation_outfile = args.outfile.replace('.csv', '_mutations.csv')
        with open(mutation_outfile, 'w', newline='') as mut_f:
            writer = csv.DictWriter(mut_f, fieldnames=[
                'seq_name', 'mutation', 'wildtype_tmd_seq', 'len_wildtype_tmd', 'wildtype_tmd_start',
                'wildtype_tmd_end', 'wildtype_tmd_topology', 'wildtype_tmd_upstream_loop_length',
                'wildtype_tmd_downstream_loop_length', 'wildtype_tmd_delta_G',
                'mutated_tmd_seq', 'len_mutated_tmd', 'mutated_tmd_start', 'mutated_tmd_end',
                'mutated_tmd_topology', 'mutated_tmd_upstream_loop_length',
                'mutated_tmd_downstream_loop_length', 'mutated_tmd_delta_G',
                'delta_delta_G'
            ])
            writer.writeheader()
            for entry in mutation_results:
                writer.writerow(entry)
        print("Finished writing mutation output CSV")


if __name__ == '__main__':
    main()
