def calc_energy(pairs, nn_pairs):
    # Values taken from the Nearest Neighbor Database
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2808915/
    # http://rna.urmc.rochester.edu/NNDB/
    # Turner2004 values
    watson_crick_nn_pair = nn_pairs[0]
    wobble_nn_pair = nn_pairs[1]
    matches = ('AU', 'UA', 'GC', 'CG', 'GU', 'UG')
    mismatches = ('AA', 'GG', 'GA', 'AG', 'UU', 'CC', 'CU', 'UC', 'AC', 'CA')
    # Initialize energy to 0
    # Per codon-anticodon interaction.
    energy = 0
    # Apply intermolecular initiation energy:
    intermolecular_init = 4.09
    energy += intermolecular_init
    # Apply AU end penalty if applicable (only if there are watson-cricks in previous position)
    if pairs[2] in ('GU', 'UG', 'AU', 'UA') and pairs[1] in matches:
        energy += 0.45
    if all([pairs[1] in ('GU', 'UG', 'AU', 'UA'), pairs[0] in matches, pairs[2] in mismatches]):
        energy += 0.45
    # If all pairs are mismatches, use +3 bulge penalty
    if all([pairs[0] in mismatches, pairs[1] in mismatches, pairs[2] in mismatches]):
        energy += 3.2

    # Watson-Crick pairs - search both pairs in codons
    WC_pairs_dict =       {('AA', 'UU'):-0.93,
                           ('UU', 'AA'):-0.93,
                           ('AU', 'UA'):-1.10,
                           ('UA', 'AU'):-1.33,
                           ('CU', 'GA'):-2.08,
                           ('AG', 'UC'):-2.08,
                           ('CA', 'GU'):-2.11,
                           ('UG', 'AC'):-2.11,
                           ('GU', 'CA'):-2.24,
                           ('AC', 'UG'):-2.24,
                           ('GA', 'CU'):-2.35,
                           ('UC', 'AG'):-2.35,
                           ('CG', 'GC'):-2.36,
                           ('GG', 'CC'):-3.26,
                           ('CC', 'GG'):-3.26,
                           ('GC', 'CG'):-3.42}

    WC_pairs_list = list(WC_pairs_dict.keys())
    if watson_crick_nn_pair in WC_pairs_list:
        energy += WC_pairs_dict[watson_crick_nn_pair]
    if wobble_nn_pair in WC_pairs_list:
        energy += WC_pairs_dict[wobble_nn_pair]


    # GU pairs - search both pairs in codons
    GU_pairs_dict =       {('AG', 'UU'):-0.55,
                           ('UU', 'GA'):-0.55,
                           ('AU', 'UG'):-1.36,
                           ('GU', 'UA'):-1.36,
                           ('CG', 'GU'):-1.41,
                           ('UG', 'GC'):-1.41,
                           ('CU', 'GG'):-1.53,
                           ('GG', 'UC'):-1.53,
                           ('GU', 'CG'):-2.51,
                           ('GC', 'UG'):-2.51,
                           ('GA', 'UU'):-1.27,
                           ('UU', 'AG'):-1.27,
                           ('GG', 'UU'):+0.50,
                           ('UU', 'GG'):+0.50,
                           ('GU', 'UG'):+1.29,
                           ('UG', 'AU'):-1.00,
                           ('UA', 'GU'):-1.00,
                           ('UG', 'GU'):+0.30}

    GU_pairs_list =  list(GU_pairs_dict.keys())
    if watson_crick_nn_pair in GU_pairs_list:
        energy += GU_pairs_dict[watson_crick_nn_pair]
    if wobble_nn_pair in GU_pairs_list:
        energy += GU_pairs_dict[wobble_nn_pair]

    # Terminal mismatches - only search last two pairs in codon
    term_mismatch_dict =  {('AA', 'UA'):-0.80,
                           ('AC', 'UA'):-0.60,
                           ('AG', 'UA'):-0.80,
                           ('AA', 'UC'):-1.00,
                           ('AC', 'UC'):-0.70,
                           ('AU', 'UC'):-0.80,
                           ('AA', 'UG'):-0.80,
                           ('AG', 'UG'):-0.80,
                           ('AC', 'UU'):-0.70,
                           ('AU', 'UU'):-0.80,
                           ('CA', 'GA'):-1.50,
                           ('CC', 'GA'):-1.00,
                           ('CG', 'GA'):-1.40,
                           ('CA', 'GC'):-1.50,
                           ('CC', 'GC'):-1.10,
                           ('CU', 'GC'):-1.40,
                           ('CA', 'GG'):-1.40,
                           ('CG', 'GG'):-1.60,
                           ('CC', 'GU'):-0.80,
                           ('CU', 'GU'):-1.20,
                           ('GA', 'CA'):-1.10,
                           ('GC', 'CA'):-1.10,
                           ('GG', 'CA'):-1.60,
                           ('GA', 'CC'):-1.50,
                           ('GC', 'CC'):-0.70,
                           ('GU', 'CC'):-1.00,
                           ('GA', 'CG'):-1.30,
                           ('GG', 'CG'):-1.40,
                           ('GC', 'CU'):-0.50,
                           ('GG', 'CU'):-1.50,
                           ('GU', 'CU'):-0.70,
                           ('GA', 'UA'):-0.30,
                           ('GC', 'UA'):-0.60,
                           ('GG', 'UA'):-0.60,
                           ('GA', 'UC'):-1.00,
                           ('GC', 'UC'):-0.70,
                           ('GU', 'UC'):-0.80,
                           ('GA', 'UG'):-0.80,
                           ('GG', 'UG'):-0.80,
                           ('GC', 'UU'):-0.70,
                           ('GU', 'UU'):-0.60,
                           ('UA', 'AA'):-1.00,
                           ('UC', 'AA'):-0.70,
                           ('UG', 'AA'):-1.10,
                           ('UA', 'AC'):-0.80,
                           ('UC', 'AC'):-0.60,
                           ('UU', 'AC'):-0.60,
                           ('UA', 'AG'):-1.10,
                           ('UG', 'AG'):-1.20,
                           ('UC', 'AU'):-0.50,
                           ('UU', 'AU'):-0.50,
                           ('UA', 'GA'):-1.00,
                           ('UC', 'GA'):-0.70,
                           ('UG', 'GA'):-0.50,
                           ('UA', 'GC'):-0.80,
                           ('UC', 'GC'):-0.60,
                           ('UU', 'GC'):-0.60,
                           ('UA', 'GG'):-1.10,
                           ('UC', 'GG'):-0.70,
                           ('UG', 'GG'):-0.80,
                           ('UC', 'GU'):-0.50,
                           ('UU', 'GU'):-0.50}

    term_mismatch_list = list(term_mismatch_dict.keys())
    # Only check Watson-Crick NN pair if second AND third positions are mismatches
    if pairs[1] in mismatches and pairs[2] in mismatches:
        if pairs[0] in matches:
            if watson_crick_nn_pair in term_mismatch_list:
                energy += term_mismatch_dict[watson_crick_nn_pair]

    # Check wobble position NN pair only if third position is a mismatch and other positions are matches
    if pairs[2] in mismatches:
        if pairs[0] in matches and pairs[1] in matches:
            if wobble_nn_pair in term_mismatch_list:
                energy += term_mismatch_dict[wobble_nn_pair]

    # If there are mismatches in the first two positions, reverse the sequence of the Wobble position pairs
    if pairs[0] in mismatches and pairs[1] in mismatches:
        wobble_nn_pair_rev = (wobble_nn_pair[1][::-1], wobble_nn_pair[0][::-1])
        # Check GU and terminal mismatch pairs now
        if wobble_nn_pair_rev in term_mismatch_list:
            energy += term_mismatch_dict[wobble_nn_pair_rev]


    #Simulate central mismatches as 1x1 internal loops?
    # Row 1
    if pairs[0] == 'AU' and pairs[2] == 'AU':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.5
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9

    if pairs[0] == 'AU' and pairs[2] == 'CG':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 0.8
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2
    
    if pairs[0] == 'AU' and pairs[2] == 'GC':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 0.8
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2

    if pairs[0] == 'AU' and pairs[2] == 'UA':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.2
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9
    
    if pairs[0] == 'AU' and pairs[2] == 'GU':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.6
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9
    
    if pairs[0] == 'AU' and pairs[2] == 'UG':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.2
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9
    # Row 2
    if pairs[0] == 'CG' and pairs[2] == 'AU':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2
    
    if pairs[0] == 'CG' and pairs[2] == 'CG':
        if pairs[1] == 'AA':
            energy += 0.9
        elif pairs[1] == 'CA':
            energy += 0.3
        elif pairs[1] == 'GA':
            energy += -0.1
        elif pairs[1] == 'AC':
            energy += -0.4
        elif pairs[1] == 'UC':
            energy += 0
        elif pairs[1] == 'GG':
            energy += -2.2
        elif pairs[1] == 'UU':
            energy += -0.1
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 0.5
    
    if pairs[0] == 'CG' and pairs[2] == 'GC':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 0.4
        elif pairs[1] == 'AA':
            energy += 0.9
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 0.5
    
    if pairs[0] == 'CG' and pairs[2] == 'UA':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 0.8
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2

    if pairs[0] == 'CG' and pairs[2] == 'GU':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 1.1
        elif pairs[1] == 'AA':
            energy += 2.2
        elif pairs[1] == 'CC':
            energy += 1.7
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2
    
    if pairs[0] == 'CG' and pairs[2] == 'UG':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 1.1
        elif pairs[1] == 'AA':
            energy += 0.6
        elif pairs[1] == 'GA':
            energy += -0.2
        elif pairs[1] == 'AC':
            energy += 0.5
        elif pairs[1] == 'UC':
            energy += 1.0
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2
    
    # Row 3
    if pairs[0] == 'GC' and pairs[2] == 'AU':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2
    
    if pairs[0] == 'GC' and pairs[2] == 'CG':
        if pairs[1] == 'GG':
            energy += -2.3
        elif pairs[1] == 'UU':
            energy += -0.6
        elif pairs[1] == 'AA':
            energy += 0.8
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 0.5
    
    if pairs[0] == 'GC' and pairs[2] == 'GC':
        if pairs[1] == 'GG':
            energy += -2.2
        elif pairs[1] == 'UU':
            energy += -0.1
        elif pairs[1] == 'AA':
            energy += 0.9
        elif pairs[1] == 'CA':
            energy += -0.4
        elif pairs[1] == 'AC':
            energy += 0.3
        elif pairs[1] == 'UC':
            energy += 0.6
        elif pairs[1] == 'AG':
            energy += -0.1
        elif pairs[1] == 'CU':
            energy += 0
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 0.5
    
    if pairs[0] == 'GC' and pairs[2] == 'UA':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += -0.8
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2
    
    if pairs[0] == 'GC' and pairs[2] == 'GU':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 0.7
        elif pairs[1] == 'AA':
            energy += 1.6
        elif pairs[1] == 'AG':
            energy += 1.0
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2
    
    if pairs[0] == 'GC' and pairs[2] == 'UG':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 1.5
        elif pairs[1] == 'AA':
            energy += 1.9
        elif pairs[1] == 'AG':
            energy += 1.5
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2

    # Row 4
    if pairs[0] == 'UA' and pairs[2] == 'AU':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.7
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9
    
    if pairs[0] == 'UA' and pairs[2] == 'CG':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2

    if pairs[0] == 'UA' and pairs[2] == 'GC':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2
    
    if pairs[0] == 'UA' and pairs[2] == 'UA':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.5
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9

    if pairs[0] == 'UA' and pairs[2] == 'GU':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9

    if pairs[0] == 'UA' and pairs[2] == 'UG':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.6
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9
    # Row 5
    if pairs[0] == 'GU' and pairs[2] == 'AU':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.6
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9

    if pairs[0] == 'GU' and pairs[2] == 'CG':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 1.5
        elif pairs[1] == 'AA':
            energy += 1.9
        elif pairs[1] == 'GA':
            energy += 1.5
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2

    if pairs[0] == 'GU' and pairs[2] == 'GC':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 1.1
        elif pairs[1] == 'AA':
            energy += 0.6
        elif pairs[1] == 'CA':
            energy += 0.5
        elif pairs[1] == 'AG':
            energy += -0.2
        elif pairs[1] == 'CU':
            energy += 1.0
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2

    if pairs[0] == 'GU' and pairs[2] == 'UA':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.2
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9

    if pairs[0] == 'GU' and pairs[2] == 'GU':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.6
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9

    if pairs[0] == 'GU' and pairs[2] == 'UG':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.2
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9
    # Row 5
    if pairs[0] == 'UG' and pairs[2] == 'AU':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9

    if pairs[0] == 'UG' and pairs[2] == 'CG':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 0.7
        elif pairs[1] == 'AA':
            energy += 1.6
        elif pairs[1] == 'GA':
            energy += 1.0
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2

    if pairs[0] == 'UG' and pairs[2] == 'GC':
        if pairs[1] == 'GG':
            energy += -1.4
        elif pairs[1] == 'UU':
            energy += 1.1
        elif pairs[1] == 'AA':
            energy += 2.2
        elif pairs[1] == 'CA':
            energy += 1.3
        elif pairs[1] == 'CC':
            energy += 1.7
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.2

    if pairs[0] == 'UG' and pairs[2] == 'UA':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.6
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9

    if pairs[0] == 'UG' and pairs[2] == 'GU':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9

    if pairs[0] == 'UG' and pairs[2] == 'UG':
        if pairs[1] == 'GG':
            energy += -0.7
        elif pairs[1] == 'UU':
            energy += 1.6
        elif pairs[1] in ('AU', 'UA', 'GU', 'UG', 'GC', 'CG'):
            energy += 0
        else:
            energy += 1.9

    return energy

