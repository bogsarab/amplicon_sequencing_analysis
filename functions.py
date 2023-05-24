from Bio import AlignIO
import pandas as pd
from typing import List
import numpy as np
import pandas as pd
import logomaker
import matplotlib.pyplot as plt


def reverse_complement(sequence: str):
    rc_dict = {'A': 'T',
               'C': 'G',
               'G': 'C',
               'T': 'A'}
    return ''.join([rc_dict[base] for base in sequence[::-1]])


def aln_site_composition_df(aln, characters="ACGTUWSMKRYBDHVN"):
    aln_rows = aln.get_alignment_length()
    comp_dict = {char: [0] * aln_rows for char in characters}
    for record in aln:
        header = record.id
        seq = record.seq

        for base_pos in range(len(seq)):
            base = seq[base_pos]
            if base in characters:
                comp_dict[base][base_pos] += 1
    return pd.DataFrame.from_dict(comp_dict)


def make_frequency_table(sequences: List[str], length_of_seq: int, possible_bases: str, path: str): # PATH A VÉGEREDMÉNYT IDE ÍRJA KI
    print('\tMaking frequency table')
    base_dict = {base: np.repeat(0, length_of_seq) for base in possible_bases}
    counting_df = pd.DataFrame(base_dict, index=range(0, length_of_seq))

    for seq in sequences:
        for idx, base in enumerate(seq): # oszl név bázisok és a píció az index a szekvencián belül
            counting_df.loc[idx, base] += 1 # pl 1. oszlop A akkor ahhoz +1

    frequency_df = counting_df.apply(lambda x: x/sum(x), 1) # soronként frekv. számol
    print('\tTable made.')
    frequency_df.to_csv(path)
    return frequency_df


def save_seqlogo(frequency_df: pd.DataFrame, sample: str, path: str): #csomag neve a logomaker ez csinálta az ábrát
    print('\tStarts making seqlogo. ')
    color_dict = {
        'A': 'green',
        'C': 'blue',
        'G': 'orange',
        'T': 'red'}
    for letter in 'RYSWKMBN':
        color_dict[letter] = 'grey' # color shame összes többi legyen szürke

    logo = logomaker.Logo(df=frequency_df,
                          color_scheme=color_dict,
                          alpha=0.7) # alpha - mennyire legyen áttetsző
    plt.title(sample)
    plt.savefig(path)


