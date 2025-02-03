# Goal: rough cut of somatic mutation library preparation simulation 
import numpy as np 
import pandas as pd 


# step 1: ligate on adapter molecules 
# unfortunately CS adapters are "propriatary" so here is my best guess
i5 = 'ACACTCTTTCCCTACACGACATGATGATGATGATGATGC'
i7 = 'GTGACTGGAGTTCAGACGTGTATGATGATGATGATGATGT'

#i have left out the 3bp UMI, which I will randomly sample and tack on with each read 
def rev_comp(dna_sequence): #get the reverse complement 
    table = str.maketrans({
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G"})
    rec_comp=dna_sequence.translate(table)[::-1]
    return rec_comp

def add_umi(fragment):
    pool = ['C', 'A', 'T', 'G']
    umi_i5 = sample(pool,3)
    umi_i7 = sample(pool,3)
    c_i5 = ''.join(umi_5)
    c_i7 = ''.join(umi_7)
    adapted = i5+c_i5+fragment+rev_comp(i7+c_i7)
    return adapted 

