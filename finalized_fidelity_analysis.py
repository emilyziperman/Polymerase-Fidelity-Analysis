#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 16:11:17 2024

@author: emilyziperman
"""

import pandas as pd
from Bio import SeqIO
import numpy as np

### forward and reverse reference sequences
total_ref1 = "ACAGCTAAATGGCGTTGTTCAAGCCCTACCCAAAGATTGGCGATATTCGTAAGGCGCGCTGCATGTTGCAGCACACCTTGCACCACCGGACCAACAAGCAGCCCAGCTACCGCAGGAGGTTGAAGACCCTCATCCCCCTCTTCAGGCGGTGCATGCTCGGCTCCGGTT"

list_of_characters_1 = [*total_ref1]

df = pd.read_csv('therminator_subs.csv')

### creating our first dataframe for substitution analysis
forward_reads_df_sub = df
forward_reads_df_sub = forward_reads_df_sub.drop(columns=['Reads'])

forward_reads_df_sub[[f'base_{i}' for i in range(len(forward_reads_df_sub['TargetSequence'][0]))]] = forward_reads_df_sub.TargetSequence.str.split('', expand=True).iloc[:, 1:-1]
forward_reads_df_sub = forward_reads_df_sub.drop(columns=['TargetSequence'])

### creating a reference data frame for the forward reads
ref_df1 = pd.DataFrame([x.split(';') for x in total_ref1.split()])
ref_df1['TargetSequence'] = ref_df1[0].astype(str)
ref_df1[[f'base_{i}' for i in range(len(ref_df1['TargetSequence'][0]))]] = ref_df1.TargetSequence.str.split('', expand=True).iloc[:, 1:-1]

ref_df1 = ref_df1.drop(columns=[0, 'TargetSequence'])
ref_len_incr_1 = len(forward_reads_df_sub.index)

ref_df1 = ref_df1.loc[ref_df1.index.repeat(ref_len_incr_1)]
ref_df1 = ref_df1.reset_index()
ref_df1 = ref_df1.drop(columns=['index'])

### identifying forward substitutions
for_compare = forward_reads_df_sub!=ref_df1 #Reports false if the values match and true if they are a mismatch
forward_reads_df_sub_1 = forward_reads_df_sub.mask(~for_compare)
forward_reads_df_sub_2 = forward_reads_df_sub_1.fillna('X')

### create a list of sequence weights
factor_list = df['Reads'].tolist()

### prepping the eventual mutations master dataframe

mutations = []
mutations = pd.DataFrame(mutations)

### A substitutions

A_subs = forward_reads_df_sub_2.replace({'C': 0, 'T': 0, 'G': 0, 'X': 0, 'A': 1})
A_subs = A_subs.T
A_subs.iloc[1:,:] = A_subs.iloc[1:,:]*factor_list
mutations['A'] = A_subs.sum(axis=1)

### T substitutions

T_subs = forward_reads_df_sub_2.replace({'C': 0, 'T': 1, 'G': 0, 'X': 0, 'A': 0})
T_subs = T_subs.T
T_subs.iloc[1:,:] = T_subs.iloc[1:,:]*factor_list
mutations['T'] = T_subs.sum(axis=1)

### G substitutions

G_subs = forward_reads_df_sub_2.replace({'C': 0, 'T': 0, 'G': 1, 'X': 0, 'A': 0})
G_subs = G_subs.T
G_subs.iloc[1:,:] = G_subs.iloc[1:,:]*factor_list
mutations['G'] = G_subs.sum(axis=1)

### C substitutions

C_subs = forward_reads_df_sub_2.replace({'C': 1, 'T': 0, 'G': 0, 'X': 0, 'A': 0})
C_subs = C_subs.T
C_subs.iloc[1:,:] = C_subs.iloc[1:,:]*factor_list
mutations['C'] = C_subs.sum(axis=1)

### total bases analyzed

number_of_reads = sum(factor_list)

### gotta calculate the total number of C's T's A's and G's

number_of_As = total_ref1.count('A')
number_of_Ts = total_ref1.count('T')
number_of_Gs = total_ref1.count('G')
number_of_Cs = total_ref1.count('C')

total_As = number_of_As*number_of_reads
total_Ts = number_of_Ts*number_of_reads
total_Gs = number_of_Gs*number_of_reads
total_Cs = number_of_Cs*number_of_reads

### calculating error rates

mutations['new_index'] = list_of_characters_1
mutations = mutations.groupby(["new_index"]).sum()

mutations = mutations.T

mutation_rates = []
mutation_rates = pd.DataFrame(mutation_rates)

mutation_rates['a_rate'] = (mutations['A']/total_As)/4
mutation_rates['c_rate'] = (mutations['C']/total_Cs)/4
mutation_rates['g_rate'] = (mutations['G']/total_Gs)/4
mutation_rates['t_rate'] = (mutations['T']/total_Ts)/4

mutation_rates = mutation_rates.T

















