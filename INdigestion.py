#!/usr/bin/python3
# -*- coding: utf-8 -*-

# INdigestion
# "I Need a digestion" allows the user to specify a list of restriction enzymes
# and a plasmid sequence and have the program automatically select the best
# one or two enzymes to use to get the best gel

import re
import Bio
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
import csv
import numpy as np
from itertools import combinations
import argparse
import sys

# location of enzyme dictionary file
print("===============================================")

# Initiate the parser
parser = argparse.ArgumentParser()
# Add long and short argument
parser.add_argument("--plasmid", "-p", help="(required) plasmid file, .gb or .ape, insert = \"gene\"")
parser.add_argument("--enzymes", "-e", help="file with list of your enzymes")
parser.add_argument("--min_bands", "-n", help="minimum number of desired bands")
parser.add_argument("--max_bands", "-x", help="maximum number of desired bands")
parser.add_argument("--small_bands", "-s", help="smallest size of desired bands")
parser.add_argument("--big_bands", "-b", help="biggest size of desired bands")
parser.add_argument("--band_gap", "-g", help="minimum gap size between bands")
parser.add_argument("--no_insert", "-i", help="no insert in plasmid", action='store_true')
parser.add_argument("--fasta", "-f", help="plasmid in fasta format (no insert)", action='store_true')
parser.add_argument("--test", "-t", help="tests program with .gb test file", action='store_true')

# Read arguments from the command line
args = parser.parse_args()
if args.plasmid:
    plasmid_file = args.plasmid
# ---------------------
if args.enzymes:
    my_enzymes_file = args.enzymes
else:
    my_enzymes_file = str(sys.path[0] + "/my_enzymes.csv") 
# ---------------------
if args.min_bands:  # minimum number of bands on gel
    min_frags = int(args.min_bands)
else:
    min_frags = 2  
# ---------------------
if args.max_bands:  # maximum number of bands on gel
    max_frags = int(args.max_bands)
else:
    max_frags = 6                  
# ---------------------
if args.small_bands:    # minimum size of bands on gel
    min_frag_size = int(args.small_bands)
else:
    min_frag_size = 500
# ---------------------    
if args.big_bands:  # maximum size of bands on gel
    max_frag_size = int(args.big_bands)
else:
    max_frag_size = 4000 
# ---------------------
if args.band_gap:   # minimum size difference between bands on gel
    min_frag_size_dif = int(args.band_gap)
else:
    min_frag_size_dif = 100 
# ---------------------
if args.no_insert:   # if there is no insert, this option is true
    insert = False
else: 
    insert = True
# ---------------------
if args.fasta:   # if this is a fasta file, this option is true
    fasta = True
    insert = False
else:
    fasta = False
# --------------------- 
if args.test:   # location of test plasmid file
    plasmid_file = str(sys.path[0] + "/test_plasmid.gb")  
# --------------------- 
    
    
# -------------------------------- USER-DEFINED FUNCTIONS --------------------------------
def cuts(enzyme_cut_pattern_list, cut_from_5, plasmid_seq):
    cut_sites = []
    for unique in enzyme_cut_pattern_list:
        rc_unique = str(Seq(unique).reverse_complement())
        index = 0
        while index < len(plasmid_seq):
            index = plasmid_seq.find(unique, index)
            if index != -1:
                cut_pos = index+cut_from_5
                cut_sites.append(cut_pos)
            if index == -1:
                index = plasmid_seq.find(rc_unique, index)
                if index != -1:
                    cut_pos = index+(len(unique)-cut_from_5)
                    cut_sites.append(cut_pos)
                if index == -1:
                    break
            index += 1
    cut_sites.sort()
    return cut_sites

def frag_calc(cut_list, plasmid_seq):
    list_len = len(cut_list)
    fragments = []
    for i in range(0, list_len):
        if i+1 == list_len:
            frag_size = (len(plasmid_seq)-cut_list[i]) + cut_list[0]
            fragments.append(frag_size)
        else:
            frag_size = cut_list[i+1] - cut_list[i]
            fragments.append(frag_size)
    fragments.sort()
    return fragments

def code_variants(input_list):
    output_list = []
    for sequence in input_list:
        match_obj = re.findall("[^ATGC]",sequence)
        if match_obj:
            for x in range (0,len(match_obj)):
                WeirdDict = {'B': ['C','G','T'], 'D': ['A','G','T'], 'H': ['A','C','T'],\
                             'K': ['G','T'], 'M': ['A','C'],'N': ['A','C','G','T'],'R': ['A','G'],\
                             'S': ['C','G'],'V': ['A','C','G'],'W': ['A','T'],'Y': ['C','T']}
                output_list = []
                for c in input_list: # loop for each code within input list
                    match_obj = re.search("[^ATGC]",c)
                    if match_obj != []:
                        key = match_obj[0]
                        substitution_list = WeirdDict[key]
                        for sub in substitution_list: # loop for each nt sub in dictionary
                            output_list.append(c.replace(key, sub, 1))
                input_list = output_list
                match_obj = re.search("[^ATGC]",sequence)
    return output_list

def cuts2(ep, plasmid_seq):
    #generate fragment list
    enzyme_name = ep[0]       # enzyme 1 info
    loc_tuple = np.where(enzyme_dict == enzyme_name) # stores list w/enzyme cut pattern location
    row_tuple = int(loc_tuple[0]) # this is the row that the enzyme associated w/pattern will be
    enzyme_cut_pattern = enzyme_dict[row_tuple][1]
    cut_from_5 = int(enzyme_dict[row_tuple][2])
    if re.search('[^ATGC]',enzyme_cut_pattern):     #if there's a non-ATGC nt in the pattern list
        enzyme_cut_pattern_list = code_variants([enzyme_cut_pattern])
    else:                                           #if pattern list is all ATGC
        enzyme_cut_pattern_list = [enzyme_cut_pattern]       
    cut_list1 = cuts(enzyme_cut_pattern_list, cut_from_5, plasmid_seq)
    
    enzyme_name = ep[1]        # enzyme 2 info
    loc_tuple = np.where(enzyme_dict == enzyme_name) # stores list w/barcode seq location
    row_tuple = int(loc_tuple[0]) # this is the row that the s_name associated w/seq will be
    enzyme_cut_pattern = enzyme_dict[row_tuple][1]
    cut_from_5 = int(enzyme_dict[row_tuple][2])
    if re.search('[^ATGC]',enzyme_cut_pattern):     #if there's a non-ATGC nt in the pattern list
        enzyme_cut_pattern_list = code_variants([enzyme_cut_pattern])
    else:                                           #if pattern list is all ATGC
        enzyme_cut_pattern_list = [enzyme_cut_pattern]   
    cut_list2 = cuts(enzyme_cut_pattern_list, cut_from_5, plasmid_seq)

    cut_list1.extend(cut_list2)
    cut_list1.sort()
    return cut_list1

def gene_cut(cut_list, gb_record):
    for feature in gb_record.features:
        floc = feature.location
        g_start = int(floc.start)
        g_end = int(floc.end)
        dxflag = 0
        if feature.type == "gene":
            #if feature.strand == 1:     #forward orientation of insert
                #you could do something with sidedness here in the future
            for cut in cut_list:
                if cut in range(g_start,g_end):
                    dxflag = 1
    return dxflag
        
def frag_length_violation(frag_list):
    f_flag = 0
    for frag in frag_list:
        if frag < min_frag_size or frag > max_frag_size: # restricting band sizes
            f_flag = 1
    for i in range(1, len(frag_list)):
        difference = frag_list[i] - frag_list[i-1]
        if difference < min_frag_size_dif: # difference between band sizes
            f_flag = 1
    return f_flag

# -------------------------------------- MAIN MODULE --------------------------------------

#print parameters
print("\n\tDisplaying enzymes producing between",min_frags,"and",max_frags,"bands,")
print("\twith sizes between",min_frag_size,"nt and",max_frag_size,"nt,")
print("\twith at least",min_frag_size_dif,"nt gap between bands.\n")
#print("Mode:",mode)
#print("   d = check direction of gene insert"), 
#print("   i = check for gene insertion\n")

# import file containing the list of on-hand enzymes and convert array to a simple list
with open(my_enzymes_file, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    # get all the rows as a list
    my_enzymes = list(reader)
    # transform data into numpy array and save dimensions for easy use
    my_enzyme_array = np.array(my_enzymes).astype(str)
    my_enzyme_rows = my_enzyme_array.shape[0]
    my_enzyme_cols = my_enzyme_array.shape[1]   
my_enzyme_list = []
for i in range(0,my_enzyme_rows):
    for j in range(0,my_enzyme_cols):
        if my_enzyme_array[i,j] != "":
            my_enzyme_list.append(my_enzyme_array[i,j])
            
# import file containing enzyme dictionary
with open(enzyme_dict_file, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    # get header from first row
    headers = next(reader)
    # get all the rows as a list
    enzyme_data = list(reader)
    # transform data into numpy array and save dimensions for easy use
    enzyme_dict = np.array(enzyme_data).astype(str)
    enzyme_dict_rows = enzyme_dict.shape[0]
    enzyme_dict_cols = enzyme_dict.shape[1]

#check if all the enzymes in your freezer are in the dictionary
my_final_enzyme_list = []    
for e in my_enzyme_list:
    if np.count_nonzero(enzyme_dict == e) == 0:
        print("\tThere is no match for this enzyme in the dictionary:",e)
        print("\n\tCheck enzyme spelling/capitalization, then check dictionary file located here:")
        print("\t",enzyme_dict_file)
    else:                       # if the enzymes are in the dictionary, add to polished list
        my_final_enzyme_list.append(e)
   
# ------------------------------------- CODE -------------------------------------
# ------ if there is no insert ----------
if (insert == False):
    if (fasta == True):
        plasmid_record = SeqIO.read(plasmid_file, "fasta")
    else: 
        plasmid_record = SeqIO.read(plasmid_file, "genbank")
    plasmid_seq = (plasmid_record.seq).upper()
    print("----------------------------------------------------------")
    print("Your plasmid:",plasmid_record.id,"\n")
    short_list = []
    x_flag = 0          # flag for if there's a cut site in the sequence for the parameters
    cutting_enzymes_list = [] # list of enzymes that cut the plasmid, potential for use in 2-combo enzymes
    
    #generate fragment list for 1 enzyme
    for e in my_final_enzyme_list:
        enzyme_name = e
        loc_tuple = np.where(enzyme_dict == e) # stores list w/barcode seq location
        row_tuple = int(loc_tuple[0]) # this is the row that the s_name associated w/seq will be 
        enzyme_cut_pattern = enzyme_dict[row_tuple][1]
        
        if re.search('[^ATGC]',enzyme_cut_pattern):     #if there's a non-ATGC nt in the pattern list
            enzyme_cut_pattern_list = code_variants([enzyme_cut_pattern])
        else:                                           #if pattern list is all ATGC
            enzyme_cut_pattern_list = [enzyme_cut_pattern]
            
        cut_from_5 = int(enzyme_dict[row_tuple][2])
        cut_list = cuts(enzyme_cut_pattern_list, cut_from_5, plasmid_seq)
        if cut_list != []:  # (this is to prevent a 0-cuts enzyme from showing up as a 2-combo result)
            cutting_enzymes_list.append(enzyme_name) 
            
        frag_list = frag_calc(cut_list, plasmid_seq)
        if frag_length_violation(frag_list) == 1: #if there are bad fragment lengths
            frag_list = []
        if len(frag_list) <= max_frags and len(frag_list) > 0:
            if len(frag_list) >= min_frags:
                print("Bands from",enzyme_name,":\t",frag_list)
                x_flag = 1

    #generate fragment list for 2 enyzmes
    combo_list = combinations(cutting_enzymes_list, 2)
    for ep in combo_list:
        cut_list = cuts2(ep, plasmid_seq)
        if cut_list:
            frag_list = frag_calc(cut_list, plasmid_seq)
            if frag_length_violation(frag_list) == 1: #if there are bad fragment lengths
                frag_list = []
            else: 
                if len(frag_list) <= max_frags and len(frag_list) >= min_frags:   
                    print("Bands from",ep[0],"+",ep[1],":\t",frag_list)
                    x_flag = 1
    
    if x_flag == 0:
        print("\nNo eligible cut sites in this sequence. Relax parameters and try again.")
    
# ------ if there is an insert and you need to have a cut site in it ----------
else: 
    plasmid_record = SeqIO.read(plasmid_file, "genbank")
    plasmid_seq = plasmid_record.seq
    print("----------------------------------------------------------")
    print("Your plasmid:",plasmid_record.id)
    short_list = []
    g_flag = 0          # flag for if there's a gene in the sequence
    x_flag = 0          # flag for if there's a cut site in the sequence for the parameters
    cutting_enzymes_list = [] # list of enzymes that cut the plasmid, potential for use in 2-combo enzymes
    
    # display coordinates of insert
    for feature in plasmid_record.features:
        if feature.type == "gene":
            g_flag = 1
            floc = feature.location
            g_start = int(floc.start)
            g_end = int(floc.end)
            print("Your gene insert location and strand:",floc,"\n")
    if g_flag == 0:
        print("Error: no genes denoted in this sequence. Annotate gene feature in Ape and try again.")
         
    #generate fragment list for 1 enzyme
    for e in my_final_enzyme_list:
        enzyme_name = e
        loc_tuple = np.where(enzyme_dict == e) # stores list w/barcode seq location
        row_tuple = int(loc_tuple[0]) # this is the row that the s_name associated w/seq will be 
        enzyme_cut_pattern = enzyme_dict[row_tuple][1]
        
        if re.search('[^ATGC]',enzyme_cut_pattern):     #if there's a non-ATGC nt in the pattern list
            enzyme_cut_pattern_list = code_variants([enzyme_cut_pattern])
        else:                                           #if pattern list is all ATGC
            enzyme_cut_pattern_list = [enzyme_cut_pattern]
            
        cut_from_5 = int(enzyme_dict[row_tuple][2])
        cut_list = cuts(enzyme_cut_pattern_list, cut_from_5, plasmid_seq)
        if cut_list != []:  # (this is to prevent an enzyme that doesn't cut from showing up as a 2-combo result)
            cutting_enzymes_list.append(enzyme_name) 
        frag_list = frag_calc(cut_list, plasmid_seq)
        if frag_length_violation(frag_list) == 1: #if there are bad fragment lengths
            frag_list = []
        if len(frag_list) <= max_frags and len(frag_list) > 0:
            gene_dx = gene_cut(cut_list, plasmid_record)
            if gene_dx == 1: # if there are cuts within the gene
                if len(frag_list) >= min_frags:
                    print("Bands from",enzyme_name,":\t",frag_list)
                    x_flag = 1
    
    #generate fragment list for 2 enyzmes
    combo_list = combinations(cutting_enzymes_list, 2)
    for ep in combo_list:
        cut_list = cuts2(ep, plasmid_seq)
        gene_dx = gene_cut(cut_list, plasmid_record)
        if gene_dx == 1: # if there are cuts within the gene
            if cut_list:
                frag_list = frag_calc(cut_list, plasmid_seq)
                if frag_length_violation(frag_list) == 1: #if there are bad fragment lengths
                    frag_list = []
                else: 
                    if len(frag_list) <= max_frags and len(frag_list) >= min_frags:   
                        print("Bands from",ep[0],"+",ep[1],":\t",frag_list)
                        x_flag = 1
    
    if x_flag == 0:
        print("\nNo eligible cut sites in this sequence. Relax parameters and try again.")

print("----------------------------------------------------------")
print("\nAll done!\n")