#!/usr/bin/python3
"""
    This script allows the user to specify a list of restriction enzymes
    and a plasmid sequence and have the program automatically select the best
    one or two enzymes to use to get the best gel.
"""
import csv
from Bio import SeqIO
from Bio.Restriction import *
from itertools import combinations
import argparse
import os
from utils import frag_calc, frag_length_violation, gene_cut

parser = argparse.ArgumentParser()
parser.add_argument("--plasmid", "-p", help="(required) plasmid file, .gb or .ape, insert = \"gene\"")
parser.add_argument("--enzymes", "-e", help="enzymes file (.txt)")
parser.add_argument("--min_bands", "-n", help="minimum number of desired bands")
parser.add_argument("--max_bands", "-x", help="maximum number of desired bands")
parser.add_argument("--small_bands", "-s", help="smallest size of desired bands")
parser.add_argument("--big_bands", "-b", help="biggest size of desired bands")
parser.add_argument("--band_gap", "-g", help="minimum gap size between bands")
parser.add_argument("--insert", "-i", help="must cut within plasmid insert", action='store_true')
parser.add_argument("--test", "-t", help="tests program with .gb test file", action='store_true')


if __name__=="__main__":
    pfile_extensions = [".fa", ".fas", ".fasta", ".gb", ".gbk", ".genbank", ".ape"]

    # Read arguments from the command line
    args = parser.parse_args()
    if args.plasmid:
        plasmid_file = args.plasmid
        plasmid_name = os.path.basename(plasmid_file).split('.')[0]
        extension = os.path.splitext(plasmid_file)[1]
        if extension not in pfile_extensions:
            print("Filetype is not supported. Choose fasta- or genbank-formatted file.")
            exit()
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
        min_seg_size = int(args.small_bands)
    else:
        min_seg_size = 500
    # ---------------------    
    if args.big_bands:  # maximum size of bands on gel
        max_seg_size = int(args.big_bands)
    else:
        max_seg_size = 4000 
    # ---------------------
    if args.band_gap:   # minimum size difference between bands on gel
        min_spacing = int(args.band_gap)
    else:
        min_spacing = 100 
    # ---------------------
    if args.insert:   # if there is a cloned insert, this option is true
        insert = True
    else: 
        insert = False
    # --------------------- 
    if args.enzymes:
        my_enzymes_file = args.enzymes
    else:
        my_enzymes_file = "./test/my_enzymes.txt"
    # --------------------- 
    if args.test:   # location of test plasmid file
        plasmid_file = "./test/test_plasmid.gb"
        plasmid_name = os.path.basename(plasmid_file).split('.')[0]
        extension = os.path.splitext(plasmid_file)[1]
    # ---------------------
    if extension == ".fa" or extension == ".fas" or extension == ".fasta":
        pfile_type = "fasta"
    elif extension == ".gb" or extension == ".gbk" or extension == ".genbank" or extension == ".ape":
        pfile_type = "genbank"
    # ---------------------


    # print parameters
    print("================================================================\n")
    print("\tDisplaying enzymes producing between",min_frags,"and",max_frags,"bands,")
    print("\twith sizes between",min_seg_size,"nt and",max_seg_size,"nt,")
    print("\twith at least",min_spacing,"nt gap between bands.\n")
    print("----------------------------------------------------------------")

    # load plasmid sequence
    print("Your plasmid:", plasmid_name)
    plasmid_record = SeqIO.read(plasmid_file, pfile_type)
    plasmid_seq = plasmid_record.seq

    # establish your enzymes on hand
    with open(my_enzymes_file, 'r') as f:
        enzyme_text = f.read()
    my_enzymes_list = enzyme_text.split("\n")
    rb = RestrictionBatch()
    for enzyme in my_enzymes_list:
        try: 
            rb += enzyme
        except ValueError:
            print(f"\nInvalid enzyme: '{enzyme}'\n")
    rb_analysis = Analysis(rb, plasmid_seq, linear=False)
    CutDict = rb_analysis.full()

    # analyze single enzyme cut patterns
    x_flag = 0 # assumes there are 0 cut sites in list until proven otherwise
    cutting_enzymes_list = [] # list to store enzymes that cut
    for enz1 in rb:
        one_enz_cut_pattern = CutDict[enz1]
        if one_enz_cut_pattern:
            digested_segments = enz1.catalyze(plasmid_seq, linear=False)
            num_of_segments = len(digested_segments)
            segment_len_list = []        
            #print("cut positions:", one_enz_cut_pattern) 
            # check if fragments satisfy parameters
            if insert == True:
                print("You indicated that there is an insert needing to be cut.")
                gene_dx = gene_cut(one_enz_cut_pattern, plasmid_record)
                if gene_dx != 1: # if there are no cuts within the gene insert, exit
                    print("There is no cut within the insert annotated as \"gene\".")
                    raise SystemExit
            if num_of_segments > 0:
                cutting_enzymes_list.append(enz1)
            for segment in digested_segments:
                segment_len = len(str(segment))
                segment_len_list.append(segment_len)
            segment_len_list.sort()
            if frag_length_violation(segment_len_list, min_seg_size, max_seg_size, min_spacing):
                segment_len_list = []
            else:
                if len(segment_len_list) <= max_frags and len(segment_len_list) >= min_frags:   
                    print("Bands from",enz1,":\t\t",segment_len_list)
                    x_flag = 1
        
    # analyze double enzyme cut patterns
    # create list of every two-enzyme combo
    combo_list = combinations(cutting_enzymes_list, 2)
    for enz_pair in combo_list:
        two_enz_cut_pattern_1of2 = CutDict[enz_pair[0]]
        two_enz_cut_pattern_2of2 = CutDict[enz_pair[1]]
        two_enz_cut_pattern = two_enz_cut_pattern_1of2 + two_enz_cut_pattern_2of2
        two_enz_cut_pattern.sort()
        segment_len_list = []
        #print("cut positions:", two_enz_cut_pattern)
        # check if fragments satisfy parameters
        if two_enz_cut_pattern:
            if insert == True:
                gene_dx = gene_cut(two_enz_cut_pattern, plasmid_record)
                if gene_dx != 1: # if there are no cuts within the gene insert, exit
                    raise SystemExit
            segment_len_list = frag_calc(two_enz_cut_pattern, plasmid_seq)
            if frag_length_violation(segment_len_list, min_seg_size, max_seg_size, min_spacing):
                segment_len_list = []
            else:
                if len(segment_len_list) <= max_frags and len(segment_len_list) >= min_frags:   
                    print("Bands from",enz_pair[0],"+",enz_pair[1],":\t",segment_len_list)
                    x_flag = 1
    if x_flag == 0:
        print("No eligible cut sites in this sequence. Relax parameters and try again.\n")

    print("----------------------------------------------------------------\n")
