import csv
from Bio import SeqIO
from Bio.Restriction import *
from itertools import combinations


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

def frag_length_violation(frag_list, min_frag_size, max_frag_size, min_frag_size_dif):
    f_flag = False
    for frag in frag_list:
        if frag < min_frag_size or frag > max_frag_size: # restricting band sizes
            f_flag = True
    for i in range(1, len(frag_list)):
        difference = frag_list[i] - frag_list[i-1]
        if difference < min_frag_size_dif: # difference between band sizes
            f_flag = True
    return f_flag

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
