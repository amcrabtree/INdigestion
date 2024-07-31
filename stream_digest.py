"""
    This script makes a streamlit app to aid with restriction digest selection.
"""
from Bio import SeqIO
from Bio.Restriction import *
from itertools import combinations
import streamlit as st
from io import StringIO
from utils import frag_calc, frag_length_violation, gene_cut
import pandas as pd
# from streamlit_feedback import streamlit_feedback

# def _submit_feedback(user_response, emoji=None):
#     st.toast(f"Feedback submitted: {user_response}", icon=emoji)
#     return user_response.update({"some metadata": 123})
    
st.set_page_config(page_title="INdigestion")

# File input
st.header("File input")

seq_file_extensions = ["fa", "fas", "fasta", "gb", "gbk", "genbank", "ape"]
plasmid_file = st.file_uploader("##### Plasmid file:", type=seq_file_extensions)
st.link_button("Example file", "https://github.com/amcrabtree/INdigestion/blob/master/test/test_plasmid.gb")

enzyme_file = st.file_uploader("##### Enzymes file:", type=["csv", "tsv", "txt"])
st.link_button("Example file", "https://github.com/amcrabtree/INdigestion/blob/master/test/my_enzymes.txt")

# Sidebar
with st.sidebar:
    num_bands = st.slider("Number of desired bands", 1, 15, value=(2, 6), step=1)
    bp_size = st.slider("Size of desired bands (bp)", 25, 6000, value=(500, 4000), step=25)
    min_spacing = st.slider("minimum gap size between bands", 25, 1000, value=100, step=25)
    insert = st.toggle("Cut inside annotated insert?")
    # feedback = streamlit_feedback(feedback_type="thumbs",
    #                               optional_text_label="Feedback details",
    #                               on_submit=_submit_feedback)

# Digest
st.header("Digests")
if plasmid_file is not None and enzyme_file is not None:
    plasmid_io = StringIO(plasmid_file.getvalue().decode("utf-8"))
    enzyme_text = enzyme_file.read().decode("utf-8")

    ext = plasmid_file.name.split(".")[-1]
    if (ext == "fa") or (ext == "fas") or (ext == "fasta"):
        pfile_type = "fasta"
    elif (ext == "gb") or (ext == "gbk") or (ext == "genbank") or (ext == "ape"):
        pfile_type = "genbank"

    # load plasmid sequence
    plasmid_name = str(plasmid_file.name).split(".")[0]
    st.write(f"##### Your plasmid:", plasmid_name)
    plasmid_record = SeqIO.read(plasmid_io, pfile_type)
    plasmid_seq = plasmid_record.seq

    # Load enzymes list into RestrictionBatch
    my_enzymes_list = enzyme_text.split("\n")
    rb = RestrictionBatch()
    for enzyme in my_enzymes_list:
        try: 
            rb += enzyme
        except ValueError:
            st.write(f"Invalid enzyme: '{enzyme}'")
    rb_analysis = Analysis(rb, plasmid_seq, linear=False)
    CutDict = rb_analysis.full()

    # Analyze single enzyme cut patterns
    digest_dict = {'enzyme1':[], 'enzyme2':[], 'bands':[]}
    x_flag = 0 # assumes there are 0 cut sites in list until proven otherwise
    cutting_enzymes_list = [] # list to store enzymes that cut
    
    for enz1 in rb:
        one_enz_cut_pattern = CutDict[enz1]
        if one_enz_cut_pattern:
            digested_segments = enz1.catalyze(plasmid_seq, linear=False)
            num_of_segments = len(digested_segments)
            segment_len_list = []        

            # check if fragments satisfy parameters
            if insert == True:
                st.write("You indicated that there is an insert needing to be cut.")
                gene_dx = gene_cut(one_enz_cut_pattern, plasmid_record)
                if gene_dx != 1: # if there are no cuts within the gene insert, exit
                    st.write("There is no cut within the insert annotated as \"gene\".")
                    exit()

            if num_of_segments > 0:
                cutting_enzymes_list.append(enz1)

            for segment in digested_segments:
                segment_len = len(str(segment))
                segment_len_list.append(segment_len)

            segment_len_list.sort()
            if frag_length_violation(segment_len_list, bp_size[0], bp_size[1], min_spacing):
                segment_len_list = []
                
            else:
                if len(segment_len_list) <= num_bands[1] and len(segment_len_list) >= num_bands[0]:   
                    digest_dict['enzyme1'].append(enz1)
                    digest_dict['enzyme2'].append("")
                    digest_dict['bands'].append(str(segment_len_list))
                    x_flag = 1
        
    # Analyze double enzyme cut patterns
    # create list of every two-enzyme combo
    combo_list = combinations(cutting_enzymes_list, 2)
    for enz_pair in combo_list:
        two_enz_cut_pattern_1of2 = CutDict[enz_pair[0]]
        two_enz_cut_pattern_2of2 = CutDict[enz_pair[1]]
        two_enz_cut_pattern = two_enz_cut_pattern_1of2 + two_enz_cut_pattern_2of2
        two_enz_cut_pattern.sort()
        segment_len_list = []

        # check if fragments satisfy parameters
        if two_enz_cut_pattern:
            if insert == True:
                gene_dx = gene_cut(two_enz_cut_pattern, plasmid_record)
                if gene_dx != 1: # if there are no cuts within the gene insert, exit
                    raise SystemExit
            segment_len_list = frag_calc(two_enz_cut_pattern, plasmid_seq)
            if frag_length_violation(segment_len_list, bp_size[0], bp_size[1], min_spacing):
                segment_len_list = []
            else:
                if len(segment_len_list) <= num_bands[1] and len(segment_len_list) >= num_bands[0]:   
                    digest_dict['enzyme1'].append(enz_pair[0])
                    digest_dict['enzyme2'].append(enz_pair[1])
                    digest_dict['bands'].append(str(segment_len_list))
                    x_flag = 1
    if x_flag == 0:
        print("No eligible cut sites in this sequence. Relax parameters and try again.\n")
    
    df = pd.DataFrame(digest_dict)
    st.dataframe(df, use_container_width=True, hide_index=True)
