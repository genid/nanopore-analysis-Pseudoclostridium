
import numpy as np
import os
from itertools import combinations
from collections import Counter


def get_seqs_motifs(dict_motifs, motif, indexes_match_list):
    
    list_seqs = []
    for value in dict_motifs.values():        
        for v in value:
            s = ""     
            c = 0            
            for i in range(len(motif)):
                #print(v[c])                
                if i in indexes_match_list:
                    s += v[c]
                    c += 1
                else:
                    s += motif[i]            
            list_seqs.append(s)
    return list_seqs            


def get_seqs_motifs_mod(dict_motifs, list_seqs, base):    
    update = []
    for value in dict_motifs.values():                
        for v in value:                                               
            v = v[0]                
            for list_seq in list_seqs:                                                
                update.append(list_seq.replace(base,v))                                
    return update


motifs_list = ["CACNNNNNNNTNGC","GCNANNNNNNNGTG","GATNNNNCTC","DGAGNNNNATC","TCABNNNNNNTARG","CYTANNNNNNVTGA"]

all_motifs_dict = {}
for motif in motifs_list:    
    list_seqs = []
    
    if len(np.where(np.array(list(motif)) == "N")[0]) > 0:            
        dict_motifs = {}
        indexes_match_list = (np.where(np.array(list(motif)) == "N")[0])        
        t = "ACGT" * len(indexes_match_list)
        comb = list(combinations(t,len(indexes_match_list)))        
        dict_motifs[motif] = list(set(comb))    
        list_seqs = get_seqs_motifs(dict_motifs, motif, indexes_match_list)    
    
    if len(np.where(np.array(list(motif)) == "B")[0]) > 0:            
        dict_motifs = {}
        indexes_match_list = (np.where(np.array(list(motif)) == "B")[0])                
        t = "CGT" * len(indexes_match_list)
        comb = list(combinations(t,len(indexes_match_list)))
        dict_motifs[motif] = list(set(comb))                
        list_seqs = get_seqs_motifs_mod(dict_motifs, list_seqs, "B")                
    
    if len(np.where(np.array(list(motif)) == "D")[0]) > 0:            
        dict_motifs = {}
        indexes_match_list = (np.where(np.array(list(motif)) == "D")[0])                
        t = "AGT" * len(indexes_match_list)
        comb = list(combinations(t,len(indexes_match_list)))
        dict_motifs[motif] = list(set(comb))        
        list_seqs = get_seqs_motifs_mod(dict_motifs, list_seqs, "D")    
    
    if len(np.where(np.array(list(motif)) == "V")[0]) > 0:            
        dict_motifs = {}
        indexes_match_list = (np.where(np.array(list(motif)) == "V")[0])        
        t = "ACG" * len(indexes_match_list)
        comb = list(combinations(t,len(indexes_match_list)))
        dict_motifs[motif] = list(set(comb))        
        list_seqs = get_seqs_motifs_mod(dict_motifs, list_seqs, "V")    
    
    if len(np.where(np.array(list(motif)) == "Y")[0]) > 0:            
        dict_motifs = {}
        indexes_match_list = (np.where(np.array(list(motif)) == "Y")[0])        
        t = "CT" * len(indexes_match_list)
        comb = list(combinations(t,len(indexes_match_list)))
        dict_motifs[motif] = list(set(comb))                
        list_seqs = get_seqs_motifs_mod(dict_motifs, list_seqs, "Y")    
    
    if len(np.where(np.array(list(motif)) == "R")[0]) > 0:            
        dict_motifs = {}
        indexes_match_list = (np.where(np.array(list(motif)) == "R")[0])        
        t = "AG" * len(indexes_match_list)
        comb = list(combinations(t,len(indexes_match_list)))
        dict_motifs[motif] = list(set(comb))                
        list_seqs = get_seqs_motifs_mod(dict_motifs, list_seqs, "R")        
    all_motifs_dict[motif] = list_seqs
    
    
for k,value in all_motifs_dict.items():        
    filename = k+".txt"
    print(filename)
    with open(filename,"w") as f:
        for seq in value:            
            header = ">"+seq            
            f = open(filename, "a")            
            f.write(header)
            f.write("\n")
            f.write(seq)
            f.write("\n")                        