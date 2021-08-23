

import multiprocessing

import pandas as pd
import numpy as np

import os

from functools import partial
from collections import Counter


def get_strand(bool_array):
    if bool_array == True:
        return "+"
    else:
        return "-"
   
    
def parse_blastfmt7(out_blast):    
    dictionary = {} #Accession as key, sequence as value
    i = 0
    for entry in out_blast:        
        if "#" not in entry:
            line = entry.split("\t")            
            dictionary[i] = line[0:]
        i+=1
    return dictionary


def process_output_blast(dictionary):
    # Fields: subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    df = pd.DataFrame(dictionary).T
    df[2] = np.array(df[2], dtype=float)
    df[3] = np.array(df[3], dtype=int)
    df[4] = np.array(df[4], dtype=int)
    df[5] = np.array(df[5], dtype=int)
    df[6] = np.array(df[6], dtype=int)
    df[7] = np.array(df[7], dtype=int)
    df[8] = np.array(df[8], dtype=int)
    df[9] = np.array(df[9], dtype=int)
    df[10] = np.array(df[10], dtype=float)
    df[11] = np.array(df[11], dtype=float)
    df["strand"] = list(map(get_strand, (df[8] < df[9]) ))
    df["motif"] = output_blast
    return df


def process_wig(wig_file):
    wig = [[line.rstrip('\n')] for line in open(wig_file)]
    wig = wig[2:]        
    wig = pd.DataFrame(wig)
    wig[['pos','mod']] = wig[0].str.split(" ",expand=True,)
    wig = wig.drop([0], axis=1)
    wig['pos'] = np.array(wig['pos'], dtype=int)
    wig['mod'] = np.array(wig['mod'], dtype=float)
    return wig


def get_filter_pos(wig_files):
    dict_barcodes = []    
    for wig_file in wig_files:            
        df_wig = process_wig(wig_file)
        dict_barcodes.append(df_wig["pos"].values)                
    counter_wig = Counter(np.concatenate(dict_barcodes))
    list_filter_pos = []
    for key, value in counter_wig.items():
        if value == 5:
            list_filter_pos.append(key)                
    return list_filter_pos


def get_modifications_plus(df_tmp_blast, mismatch, wig):                

    row = df_tmp_blast[(df_tmp_blast[8] <= (wig[0])-1) & (df_tmp_blast[9] > (wig[0])-1) ]        
    if not row.empty:              
        len_motif = len(row[0].values[0])
        len_hit = row[3].values[0]
        start_hit = row[6].values[0]                
        if start_hit == 1 and (len_motif-len_hit) <= mismatch:                
            return wig[1]    
        
        
def get_modifications_minus(df_tmp_blast, mismatch, wig):                
         
    row = df_tmp_blast[(df_tmp_blast[9] <= (wig[0])-1) & (df_tmp_blast[8] > (wig[0])-1) ]        
    if not row.empty:              
        len_motif = len(row[0].values[0])
        len_hit = row[3].values[0]        
        end_hit = row[7].values[0]                        
        if end_hit == len_motif and (len_motif-len_hit) <= mismatch:                
            return wig[1]        
        
        

if __name__ == "__main__":
    
    motifs_dirpath = os.sys.argv[1] #motifs/    
    mismatch = int(os.sys.argv[2]) #0
    n_threads = int(os.sys.argv[3]) #10               
    
    wig_files_plus = [   "barcode01.fraction_modified_reads.plus.wig",
                         "barcode02.fraction_modified_reads.plus.wig", 
                         "barcode03.fraction_modified_reads.plus.wig", 
                         "barcode04.fraction_modified_reads.plus.wig", 
                         "barcode05.fraction_modified_reads.plus.wig"]

    wig_files_minus = [  "barcode01.fraction_modified_reads.minus.wig",
                         "barcode02.fraction_modified_reads.minus.wig",
                         "barcode03.fraction_modified_reads.minus.wig",
                         "barcode04.fraction_modified_reads.minus.wig",
                         "barcode05.fraction_modified_reads.minus.wig"]
    
        
    
    ## filter positions by overlapping 5 barcodes
    print("Filtering position by overlapping from the 5 barcodes...")
    list_filter_plus = get_filter_pos(wig_files_plus)
    list_filter_minus = get_filter_pos(wig_files_minus)

    
    list_out_blast = []
    dirpath = motifs_dirpath
    dirpath = os.walk(dirpath)
    for dirpath, dirnames, filenames in dirpath:
        for filename in [f for f in filenames if f.endswith(".out")]:            
            list_out_blast.append(os.path.join(dirpath, filename))
        
                
    for wig_plus, wig_minus in zip(wig_files_plus, wig_files_minus):    
        print("--- Processing {} and {} ---".format(wig_plus,wig_minus))   
        for output_blast in list_out_blast:                

            print(output_blast)     
            name_motif = output_blast.split("/")[-1]
            nameout_plus = wig_plus+"_"+name_motif                
            nameout_minus = wig_minus+"_"+name_motif                

            lines = [line.rstrip('\n') for line in open(output_blast)]
            dictionary = parse_blastfmt7(lines)
            df_blast_out = process_output_blast(dictionary)                                            

            # plus strand (+)
            tmp_plus = df_blast_out[df_blast_out["strand"] == "+"].sort_values(by=[8, 9])                                                
            df_wig = process_wig(wig_plus)
            df_wig = df_wig[df_wig["pos"].isin(list_filter_plus)]
            print("Starting + ...")
            func = partial(get_modifications_plus,tmp_plus, mismatch)
            pool = multiprocessing.Pool(n_threads)                                        
            list_a_mod = list(pool.map(func, df_wig.values))

            list_a_mod = np.array(list_a_mod)
            list_a_mod = list_a_mod[list_a_mod != np.array(None)]
            list_a_mod = list_a_mod.astype('float64')         
            pd.DataFrame(list_a_mod).to_csv(nameout_plus, header=None, index=None)        
            print("Finished + ...")        

            # minus strand (-)
            tmp_minus = df_blast_out[df_blast_out["strand"] == "-"].sort_values(by=[8, 9])
            df_wig = process_wig(wig_minus)
            df_wig = df_wig[df_wig["pos"].isin(list_filter_minus)]
            print("Starting - ...")
            func = partial(get_modifications_minus,tmp_minus, mismatch)
            pool = multiprocessing.Pool(n_threads)                                        
            list_a_mod = list(pool.map(func, df_wig.values))
            list_a_mod = np.array(list_a_mod)
            list_a_mod = list_a_mod[list_a_mod != np.array(None)]
            list_a_mod = list_a_mod.astype('float64') 
            pd.DataFrame(list_a_mod).to_csv(nameout_minus, header=None, index=None)        
            print("Finished - ...")
            
            break
        break
        