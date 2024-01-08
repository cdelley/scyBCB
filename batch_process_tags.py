#!/usr/bin/env python

# Cyrille L. Delley 2021
#
# usage:
#python3 batch_process_tags.py -i /your_path_to_fastq_files \
#-bc /your_path/git/scyBCB/whitelists/bc1_whitelist_num.json \
#-bc /your_path/git/scyBCB/whitelists/bc2_whitelist_num.json \
#-bc /your_path/git/scyBCB/whitelists/bc2_whitelist_num.json \
#-o t --feature /your_path/git/scyBCB/whitelists/TotalSeqA_Ab_to_num.json \
#--sample_name SAMPLE_NAME --split_name _Ab_ --extract_cells 0 --anndata 0
#
#
#

import os
import sys
import argparse
import logging
import itertools
import subprocess
import pandas as pd
import json
import time
import matplotlib
import numpy as np
import CyFeat as dabm
from  config_beads import barcode_beads

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='''
    
    dab-seq: single-cell dna genotyping and antibody sequencing
    
    ''', formatter_class=argparse.RawTextHelpFormatter)
    
    parser = argparse.ArgumentParser(description='Demultiplexing of fastq files')
    parser.add_argument('--input', '-i', help='folder containing fastq files', metavar='./fastq', type=str, required=True)
    parser.add_argument('--bc_whitelist', '-bc', help = 'barcode list', metavar = 'barcode1.json barcode2,json etc.', action='append', type=str, required=True)
    parser.add_argument('--output', '-o', help = "Output directory to write individual fastq files to.", type = str, metavar = 'fastq', default = None)
    parser.add_argument('--feature', help = 'provide a feature correction dic in .json foramt', type=str, default=None)
    parser.add_argument('--sample_name', help = "A suffix to append to individual fastq files.", type = str, metavar = 'tube-1', default = '')    
    parser.add_argument('--split_name', help = "split fastq names by key word", type = str, metavar = 'tube-1', default = None)       
    parser.add_argument('--beads', help = "BHBv2_Ab (default) or BHBv2_spat, 10x_spat", type = str, metavar = 'BHBv2_Ab', default = 'BHBv2_Ab')
    parser.add_argument(
        '--extract_cells', 
        help = "if cell barcodes should be extrracted form fastq files and condensed files should be created", 
        type = int, 
        default = 1
    )
    parser.add_argument(
        '--anndata', 
        help = "if anndata output should be generated (colapsing UMIS with the 'directional' method from umi-tools)", 
        type = int, 
        default = 1
    )
    
    args = parser.parse_args()
    args.output = os.path.abspath(args.output) if args.output != None else os.path.abspath(os.getcwd())
    
    try:
        seq_pos_dict = barcode_beads[args.beads]
    except KeyError as exept:
        raise Exception("""
            unsupported --beads option used. Currently supported are: {} but {} was used \n
            you can edit these coordinates and add new ones in the config_beads file.py
            """.format(
                list(barcode_beads.keys()), args.beads
            )
        ) from e
    
    #_________________________________________________________________________________________________________
    # get the fstq files for processing
    
    fastq_R1_files = [os.path.join(args.input, f) for f in np.sort(os.listdir(args.input)) if '_R1' in f] 
    fastq_R2_files = [os.path.join(args.input, f) for f in np.sort(os.listdir(args.input)) if '_R2' in f]
    if len(args.sample_name) != len(fastq_R1_files):
        if args.split_name:
            args.sample_name = [f.split('/')[-1].split(args.split_name)[0].replace('_','') for f in fastq_R1_files]
            print(' '.join(args.sample_name))
        else:
            print('not enough sample names provided, and no split_name string given. Trying to split by the word: _R1')
            args.sample_name = [f.split('/')[-1].split('_R1')[0].replace('_','') for f in fastq_R1_files] 
            print(' '.join(args.sample_name))
    
    #_________________________________________________________________________________________________________
    # load the barcode and feature whitelists
         
    barcodes = args.bc_whitelist
    bc_dicts = []
    for bc_json in barcodes:
        with open(bc_json, 'r') as fin:
            bc_dicts.append(json.load(fin))
            
    
    with open(args.feature , 'r') as f_dict_in:
        feature_dic = json.load(f_dict_in)  
    
    #_________________________________________________________________________________________________________
    # extract reads
    cond_files = []
    
    
    for f1, f2, sample_basename in zip(fastq_R1_files, fastq_R2_files, args.sample_name):
        if args.extract_cells:
            print(sample_basename)
            start_time = time.time()
            exp = dabm.features.assign_features(file1=f1,
                file2=f2, 
                bc_correction_dic=bc_dicts,
                feature_seq_dic=feature_dic,
                seq_positions=seq_pos_dict,
                    verbose=True, down_sample=10**5)
            print('barcode extraction took: {}s\n'.format(np.round(time.time() - start_time,0)))

    #_________________________________________________________________________________________________________
    # write condensed file
                    
        if args.sample_name == '':
            name = './compacted_files/' + args.fastq.split('_R1_')[0] + '.cond.ascii'
        else:
            name = './compacted_files/' + sample_basename + '.cond.ascii'
        cond_files.append(name)
        
        if args.extract_cells:
            start_time = time.time()
            try:
                exp.write_condensed_file(filename=name)
            except OSError:
                os.mkdir('./compacted_files')
                exp.write_condensed_file(filename=name)
            print('umi collapse and file writing took: {}s\n'.format(np.round(time.time() - start_time,0)))
            del exp
            
    #_________________________________________________________________________________________________________
    # colapse umis and produce anndata table
    if args.anndata:
        sample_name = [f.split('/')[-1].split('.')[0] for f in cond_files] 
        print(' '.join(sample_name))
        
        for f1, sample_basename in zip(cond_files, sample_name):
            print('\n'+sample_basename)
            
            start_time = time.time()
            adat = dabm.count_table_feature(filename=f1)#, condensed_file=f1.split('/')[0]+'/Hcomp_'+f1.split('/')[1])#, merge_threshold=2)
            print('UMI hamming collapsing took: {}s\n'.format(np.round(time.time() - start_time,0)))
            
            #try:
            start_time = time.time()
            name = './adata_files/' + sample_basename + '.adata.h5ad'
            try:
                adat.write(filename=name, compression='gzip')
            except OSError:
                os.mkdir('./adata_files')
                adat.write(filename=name, compression='gzip')
            print('adata writing took: {}s\n'.format(np.round(time.time() - start_time,0)))
