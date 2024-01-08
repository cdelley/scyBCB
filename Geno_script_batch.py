#!/usr/bin/env python

#
#original script from Ben Demaree and Cyrille Delley, DAb-seq, Abate lab.
#expanded by Cyrille L. Delley
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
import CyGeno as dabm
import resources
import copy
from multiprocessing.pool import ThreadPool
from multiprocessing import Process

def chunker_list(seq, size):
    return (seq[i::size] for i in range(size))

def wait(processes):
    # waits for processes to finish
    return [process.communicate() for process in processes]
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='''
    
    dab-seq: single-cell dna genotyping and antibody sequencing
    
    ''', formatter_class=argparse.RawTextHelpFormatter)
    
    
    parser = argparse.ArgumentParser(description='Demultiplexing of fastq files')
    parser.add_argument('--input', '-i', help='folder containing fastq files', metavar='./fastq', type=str, required=True)
    parser.add_argument('--bc_whitelist', '-bc', help = 'barcode list', metavar = 'barcode1.json barcode2,json etc.', action='append', type=str, required=True)
    parser.add_argument('--output', '-o', help = "Output directory to write individual fastq files to.", type = str, metavar = 'fastq', default = None)
    parser.add_argument('--feature', help = 'fprovide a feature correction dic in .json foramt', action='append', type=str, default=None) # might be deducable later from imput feature dict
    parser.add_argument('--sample_name', help = "A suffix to append to individual fastq files.", action='append', type = str, metavar = 'tube-1', default = '') 
    parser.add_argument('--split_name', help = "split fastq names by key word", type = str, metavar = 'tube-1', default = None)   
    parser.add_argument('--ex', help = "bool sequence for what to perform (barcode extraction, valid cell calling, demux by cell, aligning cells, variant calling)", 
    type = str, metavar = '11111', default = 111111)
    parser.add_argument('--cfg_file', type=str, help='config filename')
    parser.add_argument('--skip-flt3', action='store_true', default=True, help='option to skip FLT3-ITD calling')
        
    args = parser.parse_args()
    args.output = os.path.abspath(args.output) if args.output != None else os.path.abspath(os.getcwd())
    cfg_f = args.cfg_file
    
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
    for i in args.sample_name:
        out_folder = os.path.join('./', i)
        if not os.path.exists(out_folder):
                os.mkdir(out_folder)
    
    seq_pos_dict = {'bc'  : [(0,0,11), (0,12,20), (0,24,32)], # three barcodes
                    'umi' : [(0,45,57)], #(0,16,28), (1,0,6) 
                    'seq' : [(0,65,150), (1,18,150)],
                    'feature' : [(0,51,69), (1,0,18)]}
             
    barcodes = args.bc_whitelist
    bc_dicts = []
    for bc_json in barcodes:
        with open(bc_json, 'r') as fin:
            bc_dicts.append(json.load(fin))
            
    feature_dic_list = []
    for feat_json in args.feature:      
        with open(feat_json , 'r') as f_dict_in:
            feature_dic_list.append(json.load(f_dict_in))  
    
    for f1, f2, sample_basename in zip(fastq_R1_files, fastq_R2_files, args.sample_name):
        
        if int(args.ex[0]):
            # Write barcode in header and identify HIV reads / human reads based on primer sequences
            start_time = time.time()
            exp = dabm.seq_exp.write_demux(file1=f1,
                                                file2=f2, 
                                                bc_correction_dic=bc_dicts,
                                                seq_positions=seq_pos_dict,
                                                feature_dic=feature_dic_list, 
                                                batch_label=sample_basename, 
                                                out_path=args.output,
                                                verbose=True)
            print('barcode extraction took: {}s\n'.format(np.round(time.time() - start_time,0)))
            
            # Export primer count statistics
            name = './adata_files/amplicon_count_' + sample_basename + '.adata.h5ad'
            try:
                exp.adata.write(filename=name, compression='gzip')
            except OSError:
                os.mkdir('./adata_files')
                exp.adata.write(filename=name, compression='gzip')
            
            # merge fastq for pipeline mapping 1
            try:
                os.makedirs('./fastq_out')
            except FileExistsError:
                # directory already exists
                pass
                
            print('merging files...')
            start = time.time()
            process = subprocess.run('cat {}/{}*_R1.fastq.gz > fastq_out/{}_R1.fastq.gz'.format(args.output, sample_basename, sample_basename), shell=True)
            process = subprocess.run('cat {}/{}*_R2.fastq.gz > fastq_out/{}_R2.fastq.gz'.format(args.output, sample_basename, sample_basename), shell=True)
            print('fastq merging took: {}'.format(time.time()-start))
        
        
        ##############################################
        ##############################################
        
        #Start the DAbseq pipeline
        print('Initializing pipeline...       ')

        # load config file variables
        # be careful about using exec
        if not os.path.isfile(cfg_f):
            print('Config file not found! Please check the file name and path.')
            raise SystemExit
        
        else:
            with open(cfg_f, 'r') as cfg:
                for line in cfg:
                    if line[0] == '#' or line[0] == '[' or line[0] == ' ':
                        continue
                    else:
                        var = line.split("#", 1)[0].strip()  # to remove inline comments
                        exec(var)
                        
        all_vars = copy.copy(globals())
        # create all other directories for this run
        # if it already exists, ignore and continue
        dirs = [all_vars[d] for d in all_vars if '_dir' in d]
        dirs.sort()
        for d in dirs:
            if not os.path.exists(d):
                os.mkdir(d)
        
        # get reads for all panel samples from read 1
#        panel_fastq_R1 = './fastq_out/{}_R1.fastq.gz'.format(sample_basename)
#        panel_r1_files = [os.path.abspath(panel_fastq_R1)]
#        all_tsv_r1 = barcode_dir + sample_basename + '_r1.all.tsv'  # alignment counts for all barcodes
        
        # get reads for all panel samples from read 2
        panel_fastq_R2 = './fastq_out/{}_R2.fastq.gz'.format(sample_basename)
        panel_r2_files = [os.path.abspath(panel_fastq_R2)]
        all_tsv_r2 = barcode_dir + sample_basename + '_r2.all.tsv'  # alignment counts for all barcodes
                           
        if int(args.ex[1]):
            print('''
            ###################################################################################
            Step 4: count read alignments to inserts
            ###################################################################################
            ''')
            
            # align r1 reads to inserts to obtain read counts across all barcodes
#            resources.count_alignments(panel_r1_files, amplicon_file, human_fasta_file, all_tsv_r1, temp_dir)
            
            # align r2 reads to inserts to obtain read counts across all barcodes
            resources.count_alignments(panel_r2_files, amplicon_file, human_fasta_file, all_tsv_r2, temp_dir)
        
        if int(args.ex[2]) or int(args.ex[4]) or int(args.ex[5]):
            print('''
            ###################################################################################
            Step 5: call valid cells using selected method
            ###################################################################################
            ''')
            
            # call valid cells using cell_caller function            
            # load alignment df from all barcodes into dataframe
            all_df = pd.read_csv(all_tsv_r2, sep='\t', header=0, index_col=0)

            # extract barcodes and read totals
            barcodes = list(all_df.index)
            reads_per_cell = [int(i) for i in list(all_df.sum(axis=1))]
            reads_per_cell, barcodes = (list(t) for t in zip(*sorted(zip(reads_per_cell, barcodes), reverse=True)))
            
            valid_cells = dabm.call_cells([i for i in zip(reads_per_cell, barcodes)], plot=True, low_count_filter=100, batch_label=sample_basename)
            

            # check that valid cell barcodes were found
            if len(valid_cells) == 0:
                print('No valid cells found! Please check input files.')
                raise SystemExit

            # create SingleCell objects for each valid cell
            cells = [resources.SingleCell(barcode,
                                          by_cell_fastq_dir,
                                          by_cell_bam_dir,
                                          by_cell_gvcf_dir,
                                          by_cell_flt3_dir)
                     for barcode in valid_cells]

            print('%s valid cells found!\n' % len(cells))
            
            # check that all barcodes have the same length for demultiplexing
            bar_length = list(set([len(k) for k in valid_cells]))
            assert len(bar_length) == 1, 'All barcodes must have the same length!'

            # split cell fastq files for panel reads associated with valid cells
            barcode_names = temp_dir + 'barcodes.txt'
            with open(barcode_names, 'w') as f:
                for b in valid_cells:
                    f.write(b + '\n')
            



        if int(args.ex[3]):
        
            print('''
            ###################################################################################
            # write valid cells from panel reads to separate fastq files
            ###################################################################################
            ''')
            
            if human_concordant:
                R1_file = os.path.join(args.output, '{}_Hu_R1.fastq.gz'.format(sample_basename)) 
                R2_file = os.path.join(args.output, '{}_Hu_R2.fastq.gz'.format(sample_basename))
            else:
                R1_file = './fastq_out/{}_R1.fastq.gz'.format(sample_basename)
                R2_file = './fastq_out/{}_R2.fastq.gz'.format(sample_basename)
            
            # split files by cell barcode using bbmap demuxbyname.sh
            demux_cmd = 'demuxbyname.sh prefixmode=f -Xmx10g prefixmode=f delimiter=_ in={} in2={} out={} names={}'.format( 
                R1_file, R2_file,
                by_cell_fastq_dir + '%.fastq.gz',
                barcode_names) # Xmx10g allocates 10 gigabyte memory
            p = subprocess.Popen(demux_cmd, shell=True)
            wait([p])
        
        if int(args.ex[4]):
        
            print('''
            ####################################################################################
            # Step 7: align, convert, sort, index panel reads (optional: call FLT3-ITDs)
            ####################################################################################
            ''')

            # limit number of cells to preprocess at a time (based on hardware limitations)
            n_preprocess = 24

            # create pool of workers and run through all samples
            preprocess_pool = ThreadPool(processes=n_preprocess)

            # align and index cells
            for c in cells:
                preprocess_pool.apply_async(resources.SingleCell.align_and_index, args=(c, bt2_ref_2,))

            preprocess_pool.close()
            preprocess_pool.join()
            
            # optionally, call FLT3-ITDs using ITDseek
            if not args.skip_flt3:

                preprocess_pool = ThreadPool(processes=n_preprocess)

                for c in cells:
                    preprocess_pool.apply_async(resources.SingleCell.call_flt3, args=(c, human_fasta_file_2,))

                preprocess_pool.close()
                preprocess_pool.join()

            else:
                os.rmdir(by_cell_flt3_dir)
        
        if int(args.ex[5]):
        
            print('''
            ####################################################################################
            # Step 8: perform variant calling for all cells
            ####################################################################################''')

            # limit number of cells to call variants at a time (based on hardware limitations)
            n_call_variants = 24

            # create pool of workers and run through all samples
            call_variants_pool = ThreadPool(processes=n_call_variants)

            for c in cells:
                call_variants_pool.apply_async(resources.SingleCell.call_variants, args=(c,
                                                                                         human_fasta_file_2,
                                                                                         interval_file_2,))

            call_variants_pool.close()
            call_variants_pool.join()

            ################################################################################
        
    print('Pipeline complete!')
        
