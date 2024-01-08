"""
Created on Tue Nov 28 18:12:22 2017

@author: cyrille
"""
import random
import numpy as np
import os
import gzip
import json

from collections import Counter
import itertools
import anndata as ad
import pandas as pd
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

from CyBCB import seq_exp
from xopen import xopen


class sequences(seq_exp):
    """Class contains methods to evaluate sequencing experiment"""
    
    def __init__(self, files:list, bc_correction_dic:list, seq_positions:dict):
        """         
        Initializes the features class with given files, barcode correction dictionaries, and sequence positions.
        Inherits from seq_exp class and uses its initialization method. 
        """
        super().__init__(files, bc_correction_dic, seq_positions)
    
    @classmethod
    def get_cells(cls, file1:str, file2:str, bc_correction_dic:list, seq_positions:dict, feature_dic:list=[], verbose:bool=False):

        inst = cls((file1, file2), bc_correction_dic, seq_positions, feature_dic) # create the instance
        
        with gzip.open(file1, 'rt') as f1:
            with gzip.open(file2, 'rt') as f2:
                for read in inst.extract_barcodes(f1, f2):
                
                    try:
                        inst.bc_groups[read['BC']][read['ID']] = read
                    except KeyError:
                        inst.bc_groups[read['BC']] = {read['ID'] : read}

        if verbose:
            print('total reads processed: {}'.format(inst.output_stats[0]))
            print('{} good barcodes discovered'.format(inst.output_stats[1]))
            print('{} reads failed  barcode error threshold '.format(inst.output_stats[2]))
        return inst
        
        
    @classmethod
    def write_demux(
        cls, 
        file1:str, 
        file2:str, 
        bc_correction_dic:list, 
        seq_positions:dict,
        feature_dic:list, 
        batch_label:str='', 
        out_path:str='./', 
        verbose:bool=False
    ):

        inst = cls((file1, file2), bc_correction_dic, seq_positions, feature_dic) # create the instance
        
        feat_names1 = np.sort(list(set([name[1] for name in feature_dic[0].values()]))+['0_unk'])
        feat_idx1 = {fn:i for (i, fn) in enumerate(feat_names1)}
        feat_names2 = np.sort(list(set([name[1] for name in feature_dic[1].values()]))+['1_unk'])
        feat_idx2 = {fn:i for (i, fn) in enumerate(feat_names2)}
        
        # separate the reads human from HIV and nonesense 
        HIV_handles = (xopen(os.path.join(out_path, '{}_HIV_R1.fastq.gz'.format(batch_label)), 'w'), xopen(os.path.join(out_path, '{}_HIV_R2.fastq.gz'.format(batch_label)), 'w'))
        human_handles = (xopen(os.path.join(out_path, '{}_Hu_R1.fastq.gz'.format(batch_label)), 'w'), xopen(os.path.join(out_path, '{}_Hu_R2.fastq.gz'.format(batch_label)), 'w'))
        crap_handles = (xopen(os.path.join(out_path, '{}_unk_R1.fastq.gz'.format(batch_label)), 'w'), xopen(os.path.join(out_path, '{}_unk_R2.fastq.gz'.format(batch_label)), 'w'))
        
        inst.amplicon_count = np.zeros([len(feat_names1), len(feat_names2)])
        
        with xopen(file1, 'rt') as f1:
            with xopen(file2, 'rt') as f2:
                for read in inst.extract_barcodes(iter(f1), iter(f2)):
                    #amplicon count
                    inst.amplicon_count[feat_idx1[read['feature'][0]], feat_idx2[read['feature'][1]]] += 1
                    
                    #fastq out
                    r1 = '@'+read['ID']+'-R1-'+read['feature'][0][:-4]+'_'+read['BC']+batch_label+'\n'+read['seq'][0]+'\n'+'+\n'+read['qual'][0]+'\n'
                    r2 = '@'+read['ID']+'-R2-'+read['feature'][1][:-4]+'_'+read['BC']+batch_label+'\n'+read['seq'][1]+'\n'+'+\n'+read['qual'][1]+'\n'
                    
                    # sort by being not bonafide human
                    if read['feature'][0][:-4] == read['feature'][1][:-4]:
                        fout1, fout2 = human_handles
                    elif read['feature'][0][:3] == 'HIV' or read['feature'][1][:3] == 'HIV':
                        fout1, fout2 = HIV_handles
                    else:
                        fout1, fout2 = crap_handles
                        
                    fout1.write(r1)
                    fout2.write(r2)
                    
        # close the file handles
        for f in HIV_handles:
            f.close()
        for f in human_handles:
            f.close()
        for f in crap_handles:
            f.close()
            
        inst.adata = ad.AnnData(inst.amplicon_count, obs=pd.DataFrame(index=feat_names1), var=pd.DataFrame(index=feat_names2))
        
        if verbose:
            print('total reads processed: {}'.format(inst.output_stats[0]))
            print('{} good barcodes discovered'.format(inst.output_stats[1]))
            print('{} reads failed  barcode error threshold '.format(inst.output_stats[2]))
            print('{}'.format(inst.output_stats))
        return inst 
               
def call_cells(barcode_count_touples, plot=True, low_count_filter=100, batch_label=''):
    """call cells using the knee method, returns list of valid cell barcodes."""
    # count reads per barcode
    bcs = []
    reads = []
    for i in barcode_count_touples:
        bcs.append(i[1])
        reads.append(i[0])
       
    reads, bcs = (list(t) for t in zip(*sorted(zip(reads, bcs), reverse=True)))

    # second derivative method (inflection point) of the knee plot to identify cells
    # 1 - first derivative of cell rank plot
    # exclude barcodes with low numbers of reads
    rpc_thresh = [x for x in reads if x >= low_count_filter]
    
    x = np.log10(range(1, len(rpc_thresh)+1))
    y = np.log10(np.array(rpc_thresh))

    dy = np.zeros(y.shape, np.float)
    dy[0:-1] = np.diff(y) / np.diff(x)
    dy[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    dy = -dy  # invert for positive graph

    # smooth the data by savgol filtering and call first peak
    try:
        yhat = savgol_filter(dy, int(len(dy)/5), 3)  # window size, polynomial order
    except ValueError:
        yhat = savgol_filter(dy, int(len(dy)/5)+1, 3)  # window size, polynomial order

    # prominence of peak (0.1 should be adequate for most mammalian cell panels)
    prominence = 0.25
    height = 0.5

    peaks = find_peaks(yhat, height=height, prominence=prominence)
    max_peak_i = np.argmax(peaks[1]['prominences'])
    max_peak_old = peaks[0][max_peak_i]
    max_peak = peaks[0][0]

    # first n cell barcodes are valid
    n_cells = max_peak    #n_cells =  int(bin_centers[max_peak])
    cell_barcodes = np.sort(bcs[:n_cells])

    if plot:
        cells=[]

        for i, j in barcode_count_touples:
            cells.append(i)
        cells = np.array(cells)
        
        fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8,8), sharex=True)
        l1 = ax1.plot(np.sort(cells[cells>5])[::-1])
        ax1.plot([n_cells, n_cells],[300, np.max(cells)], 'k:')
        ax1.plot([max_peak_old, max_peak_old],[300, np.max(cells)], 'r:')
        ax1.set_yscale('symlog')
        ax1.set_xscale('symlog')
        ax1.set_title('read per barcode distribution, {} cells'.format(n_cells))
        ax1.set_xlabel('barcode goup rank')
        ax1.set_ylabel('log read per barcode count')

        ax2.plot(range(len(yhat)),yhat)
        ax2.set_title('Savitzky-Golay filter smoothed read counts')
        ax2.set_xlabel('barcode goup rank')
        ax2.set_ylabel('dy/dx')
        ax2.set_xscale('symlog')
        CHECK_FOLDER = os.path.isdir('./qc_output')
        if not CHECK_FOLDER:
            os.makedirs('./qc_output')
        plt.savefig('./qc_output/Cell_calling_{}.png'.format(batch_label), dpi=600)
        plt.show()
    return cell_barcodes
