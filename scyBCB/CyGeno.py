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

from scyBCB.CyBCB import seq_exp
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
    def get_cells(
        cls, 
        file1:str, 
        file2:str, 
        bc_correction_dic:list, 
        seq_positions:dict, 
        feature_dic:list=[], 
        verbose:bool=False,
        down_sample=-1,
    ):

        inst = cls((file1, file2), bc_correction_dic, seq_positions) # create the instance
        counter = 0
        
        with gzip.open(file1, 'rt') as f1:
            with gzip.open(file2, 'rt') as f2:
                for read in inst.extract_barcodes(f1, f2):
                    if counter == down_sample:
                        break
                    counter += 1
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
    def write_good_reads(
        cls, 
        file1:str, 
        file2:str, 
        bc_correction_dic:list, 
        seq_positions:dict,
        primer_dic:list = None, 
        batch_label:str = '', 
        out_path:str = './', 
        verbose:bool=False,
        down_sample=-1,
    ):

        inst = cls((file1, file2), bc_correction_dic, seq_positions) # create the instance
        
        if primer_dic is not None:
            feat_names1 = np.sort(list(set([name for name in primer_dic[0].values()]))+['0_unk'])
            feat_idx1 = {fn:i for (i, fn) in enumerate(feat_names1)}
            feat_names2 = np.sort(list(set([name for name in primer_dic[1].values()]))+['1_unk'])
            feat_idx2 = {fn:i for (i, fn) in enumerate(feat_names2)}
        
            inst.amplicon_count = np.zeros([len(feat_names1), len(feat_names2)])
        
        fout1, fout2 = (xopen(os.path.join(out_path, '{}R1.fastq.gz'.format(batch_label)), 'w'), xopen(os.path.join(out_path, '{}R2.fastq.gz'.format(batch_label)), 'w'))
        
        counter = 0
        with xopen(file1, 'rt') as f1:
            with xopen(file2, 'rt') as f2:
                for read in inst.extract_barcodes(iter(f1), iter(f2)):
                    
                    if counter == down_sample:
                        break
                    counter += 1
                    
                    if primer_dic is not None:
                        #amplicon count
                        try:
                            _p1 = primer_dic[0][read['feature'][0]]
                        except KeyError:
                            _p1 = '0_unk'
                        try:
                            _p2 = primer_dic[1][read['feature'][1]]
                        except KeyError:
                            _p2 = '1_unk'
                        inst.amplicon_count[feat_idx1[_p1], feat_idx2[_p2]] += 1
                        
                    try:
                        inst.bc_groups[read['BC']] += 1
                    except KeyError:
                        inst.bc_groups[read['BC']] = 1
                    
                    #fastq out
                    r1 = '@'+read['ID']+'-R1-'+read['feature'][0][:-4]+'_'+str(read['BC'])+batch_label+'\n'+read['seq'][0]+'\n'+'+\n'+read['qual'][0]+'\n'
                    r2 = '@'+read['ID']+'-R2-'+read['feature'][1][:-4]+'_'+str(read['BC'])+batch_label+'\n'+read['seq'][1]+'\n'+'+\n'+read['qual'][1]+'\n'
                        
                    fout1.write(r1)
                    fout2.write(r2)
        fout1.close()
        fout2.close()
        
        if primer_dic is not None:           
            inst.adata = ad.AnnData(inst.amplicon_count, obs=pd.DataFrame(index=feat_names1), var=pd.DataFrame(index=feat_names2))
        
        if verbose:
            print('total reads processed: {}'.format(inst.output_stats[0]))
            print('{} good barcodes discovered'.format(inst.output_stats[1]))
            print('{} reads failed  barcode error threshold '.format(inst.output_stats[2]))
            print('{}'.format(inst.output_stats))
        return inst 
               
def call_cells(barcode_count_dict, plot=True, low_count_filter=100, batch_label=''):
    """call cells using the knee method, returns list of valid cell barcodes."""
    # count reads per barcode
    reads = np.array(list(barcode_count_dict.values()))
    _thresh_idx = reads >= low_count_filter

    reads = reads[_thresh_idx]
    bcs = np.array(list(barcode_count_dict.keys()))[_thresh_idx]

    _srt_idx = np.argsort(reads)[::-1]
    reads = reads[_srt_idx]
    bcs = bcs[_srt_idx]  
    
    # second derivative method (inflection point) of the knee plot to identify cells
    # 1 - first derivative of cell rank plot
    x = np.log10(range(1, len(reads)+1))
    y = np.log10(reads)

    dy = np.zeros(y.shape, float)
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

    barcode_count_touples = list(barcode_count_dict.items())

    if plot:
        cells=[]

        for i, j in barcode_count_touples:
            cells.append(i)
        cells = np.array(cells)
        
        fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6,8), sharex=True)
        l1 = ax1.plot(np.sort(list(barcode_count_dict.values()))[::-1])
        ax1.plot([n_cells, n_cells],[300, np.max(reads)], 'k:')
        ax1.plot([max_peak_old, max_peak_old],[300, np.max(reads)], 'r:')
        ax1.set_title('read per barcode distribution, {} cells'.format(n_cells))
        ax1.set_xlabel('barcode goup rank')
        ax1.set_ylabel('log read per barcode count')
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        
        ax2.plot(range(len(yhat)),yhat)
        ax2.set_title('Savitzky-Golay filter smoothed read counts')
        ax2.set_xlabel('barcode goup rank')
        ax2.set_ylabel('dy/dx')
        ax2.set_xscale('symlog')
        os.makedirs('./qc_output', exist_ok=True) # ensure qc folder exists
        plt.savefig('./qc_output/Cell_calling_{}.png'.format(batch_label), dpi=600)
        plt.show()
    return cell_barcodes
