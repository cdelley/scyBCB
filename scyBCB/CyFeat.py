"""
Created on Tue Nov 28 18:12:22 2017

@author: Cyrille L. Delley
"""
#from typing import Iterator
from collections import Counter
from xopen import xopen
import numpy as np
import time
import scipy.sparse as sps
import anndata as ad
import pandas as pd

import scyBCB.pyscBUS.BUSlib as bus
from scyBCB.CyBCB import seq_exp
from scyBCB.utilities.network import UMIClusterer

def decode_umi(bin_code):
    return str(int(bin(bin_code)[2::2]) + int(bin(bin_code)[3::2])*2).replace('0','A').replace('1','G').replace('2','C').replace('3','T')[1:]

class features(seq_exp):
    """Class for processing and assigning features in single-cell sequencing data, particularly for antibody barcodes."""
    
    def __init__(self, files:list, bc_correction_dic:list, seq_positions:dict):
        """         
        Initializes the features class with given files, barcode correction dictionaries, and sequence positions.
        Inherits from seq_exp class and uses its initialization method. 
        """
        super().__init__(files, bc_correction_dic, seq_positions)
        self._is_sorted = False
    
    @classmethod
    def assign_features(
            cls,
            file1:str,
            file2:str, 
            bc_correction_dic:list,
            feature_seq_dic:dict,
            seq_positions:dict,
            verbose=True, 
            down_sample=-1,
        ):
                        
        """
        Class method to assign features from sequencing data.

        Processes sequencing data from two files to extract and assign features based on provided dictionaries and sequence positions.
        It tracks various statistics and unknown features during processing.

        Args:
            file1 (str): Path to the first sequencing data file.
            file2 (str): Path to the second sequencing data file.
            bc_correction_dic (list): List of barcode correction dictionaries.
            feature_seq_dic (dict): Dictionary mapping feature sequences.
            seq_positions (dict): Dictionary specifying positions of sequences.
            verbose (bool): Flag to control printing of process statistics.
            down_sample (int): Number of reads to process, -1 for all.

        Returns:
            features: An instance of the features class with processed data.
        """
        
        inst = cls((file1, file2), bc_correction_dic, seq_positions) # create the instance
        
        unknown_features = []
        counter = 0
        inst.umis = []
        inst.bc_feat = []
        with xopen(file1, 'rt') as f1:
            with xopen(file2, 'rt') as f2:
                for read in inst.extract_barcodes(iter(f1), iter(f2)):
                    if counter == down_sample:
                        break
                    counter += 1
                    try:
                        p_name  = feature_seq_dic[read['feature'][0]]
                        inst.output_stats[3] += 1
                    except KeyError:
                        unknown_features.append(read['feature'][0])
                        inst.output_stats[4] += 1
                        continue
                  
                    inst.bc_feat.append(read['BC']*10000 + p_name)
                    inst.umis.append(read['UMI'])
        
        inst.bc_feat = np.array(inst.bc_feat)                             
        inst.unknown_features = Counter(unknown_features)
        
        if verbose:
            print('total reads processed: {}'.format(inst.output_stats[0]))
            print('{} good barcodes discovered'.format(inst.output_stats[1]))
            print('{} reads failed  barcode error threshold '.format(inst.output_stats[2]))
            print('{} reads failed  UMI with N '.format(inst.output_stats[5]))
            print('{} good feature sequences identified'.format(inst.output_stats[3]))
            print('{} feature sequences failed error threshold '.format(inst.output_stats[4]))
            
        return inst
    
    def sort_bc(self) -> None:
        """Sorts the barcode group list and the UMIs for easy access."""
        self.srt_idx = self.bc_feat.argsort()
        self.umis = np.array(self.umis)[self.srt_idx].astype(int)
        self.bc_feat = np.array(self.bc_feat)[self.srt_idx]
        self._is_sorted = True
        
    def prepare_BUS(self) -> None:
        """Prepare the output for a BUS file format"""
        
        if not self._is_sorted:
            self.sort_bc()
            
        umis = [self.umis[0]]
        bc_feat0 = self.bc_feat[0]

        header = '#01 '+ self.files[0]+' \n'
        header += '#02 '+time.strftime("%Y-%m-%d %H:%M")+'\n'
        header += '#03 {:> 11} TOTAL_READS \n#04 {:> 11} GOOD_BARCODES \n#05 {:> 11} GOOD_FEATURES \n'.format(
            self.output_stats[0], self.output_stats[1], self.output_stats[3])
        header += '#06 {:> 11} FAILED_BARCODE \n#07 {:> 11} FAILED_UMI_N \n#08 {:> 11} FAILED_FEATURE \n#09\n'.format(
            self.output_stats[2], self.output_stats[5], self.output_stats[4])

        self.bus = bus.BUSFIle(
            text = header,
            bcd = dict(),
            bclen = self.bc_len,
            umilen = self.umi_len,
        )
        for i in range(self.output_stats[3]-1):
            if self.bc_feat[i+1] != bc_feat0:
                # compress UMIs
                umi_dic = np.array(list(Counter(umis).items())) 
                # prepare output
                bc = bc_feat0 // 10**4
                feat = bc_feat0 % 10**4
                flag = 0
                for umi, num in umi_dic[umi_dic[:, 1].argsort()][::-1]:
                    try:
                        self.bus.bcd[bc][feat].append((umi, num, flag))
                    except KeyError:
                        try:
                            self.bus.bcd[bc][feat] = [(umi, num, flag)]
                        except KeyError:
                            self.bus.bcd[bc] = {feat : [(umi, num, flag)]}
                    self.bus.n_entries += 1
                    
                # prepare next cycle
                bc_feat0 = self.bc_feat[i+1]
                umis = [self.umis[i+1]]
                
            else:
                umis.append(self.umis[i+1])
    
    def write_condensed_file(self, filename) -> None:
        """ write a condensed feature file"""
        if not self._is_sorted:
            self.sort_bc()
            
        with open(filename, 'w') as fout:
            fout.write('#01 '+ self.files[0]+'\n')
            fout.write('#02 '+time.strftime("%Y-%m-%d %H:%M")+'\n')
            fout.write('#03 {:> 11} TOTAL_READS\n#04 {:> 11} GOOD_BARCODES\n#05 {:> 11} GOOD_FEATURES\n'.format(
                self.output_stats[0], self.output_stats[1], self.output_stats[3]))
            fout.write('#06 {:> 11} FAILED_BARCODE\n#07 {:> 11} FAILED_UMI_N\n#08 {:> 11} FAILED_FEATURE\n#09\n'.format(
                self.output_stats[2], self.output_stats[5], self.output_stats[4]))
            
            umis = [self.umis[0]]
            bc_feat0 = self.bc_feat[0]
            for i in range(self.output_stats[3]-1):
                if self.bc_feat[i+1] != bc_feat0:
                    umi_dic = np.array(list(Counter(umis).items()))
                    bc = bc_feat0 // 10**4
                    feat = bc_feat0 % 10**4
                    # write to file
                    for umi, c in umi_dic[umi_dic[:, 1].argsort()][::-1]:
                        fout.write('{:>8} {:>4} {} {}\n'.format(bc, feat, umi, c))
                    bc_feat0 = self.bc_feat[i+1]
                    umis = [self.umis[i+1]]
                else:
                    umis.append(self.umis[i+1])    

def cluster_umis(umi_dict, method='directional'):
    # clusters the umis using the specified method (or all)
    # uses functions from umi-tools paper (Genome Research, 2017)

    # split umi dict into umis (keys) and counts
    umis = list(umi_dict.keys())
    counts = umi_dict

    # set up UMIClusterer functor with parameters specific to specified method
    # choose method = 'all' for all available methods
    # otherwise provide methods as a list of methods
    processor = UMIClusterer()  # initialize UMIclusterer

    # cluster the umis
    clusters = processor(
        umis,
        counts,
        threshold=1,
        cluster_method=method)
    return clusters

def count_table_feature(filename, method='directional'):
    
    barcode = []
    feature = []
    count = []
    count_raw = []

    hamming_merged = 0
    with open(filename, 'r') as fin:
        for i in range(9):
            fin.readline()
    
        bc0, feat0, umi, c = fin.readline().strip().split()
        umis = [decode_umi(int(umi))]
        counts = [int(c)]

        for i, line in enumerate(fin):
            bc, feat, umi, c = line.strip().split()
            if (bc, feat) != (bc0, feat0): 
                # row, column index of coo sparse matrix
                barcode.append(int(bc0))
                feature.append(int(feat0))
                # values of coo sparse matrix
                umi_dic = {u:c for (u,c) in zip(umis, counts)}
                count_raw.append(sum(counts))
                clusters = cluster_umis(umi_dic, method)
                corr_counts = len(clusters[method])
                count.append(corr_counts)
                hamming_merged += len(umis) - corr_counts
                
                #reset
                bc0 = bc
                feat0 = feat
                umis = [decode_umi(int(umi))]
                counts = [int(c)]                  
            else:
                umis.append(decode_umi(int(umi)))
                counts.append(int(c))
        # row, column index of coo sparse matrix
        barcode.append(int(bc0))
        feature.append(int(feat0))
        # values of coo sparse matrix
        umi_dic = {u:c for (u,c) in zip(umis, counts)}
        count_raw.append(sum(counts))
        clusters = cluster_umis(umi_dic, method)
        corr_counts = len(clusters[method])
        count.append(corr_counts)
        hamming_merged += len(umis) - corr_counts
            
    print('Hamming merged a total of {} umis out of {} total'.format(hamming_merged, i))
    
    # make the raw array
    matr = sps.coo_matrix(
        (np.array(count_raw), (np.array(barcode), np.array(feature)))
    )
    matr = sps.csr_matrix(matr)
    non_zero = matr.getnnz(axis=1) != 0
    adata = ad.AnnData(matr[non_zero], 
    obs=pd.DataFrame(index=np.arange(matr.shape[0])[non_zero]), 
    var=pd.DataFrame(index=np.arange(matr.shape[1])))
    adata.obs['raw_counts'] = np.array(adata.X.sum(axis=1)).squeeze()
    adata.layers['raw_counts'] = adata.X
    
    # make the corrected array
    matr = sps.coo_matrix(
        (np.array(count), (np.array(barcode), np.array(feature)))
    )
    matr = sps.csr_matrix(matr)
    adata.X = matr[non_zero]
    adata.obs['umi_counts'] = np.array(adata.X.sum(axis=1)).squeeze()
    print(adata)        
    return adata
