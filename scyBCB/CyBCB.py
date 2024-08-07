"""
Created on Tue Nov 28 18:12:22 2017

@author: Cyrille L Delley
"""
from typing import Iterator
import numpy as np

class seq_exp(object):
    """
    A class for evaluating single-cell sequencing experiments.
    
    This class provides methods for processing and analyzing sequencing data, 
    particularly focusing on barcode extraction and feature assignment in 
    single-cell sequencing experiments.

    Attributes:
        files (list): List of file paths for the sequencing data.
        bc_correction_dic (list): A list of dictionaries for barcode correction.
        seq_positions (dict): Dictionary specifying the positions of barcodes, UMIs, and features in the sequence.
        bc_groups (dict): Dictionary to store barcode groups (initialized empty).
        output_stats (numpy.ndarray): Array to store statistics of the processing steps.
    """
    
    def __init__(self, files:list, bc_correction_dic:list, seq_positions:dict):
        """
        Initializes the seq_exp object with provided files, barcode correction dictionaries, and sequence positions.

        Args:
            files (list): List of file paths for the sequencing data.
            bc_correction_dic (list): A list of dictionaries for barcode correction.
            seq_positions (dict): Dictionary specifying the positions of barcodes, UMIs, and features in the sequence.
            
        Example of `self.seq_positions` for reference:
            self.seq_positions = {
                'bc': [(0, 11, 21), (1, 2, 12)],
                'umi': [(0, 23, 30)],
                'feature': [(0, 30, 50)]
            }
        """
        self.files = files
        self.bc_correction_dic = bc_correction_dic
        self.seq_positions = seq_positions
        self.bc_groups = {}
        
        # add non specified keywords
        try:
            self.seq_positions['umi']
        except KeyError:
            self.seq_positions['umi'] = []
        try:
            self.seq_positions['seq']
        except KeyError:
            self.seq_positions['seq'] = []
        try:
            self.seq_positions['feature']
        except KeyError:
            self.seq_positions['feature'] = []
            
        self.umi_len = np.sum([i[-1] - i[-2] for i in seq_positions['umi']])
        self.bc_len = np.sum([i[-1] - i[-2] for i in seq_positions['bc']])
        self.feature_len = np.sum([i[-1] - i[-2] for i in seq_positions['feature']])

        # processed reads, correct bc, failed bc, good feature, corrected feature, ambigous feature, error to high
        self.output_stats = np.zeros(7).astype('int32') 

    
    def extract_barcodes(self, handle1:object, handle2:object) -> Iterator[dict]:
        """
        Extracts barcodes, UMIs, and feature sequences from two FASTQ file handles.

        This method iterates over two open FASTQ file handles, processing four lines at a time. It extracts 
        barcodes (BC), UMIs, and feature sequences based on positions defined in `self.seq_positions`.

        Each element in 'self.seq_positions' is a tuple representing a sequence section. The first element
        of the tuple indicates whether it's from read 1 or read 2 file (0 or 1), the second is the start basepair
        position, and the third is the end basepair position. Multiple tuples can be used for each sequence type.

        For barcode ('bc') sequences, each sequence is matched against a barcode whitelist and encoded with a numeric
        ID. For UMI sequences, the fragments are concatenated and translated into a binary sequence for storage; reads
        with ambiguous bases are discarded. Feature sequences are extracted and returned as is.

        Args:
            handle1 (object): An open file handle to the first FASTQ file.
            handle2 (object): An open file handle to the second FASTQ file.

        Yields:
            dict: A dictionary containing the processed read information with keys 'BC', 'UMI', and 'feature'.

        Example of returned 'read' dictionary:
            read = {'BC': 9696, 'UMI': 1236, 'feature': 'GATGGAT...'}
        """
        while True:
            try:
                ID1 = next(handle1).rstrip()[1:]
                seq1 = next(handle1).rstrip()
                next(handle1)
                qual1 = next(handle1).rstrip()
                
                ID2 = next(handle2).rstrip()[1:]
                seq2 = next(handle2).rstrip()
                next(handle2)
                qual2 = next(handle2).rstrip()
                if not ID1:
                    break
                
                # assert read ids match
                ID = ID1.split()[0]     
                assert ID == ID2.split()[0]
                self.output_stats[0] += 1
                
                _seq = {0 : seq1, 1 : seq2}
                _qual = {0 : qual1, 1 : qual2}
                _off = {0 : 0, 1 : 0}
                    
                # identify barcode from whitelist
                read = {'BC':0, 'UMI':'', 'feature':[], 'ID':ID, 'seq':[], 'qual':[]}
                try:
                    for se, bc_di in zip(self.seq_positions['bc'], self.bc_correction_dic):
                        _bc_of = bc_di[_seq[se[0]][se[1] + _off[se[0]] : se[2] + _off[se[0]]]]
                        
                        # for each barcode group we right shift the previous barcodes
                        read['BC'] = read['BC']*100 + _bc_of[0]
                        _off[se[0]] += _bc_of[1]
                        
                    # to ensure all barcodes have the same number of digit we add a leading 1
                    read['BC'] += 100**len(self.seq_positions['bc'])
                    
                except KeyError:
                    self.output_stats[2] += 1
                    continue
                self.output_stats[1] += 1

                # get umi sequence if present                
                for i in self.seq_positions['umi']:
                    read['UMI'] += _seq[i[0]][i[1] + _off[i[0]] : i[2] + _off[i[0]]]
                    
                # translate umi to binary encoding. This will fail if 'N's are present
                try:
                    read['UMI'] = int(
                        '11' + read['UMI'].replace('A','00').replace('C','01').replace('G','10').replace('T','11'),
                        2
                    )
                except ValueError:
                    self.output_stats[5] += 1
                    continue
                
                # get feature sequences if present
                for i in self.seq_positions['feature']:
                    read['feature'].append(_seq[i[0]][i[1] + _off[i[0]] : i[2] + _off[i[0]]])
                    
                # get sequence fragments of interest
                for i in self.seq_positions['seq']:
                    read['seq'].append(_seq[i[0]][i[1] + _off[i[0]] : i[2] + _off[i[0]]])
                    read['qual'].append(_qual[i[0]][i[1] + _off[i[0]] : i[2] + _off[i[0]]])
                
                yield read
            except StopIteration:
                return

