#
#   Author: Cyrille L. Delley
#   
#   Functions to read and write uncompressed BUS files with python.
#    
#   The file format was introduced by Páll Melsted, Sina Booeshaghi, Lior Pachtr
#   and coworkers for kallisto.
#
#   Melsted, Páll, Booeshaghi, A. Sina et al. Modular and efficient pre-processing
#   of single-cell RNA-seq. BioRxiv (2019): 673285, doi.org/10.1101/673285.   
#   See file specifications at:
#
#   https://github.com/BUStools/BUS-format
#   https://github.com/BUStools/bustools
#

import csv
import math
import struct
from typing import List, Dict, Tuple, Union, BinaryIO
import warnings
from xopen import xopen

BUSFORMAT_VERSION = 1

def stringToBinary(s: str, length: int = 32, formatting: str = '064b')  -> int:
    """
    Converts DNA string to binary number representation with:
    A   -> 00
    C   -> 01
    G   -> 10
    T   -> 11
    
    Input
        s: DNA string
        length: length of s, if not provided uses len(s) (default None)
        formatting: binary format of the output, (default '032b')
    
    Returns
        int (formatted accoding to 'formatting')
    """
    bs = s[:length].replace(
        'A','00').replace(
        'C','01').replace(
        'G','10').replace(
        'T','11')
    r = int(bs, 2)
    #return format(r, formatting)
    return int(r)
    
def binaryToString(x: int, length: int)  -> str:
    """
    Reverts an binary integer number into a DNA string with:
    A   <- 00
    C   <- 01
    G   <- 10
    T   <- 11
    
    string length needs to be specified to get the appropriate number
    of leading 'A's
    """
    s = ["N"] * length
    sh = length - 1
    for i in range(length):
        c = "N"
        if ((x >> (2 * sh)) & 0x03) == 0x00:
            c = "A"
        elif ((x >> (2 * sh)) & 0x03) == 0x01:
            c = "C"
        elif ((x >> (2 * sh)) & 0x03) == 0x02:
            c = "G"
        elif ((x >> (2 * sh)) & 0x03) == 0x03:
            c = "T"
        sh -= 1
        s[i] = c
    return "".join(s)

class BUSFIle(object):

    def __init__(
        self,
        text = "",
        bcd = dict(),
        bclen = 0,
        umilen = 0,
        textlen = 0,
        n_entries = 0
    ):
        self.text = text
        self.bcd = bcd     # type: dict[BC : [(UMI1, ecs1, count, flag), 
                           #                  (UMI2, ecs1, count, flag), 
                           #                  (UMI1, ecs2, count, flag), ....]
        self.version = BUSFORMAT_VERSION
        self.bclen = bclen
        self.umilen = umilen
        self.textlen = textlen
        self.n_entries = n_entries
        self.format = ''
    
    @classmethod
    def read(cls, in_file: str) -> "BUSFIle":
        """
        read BUS uncompressed files and return a BUSFile instance
        """
        instance = cls()
        instance.in_file_path = in_file
        with xopen(in_file, 'rb', threads=2, compresslevel=3) as inf:
            instance.parseHeader(inf)
            if instance.format == 'BUS':
                instance.parseBUS(inf)
            elif instance.format == 'CBIN':
                instance.parseCBIN(inf)
            else:
                print('file format {} not implemented'.format(instance.format))
        return instance

    def parseBUS(self, inf: BinaryIO) -> None:
        """
        parse the data blocks of an uncompressed BUS file
        """
        b = inf.read(32)
        self.n_entries += 1
        while b:
            bc, umi, ec, num, flag, pad = struct.unpack('<QQIIII', b)
            b = inf.read(32)
            self.n_entries += 1
            try:
                self.bcd[bc][ec].append((umi, num, flag))
            except KeyError:
                try:
                    self.bcd[bc][ec] = [(umi, num, flag)]
                except KeyError:
                    self.bcd[bc] = {ec : [(umi, num, flag)]}
    
    def BUS_to_txt(self) -> None:
        # binaryToString(bc, self.bclen), binaryToString(umi, self.umilen), ec, num
        pass
        
    def parseHeader(self, inf: BinaryIO) -> None:
        """
        parse the header of an uncompressed BUS file
        """
        magic = struct.unpack("<4s", inf.read(4))[0]
        if magic == b"BUS\x00":
            self.format = 'BUS'
        elif magic == b"CBIN":
            self.format = magic.decode('utf-8')
        else:
            self.format = magic    
            return       
        self.version, self.bclen, self.umilen, self.textlen = struct.unpack('<IIII', inf.read(16))
        if self.version != BUSFORMAT_VERSION:
            warnings.warn("The BUS file version is different from the PyBUS decoder version")
        self.text = struct.unpack("<{}s".format(self.textlen), inf.read(self.textlen))[0]
    
    def write(self, outf: str, formatting: str = 'BUS') -> None:
        """
        write an uncompressed BUS file
        """
        self.out_file_path = outf
        with xopen(outf, 'wb', threads=2, compresslevel=3) as fout:
            self.writeHeader(fout, formatting)
            if formatting == 'BUS':
                self.writeBUS(fout)
            elif formatting == 'CBIN':
                self.writeCBIN(fout)
    
    def writeHeader(self, outf: BinaryIO, formatting: str) -> None:
        """
        write the header of an uncompressed BUS file
        """
        try:
            text_bytes = self.text.encode()
        except AttributeError:
            text_bytes = self.text
        outf.write(struct.pack("<4s", formatting.encode()))
        outf.write(struct.pack('<IIII', self.version, self.bclen, self.umilen, len(text_bytes)))
        outf.write(text_bytes)
        
    def writeBUS(self, outf: BinaryIO, sort_key: str = None) -> None:
        """
        write the data blocks of an uncompressed BUS file
        """
        if not sort_key:
            for bc, ec_dict in self.bcd.items():
                for ec, value in ec_dict.items():
                    for tup in value:
                        block = struct.pack('<QQIIII', bc, tup[0], ec, tup[1], tup[2], 0)
                        outf.write(block)

    def writeCBIN(self, outf: BinaryIO, sort_key: str = None) -> None:
        """
        write the data blocks formated in cy binary
        """
        if not sort_key:
            bc0 = -1
            for bc, ec_dict in self.bcd.items():
                block1 = struct.pack('<QI', bc, len(ec_dict))
                outf.write(block1)
                for ec, value in ec_dict.items():
                    block2 = struct.pack('<II', ec, len(value))
                    outf.write(block2)
                    for tup in value:
                        block3 = struct.pack('<QI', tup[0], tup[1])
                        outf.write(block3)
                    
    def parseCBIN(self, inf: BinaryIO) -> None:
        """
        parse the data blocks of a cy binary file
        """
        b = inf.read(12)
        while b:
            bc, n_ec = struct.unpack('<QI', b)
            self.bcd[bc] = dict()
            for i in range(n_ec):
                e = inf.read(8)
                ec, n_umi = struct.unpack('<II', e)
                u = inf.read(int(n_umi*12))
                umi_c = iter(struct.unpack('<'+'QI'*n_umi, u))
                self.bcd[bc][ec] = list(zip(umi_c, umi_c))
            b = inf.read(12)
            self.n_entries += 1
                    
