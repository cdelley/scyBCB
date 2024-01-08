# Cyrille L. Delley 2023
#
# Each dictionary input represents a different barcode bead structure.
# each entry consists of a dictionary with three keys:
# 'bc'     : the position of the split pool barcode blocks
# 'umi'    : the position of the umi
# 'feature': the read segment that contains biological information
#
# all positions are encoded as touple (X, Y, Z) where X denotes
# the read file (0 or 1) from a pired end experiment. Y and Z are the
# the start and end positions of the DNA string. All counting is zero based.
# each entry can contain arbitrary many touples for instance to encode split-pool
# barcoding shemes with two blocks like InDrops [(x1, y1, z1), (x2, y2, z2)] or four
# for instance for Fluent_v3 beads [(), (), (), ()]. 'umi' and 'feature' can be empty lists.


barcode_beads = {
    'BHBv2_spat' : {
        'bc'  : [(0,0,11),(0,12,20),(0,24,32)],
        'umi' : [(1,0,6)],# no UMI on primer
        'feature' : [(1,6,15)]
    },
    'BHBv2_spat1' : {
        'bc'  : [(0,0,11),(0,12,20),(0,24,32)],
        'umi' : [(1,0,6)],# no UMI on primer
        'feature' : [(1,6,14)]
    },
    'BHBv2' : {
        'bc'  : [(0,0,11),(0,12,20),(0,24,32)],
        'umi' : [(0,32,38)],
        'feature' : [(0,24,32)]
    },
    '10x_spat' : {
        'bc'  : [(0,0,16)],
        'umi' : [(0,16,28)], #(0,16,28) on bead, (1,0,6) on tag
        'feature' : [(1,6,15)]
    },
    '10x_spat1' : {
        'bc'  : [(0,0,16)],
        'umi' : [(0,16,28)], #(0,16,28), (1,0,6) 
        'feature' : [(1,6,14)]
    },
    '10x_spat2' : {
        'bc'  : [(0,0,16)],
        'umi' : [(0,16,28)], #(0,16,28), (1,0,6) 
        'feature' : [(1,5,14)]
    },
    '10x_spat3' : {
        'bc'  : [(0,0,16)],
        'umi' : [(0,16,28)], #(0,16,28), (1,0,6) 
        'feature' : [(1,7,16)]
    },
    'BHBv2_Ab' : {
        'bc'  : [(0,0,11), (0,12,20), (0,24,32)], # three barcodes
        'umi' : [(0,47,57)], 
        'feature' : [(1,0,15)]
    },
    'fluentv3_spat' : {
        'bc'  : [(0,0,8), (0,11,17), (0,20,26), (0,31,39)], # four barcodes
        'umi' : [(0,39,51)], 
        'feature' : [(1,0,9)]
    },
    'BHBv2_spat_fwd' : {
        'bc'  : [(0,0,11),(0,12,20),(0,24,32)],
        'umi' : [(0,61,67)],
        'feature' : [(0,67,76)]
    },
}

