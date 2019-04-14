# compute_rnaQuast.py
#
# Laura Tung
#
# Usage: compute_rnaQuast.py <result_dir>
# 
# <result_dir> is the directory containing rnaQUAST_output and rnaQUAST_output_1

import sys
import numpy as np


def load_data(dataset):

    loaded_isoform = np.loadtxt(dataset+"/isoform_data", dtype='int', usecols=(2, 3, 4))

    loaded_matched =  np.loadtxt(dataset+"/matched_data", dtype='int', usecols=(1, 2, 3))

    return loaded_isoform, loaded_matched


if __name__ == "__main__":
    
    # load rnaQUAST_output (75% and 95%)
    loaded_isoform, loaded_matched = load_data(sys.argv[1]+"/rnaQUAST_output")

    # load rnaQUAST_output_1 (0% and 50%)
    loaded_isoform_1, loaded_matched_1 = load_data(sys.argv[1]+"/rnaQUAST_output_1")

    print("\t\t\t\tScallop\t\tStringTie\tIsoseq")

    # assembled isoforms
    print("0-50%-assembled-isoforms\t"+str(loaded_isoform_1[0][0]-loaded_isoform_1[1][0])+"\t\t"+str(loaded_isoform_1[0][1]-loaded_isoform_1[1][1])+"\t\t"+str(loaded_isoform_1[0][2]-loaded_isoform_1[1][2]))
    print("50-75%-assembled-isoforms\t"+str(loaded_isoform_1[1][0]-loaded_isoform[0][0])+"\t\t"+str(loaded_isoform_1[1][1]-loaded_isoform[0][1])+"\t\t"+str(loaded_isoform_1[1][2]-loaded_isoform[0][2]))
    print("75-95%-assembled-isoforms\t"+str(loaded_isoform[0][0]-loaded_isoform[1][0])+"\t\t"+str(loaded_isoform[0][1]-loaded_isoform[1][1])+"\t\t"+str(loaded_isoform[0][2]-loaded_isoform[1][2]))
    print("95-100%-assembled-isoforms\t"+str(loaded_isoform[1][0])+"\t\t"+str(loaded_isoform[1][1])+"\t\t"+str(loaded_isoform[1][2]))

    # matched
    print("0-50%-matched-transcripts\t"+str(loaded_matched_1[0][0]-loaded_matched_1[1][0])+"\t\t"+str(loaded_matched_1[0][1]-loaded_matched_1[1][1])+"\t\t"+str(loaded_matched_1[0][2]-loaded_matched_1[1][2]))
    print("50-75%-matched-transcripts\t"+str(loaded_matched_1[1][0]-loaded_matched[0][0])+"\t\t"+str(loaded_matched_1[1][1]-loaded_matched[0][1])+"\t\t"+str(loaded_matched_1[1][2]-loaded_matched[0][2]))
    print("75-95%-matched-transcripts\t"+str(loaded_matched[0][0]-loaded_matched[1][0])+"\t\t"+str(loaded_matched[0][1]-loaded_matched[1][1])+"\t\t"+str(loaded_matched[0][2]-loaded_matched[1][2]))
    print("95-100%-matched-transcripts\t"+str(loaded_matched[1][0])+"\t\t"+str(loaded_matched[1][1])+"\t\t"+str(loaded_matched[1][2]))

