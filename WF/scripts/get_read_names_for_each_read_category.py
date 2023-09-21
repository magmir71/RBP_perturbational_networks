#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# Author : Aleksei Mironov
# Company: Mihaela Zavolan, Biozentrum, Basel
# This script is part of the Zavolan lab ZARP pipeline.
# -----------------------------------------------------------------------------

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess
import pandas as pd
import itertools

def main():
    """For each read category (duplication level x (mm or um)), get two lists: a list of read names and their weights"""

    __doc__ = "For each read category (duplication level x (mm or um)), get two lists: a list of read names and their weights"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--grouped_bed",
                        dest="grouped_bed_file_path",
                        help="Path to the input grouped bed file produced with groupBamToBedFile",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--bed",
                        dest="bed_file_path",
                        help="path for the input bed file",
                        required=True,
                        metavar="FILE")
    
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    grouped_bed_file_path = options.grouped_bed_file_path
    bed_file_path = options.bed_file_path
    
    outdir = os.path.dirname(grouped_bed_file_path)+'/read_categories/'
    sample_id = os.path.basename(grouped_bed_file_path).replace('.dedup.sorted.indexed.grouped.bed','')
    
    d_cats = [0,1,2,3,4,5]
    mm_modes = ['um','mm']
    
    read_categories = list(itertools.product(*[d_cats,mm_modes])) # cartesian product
    
    grouped_bed_file = pd.read_csv(grouped_bed_file_path,delimiter="\t",index_col=None,header=None)
    grouped_bed_file.columns = [0,1,2,'name','w',5,'d','wd','d_cat']

    bed_file = pd.read_csv(bed_file_path,delimiter="\t",index_col=None,header=None)
    bed_file.columns = [0,1,2,'name',4,5]
    
    for read_category in read_categories:
        cat_name = '_'.join(list(pd.Series(read_category).astype('str')))
        out_subdir = outdir+cat_name+'/'
        command = 'mkdir -p '+out_subdir
        out = subprocess.check_output(command, shell=True)    
    
        output_read_names_path = out_subdir+sample_id+'.'+cat_name+'.read_names.txt'
        
        d_cat = read_category[0]
        if read_category[1]=='um':
            cur_reads = grouped_bed_file.loc[(grouped_bed_file['d_cat']==d_cat)&(grouped_bed_file['w']==1)]
        else:
            cur_reads = grouped_bed_file.loc[(grouped_bed_file['d_cat']==d_cat)&(grouped_bed_file['w']<1)]
        if len(cur_reads)==0:
            command = 'echo "empty" > '+output_read_names_path
            out = subprocess.check_output(command, shell=True)
            continue
        cur_reads = pd.merge(bed_file,cur_reads[[0,1,2,5]],how='inner',on=[0,1,2,5])
        cur_reads[['name']].drop_duplicates().to_csv(output_read_names_path, sep=str('\t'),header=False,index=None)
            
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)