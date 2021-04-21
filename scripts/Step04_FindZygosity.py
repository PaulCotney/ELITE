import gzip
import time
import pandas as pd
import numpy as np
import subprocess
import csv
import sys
from string import maketrans
import pysam
import MUSCython
import MUSCython.MultiStringBWTCython as MultiStringBWT
import os
import glob
import Levenshtein as lv
import subprocess
from string import maketrans
import argparse
from datetime import datetime
import glob

def getSamples(sample_file):
    df = pd.read_csv(filepath_or_buffer=sample_file, sep=',')
    df = df.replace(np.nan, '', regex=True)
    sample_list = df.loc[:,['sample', 'bwtfile', 'threshold']].values
    return sample_list

def get_kmers(sample, ref_te):
    file_dir = "tmp/FinalOutput"
    filename = "%s/%s_%d.csv" %(file_dir, sample, ref_te)
    kmer_list = set()

    df = pd.read_csv(filepath_or_buffer=filename, dtype= {'chromo': str, 'pos': int}, sep=',')
    df = df.replace(np.nan, '', regex=True)
    header = list(df)
    data = df.iloc[:,].values
    for TEi_id,my_id,side,chromo,pos,strand,context,TE,TSD in data:
        kmer_list.add((TEi_id,my_id,side,chromo,pos,strand,ref_te,context,TE,TSD))
    return list(kmer_list)

def get_kmer_count(msbwt, kmer):
    c1 = msbwt.countOccurrencesOfSeq(kmer)
    c2 = msbwt.countOccurrencesOfSeq(MultiStringBWT.reverseComplement(kmer))
    return c1+c2

def growContext(msbwt, low, high, prefix, expected_length, context_List):
    if len(prefix) == expected_length:
        context_List.add(prefix)
        return context_List
    for base in "ACGT":
        lo, hi = msbwt.findIndicesOfStr(base, (low, high))
        new_seq = base + prefix
        if (hi - lo) > 0:
            context_List = context_List | growContext(msbwt, lo, hi, new_seq, expected_length, context_List)
    return context_List

def get_zygosity(sample, bwtfile, threshold, expected_length, ref_te):
    msbwt = MultiStringBWT.loadBWT(bwtfile, useMemmap = False)
    oc_length = 25
    ed_th = .2*oc_length
    kmer_list = get_kmers(sample,ref_te)
    zygosity_data = []
    zygo_dict = {}
    for TEi_id,my_id,side,chromo,pos,strand,ref_te,context,TE,TSD in kmer_list:
        #if (TEi_id, ref_te) != (47,0):
        #    continue
        #print TEi_id,my_id,side,chromo,pos,strand,ref_te,context,TE,TSD
        new_context = MultiStringBWT.reverseComplement(context) if side == 'start' else context
        lo,hi = msbwt.findIndicesOfStr(new_context)
        context_List = set()
        context_List = growContext(msbwt, lo, hi, '', oc_length, context_List)
        zygo_dict[TEi_id] = zygo_dict.get(TEi_id,set([1]))
        for oc in context_List:
            oc = MultiStringBWT.reverseComplement(oc) if side == 'start' else oc
            ed = lv.distance(oc,TE[:oc_length]) if side == 'start' else lv.distance(oc,TE[-oc_length:])
            if ed > ed_th:
                zygo_dict[TEi_id] = zygo_dict[TEi_id] | set([0])
        zygosity_data.append([TEi_id,my_id,side,chromo,pos,strand,ref_te,context,TE,TSD])
    
    individual_file = "IndividualAnalysis/%s_%s.csv" %(sample,ref_te)
    zygosity_data.sort(key=lambda x:(x[3],x[4]))
    header = ['TEi_id', 'my_id', 'side', 'chromo', 'pos', 'strand', 'ref_te', 'context','TE','TSD', 'zygosity']
    with open(individual_file,'wb') as fp:
        a = csv.writer(fp,delimiter=',')
        a.writerows([header])
        for d in zygosity_data:
            zygosity = 'heterozygous' if len(zygo_dict[d[0]]) == 2 else 'homozygous'
            d.append(zygosity)
            a.writerows([d])
    print "Wrote file: %s [%d lines]" %(individual_file, len(zygosity_data))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_file', help='List of Sample and corresponding bwtdir')
    parser.add_argument('te_file', help='Path to TE Template')
    parser.add_argument('C', help='Expected length of context')
    
    args = parser.parse_args()
    
    sample_file = args.sample_file
    te_file = args.te_file
    C = int(args.C)
    
    if not(os.path.isdir("./IndividualAnalysis")):
        command = "mkdir ./IndividualAnalysis"
        rval = os.system(command)
        
    sample_list = getSamples(sample_file)
    for sample, bwtdir, threshold in sample_list:
        bwtfile = "%s/%s" %(bwtdir,sample)
        print "Started for %s at: %s" %(sample, datetime.now())
        for ref_te in [0,1]:
            get_zygosity(sample, bwtdir, threshold, C, ref_te)