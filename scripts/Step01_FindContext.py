#!/bin/bash
import numpy as np
import pandas as pd
import sys
import math
import sqlite3
import os
import time
import glob
import Levenshtein as lv
import pysam
import MUSCython
import MUSCython.MultiStringBWTCython as MultiStringBWT
from datetime import datetime
import alignment as AL
import itertools
from multiprocessing import Pool
import argparse

dirs = ['TE', 'context', 'bowtie_data']
context_List = set()

if not(os.path.isdir("./tmp")):
    command = "mkdir ./tmp"
    rval = os.system(command)

for directory in dirs:
    new_dir = './tmp/%s' %directory
    if (os.path.isdir(new_dir)):
        continue
    command = "mkdir %s" %(new_dir)
    rval = os.system(command)

def getSamples(sample_file):
    df = pd.read_csv(filepath_or_buffer=sample_file, sep=',')
    df = df.replace(np.nan, '', regex=True)
    sample_list = df.loc[:,['sample', 'bwtfile', 'threshold']].values
    return sample_list

def get_te_info(te_file):
    df = pd.read_csv(filepath_or_buffer=te_file, sep=',')
    df = df.replace(np.nan, '', regex=True)
    data = df.loc[:,['my_id', 'TE', 'start_seed', 'end_seed']].values
    te_data = []
    for my_id,TE,start_seed,end_seed in data:
        te_data.append((my_id,TE,start_seed,end_seed))
    return te_data

def get_kmer_count(msbwt,kmer):
    cnt1 = msbwt.countOccurrencesOfSeq(kmer)
    cnt2 = msbwt.countOccurrencesOfSeq(MultiStringBWT.reverseComplement(kmer))
    return cnt1, cnt2

def growPrefix(msbwt, low, high, prefix, orig_kmer, startTE, current_ind, K, TEprefix_List):
    ed = int(math.ceil((len(startTE)*10.)/80))
    if len(prefix) > current_ind + 20:
        _,TEseq = AL.align_TE(prefix+orig_kmer,startTE[:current_ind+K])
        seq_partial = startTE[:current_ind+K]
        #print new_seq, seq_partial
        new_edit_distance = lv.distance(seq_partial, TEseq)
        #print len(TE), current_ind+K-5, new_edit_distance
        if new_edit_distance < ed and len(TEseq) >= current_ind-5:
            TEprefix_List.add(TEseq)
        return TEprefix_List
    
    for base in "ACGT":
        lo, hi = msbwt.findIndicesOfStr(base, (low, high))
        new_seq = base + prefix
        
        if len(prefix) < current_ind:
            seq_partial = startTE[:current_ind][-len(new_seq):]
            #print new_seq, seq_partial
            new_edit_distance = lv.distance(seq_partial, new_seq)
        else:
            new_edit_distance = 0

        if (hi - lo) > 0 and new_edit_distance < ed:
            TEprefix_List = TEprefix_List | growPrefix(msbwt, lo, hi, new_seq, orig_kmer, startTE, current_ind, K, TEprefix_List)
    return TEprefix_List

def growSuffix(msbwt, low, high, prefix, startTE, current_ind, K, TE_List):
    #startTE_rev = MultiStringBWT.reverseComplement(startTE)
    ed = int(math.ceil((len(startTE)*10.)/80))
    if len(prefix) == len(startTE):
        prefix_rev = MultiStringBWT.reverseComplement(prefix)
        new_edit_distance = lv.distance(startTE, prefix_rev)
        if new_edit_distance < ed:
            TE_List.add(prefix_rev)
        #print prefix_rev
        return TE_List
    
    for base in "ACGT":
        lo, hi = msbwt.findIndicesOfStr(base, (low, high))
        new_seq = base + prefix
        seq_partial = startTE[:len(prefix)]
        new_seq_rev = MultiStringBWT.reverseComplement(new_seq)
        
        new_edit_distance = lv.distance(seq_partial, new_seq_rev)
        #print new_seq, new_seq_rev, new_edit_distance
        
        if (hi - lo) > 0 and new_edit_distance < ed:
            TE_List = TE_List | growSuffix(msbwt, lo, hi, new_seq, startTE, current_ind, K, TE_List)
    return TE_List

def growContext(msbwt, low, high, prefix, expected_length, threshold):
    global context_List
    for base in "ACGT":
        lo, hi = msbwt.findIndicesOfStr(base, (low, high))
        new_seq = base + prefix
        if (hi - lo) >= threshold:
            growContext(msbwt, lo, hi, new_seq, expected_length, threshold)
        else:
            if len(prefix) == expected_length:
                context_List.add(prefix)
                return
                
def getSeed(sample, msbwt, TE, T, K):
    start_seeds = []
    for i in range(20,T-K+1,1):
        kmer = TE[i:i+K]
        c1,c2 = get_kmer_count(msbwt,kmer)
        start_seeds.append((kmer,c1+c2))
    start_seeds.sort(key=lambda x: (-x[1]))
    start_seed = start_seeds[0][0]
    
    te_len = len(TE)
    end_seeds = []
    for i in range(20,T-K+1,1):
        start_ind = te_len-K-i
        kmer = TE[start_ind:start_ind+K]
        c1,c2 = get_kmer_count(msbwt,kmer)
        end_seeds.append((kmer,c1+c2))
    end_seeds.sort(key=lambda x: (-x[1]))
    end_seed = end_seeds[0][0]
    
    return start_seed, TE[:T], end_seed, TE[-T:]

def getContext(sample, msbwt, start_seed, startTE, end_seed, endTE, T, expected_length, side, threshold, my_id):
    global context_List
    if side == 'start':
        seq = startTE
        seed = start_seed
        current_ind = seq.find(seed)
        K = len(seed)
    else:
        seq = MultiStringBWT.reverseComplement(endTE)
        seed = MultiStringBWT.reverseComplement(end_seed)
        current_ind = seq.find(seed)
        K = len(seed)
        
    te_len = len(seq)
    TEprefix_List = set()
    TE_List = set()
    context_List = set()
    
    lo, hi = msbwt.findIndicesOfStr(seed)
    
    TEprefix_List = growPrefix(msbwt, lo, hi, '', seed, seq, current_ind, K, TEprefix_List)
    
    for TEprefix in TEprefix_List:
        TE_without_seed = TEprefix[:current_ind]
        TE_without_seed_rev = MultiStringBWT.reverseComplement(TE_without_seed)
        lo, hi = msbwt.findIndicesOfStr(TE_without_seed_rev)
        newTE_List = set()
        newTE_List = growSuffix(msbwt, lo, hi, TE_without_seed_rev, seq, current_ind, K, newTE_List)
        TE_List = TE_List | newTE_List
    TEprefix_List = set()
    
    write_type = 'w'
    context_file = "tmp/context/%s_%s_%s.csv" %(sample,my_id,side)
    fp_context = open(context_file, write_type)
    fp_context.write("%s,%s,%s,%s,%s\n" %("my_id","context","TE","Fcnt","Rcnt"))
    
    TE_file = "tmp/TE/%s_%s_%s.csv" %(sample,my_id,side)
    fp_TE = open(TE_file, write_type)
    fp_TE.write("%s,%s,%s,%s\n" %("my_id","TE","Fcnt","Rcnt"))
    
    for TEseq in TE_List:
        lo, hi = msbwt.findIndicesOfStr(TEseq)
        context_List = set()
        growContext(msbwt, lo, hi, '', expected_length, threshold)
        for context in context_List:
            if side == 'end':
                new_context = MultiStringBWT.reverseComplement(context)
                new_TE = MultiStringBWT.reverseComplement(TEseq)
                kmer = new_TE + new_context
            if side == 'start':
                new_context = context
                new_TE = TEseq
                kmer = new_context + new_TE
            c1,c2 = get_kmer_count(msbwt, kmer)
            if c1+c2 > threshold:
                fp_context.write("%s,%s,%s,%d,%d\n" %(my_id,new_context,new_TE,c1,c2))
                c11,c22 = get_kmer_count(msbwt, new_TE)
                fp_TE.write("%s,%s,%d,%d\n" %(my_id,new_TE,c11,c22))
        context_List = set()
    
    fp_context.close()
    fp_TE.close()
    
def solve(sample,bwtfile,te_list_part,T,K,expected_length,threshold):
    start_time = time.time()
    print sample, bwtfile
    msbwt = MultiStringBWT.loadBWT(bwtfile, useMemmap=False)
    print "loaded bwt in %s" %(time.time()-start_time)
    for (my_id,TE,start_seed,end_seed) in te_list_part:
        if start_seed == '' or end_seed == '':
            start_seed, startTE, end_seed, endTE = getSeed(sample, msbwt, TE, T, K)
        else:
            startTE, endTE = TE[:T], TE[-T:]
        print my_id,TE,start_seed,end_seed, startTE, endTE
        for side in ["start", "end"]:
            getContext(sample,msbwt,start_seed,startTE,end_seed,endTE,T,expected_length,side,threshold,my_id)
            
def run(sample_file,te_file,expected_length,T,K):
    sample_list = getSamples(sample_file)
    te_data = get_te_info(te_file)
    for sample, bwtfile, threshold in sample_list:
        start_time = time.time()
        print "Started for %s at: %s" %(sample, datetime.now())
        solve(sample,bwtfile,te_data,T,K,expected_length,threshold)
        end_time = time.time()
        print "Total time for %s: %s\n" %(sample, time.time()-start_time)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_file', help='List of Sample and corresponding bwtdir')
    parser.add_argument('te_file', help='Path to TE Template')
    parser.add_argument('C', help='Expected length of context')
    parser.add_argument('T', help='Length of proximal and distal TE')
    parser.add_argument('K', help='Length of proximal and distal seed')
    args = parser.parse_args()
    
    sample_file = args.sample_file
    te_file = args.te_file
    C = int(args.C)
    T = int(args.T)
    K = int(args.K)

    print "Path to SampleList: %s" %sample_file
    print "Path to TE template: %s" %te_file
    print "expected length of context: %d" %C
    print "Length of proximal and distal TE: %d" %T
    print "Length of proximal and distal seed: %d" %K
    run(sample_file,te_file,C,T,K)
