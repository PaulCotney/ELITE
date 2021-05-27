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
import argparse
import csv

extension, extension2, extension3, extension4 = set(), set(), set(), set()
   
def getSamples(sample_file):
    df = pd.read_csv(filepath_or_buffer=sample_file, sep=',')
    df = df.replace(np.nan, '', regex=True)
    sample_list = df.loc[:,['sample', 'bwtfile', 'threshold']].values
    return sample_list

def merge_files(sample):
    filename1 = "tmp/bowtie_data/%s_start_UNIQUE.csv" %sample
    df1 = pd.read_csv(filepath_or_buffer=filename1, dtype= {'chromo': str, 'pos': int}, sep=',')
    df1['side'] = 'start'
    
    filename2 = "tmp/bowtie_data/%s_end_UNIQUE.csv" %sample
    df2 = pd.read_csv(filepath_or_buffer=filename2, dtype= {'chromo': str, 'pos': int}, sep=',')
    df2['side'] = 'end'
    
    merged = pd.concat([df1,df2])
    merged.sort_values(by=['chromo','pos'], inplace=True)
    #merged['distance'] = merged['pos'].diff()
    
    merged.to_csv("tmp/Tables/%s_merged.csv" %sample, index = False)
    print "Wrote file for %s" %sample 
    
def remove_duplicate(sample):
    filename = "tmp/Tables/%s_merged.csv" %sample
    df = pd.read_csv(filepath_or_buffer=filename, dtype= {'chromo': str, 'pos': int}, sep=',')
    df = df.drop_duplicates(subset = ['chromo', 'pos', 'strand'], keep = 'first')
    df.to_csv("tmp/Tables/%s_mapped.csv" %sample, index = False)
    
def get_kmer_count(msbwt, kmer):
    c1 = msbwt.countOccurrencesOfSeq(kmer)
    c2 = msbwt.countOccurrencesOfSeq(MultiStringBWT.reverseComplement(kmer))
    return c1+c2

def get_te_dict(te_file):
    df = pd.read_csv(filepath_or_buffer=te_file, sep=',')
    data = df.iloc[:,].values
    te_dict = {}
    for d in data:
        [my_id, TEseq] = d[0:2]
        te_dict[my_id] = TEseq
    return te_dict

def find_fdr(a,b):
    la = len(a)
    lb = len(b)
    FDR = ""
    for i in range(2,la,1):
        current_fdr = a[-i:]
        if b.startswith(current_fdr) and len(current_fdr) > len(FDR):
            FDR = current_fdr
    return FDR
    
def merge_sides_v1(sample, ref_te):
    filename = "tmp/Tables/%s_mapped.csv" %(sample)
    df = pd.read_csv(filepath_or_buffer=filename, dtype= {'chromo': str, 'pos': int}, sep=',')
    df = df.replace(np.nan, '', regex=True)
    df = df[(df['ref_te'] == ref_te)]
    chromosomes = sorted(set(df.loc[:,'chromo'].values))

    header = list(df)
    count_dict = {}
    all_differences = []
    sides = set(['start', 'end'])
    seq_dict = {}
    for c in chromosomes:
        total = set()
        pos_list = set()
        new_df = df[(df['chromo'] == c) & (df['strand'] == '+')]
        new_df = new_df.sort_values(by=['chromo','pos','strand'])
        data = new_df.iloc[:,].values
        for d in data:
            my_id = d[0]
            [chromo,pos,strand,side] = d[6:10]
            context, te_seq = d[1:3]
            total.add((chromo,pos,strand))
            pos_list.add((pos,strand,side,context,te_seq,my_id))
            #pos_list.add(int(pos.replace(',','')))
        pos_list = sorted(pos_list)
        diff = [(pos_list[i+1][0]-pos_list[i][0], (pos_list[i],pos_list[i+1])) for i in range(len(pos_list)-1) if pos_list[i][1]==pos_list[i+1][1]]
        #print "1. len of all diff: %d" %len(diff)
        
        for (d,((p1,s1,side1,context1,te_seq1,my_id1),(p2,s2,side2,context2,te_seq2,my_id2))) in diff:
            all_differences.append((d, (side1,c,p1,s1,context1,te_seq1,my_id1), (side2,c,p2,s2,context2,te_seq2,my_id2)))

        pos_list = set()
        new_df = df[(df['chromo'] == c) & (df['strand'] == '-')]
        new_df = new_df.sort_values(by=['chromo','pos','strand'])
        data = new_df.iloc[:,].values
        for d in data:
            my_id = d[0]
            [chromo,pos,strand,side] = d[6:10]
            context,te_seq = d[1:3]
            total.add((chromo,pos,strand))
            pos_list.add((pos,strand,side,context,te_seq,my_id))
        pos_list = sorted(pos_list)
        diff = [(pos_list[i+1][0]-pos_list[i][0], (pos_list[i],pos_list[i+1])) for i in range(len(pos_list)-1) if pos_list[i][1]==pos_list[i+1][1]]
        #print "2. len of all diff: %d" %len(diff)
        for (d,((p1,s1,side1,context1,te_seq1,my_id1),(p2,s2,side2,context2,te_seq2,my_id2))) in diff:
            all_differences.append((d, (side1,c,p1,s1,context1,te_seq1,my_id1), (side2,c,p2,s2,context2,te_seq2,my_id2)))
        count_dict[c] = len(total)
    
    all_differences = sorted(all_differences)
    
    mapped_both_ends = set()
    skip_pos = set()
    found_both_sides = set()
    fdr_dict = {}
    te_dict = get_te_dict(te_file)
    for d, (side1,c1,p1,s1,context1,te_seq1,my_id1), (side2,c2,p2,s2,context2,te_seq2,my_id2) in all_differences:
        expected_d = max(len(te_dict[my_id1]), len(te_dict[my_id2])) + 100 if ref_te == 1 else 20
        if d > expected_d or side1 == side2:
            continue
        
        if side1 == 'start':
            a,b = context1, context2
            te_start, te_end = te_seq1, te_seq2
            skip_pos.add((c2,p2,s2))
        if side2 == 'start':
            a,b = context2, context1
            te_start, te_end = te_seq2, te_seq1
            skip_pos.add((c1,p1,s1))
        #fdr1 = find_fdr(a,b)
        fdr = b[:d] if ref_te == 0 else find_fdr(a,b)
        fdr_dict[(side1,c1,p1,s1)] = fdr
        fdr_dict[(side2,c2,p2,s2)] = fdr

        mapped_both_ends.add((side1,c1,p1,s1))
        found_both_sides.add((side1,c1,p1,s1))
        found_both_sides.add((side2,c2,p2,s2))
        seq_dict[(c1,p1,s1)] = (a,te_start,te_end,b)
        seq_dict[(c2,p2,s2)] = (a,te_start,te_end,b)
    
    #print skip_pos
    filename = "tmp/Tables/%s_mapped.csv" %(sample)
    df = pd.read_csv(filepath_or_buffer=filename, dtype= {'chromo': str, 'pos': int}, sep=',')
    df = df.replace(np.nan, '', regex=True)
    df = df[(df['ref_te'] == ref_te)]
    data = df.iloc[:,].values
    new_data = []
    need_other_end = []
    #my_id,chromo,pos,strand,prefix,TEstart,TEend,suffix,FDR
    for d in data:
        my_id = d[0]
        side = d[9]
        a,te_start,te_end,b = '', '', '', ''
        if side == 'start':
            a,te_start = d[1:3]
        if side == 'end':
            te_end,b = d[1:3]
        [chromo,pos,strand,side] = d[6:10]
        
        row_data = [my_id]
        FDR = ""
        if (side,chromo,pos,strand) in found_both_sides:
            FDR = fdr_dict[(side,chromo,pos,strand)]
            (a,te_start,te_end,b) = seq_dict[(chromo,pos,strand)]
            
        if len(FDR) == 0 or (chromo,pos,strand) in skip_pos:
            continue
        
        row_data = [my_id, 'start',chromo,pos,strand,a,te_start,FDR]
        new_data.append(row_data)
        row_data = [my_id, 'end',chromo,pos,strand,b,te_end,FDR]
        new_data.append(row_data)
    
    #print sorted([len(fdr) for fdr in fdr_dict.values()])
    header = ['my_id', 'side', 'chromo', 'pos', 'strand', 'context','TE', 'TSD']
    filename_both_sides = "tmp/FoundTSD/%s_%d.csv" %(sample,ref_te)
    with open(filename_both_sides,'wb') as fp:
        a = csv.writer(fp,delimiter=',')
        a.writerows([header])
        for d in new_data:
            a.writerows([d])
    print "Wrote file: %s [%d lines]" %(filename_both_sides, len(new_data))

def get_ref_dict(sample,side,expected_length):
    filename = 'tmp/bowtie_data/%s_%s_UNIQUE.csv' %(sample,side)
    df = pd.read_csv(filepath_or_buffer=filename, sep=',',dtype= {'chromo': str, 'pos': int})
    data = df.loc[:,].values
    ref_dict = {}
    for d in data:
        [my_id, context, te, ref_te, ref_prefix, ref_suffix,chromo, pos, strand] = d[:9]
        ref_dict[(chromo,pos,strand)] = [ref_prefix[-expected_length:],ref_suffix[:expected_length]]
    return ref_dict

def still_unknown_FDR(sample, ref_te):  
    found_TE = set()
    filename = 'tmp/FoundTSD/%s_%d.csv' %(sample,ref_te)
    df1 = pd.read_csv(filepath_or_buffer=filename, sep=',',dtype= {'chromo': str, 'pos': int})
    df1 = df1.replace(np.nan, '', regex=True)

    data1 = df1.iloc[:,].values
    for my_id,side,chromo,pos,strand,context,TE,FDR in data1:
        found_TE.add((context,TE))

    filename = 'tmp/Tables/%s_mapped.csv' %(sample)
    df2 = pd.read_csv(filepath_or_buffer=filename, sep=',',dtype= {'chromo': str, 'pos': int})
    df2 = df2.replace(np.nan, '', regex=True)
    #df2 = df2[(df2['ref_te'] == 0) & (df2['side'] == side)]
    df2 = df2[(df2['ref_te'] == ref_te)]
    
    #print found_TE
    data2 = df2.iloc[:,].values
    need_to_find = set()
    
    for my_id,context,TE,ref_te,ref_prefix,ref_suffix,chromo,pos,strand,side in data2:
        if (context,TE) in found_TE:
            continue
        need_to_find.add((my_id,side,chromo,pos,strand,context,TE))
    return need_to_find
    
def merge_all_results(sample, ref_te):
    filename = 'tmp/FoundTSD/%s_%d.csv' %(sample,ref_te)
    df = pd.read_csv(filepath_or_buffer=filename, sep=',',dtype= {'chromo': str, 'pos': int})
    df = df.replace(np.nan, '', regex=True)
    
    noTSD = still_unknown_FDR(sample,ref_te)
    new_data = []
    data = df.iloc[:,].values
    for d in data:
        new_data.append(list(d))
    for my_id,side,chromo,pos,strand,context,TE in noTSD:
        new_data.append([my_id,side,chromo,pos,strand,context,TE,''])
        
    header = ['my_id','side','chromo','pos','strand','context','TE','TSD']
    filename = 'tmp/FinalOutput/%s_%d.csv' %(sample,ref_te)
    with open(filename,'wb') as fp:
        a = csv.writer(fp,delimiter=',')
        a.writerows([header])
        for row_data in new_data:
            a.writerows([row_data])
    print "Wrote file %s [%d lines]" %(filename, len(new_data))

def merge_close_pos(sample, d, ref_te):
    filename = 'tmp/FinalOutput/%s_%d.csv' %(sample,ref_te)
    df = pd.read_csv(filepath_or_buffer=filename, sep=',',dtype= {'chromo': str, 'pos': int})
    df = df.replace(np.nan, '', regex=True)
    chromosomes = sorted(set(df.loc[:,'chromo'].values))
    cluster_by_chromo = []
    pos_to_keep = set()
    for chromo in chromosomes:
        filename = 'tmp/FinalOutput/%s_%d.csv' %(sample,ref_te)
        df = pd.read_csv(filepath_or_buffer=filename, sep=',',dtype= {'chromo': str, 'pos': int})
        df = df.replace(np.nan, '', regex=True)
        df = df.sort_values(by = ['chromo', 'pos'])
        df = df[(df['chromo'] == chromo)]
        if df.shape[0] == 0:
            continue
        pos_list = df.loc[:,'pos'].values
        cluster_rep = pos_list[0]
        cluster_mem  = [cluster_rep]
        clusters = [(cluster_rep, cluster_mem)]
        for pos in pos_list[1:]:
            prev = clusters[-1]
            cluster_rep = prev[0]
            cluster_mem = prev[1]
            new_cluster_rep = pos
            if new_cluster_rep - cluster_rep < d:
                cluster_mem.append(pos)
                clusters[-1] = ((new_cluster_rep, cluster_mem))
            else:
                clusters.append((new_cluster_rep, [new_cluster_rep]))
        cluster_by_chromo = cluster_by_chromo + [(chromo,i[0],i[1]) for i in clusters]
    
    pos_to_TEi_id = {}
    for TEi_id,(chromo,cluster_rep,cluster_mem) in enumerate(cluster_by_chromo):
        for mem in cluster_mem:
            pos_to_TEi_id[(chromo,mem)] = TEi_id
            
    print len(pos_to_TEi_id)
    filename = 'tmp/FinalOutput/%s_%d.csv' %(sample,ref_te)
    df = pd.read_csv(filepath_or_buffer=filename, sep=',',dtype= {'chromo': str, 'pos': int})
    df = df.replace(np.nan, '', regex=True)
    new_data = []
    header = list(df)
    data = df.iloc[:,].values
    for my_id,side,chromo,pos,strand,context,TE,FDR in data:
        if (chromo,pos) not in pos_to_TEi_id.keys():
            continue
        TEi_id = pos_to_TEi_id[(chromo,pos)]
        new_data.append([TEi_id,my_id,side,chromo,pos,strand,context,TE,FDR])
    
    new_header = ['TEi_id'] + header
    new_data.sort(key=lambda x:(x[3],x[4]))
    with open(filename,'wb') as fp:
        a = csv.writer(fp,delimiter=',')
        a.writerows([new_header])
        for row_data in new_data:
            a.writerows([row_data])
    print "Wrote file %s [%d lines]" %(filename, len(new_data))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_file', help='List of Sample and corresponding bwtdir')
    parser.add_argument('te_file', help='Path to TE Template')
    args = parser.parse_args()
    
    sample_file = args.sample_file
    te_file = args.te_file
    
    print "Path to SampleList: %s" %sample_file
    print "Path to TE template: %s" %te_file
    
    dirs = ['FoundTSD', 'Tables', 'FinalOutput']
    if not(os.path.isdir("./tmp")):
        command = "mkdir ./tmp"
        rval = os.system(command)

    for directory in dirs:
        new_dir = './tmp/%s' %directory
        if (os.path.isdir(new_dir)):
            continue
        command = "mkdir %s" %(new_dir)
        rval = os.system(command)

    sample_list = getSamples(sample_file)
    for sample, bwtdir, threshold in sample_list:
        bwtfile = "%s/%s" %(bwtdir,sample)
        print "Started for %s at: %s" %(sample, datetime.now())
        start_time = time.time()
        #msbwt = MultiStringBWT.loadBWT(bwtfile, useMemmap = False)
        print "loaded bwt in %s" %(time.time()-start_time)
        merge_files(sample)
        remove_duplicate(sample)
        for ref_te in [0,1]:
            merge_sides_v1(sample,ref_te)
            merge_all_results(sample,ref_te)
            merge_close_pos(sample, 20, ref_te)
