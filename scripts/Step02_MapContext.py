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

dnatab = maketrans("ACGTN+$", "TGCAN+$")
# MAPQ: 0: read unmapped, 1: read mapped to multiple locations >3: read map quality score 

qToCopy = {0:10, 1:5, 2:3, 3:2}

def revcomp(seq):
    return ''.join(reversed(seq.translate(dnatab)))

def loadFasta(filename):
    """ Parses a classically formatted and possibly 
        compressed FASTA file into a list of headers 
        and fragment sequences for each sequence contained"""
    if (filename.endswith(".gz")):
        fp = gzip.open(filename, 'rb')
    else:
        fp = open(filename, 'rb')
    # split at headers
    data = fp.read().split('>')
    fp.close()
    # ignore whatever appears before the 1st header
    data.pop(0)     
    headers = []
    sequences = []
    for sequence in data:
        lines = sequence.split('\n')
        headers.append(lines.pop(0).split()[0])
        # add an extra "+" to make string "1-referenced"
        sequences.append('+' + ''.join(lines))
    return (headers, sequences)

def getSamples(sample_file):
    df = pd.read_csv(filepath_or_buffer=sample_file, sep=',')
    df = df.replace(np.nan, '', regex=True)
    sample_list = df.loc[:,['sample', 'bwtfile', 'threshold']].values
    return sample_list

def merge_contexts(sample, side):
    filelist = glob.glob("./tmp/context/%s*%s.csv" %(sample,side))
    context_set = set()
    
    for filename in filelist:
        print filename
        try:
            df = pd.read_csv(filepath_or_buffer=filename, sep=',')
            data = df.iloc[:,].values
            for my_id,context,te,c1,c2 in data:
                context_set.add((my_id,context,te))
        except:
            pass

    infile = "tmp/bowtie_data/%s_%s.csv" %(sample,side)
    with open(infile,'wb') as fp:
        a = csv.writer(fp,delimiter=',')
        header = ['my_id', 'context', 'te']
        a.writerows([header])
        for (my_id,context,te) in context_set:
            a.writerows([[my_id,context,te]])
    print "Wrote file: %s [%d lines]\n" %(infile, len(context_set)) 
    
def TestForUnique(sample,side,bowtie_dir,species,expected_length):
    global genome
    bt2cmd = "module load python/3.5.1;bowtie2 -x %s/%s --no-head -r --end-to-end -k 4 %s.seq > %s.sam"
    designfile = "tmp/bowtie_data/%s_%s.csv" %(sample,side)
    
    t = designfile.rfind('.')
    root = designfile[:t]
    outfile = root + ".seq"
    
    fp = open(outfile, 'wb')
    design = pd.read_csv(filepath_or_buffer=designfile, sep=',')
    N = design.shape[0]
    probes = {}
    distinct_context= set()
    for index, row in enumerate(design.values):
        seq = row[1]
        fp.write(seq+"\n")
        distinct_context.add(seq)
    print "TOtal distinct contexts: %d" %(len(distinct_context))
    fp.close()
    print "Wrote %s (%d lines)" % (outfile, N)
    sys.stdout.flush()
    code = subprocess.call(bt2cmd % (bowtie_dir, species, root, root), shell=True)
    if (code == 0):
        print "Alignment completed"
    else:
        print "Alignment failed:" + (bt2cmd % (bowtie_dir, species, root, root))
        return
    
    samfile = outfile.replace('.seq', '.sam')
    columns = [str(i) for i in range(20)]
    df = pd.read_csv(filepath_or_buffer=samfile, names=columns, sep='\t', header=None)
    df = df.drop_duplicates(subset=['0'], keep=False)
    #df.to_csv("tmp/bowtie_data/test.csv")
    data = df.iloc[:,].values
    unique_locations = {}
    unmapped = 0
    pos_set = set()
    new_data = []
    unique = set()
    for fields in data:
        index = fields[0]
        chromo = fields[2]
        pos = fields[3]
        flags = int(fields[1])
        if chromo == '*':
            unmapped += 1
            continue
        alignment_score = -100
        if fields[11].find("AS:i:") == 0:
            alignment_score = int(fields[11].split(":")[-1])
        if alignment_score < 0:
            continue
        
        strand = '-' if flags & 16 else '+'
        new_seq = revcomp(fields[9]) if flags & 16 else fields[9]
        new_data.append([new_seq,chromo,pos,strand])
        unique.add(new_seq)
        
    locationfile = outfile.replace('.seq', '_location.csv')
    with open(locationfile,'wb') as fp:
        a = csv.writer(fp,delimiter=',')
        header = ['context', 'chromo', 'pos', 'strand']
        a.writerows([header])
        for d in new_data:
            a.writerows([d])
    
    print "Wrote file: %s [%d lines]" %(locationfile, len(new_data)) 
    print "unmapped: %d" %unmapped
    print "unique: %d" %len(unique)
    ########
    
    df1 = pd.read_csv(filepath_or_buffer=designfile, sep=',')
    df2 = pd.read_csv(filepath_or_buffer=locationfile, sep=',')
    
    result = pd.merge(df1, df2, how='right', on=['context'])
    data = result.iloc[:,].values
    new_data = []
    
    # ['id', 'side', 'context', 'chromo', 'pos', 'strand']
    for d in data:
        [my_id,context,te,chromo,pos,strand] = d[0:6]
        # Remove the prepedended "CHR"
        if chromo[0:2].lower() == 'chr':
            chromo = chromo[3:]
        if chromo not in genome.keys():
            continue
        plen = len(context)
        
        if strand == '+' and side == 'start':
            ref_prefix = genome[chromo][pos:pos+plen]
            ref_suffix = genome[chromo][pos+plen:pos+plen+25]
            other_context = genome[chromo][pos+plen:pos+plen+25]
            pos = pos + plen
        if strand == '+' and side == 'end':
            ref_prefix = genome[chromo][pos-25:pos]
            ref_suffix = genome[chromo][pos:pos+plen]
            other_context = genome[chromo][pos-25:pos]
        if strand == '-' and side == 'start':
            ref_prefix = genome[chromo][pos:pos+25]
            ref_suffix = genome[chromo][pos-plen:pos]
            other_context = genome[chromo][pos-25:pos]
        if strand == '-' and side == 'end':
            ref_prefix = genome[chromo][pos+plen:pos+plen+25]
            ref_suffix = genome[chromo][pos:pos+plen]
            other_context = genome[chromo][pos+plen:pos+plen+25]
            pos = pos + plen
        other_context = MultiStringBWT.reverseComplement(other_context) if strand == '-' else other_context
        if strand == '-':
            ref_prefix = MultiStringBWT.reverseComplement(ref_prefix)
            ref_suffix = MultiStringBWT.reverseComplement(ref_suffix)
            
        ed_th = .2*25
        if side == 'start':
            ed = lv.distance(ref_suffix[:25],te[:25])
            ref_te = 1 if ed <= ed_th else 0
        if side == 'end':
            ed = lv.distance(ref_prefix[-25:],te[-25:])
            ref_te = 1 if ed <= ed_th else 0
        new_data.append([my_id,context,te,ref_te,ref_prefix,ref_suffix,chromo,pos,strand,len(context)])
        #new_data.append([my_id,context,te,ref_prefix,ref_suffix,chromo,pos,strand])
    
    finalfile = locationfile.replace('location', 'UNIQUE')
    with open(finalfile,'wb') as fp:
        a = csv.writer(fp,delimiter=',')
        header = ['my_id', 'context', 'TE', 'ref_te', 'ref_prefix', 'ref_suffix', 'chromo', 'pos', 'strand', 'clen']
        a.writerows([header])
        for d in new_data:
            a.writerows([d])
    remove_duplicates(finalfile)
    return
    command = "rm ./tmp/bowtie_data/*.seq"
    rval = os.system(command)
    command = "rm ./tmp/bowtie_data/*.sam"
    rval = os.system(command)
    command = "rm ./tmp/bowtie_data/*_location.csv"
    rval = os.system(command)
    
def remove_duplicates(finalfile):
    df = pd.read_csv(filepath_or_buffer=finalfile, dtype= {'chromo': str, 'pos': int}, sep=',')
    df = df.replace(np.nan, '', regex=True)
    df.sort_values(by = ['clen'], inplace=True)
    df = df.drop_duplicates(subset=['chromo', 'pos'], keep = 'first')
    df.drop(columns=['clen'], inplace=True)
    df.sort_values(by = ['chromo', 'pos'], inplace=True)
    df.to_csv(finalfile, index=False)
    print "Wrote file: %s [%d lines]\n" %(finalfile, df.shape[0]) 
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_file', help='List of Sample and corresponding bwtdir')
    parser.add_argument('reference', help='Path to reference sequence')
    parser.add_argument('bowtie_dir', help='directory of bowtie index of reference genome')
    parser.add_argument('species', help='Name of species')
    parser.add_argument('C', help='Expected length of context')
    args = parser.parse_args()
    
    sample_file = args.sample_file
    reference = args.reference
    bowtie_dir = args.bowtie_dir
    species = args.species
    C = int(args.C)
    
    print "Path to SampleList: %s" %sample_file
    print "Path to REF: %s" %reference
    print "Name of species: %s" %species
    print "Path to the bowtie index: %s" %bowtie_dir
    print "expected length of context: %d" %C
    
    global genome
    genome = {}
    start_time = time.time()
    hdr, seq = loadFasta(reference)
    for i in xrange(len(hdr)):
        genome[hdr[i]] = seq[i]
    print "Reference sequence loaded in %s sec" %(time.time()-start_time)
    print genome.keys()
    
    if not os.path.exists(bowtie_dir):
        print "Need to build bowtie index"
        command = "mkdir %s" %bowtie_dir
        rval = os.system(command)
        build_command = "module load python/3.5.1;bowtie2-build %s ./%s/%s" %(reference,bowtie_dir,species)
        rval = os.system(build_command)

    command = "mkdir ./tmp/bowtie_data"
    rval = os.system(command)
        
    sample_list = getSamples(sample_file)
    for sample, bwtdir, threshold in sample_list:
        print "Started Mapping step for %s at: %s\n" %(sample, datetime.now())

        merge_contexts(sample, 'start')
        merge_contexts(sample, 'end')

        TestForUnique(sample, 'start', bowtie_dir, species, C)
        TestForUnique(sample, 'end', bowtie_dir, species, C)
