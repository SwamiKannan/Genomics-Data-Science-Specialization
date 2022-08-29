#!/usr/bin/env python
# coding: utf-8
__author__ = "Swaminathan Kannan"

def readf(textf): #Reads the fasta file
    seq_dict={}
    with open(textf,'r') as f:
        for r in f.readlines():
            if r[0]=='>':
                name=r
                seq=''
            else:
                seq=seq+r.rstrip()
        seq_dict[name]=seq
    return seq_dict
    
def reverse_complement(seq):
    complementary_seq=''
    complement={'A':'T','T':'A','G':'C','C':'G'}
    for s in seq:
        complementary_seq=complement[s]+complementary_seq
    return complementary_seq

__author__ = "Swaminathan Kannan"

# # <center> Function compilation for algorithms in Week 1 </center>
def naive_comparison(pattern,text):
    occurances=[]
    for t in range(len(text)-len(pattern)+1):
        match=True
        for p in range(len(pattern)):
            if text[p+t]!=pattern[p]:
                match=False
                break
        if match==True:
            occurances.append(t) # record the index in t where p matches
    return list(set(occurances)) #also return whether a match was present or not
    
def naive_2mm(pattern, text):
    count_m=[]
    for t in range(len(text)-len(pattern)+1):
        match=True
        miscount=0
        for p in range(len(pattern)):
            if text[t+p]!=pattern[p]:
                miscount+=1
                if miscount>2:
                    match=False
                    break
        if match:
            count_m.append(t)
    return count_m

def read_q(text_q): #Read a fastq file (with quality scores)
    seqa=[]
    quala=[]
    with open(text_q,'r') as f:
        lines=f.readlines()
        for line in range(len(lines)):
            if lines[line][0]=='@':
                seqa.append(lines[line+1].rstrip())
                quala.append(lines[line+3].rstrip())
            if line==len(lines)-4:
                break
    return seqa,quala

def get_score(chara): # convert the ASCII value of a score to the actual Quality score
    score = ord(chara)-33
    return score
   
