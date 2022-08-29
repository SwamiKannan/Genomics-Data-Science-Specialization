#!/usr/bin/env python
# coding: utf-8
__author__ = "Swaminathan Kannan"

# # <center> Function compilation for algorithms in Week 2 </center>


import bm_preproc
from bm_preproc import BoyerMoore
import functions
from functions import readf
import matplotlib.pyplot as plt
import bisect
import numpy as np
from scipy.stats import norm



alphabet='ATCG' #dictionary of letters for all BoyerMoore's object


# ### Processing of pattern 'p' and not processing text 't' (Online processing)


def boyer_moore(pattern,text,bm_obj):
    occurances=[]
    for t in range(len(text)-len(pattern)+1):
        mismatch=False
        for p in reversed(range(len(pattern))):
            if text[t+p]!=pattern[p]:
                shift_bad_char=bm_obj.bad_character_rule(p, pattern[p])
                shift_good_suffix=bm_obj.good_suffix_rule(p)
                shift=max(shift_bad_char,shift_good_suffix,1)
                mismatch=True
                break
        if not mismatch:
            occurances.append(t)
            shift=max(1,bm_obj.match_skip())
        t+=shift
    return list(set(occurances))


# ### Processing of processing text 't' (Offline processing)

# #### Indexing


#Understanding the bisect function. Basically, bisect_left tells us the left most position in text where we can insert a specific string 
import bisect
def find_insertion(text, pattern):
    check=bisect.bisect_left(sorted(text),pattern) #Sorted is required to ensure that the correct position is chosen. Analogy to a dictionary: if a dictionary wasn't sorted, we wouldn't be able to find the word
    return check,sorted(text)
#We use pattern[0] because later, we need to use this to find the position in an index which is an array of (kmer,position of that kmer in text)


# ### K-mer lookup in the text index

# <ol><li> We create an index of all combination of k substrings (integer k < len(pattern)) of the text and create a lookup table of the same (containing the k-mer and its index/ offset. ----> <b>create_index() </b></li>
# <li> We then create similar k-mers of the pattern and look up the index of the pattern k-mer to see if it exists. ----> <b>find_all_indexes()</b></li>
# <li> Once we confirm that the kmer exists, we use the same to VERIFY if the rest of the pattern is also present at the offset ---> <b>verify_pattern()</b></li></ol>


def create_index(text, kmer): #Create the index of all kmers from text 't'
    index={}
    for i in range(len(text)-kmer+1):
        if text[i:i+kmer] in index:
            index[text[i:i+kmer]]+=[i]
        else:
            index[text[i:i+kmer]]=[i]
    return index


# find_insertion_index() only provides us with the first insertion point, if any, for the pattern <b> k-mer </b>. We still need to perform verification for the rest of the pattern

def find_all_indexes(index, pattern,kmer): #does not use bisect
    sample=pattern[:kmer]
    if sample in index.keys():
        return index[sample]


def verify_pattern(index, text, pattern,kmer):
    hits=find_all_indexes(index, pattern, kmer)
    print('hits:',hits)
    matches=True
    occurances=[]
    for h in hits:
        for i in range(len(pattern)-kmer+1):
            if text[h+i]!=pattern[i]:
                matches=False
                break
        if matches:
            occurances.append(h)
    return occurances


# ### Approximate Matching

# ### Sub-sequence Matching

# Subsequence matching would pretty much follow the K-mer lookup process except we need an additional function that creates the logic for the substrings

def create_subsequence(text, leng,interval):#creates subsequences of length "leng" whose initial chars are "interval" values apart
    subsequence=[]
    for i in range(len(text)):
        str1=''
        for r in range(i,len(text),interval):
            str1+=text[r]
            if len(str1)>=leng:
                break
        if len(str1)==leng:
            subsequence.append((str1,i))
        i+=1
    return subsequence

def create_subsequence_index(text, leng, interval):
    index={}
    subsequences=create_subsequence(text, leng, interval)
    subsequence_key,subsequence_value=[i for i in zip(*subsequences)]
    for i in range(len(subsequence_key)):
        if subsequence_key[i] in index:
            index[subsequence_key[i]].append(subsequence_value[i])
        else:
            index[subsequence_key[i]]=[subsequence_value[i]]
    return index,subsequences        

def find_subsequence_index(index, pattern, leng, interval):
    occurances=[]
    subsequences_pattern=create_subsequence(pattern, leng, interval)
    subsequence_key_pattern,_=[i for i in zip(*subsequences_pattern)]
    print('pattern_sub',subsequence_key_pattern)
    count=0
    for sub in subsequence_key_pattern:
        if sub in index[0]:
            occurances.append(index[0][sub])
        count+=1
    return occurances,len(subsequence_key_pattern),count


def verify_subsequence_pattern(index, text,pattern, leng, interval,mismatch_threshold=1):
    hits, pattern_key_length, total_key_matches=find_subsequence_index(index, pattern, leng, interval)
    print('hits',hits)
    matches=True
    occurances=[]
    index_hits=0
    for h in hits:
        h=list(h)
        for i in h:
            index_hits+=1
            mismatch_count=0
            for j in range(interval+1):
                if text[i+j]!=pattern[j]:
                    mismatch_count+=1
                    if mismatch_count>mismatch_threshold:
                        matches=False
                        break
        if matches:
            occurances.append(h)
    return np.squeeze(occurances).tolist(),index_hits


# ### Pigeonhole

# Basically, we need to divide the pattern into (n+1) non-overlapping but mutually exclusive lists where n is the maximum number of mismatches allowed. Then, we check every one of these lists to the text and check if atleast one of these lists are an exact match (probably using Boyle Moore)

# ##### Use BoyerMoore naive lookup

def get_BM_model(p,alphabet='ATGC'):
    return BoyerMoore(p,alphabet)


def get_patterns(p, mismatches):
    batch_number=mismatches+1
    batch_size=int(round(len(p)//batch_number))
    print(batch_size)
    batches=[]
    for i in range(0,len(p),batch_size):
        batches.append(p[i:i+batch_size])
    return batches
                       

def compare_batches_bm(batches, text, pattern, mismatch_th=2):
    complete_matches=set()
    batch_length=len(batches[0])
    for i,batch in enumerate(batches):
        bm_obj=BoyerMoore(batch,alphabet='ATCG')
        occurance_batch=boyer_moore(batch,text,bm_obj)
        for m in occurance_batch:
            if m<i*batch_length or m-(i*batch_length)+len(pattern)>len(text):                    
                continue
            mismatches=0
            for j in range(0,i*batch_length):
                if text[m-i*batch_length+j]!=pattern[j]:
                    mismatches+=1
                    if mismatches>mismatch_th:
                        break
            for j in range(min((i+1)*batch_length,len(p)),len(p)):
                if text[m-i*batch_length+j]!=pattern[j]:
                    mismatches+=1
                    if mismatches>mismatch_th:
                        break
            if mismatches<=mismatch_th:
                complete_matches.append(m-i*batch_length)
    return list(complete_matches)      


# ##### Use Index Matching

def create_index1(text, kmer): #Create the index of all kmers from text 't' #### This is just a copy-paste from above for easy readability
    index={}
    for i in range(len(text)-kmer+1):
        if text[i:i+kmer] in index:
            index[text[i:i+kmer]]+=[i]
        else:
            index[text[i:i+kmer]]=[i]
    return index

def compare_batches_index(batches, text, pattern, index,kmer, mismatch_th=2):
    complete_matches=[]
    batch_length=len(batches[0])
    for i,batch in enumerate(batches):
        occurance_batch=verify_pattern(index, text, batch,kmer)
        for m in occurance_batch:
            if m<i*batch_length or m-(i*batch_length)+len(pattern)>len(text):                    
                continue
            mismatches=0
            for j in range(0,i*batch_length):
                if text[m-i*batch_length+j]!=pattern[j]:
                    mismatches+=1
                    if mismatches>mismatch_th:
                        break
            for j in range(min((i+1)*batch_length,len(p)),len(p)):
                if text[m-i*batch_length+j]!=pattern[j]:
                    mismatches+=1
                    if mismatches>mismatch_th:
                        break
            if mismatches<=mismatch_th:
                complete_matches.append(m-i*batch_length)
    return complete_matches 

