from functions_week3 import overlap
from functions import read_q
import itertools
from itertools import permutations
from collections import defaultdict


def scs(str_list,k):
    reference_string=''.join(str_list)
    ss_set=[]
    count=0
    for strings in permutations(str_list):
        count_overlap=0
        final_str=strings[0]  #Basically, at the end, we are appending the overlap value of the second string to the final_str. Hence, if any 'strings' is not considered as the second string in overlap() e.g. first string in the permutation set, it will not be updated in the final string 
        for i in range(len(str_list)-k):
            oval=overlap(strings[i],strings[i+1],k)
            final_str+=strings[i+1][oval:]
            count_overlap+=1
        ss_set.append(final_str)
        if len(final_str)<len(reference_string):
            reference_string=final_str
        count+=1
    count_incidence=0
    for ss in ss_set:
        if len(ss)==len(reference_string):
            count_incidence+=1
    return reference_string,count_incidence,ss_set
    
def get_maximal_overlap(str_list,k): #implementation based on lecture. Need to try dictionary approach
    best_olen=0
    best_a,best_b='',''
    for a,b in permutations(str_list,2):
        oval=overlap(a,b,k)
        if oval>best_olen:
            best_a, best_b, best_olen= a,b,oval
    return best_a, best_b, best_olen
    
def greedy_scs(str_list,k):
    a,b,olen=get_maximal_overlap(str_list,k)
    while olen>0:
        str_list.remove(a)
        str_list.remove(b)
        str_list.append(a+b[olen:])
        a,b,olen=get_maximal_overlap(str_list,k)
    if len(str_list)>1:
        return str_list
    else:
        return ''.join(str_list)
        
def visualize_de_b(st,k):
    nodes,edges=get_db(st,k) #Need to implement this function separately
    dot_str = 'digraph "DeBruijn graph" {\n'
    for node in nodes:
        dot_str+=' % s [label="%s"] ; \n' % (node, node)
    for src,dst in edges:
        dot_str+=' %s->%s; \n' % (src,dst)
    return dot_str + '}\n' 
    
    
def get_db(st,k):
    edges=[]
    nodes=set()
    for i in range(len(st)-k+1):
        edges.append((st[i:i+k-1],st[i+1:i+k]))
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1:i+k])
    return nodes,edges
    
def fast_overlap(reads,kmer=6,minlength=3):
    count=0
    max_olen=0
    max_a,max_b='',''
    dict_kmers=defaultdict(set)
    for read in reads:
        for t in range(len(read)-kmer+1):
            dict_kmers[read[t:t+kmer]].add(read)    #This is a dictionary with key=kmer and value=entire read
    overlap_dict=defaultdict(set)
    for read in reads:
        suffs=read[-kmer:] #Get the last kmer/ suffix of length k from the read
        read_suffs=dict_kmers[suffs] #Get all the reads that have this suffix somewhere in their list
        for suff in read_suffs: 
            if read!=suff:#If the original read is not equal to the current read selected:
                olen=overlap(read,suff,minlength)
                if olen>max_olen:
                    max_a=read
                    max_b=suff
                    max_olen=olen
    return max_a, max_b, max_olen
    
def fast_greedy(str_list,k=3):
    a,b,olen=fast_overlap(str_list,k)
    i=0
    while olen>0:
        print('\ri {}'.format(i), end="")
        str_list.remove(a)
        str_list.remove(b)
        str_list.append(a+b[olen:])
        a,b,olen=fast_overlap(str_list,k)
        i+=1
    if len(str_list)==1:
        return str_list
    else:
        return ''.join(str_list)