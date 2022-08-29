import numpy as np

def edit_distance(a,b): #Calculate edit distance between string 'a' and string 'b'
    if len(a)==0:
        return len(b)
    if len(b)==0:
        return len(a)
    else:
        x=a[-1]
        y=b[-1]
        alpha=a[:-1]
        beta=b[:-1]
        dist_last=0 if x==y else 1
        dist=min(edit_distance(alpha,beta)+dist_last,edit_distance(alpha,b)+1,edit_distance(beta,a)+1)
    return dist        
    
    
def edit_recursive(x,y): #Calculate edit distance using dynamic programming
    matrix=np.zeros((len(x)+1,len(y)+1))
    matrix[0,0]=0 #Setting empty coordinates to 0
    for i in range(1,len(x)+1):
        matrix[i,0]=matrix[i-1,0]+1
    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            dist_last=0 if x[i-1]==y[j-1] else 1
            matrix[i,j]=min(matrix[i-1,j-1]+dist_last,matrix[i-1,j]+1, matrix[i,j-1]+1)
    return min(matrix[-1])
    
    
    
def global_alignment(pattern, text): #Expanding recursive edit distance alignment for finding patterns in text
    matrix=np.zeros((len(pattern)+1, len(text)+1))
    matrix[0,:]=0 #keeping the first row pertaining to epsilon completely 0 because we don't know where in the text, the pattern matching start
    for i in range(1,len(pattern)+1):
        matrix[i,0]=matrix[i-1,0]+1
    for i in range(1,len(pattern)+1):
        for j in range(1,len(text)+1):
            dist_last=0 if pattern[i-1]==text[j-1] else 1
            matrix[i,j]=min(matrix[i-1,j-1]+dist_last,matrix[i-1,j]+1, matrix[i,j-1]+1)
    return min(matrix[-1,:]),matrix


   
def find_index(matrix,pattern,text): #matrix = completed DP matrix; find the offset of pattern with respect to text given the the global alignment matrix
    vals=min(matrix[-1,:]) #Find the minimum value in the last row
    x_max=np.squeeze(np.where(matrix[-1,:]==vals)) #Find the column index where 'vals' is located. Basically, anything to the right of x_max in the matrix is ignored in this operation
    t1=text[:x_max]
    x=x_max
    y=len(pattern) # Basically our area of focus is the slice of the matrix matrix[:,:x_max]
    while x>=0 and y>0:
        diff=0 if t1[x-1]==pattern[y-1] else 1
        check1=matrix[y-1,x-1]-diff
        check2=matrix[y-1,x]-1
        check3=matrix[y,x-1]-1
        r=np.argmin([check1,check2,check3])
        if r==0:
            x-=1
            y-=1
        if r==1:
            y-=1
        if r==2:
            x-=1   
    return text[x-1:],x-1

#**************************************************************************
#Sample penalty matrix for the global alignment algorithm using variable penalties

bases=['A','C','G','T','s'] #This is the order in which we are making the penalty table. 's' stands for a skip
penalty=np.empty((len(bases),len(bases)))
penalty[0]=np.array([0,4,2,4,8])
penalty[1]=np.array([4,0,4,2,8])
penalty[2]=np.array([2,4,0,4,8])
penalty[3]=np.array([4,2,4,0,8])
penalty[4]=np.array([8,8,8,8,8])
penalty=np.array(penalty).astype('int32')
print(penalty)
#**************************************************************************

def global_alignment_penalty(pattern, text,penalty):
    matrix=np.zeros((len(pattern)+1, len(text)+1))
    for i in range(1,len(pattern)+1):
        matrix[i,0]=matrix[i-1,0]+penalty[bases.index(pattern[i-1]),-1] #populating the first row .This will be prev cell + penalty for a skip (since we are mapping a base to an empty string)
    for i in range(1,len(text)+1):
        matrix[0,i]=matrix[0,i-1]+penalty[-1,bases.index(text[i-1])] #populating the first column .This will be prev cell + penalty for skipping the previous letter (since we are mapping a base to an empty string)
    for i in range(1,len(pattern)+1):
        for j in range(1,len(text)+1):
            dist_last=penalty[bases.index(pattern[i-1]),bases.index(text[j-1])] #substitution
            dist_hor=penalty[-1,bases.index(text[j-1])]  # So if the previous letter in the
            dist_ver=penalty[bases.index(pattern[i-1]),-1]
            matrix[i,j]=min(matrix[i-1,j-1]+dist_last,matrix[i-1,j]+dist_ver, matrix[i,j-1]+dist_hor)
    return matrix[-1,-1],matrix
    


def overlap(a,b,minlength=3): #Find if prefix of 'b' is in suffix of 'a'
    start=0
    while True:
        start = a.find(b[:minlength],start)
        if start==-1:
            return 0
        if b.startswith(a[start:]): ###Only need the space where b starts
            return len(a)-start
        start+=1
        
from itertools import permutations #<---- new library I learnt this time :)
def naive_overlap(reads_list, min_length=3):
    overlaps_dict={}
    for a,b in permutations(reads_list, 2):
        print(a,b)
        overlap_length=overlap(a,b,min_length)
        if overlap_length>0:
            overlaps_dict[(a,b)]=overlap_length
    return overlaps_dict