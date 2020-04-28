#!/usr/bin/env python
# coding: utf-8

def naive_match(pattern,text,st_aw=False,mismatch=0):
    '''function to return the occurances of pattern in given text, 
    last argument should be passed as True for searching for the rev compliment of the pattern as well'''
    occurances = []
    num_aligns = 0 # counting number of 
    num_char_comp = 0
    len_p = len(pattern)
    len_t = len(text)
    
    rev_same=False
    if st_aw: # calculate the reverse compliment if argument is true
        revP = reverseCompliment(pattern)
        if pattern == revP: # reverse compliment is the same
            rev_same = True 
        
    for i in range(len_t-len_p+1): # iterate over all possible indices
        notMatch = 0
        num_aligns += 1
        
        for j in range(len_p): # for every character in current possible alignment
            num_char_comp += 1
            if pattern[j] != text[i+j]: # if any character not equal, set match False and take a break
                notMatch+=1
                if notMatch > mismatch: # if number of mismtaches are more than what we allow
                    break
        if notMatch <= mismatch: # if all characters are same
            occurances.append(i)

        if st_aw and not rev_same: # if the reverse compliment also needs to be checked
            notMatch = 0
            for j in range(len_p): # for every character in current possible alignment
                if revP[j] != text[i+j]: # if any character not equal, set match False and take a break
                    notMatch+=1
                    if notMatch > mismatch: # if number of mismtaches are more than what we allow
                        break
            if notMatch <= mismatch: # if all characters are same
                occurances.append(i)
            
    return occurances, num_aligns, num_char_comp


def FrequentWords(text,k):
    
    '''function to return the most frequent kmers in the text'''
    '''This can be implemented using the PatternToNumber function as well'''
    kmer_count=dict() # keep count of kmers
    
    for index in range(len(text)-k+1): # for every kmer
        
        kmer=text[index:index+k] # extract the kmer
        if kmer in kmer_count.keys(): #increment the count associated with that kmer
            kmer_count[kmer]+=1
        else:
            kmer_count[kmer]=1
            
    max_count = 0
    max_kmers = ''
    
    for key in kmer_count.keys():
        
        if kmer_count[key] > max_count: # we have a new most frequent kmer
            max_count=kmer_count[key]
            max_kmers=''
            max_kmers=key
        elif kmer_count[key] == max_count: # same frequent kmer
            max_kmers = max_kmers + ' ' + key
            
    print(max_kmers) # print the concatenated list of all most frequent kmers
            

def PatternToNumber(Pattern):
    
    '''function to convert a DNA pattern to a number - can be used for indexing kmers'''
    
    sym_to_num={'A': 0, 'C': 1, 'G': 2, 'T': 3} # each symbol representing the number
    number = 0
    length = len(Pattern)
    
    for i in range(length-1,-1,-1): #iterating the pattern backwards
        
        number += (4**(length-i-1))*(sym_to_num[Pattern[i]]) # the symbol needs to be multiplied by the power of 4 
                   # that power increases in magnitude from right to left like in numbers                   
    return number 

def NumberToPattern(number,length):
    
    '''function to return the pattern of length corresponding to the number'''
    
    num_to_sym = {0:'A', 1:'C', 2:'G', 3:'T'}
    pattern = ''
    
    while number != 0:
        rem = number % 4 # getting the remainder, this is the next symbol from right to left in the pattern
        pattern = num_to_sym[rem] + pattern 
        number = number // 4 # getting the floor quotient
        
    if len(pattern) < length: # less symbols, need to concatenate A's at the begining of the pattern
        pattern = 'A'*(length-len(pattern)) + pattern
        
    return pattern
    
def ComputingFrequencies(text,len_k):
    
    '''function to return a frequency array for all occurances of len_k kmers in the text
    the indexing of the frequency array in the PatternToNumber of the kmer'''
    
    freq_array = [0]*(4**len_k) # frequency array equal to the number of kmers of size len_k
    
    for i in range(len(text)-len_k+1): # for each kmer in text
        
        kmer = text[i:i+len_k] # extract the kmer
        num = PatternToNumber(kmer) # get the index
        freq_array[num]+=1 # increment the frequency
        
    return freq_array

def reverseCompliment(seq):
    '''function to return the reverse complement of given sequence'''
    
    compliment = {'A':'T','T':'A','C':'G','G':'C','N':'N'} # dictionary of base compliments
    revComp = ''
    
    for c in seq: # for each character in the seqeuence
        revComp = compliment[c] + revComp
        
    return revComp

def ClumpFinding(genome,klen,win_len,number):
    
    '''function to return all kmers of length klen that occur at least number times in the win_length window of the genome'''
    clumps_kmers=set() #set of distnct kmers that satisfy the condition
    
    #calculating the freq array of the first window of the genome
    freq_array = ComputingFrequencies(genome[:win_len],klen)
    #adding the kmers that have at least number occurances
    for i in range(len(freq_array)):
        if freq_array[i] >= number:
            clumps_kmers.add(NumberToPattern(i,klen))
            
    #iterating over all other windows of the genome and updating the frequency array and clump_kmers
    for offset in range(1,len(genome)-win_len+1): # all windows of win_len starting from offset 1
        #print(offset)
        win_seq = genome[offset:offset+win_len] # current window
        first_kmer = win_seq[:klen]
        last_kmer = win_seq[win_len-klen:]
        
        #decreasing the freq of the first kmer and increasing the freq of the last kmer
        first_num = PatternToNumber(first_kmer)
        last_num = PatternToNumber(last_kmer)
        freq_array[first_num]-=1
        freq_array[last_num]+=1
        
        #checking if the last kmer now satisfies the condition - as only the last freq can increase for the last kmer
        if freq_array[last_num] >= number:
            clumps_kmers.add(last_kmer) # adding the kmer if it satisfies the condition
            
    return list(clumps_kmers)

