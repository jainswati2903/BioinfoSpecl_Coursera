#!/usr/bin/env python
# coding: utf-8

from Week1 import *

def readGenome(file):
    '''function to read and return a genome as one string in fasta format'''
    genome =  ''
    f=open(file)
    for line in f:
        if line[0] == ">": #ignoring the name of the genome
            pass
        else:
            genome += line.strip()
    f.close()
    return genome

def MinimumSkew(Genome,max=False):
    
    '''function to return the positions in the genome that minimize the G-C skew'''
    
    len_g=len(Genome)
    skew=[0]*(len_g+1) # as skew is calculated from 0 - len_g (so total of len_g+1)
    min_skew = 0 # store the minimum skew
    min_skew_index = [] # store the minimum skew indices
    
    for i in range(1,len_g+1): #for all indices
        
        if Genome[i-1] == 'C':
            skew[i]=skew[i-1] - 1
        elif Genome[i-1] == 'G':
            skew[i]=skew[i-1] + 1
        else:
            skew[i]=skew[i-1]
            
        if not max: # if min skew needs to be returned
            if skew[i] < min_skew: # new minimum skew
                min_skew = skew[i]
                min_skew_index = []
                min_skew_index.append(i)
            elif skew[i] == min_skew: # if it is the same minimum
                min_skew_index.append(i)
        else: # if max skew need to be returned
            if skew[i] > min_skew: # new minimum skew
                min_skew = skew[i]
                min_skew_index = []
                min_skew_index.append(i)
            elif skew[i] == min_skew: # if it is the same minimum
                min_skew_index.append(i)
            
    return min_skew_index

def HammingDistance(p, q):
    
    '''function to calculate the hamming distance between two strings'''
    ham_dist = 0
    for index in range(len(p)):
        
        if p[index] != q[index]: #increase the distance if two chars are not equal
            ham_dist+=1
            
    return ham_dist

def Neighbors(Pattern,number):
    
    '''function to return all strings that mismatch the pattern at atmost number places'''
    '''this will be a recursive function'''
    '''better implementation would be to retrun the hamming distance of each as well'''
    
    neighbors=[]
    nucs=['A','C','G','T']
    
    #termination conditions
    if number == 0: # no mismatches allowed
        neighbors.append(Pattern)
        return neighbors
    
    if len(Pattern) == 1: #the pattern has only one letter
        return nucs # return all single nucleotides as number will always be atleast 1 if it comes here
    
    suffix = Pattern[1:]
    neighbors_suffix = Neighbors(suffix,number) # recursively calculate all neighbors of the suffix
    for n in neighbors_suffix:
        
        dist=HammingDistance(n,suffix)
        if dist < number: # it can tolerate more mismatches
            for i in nucs:
                neighbors.append(i+n) # append all
        elif dist == number: # no more mismatched needed
                neighbors.append(Pattern[0]+n) # add the original nuc at position 0 to the suffix neighbor
                
    return neighbors    
    
def FrequentWordsWithMismatches(Text, kmer_len, number, rc=True):
    
    '''function to return the kmer of length kmer_len that occurs the most times in the text with at most 
    number of mismatches'''
    '''will consider reverse compliments by default'''
    '''some neighbors for the RC can be the same as the ones for the original pattern, therefore maybe we can improve it
    by making the neighbors a set'''
    
    kmer_count=dict() # creating a dictionary to count the number of kmer occurances with at most number mismatches
    for index in range(len(Text)-kmer_len+1):# for every kmer
    
        cur_kmer=Text[index:index+kmer_len] 
        neighborhood = Neighbors(cur_kmer,number) # calculate all strings with atmost number mismatches from the current kmer
        if rc: # if reverse compliment also needs to be taken
            revComp=reverseCompliment(cur_kmer)
            neighbor_rc=Neighbors(revComp,number)
            neighborhood.extend(neighbor_rc)
    
        # every kmer will also be counted as an occurance of every pattern in its neighborhood 
        for n in neighborhood:
            if n in kmer_count: # if it already exists
                kmer_count[n]+=1
            else: # if not, create a new key
                kmer_count[n]=1
                
    # sort based on values
    sorted_kmer = sorted(kmer_count.items(), key=lambda kv: kv[1])
    freq_kmer=[]
    index = len(sorted_kmer)-1
    while sorted_kmer[index][1] == sorted_kmer[len(sorted_kmer)-1][1] and index>=0: # as long as the frequency of the kmer is equal to the maximum frequeny
        freq_kmer.append(sorted_kmer[index][0])
        index-=1
        
    return freq_kmer
