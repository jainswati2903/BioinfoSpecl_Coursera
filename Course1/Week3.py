#!/usr/bin/env python
# coding: utf-8

from Week1 import *
from Week2 import *
import numpy as np
import math

def MotifEnumeration(dna_strings,kmer_len,number):
    
    '''function to return all kmers that occur in all strings in dna_strings with atmost number mismatches'''
    kmers_in_all=dict() # dictionary to contain all number mismatch kmers for all kmers in the first string 
    num_strings=len(dna_strings)
    # the value v for each kmer will represent if it can be counted as being present in all dna_strings upto index v
    
    for str_index in range(num_strings):# for every string
        print(str_index)
        
        for index in range(len(dna_strings[str_index])-kmer_len+1): # for every kmer in the current string
            
            cur_kmer=dna_strings[str_index][index:index+kmer_len] # extracting the kmer
            neighborhood=Neighbors(cur_kmer,number) # getting the neighborhood of the current kmer
            
            for n_kmer in neighborhood: # for every neighbor of the kmer, it counts as being present in this string
                
                if str_index == 0: # this is the first string, so kmers need to be inserted into the dictionary
                    kmers_in_all[n_kmer]=0
                else: # for every other string
                    
                    if n_kmer in kmers_in_all: # if kmer is in the dictionary already, then only count keep it if it was counted as being present in all previous strings as well
                        
                        if kmers_in_all[n_kmer] == str_index - 1: # if it was counted as being present in all strings uptill now
                            kmers_in_all[n_kmer] = str_index # count it as being present in this string as well
    
    #iterate over the dictionary and keep the kmers that were counted as being present in all strings
    present=[]
    for kmer in kmers_in_all.keys():
        if kmers_in_all[kmer] == num_strings - 1:
            present.append(kmer)
            
    return present
              
def DistanceBetweenPatternAndStrings(pattern, dna):
    '''function to give the sum of distane between the pattern and every string in the list dna. The distance between a pattern
    and a string is the minimum hamming distance between pattern and kmers of the string'''
    
    kmer_len=len(pattern)
    total_dist = 0 # the sum of distances to be returned
    
    for string in dna: # for every string in the list dna
        
        min_dist = kmer_len # as this is the maximum hamming distance you can have
        for index in range(len(string)-kmer_len+1): # for every kmer in the string
            
            cur_kmer = string[index:index+kmer_len]
            cur_dist = HammingDistance(pattern,cur_kmer) # calculate the hamming distance between pattern and kmer
            
            if cur_dist < min_dist:
                min_dist = cur_dist #store the new minimum
        
        # after processing the whole string, add the min dist of this string to the total distance
        total_dist += min_dist
        
    return total_dist
            
def MedianString(dna, k,getlist=False):
    '''function that returns the kmer of length k that minimizes the distance between itself and all strings in the list dna'''
    
    min_dist=k*len(dna) + 1# this is the maximum the distance can be
    min_kmer_list = []
    
    for n in range(4**k): # for each possible kmer of length k
        
        kmer = NumberToPattern(n,k) # get the kmer corresponding to the number
        cur_dist = DistanceBetweenPatternAndStrings(kmer, dna) # distance of this kmer from the set of strings
        
        if cur_dist < min_dist: # update the minimum
            min_dist=cur_dist
            min_kmer_list = []
            min_kmer_list.append(kmer)
        elif cur_dist == min_dist and getlist:
            min_kmer_list.append(kmer)
            
    if not getlist:
        return min_kmer_list[0]
    else:
        return min_kmer_list
           
def ProfileMostProbableKmer(text, k, profile):
    '''function to return the most probable kmer of length k in the text based on the profile probability matrix'''
    
    max_prob = -1
    max_kmer = ''
    prob_dict={'A':profile[0], 'C':profile[1], 'G':profile[2], 'T':profile[3]} # converting the matrix into a dictionary for easy access
    
    for index in range(len(text)-k+1): #for every kmer in text
        
        cur_kmer = text[index:index+k]
        prob = 1
        for c in range(k):
            prob *= prob_dict[cur_kmer[c]][c] # multiplying the probability of each char at each position
            
        if prob > max_prob:
            max_prob=prob
            max_kmer = cur_kmer
            
    return max_kmer

def MotifProfile(Motifs,pseudocount=False):
    
    '''function to return the profile matrix of a given set of Motifs'''
    
    motif_len=len(Motifs[0]) # the length of each motif
    number_of_motifs=len(Motifs)
    profile = [[0]*motif_len for i in range(4)] # initializing a profile matrix
    
    #counting the numbers of each nucleotide in each column and store it in the profile
    
    for col in range(motif_len):
        
        for row in range(number_of_motifs):
            
            if Motifs[row][col] == 'A':
                profile[0][col]+=1
            elif Motifs[row][col] == 'C':
                profile[1][col]+=1
            elif Motifs[row][col] == 'G':
                profile[2][col]+=1
            elif Motifs[row][col] == 'T':
                profile[3][col]+=1       
        
    # dividing every entry by the total number of motifs to get the probability
    if not pseudocount:
        profile = np.divide(profile,number_of_motifs)
    else: # need to add pseudocount
        profile= np.add(profile,1)
        profile= np.divide(profile,(number_of_motifs+4))
        
    return profile

def MotifEntropy(Motifs):
    
    '''function to calculate the entropy of the motif matrix'''
    profile = MotifProfile(Motifs) # calculate the profile matrix for the motifs
    
    entropy = 0
    
    for row in range(len(profile)):
        
        for col in range(len(profile[0])):
            
            if profile[row][col] != 0: # as entropy is not defined for probability 0
                
                entropy += (-1)*(profile[row][col])*math.log(profile[row][col],2)
                
    return entropy
    
def MotifScore(Motifs,profile=0):
    
    '''function to calculate the score of the motifs, profile to be calculated if it is not provided'''
    
    if profile.all() == 0:
        profile=MotifProfile(Motifs) # calculate the profile if not provided
        
    nuc_dict={0:'A', 1:'C', 2:'G', 3:'T'}  
    profile_trans = profile.transpose() # transposing the profile for ease of calculating the max of each row (prev the col)
    consensus_string='' # the consensus string of the list of motifs
    score = 0
    motif_len=len(Motifs[0])
    
    for i in range(motif_len):
        
        temp_list=list(profile_trans[i]) # converting it into a list
        max_index = temp_list.index(max(temp_list))
        consensus_string+=nuc_dict[max_index]
        
    for m in Motifs: # for every motif
        
        score += HammingDistance(m,consensus_string) #sum up the hamming distance between each motif and the consensus string
        
    return score

def GreedyMotifSearch(DNA, motif_len, num_strings, pseudocount=False):
    
    '''function to return list of motifs of length motif_len, one from each string in DNA that minimize the score(motifs) '''
    BestMotifs=[] # keeping track of best motifs
    min_score = (num_strings - (num_strings/4 + 1))*motif_len + 1# the maximum score that can be for a motif matrix
    
    string_len = len(DNA[0]) # length of each string in the DNA
    for index in range(string_len-motif_len+1): # for each kmer in the first DNA string
        
        CurMotifs = []
        first_motif = DNA[0][index:index+motif_len] # the first motif from the first string
        CurMotifs.append(first_motif)
        CurProfile = MotifProfile(CurMotifs,pseudocount)
        
        for str_index in range(1,num_strings): # for every other string
            
            most_prob_motif = ProfileMostProbableKmer(DNA[str_index],motif_len,CurProfile) # get the most probable motif 
            # from the current string based on the CurProfile
            #adding the motif to the collection and updating the profile
            CurMotifs.append(most_prob_motif)
            CurProfile = MotifProfile(CurMotifs,pseudocount)
        
        #check and update the minimum score
        cur_score = MotifScore(CurMotifs,CurProfile)
        if cur_score < min_score: # if we have a new minimum
            min_score = cur_score
            BestMotifs = CurMotifs
           
    return BestMotifs
