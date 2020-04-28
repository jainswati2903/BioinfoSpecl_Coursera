#!/usr/bin/env python
# coding: utf-8

from Week1 import *
from Week2 import *
from Week3 import *
import random
import numpy as np
import copy

def RandomizedMotifSearch(DNA, motif_len, num_strings, num_steps=1000):
    
    '''function to return the set of motifs that minimize the score(motifs) after num_steps number of random searches'''
    '''one random search basically starts with one random kmer from eacg string, then construct the profile and select another 
    set of kmers that mazimize the profile for each string, and then continue the process untill the score keeps decreasing'''
    
    overall_lowest_score = 1000000
    overall_lowest_motifs = []
    num_kmers_in_string = len(DNA[0])-motif_len+1
    
    for i in range(num_steps): # for each step
    
        #print("Step: %d"%i)
        step_lowest_score = 1000000
        step_lowest_motifs = []
        cur_motifs = []
        cur_score = 1000000
        
        # select a random set of kmers from each string
        for string in DNA: # for each string in the DNA
            index = random.randrange(num_kmers_in_string)
            cur_motifs.append(string[index:index+motif_len])
            
        cur_profile=MotifProfile(cur_motifs,True) # calculate the profile and the score based on the initial random selection
        cur_score=MotifScore(cur_motifs,cur_profile)
        #step_lowest_score=cur_score # designate this as the lowest score in this step
        #step_lowest_motifs=cur_motifs
        
        while step_lowest_score > cur_score: # while the new score is less than the best score
            
            step_lowest_score = cur_score # store the previous score and motif
            step_lowest_motifs = cur_motifs
            
            cur_motifs=[]
            for string in DNA: # get the most probable motifs from each string based on the current profile
                cur_motifs.append(ProfileMostProbableKmer(string, motif_len, cur_profile))
            cur_profile = MotifProfile(cur_motifs,True) # calculate the new profile and new score
            cur_score = MotifScore(cur_motifs,cur_profile)
            
        # when this step is done, the lowest from this step is stored in the step variables
        if step_lowest_score < overall_lowest_score: # if we have a new minimum from this step
            overall_lowest_score = step_lowest_score
            overall_lowest_motifs = step_lowest_motifs
            
        #print(overall_lowest_score)
            
    return overall_lowest_motifs

def BiasedRandom(Probs):
    
    '''function to return a random number between 0 and len(Probs), each with a specified probability'''
    
    sum_p = sum(Probs) # sum of all probabilities
    Probs_sum1 = np.divide(Probs,sum_p) # now all probabilities sum upto 1
    
    #creating a cumulative array of probabilities, so that each index in the original Probs array is represented by the
    #range of cum_probs[i-1] and cum_probs[i]
    cum_prob = []
    cum_prob.append(Probs_sum1[0])
    for i in range(1,len(Probs_sum1)):
        cum_prob.append(cum_prob[i-1]+Probs_sum1[i])
        
    #generate a random number between 0 and 1
    ran_num=random.random()
    #check which range cum_probs[i-1] and cum_probs[i] this random number falls, and return that i as the random number generated
    for i in range(len(cum_prob)):
        if i == 0 and ran_num <= cum_prob[0]:
            return i
        elif ran_num > cum_prob[i-1] and ran_num <= cum_prob[i]:
            return i
        
def ProfileRandomKmer(Motifs,index,text):
    
    '''function to return the a new motif from text that replaces the Motifs[index] in the motif matrix
    the new motif is randomly seleted from kmers in text, with the probability given by profile generated 
    by all Motifs except Motifs[index]'''
    
    kmer_len = len(Motifs[0]) # length of the kmer
    
    Temp_motifs = copy.deepcopy(Motifs) # as we have to change it
    Temp_motifs.pop(index) # remove the motif in index
    Temp_profile=MotifProfile(Temp_motifs,True)
    
    prob_dict={'A':Temp_profile[0], 'C':Temp_profile[1], 'G':Temp_profile[2], 'T':Temp_profile[3]} # converting the matrix into a dictionary for easy access
    prob_kmers = []
    #calculate the probability of the each kmer in text based on the Temp_profile
    for k_i in range(len(text)-kmer_len+1): # for every kmer in text
        
        cur_kmer = text[k_i:k_i+kmer_len] # current kmer
        prob=1
        for c in range(kmer_len):
            prob *= prob_dict[cur_kmer[c]][c] # multiplying the probability of each char at each position
        prob_kmers.append(prob)
        
    chosen_index = BiasedRandom(prob_kmers) # randomly choose the kmer using the biased function based on kmer probability
    return text[chosen_index:chosen_index+kmer_len] # return the actual kmer
    
def GibbsSampler(DNA,motif_len,num_strings,iterations,num_steps=20):
    
    '''function to return the motif array with min score calculated using GibbsSampling approach with inner loops iterations time
    and num_steps random starts'''
    
    overall_lowest_score = 1000000
    overall_lowest_motifs = []
    num_kmers_in_string = len(DNA[0])-motif_len+1
    
    for i in range(num_steps): # for each step
    
        print("Step: %d"%i)
        step_lowest_score = 1000000
        step_lowest_motifs = []
        cur_motifs = []
        cur_score = 1000000
        
        # select a random set of kmers from each string
        for string in DNA: # for each string in the DNA
            index = random.randrange(num_kmers_in_string)
            cur_motifs.append(string[index:index+motif_len])
            
        cur_profile=MotifProfile(cur_motifs,True) # calculate the profile and the score based on the initial random selection
        cur_score=MotifScore(cur_motifs,cur_profile)
        step_lowest_score=cur_score # designate this as the lowest score in this step
        step_lowest_motifs=cur_motifs   
        
        for j in range(iterations): # each iteration
            
            motif_i = random.randrange(num_strings) # select one string at random
            new_motif = ProfileRandomKmer(cur_motifs,motif_i,DNA[motif_i]) # will return the new motif that becomes the new motif for string motif_i
            
            cur_motifs[motif_i] = new_motif
            cur_profile=MotifProfile(cur_motifs,True)
            cur_score=MotifScore(cur_motifs,cur_profile)
            
            if cur_score < step_lowest_score:
                step_lowest_score = cur_score
                step_lowest_motif = copy.deepcopy(cur_motifs) # as we are changing the cur_motifs array every time
            
        # when this step is done, the lowest from this step is stored in the step variables
        if step_lowest_score < overall_lowest_score: # if we have a new minimum from this step
            overall_lowest_score = step_lowest_score
            overall_lowest_motifs = step_lowest_motifs
            
        #print(overall_lowest_score)
            
    return overall_lowest_motifs
