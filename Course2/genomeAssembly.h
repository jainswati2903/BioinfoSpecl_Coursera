//
//  genomeAssembly.h
//  Course2_Assembly
//
//  Created by Swati Jain on 4/23/20.
//  Copyright © 2020 Swati Jain. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <cmath>
#include <queue>

using namespace std;

/*** Utility Functions***/
void printAdjList(map<string, vector<string> >, ostream &); // function to print the adjacency list given a map of nodes and connections
map<string, vector<string> > readAdjList(istream &); // read the adjacency list from the input stream
void printPath(vector<string>,ostream &); // print the path given by the nodes in the vector
vector<string>GenBinStrings(int,string,vector<string>);// function to generate all binary strings of given length
string transcribe(string); // function to transribe the DNA string into an RNA string
string translate(string); // function to translate the RNA string into a protein string
string reverseCompliment(string); // function to return the reverse compliment for the given DNA string

/*** Main Assembly Functions***/
vector<string> Composition_substrings(ifstream &,ofstream &,bool,bool); // generate all substrings of given length and given sequence
string PathProblem(ifstream &,ofstream &); // generate a full sequence from kmers in order given
string PathToGenome(vector<string>); // generate a full sequence from kmers in order given
void Overlap_graph(ifstream &,ofstream &); // generate an overlap graph from kmers given
map<string, vector<string> > DeBruijn_graph(ifstream &, ofstream &, bool); // generate a debruijn graph for a given text and kmer-length or set of kmers but to be read from input
map<string, vector<string> > DeBruijn_graph(vector<string>,bool); // generate a debruijn graph for given kmers passed as arguments
vector<string> EulerianCycle(map<string, vector<string> >); // generate a eulerian cycle for the given graph
vector<string> EulerianPath(map<string, vector<string> >); // generate a eulerian path for the given graph
string kUniversal(int); // function to generate the k universal string for given kmer length
string StringFromGapPattern(vector<string>,int,int); // function to return a string that is spelled by the path of paired reads
vector< vector<string> > MaximalNonBranchPaths(map<string, vector<string> >); // function to find all maximal non branching paths for a given adjacency list

/*** Main Antibiotics Function ***/
vector<string> encodeProtein(string, string); // find all substrings of given DNA that encode for the given protein
vector<int> LinearSpectrum(string, bool = false, bool = false); // function to return the linear spectrum of a peptide string, i.e. mass of all subpeptides of the peptide string
bool SpectrumConsistent(vector<int>, vector<int>); // to see if spectrum 1 is consistent with spectrum 2
vector<string> CyclopeptideSequencing(vector<int>); // function to return peptide sequences (in masses) that matches the given spectrum
int SpecScore(string,vector<int>,bool = false, bool = false); // function to return the score of the cyclic peptide with respect to the given spectrum
int peptideMass(string,bool = false); // function to return the mass of the given peptide
vector<string> LeaderBoardCycloPepSeq(vector<int>, int, int = 0); // function to return the peptide that has the highest score with respect to the spectrum using leader board also where at each step you only keep given number of candidates
map<int,int> SpecConvolution(vector<int>); // return the convolution of the spectrum, with each element counted as many times as it occurs in the convolution
