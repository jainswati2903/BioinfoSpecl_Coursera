//
//  main.cpp
//  Course2_Assembly
//
//  Created by Swati Jain on 4/23/20.
//  Copyright © 2020 Swati Jain. All rights reserved.
//

#include "genomeAssembly.h"

//uncomment specific codes segments to run specific functions
int main(int argc, char ** argv){
    
    string input_file,output_file;
    ifstream input;
    ofstream output;
    
    cin >> input_file; // take the input and output file as
    cin >> output_file;
    
    input.open(input_file);
    output.open(output_file);
    
	//running the Composition_substrings function
	/*vector<string> substrings;
	substrings = Composition_substrings(input,output);*/

	//running the PathProblem function
    /*string text;
    text=PathProblem(input,output);*/
    
    //running the Overlap_graph function
    //Overlap_graph(input,output);
    
    //running the DeBruijn_graph function
    /*map<string, vector<string> > adjacency_list;
    adjacency_list=DeBruijn_graph(input,output,false);
    printAdjList(adjacency_list,output);*/
    
    //running the EulerianCycle function
    /*vector<string> cycle;
    map<string, vector<string> > adj_list;
    adj_list = readAdjList(input);
    cycle=EulerianCycle(adj_list);
    printPath(cycle,output);*/
    
    //running the EulerianPath function
    /*vector<string> path;
    map<string, vector<string> > adj_list;
    adj_list = readAdjList(input);
    path=EulerianPath(adj_list);
    printPath(path,output);*/
    
    //runnning the string reconstruction workflow
    /*int kmer_len = 0;
    vector<string> kmers, path;
    string line,sequence;
    map<string, vector<string> > adj_list;
    input >> kmer_len;
    while(input >> line && !line.empty())
        kmers.push_back(line);
    adj_list=DeBruijn_graph(kmers,false );
    path=EulerianPath(adj_list);
    sequence=PathToGenome(path);
    output << sequence << endl;*/
    
    // running the k-universal string function
    /*int kmer_len;
    string universal;
    cin >> kmer_len;
    universal = kUniversal(kmer_len);
    cout << universal << endl;*/
    
    //running the string spelled by gapped pattern function - just string
    /*int kmer_len = 0, gap = 0;
    vector<string> gappedPairs;
    string sequence, line;
    input >> kmer_len;
    input >> gap;
    while(input >> line && !line.empty())
        gappedPairs.push_back(line);
    sequence=StringFromGapPattern(gappedPairs,kmer_len,gap);
    output << sequence << endl;*/
    
    //running the string reconstruction from paired reads workflow
    /*int kmer_len = 0, gap = 0;
    vector<string> paired_kmers, path;
    string line,sequence;
    map<string, vector<string> > adj_list;
    input >> kmer_len;
    input >> gap;
    while(input >> line && !line.empty())
        paired_kmers.push_back(line);
    adj_list=DeBruijn_graph(paired_kmers,true); // last arg true as the reads are paired
    printAdjList(adj_list, output);
    path=EulerianPath(adj_list);
    printPath(path, output);
    sequence=StringFromGapPattern(path,kmer_len,gap);
    output << sequence << endl;*/
    
    //running the MaximalNonBranchPath function
    /*map<string, vector<string> > adj_list;
    vector< vector<string> > paths;
    adj_list = readAdjList(input);
    paths = MaximalNonBranchPaths(adj_list);
    for(int i=0;i<paths.size();i++)
        printPath(paths[i], output);*/
    
    //running the Contig workflow here
    /*vector<string> kmers;
    string line;
    vector< vector<string> > paths;
    map<string, vector<string> > adj_list;
    while(input >> line && !line.empty())
        kmers.push_back(line);
    adj_list=DeBruijn_graph(kmers,false);
    paths = MaximalNonBranchPaths(adj_list);
    for(int i=0; i<paths.size(); i++)
        output << PathToGenome(paths[i]) << " ";
    output << endl;*/
    
    //running the translate function
    /*string rna, protein;
    input >> rna;
    protein=translate(rna);
    output << protein;*/
    
    //running the encode protein function
    /*string DNA, peptide;
    vector<string> substrings;
    input >> DNA;
    input >> peptide;
    substrings=encodeProtein(DNA, peptide);
    for(int i=0; i<substrings.size(); i++)
        output << substrings[i] << endl;*/
    
    //running the linear spectrum function
    /*string peptide;
    vector<int> spectrum;
    input >> peptide;
    spectrum=LinearSpectrum(peptide,true);
    for(int i=0;i<spectrum.size(); i++)
        output << spectrum[i] << " ";
    output << endl;*/
    
    //running the CyclopeptideSequencing function
    /*string temp;
    vector<int> spectrum;
    vector<string> peptides;
    while(input >> temp && !temp.empty())
        spectrum.push_back(stoi(temp));
    peptides=CyclopeptideSequencing(spectrum);
    for (int i=0; i<peptides.size(); i++)
        output << peptides[i] << " ";
    output << endl;*/
    
    //running the code for scoring spectrum
    /*string peptide, spec;
    vector<int>spectrum;
    input >> peptide;
    while(input >> spec && !spec.empty())
        spectrum.push_back(stoi(spec));
    cout << SpecScore(peptide, spectrum) << endl;*/
    
    //running the leaderBoard function
    /*int keepN = 0;
    string spec;
    vector<int> spectrum;
    vector<string> peptides;
    input >> keepN;
    while(input >> spec && !spec.empty())
        spectrum.push_back(stoi(spec));
    peptides=LeaderBoardCycloPepSeq(spectrum, keepN);
    for(int i=0;i<peptides.size();i++)
        output << peptides[i] << endl;*/
    
    //running the spectral convolution function
    /*vector<int> spectrum;
    map<int,int> count;
     map<int,int>::iterator it;
    string spec;
    while(input >> spec && !spec.empty())
        spectrum.push_back(stoi(spec));
    count=SpecConvolution(spectrum);
    for(it=count.begin();it!=count.end();it++){
     
        if(it->first!=0) // add all non zero elements to the convolution
            for(int i=0;i<it->second;i++)
                output << it->first << " "
     }
     output << endl; */
    
    //running the convolution leaderBoard function
    /*int keepN = 0, topAA = 0;
    string spec;
    vector<int> spectrum;
    vector<string> peptides;
    input >> topAA;
    input >> keepN;
    while(input >> spec && !spec.empty())
        spectrum.push_back(stoi(spec));
    peptides=LeaderBoardCycloPepSeq(spectrum, keepN, topAA);
    for(int i=0;i<peptides.size();i++)
        output << peptides[i] << endl;*/
    
    input.close();
    output.close();

	return 0;
}


