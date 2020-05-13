//
//  utility.cpp
//  Course2_Assembly
//
//  Created by Swati Jain on 4/26/20.
//  Copyright © 2020 Swati Jain. All rights reserved.
//

#include "genomeAssembly.h"

/*** function to print the adjacency list of the graph ***/
void printAdjList(map<string, vector<string> > adjacency_list, ostream & output){
    
    map<string, vector<string> >::iterator it;
    vector<string> temp;
    
    for(it=adjacency_list.begin(); it!=adjacency_list.end(); it++){
        
        output << it->first << " -> ";
        temp=it->second;
        for(int i=0;i<temp.size();i++){
            
            if(i==0)
                output << temp[i];
            else
                output << "," << temp[i];
        }
        output << endl;
    }
}

/*** function to read in the adjacency list of the graph from input ***/
map<string, vector<string> > readAdjList(istream & input){
    
    map<string, vector<string> > adjacency_list;
    vector<string> temp;
    string node,trash,edges,temp_node;
    int prev_index=0,index = 0;
    
    //reading the adjacency list from input
    while (input >> node && !node.empty()) {
        
        input >> trash; // reading in the ->
        input >> edges; // reading in the list of connections
        
        temp.clear();
        prev_index = 0;
        index = 0;
        //separating the connecting nodes from the list of edges
        index = edges.find(",");
        while(index != string::npos){ //finding all commas in the string and extracting the nodes that are in between
            
            temp_node = edges.substr(prev_index,index-prev_index);
            temp.push_back(temp_node);
            prev_index = index+1;
            index = edges.find(",",prev_index);
        }
        temp_node = edges.substr(prev_index);
        temp.push_back(temp_node);
        adjacency_list[node]=temp; // storing the list of nodes
    }
    
    return adjacency_list;
}

/*** function to print the path given by the nodes in the vector ***/
void printPath(vector<string> path, ostream & output){
    
    output << path[0];
    for(int i=1; i<path.size();i++)
        output << " -> " << path[i];
    output << endl;
}


/*** function to generate all binary strings of given length ***/
vector<string> GenBinStrings(int length, string cur_string, vector<string> kmers){
    
    if(cur_string.length() == length){ // the string is of required length
        
        kmers.push_back(cur_string);
        return kmers;
    }
    
    string string1, string2;
    string1 = cur_string + "0";
    string2 = cur_string + "1";
    
    kmers=GenBinStrings(length, string1, kmers); // call the function recursively
    kmers=GenBinStrings(length, string2, kmers);
    
    return kmers;
}


/*** function to transcribe a DNA into an RNA ***/
string transcribe(string DNA){
    
    string RNA = "";
    
    for(int i=0;i<DNA.length();i++){
        
        if(DNA.at(i) == 'T')
            RNA.append("U");
        else
            RNA = RNA + DNA.at(i);
    }
    return RNA;
}


/*** function to translate an RNA into a protein ***/
string translate(string RNA){
    
    string protein = "", r, p, cur_codon;
    map<string,string> codons;
    ifstream file;
    
    // creating a codon map
    file.open("RNA_codon_table_1.txt");
    while(file >> r && !r.empty()){
        file >> p;
        codons[r]=p;
    }
    file.close();
    
    //translate the RNA string
    for(int i=0; i<RNA.size(); i+=3){ // for each codon in the RNA string
        
        cur_codon=RNA.substr(i,3);
        if(codons[cur_codon] != "Stop") // elongate the chain
            protein.append(codons[cur_codon]);
        else // encountered a stop codon
            break;
    }
    
    return protein;
}


/*** function to return the reverse compliment of the given DNA string ***/
string reverseCompliment(string DNA){
    
    map<char,char> base_pair;
    string revComp = "";
    
    base_pair['A']='T';
    base_pair['C']='G';
    base_pair['G']='C';
    base_pair['T']='A';
    
    for(int i=DNA.length()-1; i>=0; i--) // traversing the DNA string in reverse order
        revComp += base_pair[DNA.at(i)];
    
    return revComp;
    
}
