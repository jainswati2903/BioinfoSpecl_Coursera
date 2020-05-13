//
//  genomeAssembly.cpp
//  Course2_Assembly
//
//  Created by Swati Jain on 4/23/20.
//  Copyright © 2020 Swati Jain. All rights reserved.
//

#include "genomeAssembly.h"

/*** function to generate all substrings of a given length from the given text (input and output from given file streams), arranged in lexicographical order ***/
vector<string> Composition_substrings(ifstream & input, ofstream & output, bool lex_order = true, bool toprint = true){
    
    int kmer_len = 0, text_len = 0, num_subs = 0;
    string text;
    vector<string> substrings; // array to store all the substrings
    
    input >> kmer_len;
    input >> text;
 
    text_len = text.length();
    num_subs = text_len - kmer_len + 1;
    
    //extract each substring of the given text
    for(int i=0; i < num_subs; i++){
        
        substrings.push_back(text.substr(i,kmer_len)); 
    }
    
    //sort the substrings in lexicographical order if the flag says so
    if(lex_order)
        sort(substrings.begin(),substrings.end());
    
    //for testing purposes - print to cout as well as return it
    if(toprint){
        for(int i=0; i< num_subs; i++){
        
            output << substrings[i] << endl;
        }
    }
    return substrings;
}


/*** function to generate a full sequence composed of kmers overlapping by kmer_len - 1 (input and output from given file streams) in order ***/
string PathProblem(ifstream & input, ofstream & output){
    
    string kmer,full_seq;
    int kmer_len;
    
    input >> full_seq; // getting the first kmer and putting that into the combined seq
    kmer_len = full_seq.length(); // length of each kmer 
    
    while(input >> kmer && !kmer.empty()) // read each substring untill it is empty
        full_seq.append(kmer,kmer_len-1,1); // append the last character of the string to the full sequence
    
    //for testing purposes, print to cout as well as return it
    output << full_seq << endl;
    
    return full_seq;
}


/*** function to generate a full sequence composed of kmers overlapping by kmer_len - 1 in order ***/
string PathToGenome(vector<string> path){
    
    string full_seq;
    int kmer_len = 0;
    
    full_seq = path[0]; // getting the first kmer and putting that into the combined seq
    kmer_len = full_seq.length(); // length of each kmer
    
    for(int i=1;i<path.size();i++)
        full_seq.append(path[i],kmer_len-1,1); // append the last character of the string to the full sequence
    
    return full_seq;
}


/*** function to print out an adjacency list of the overlap graph between given kmers (input, output from file streams). Assuming overlap of k - 1 for k-length kmers and no kmer is duplicate ***/
void Overlap_graph(ifstream & input, ofstream & output){
 
    vector<string> kmers; // storing the given kmers
    vector<string>::iterator it;
    vector< vector<int> > adjacency_list;
    vector<int> temp;
    string line;
    int kmer_len = 0;
    
    while (input >> line && !line.empty()) { // reading all the kmers
        
        it = find(kmers.begin(),kmers.end(),line);
        if(it == kmers.end()) // this kmer is not there in the list
            kmers.push_back(line);
    }
    kmer_len = kmers[0].length();
    
    for(int i=0; i<kmers.size(); i++){
        
        temp.clear(); // clear all contents
        
        for(int j=0; j<kmers.size(); j++){
            
            if(kmers[i].substr(1) == kmers[j].substr(0,kmer_len - 1)) // if suffix of kmer i mathces prefix of kmer j
                temp.push_back(j);
        }
        
        adjacency_list.push_back(temp); // add the list of suffixes to adjacency list
    }
    
    //print the overlap graph
    for(int i=0; i<adjacency_list.size(); i++){
        
        for(int j=0; j<adjacency_list[i].size(); j++){
            
            if(j==0) // if this is the first time, then print the kmer and the arrow and the first prefix kmer
                output << kmers[i] << "->" << kmers[adjacency_list[i][j]];
            else
                output << "," << kmers[adjacency_list[i][j]];
        }
        
        if(adjacency_list[i].size() != 0)
            output << endl; // printing a new line
    }
}


/*** function to return the adjacency list for the DeBruijn graph for the passed kmers ***/
map<string, vector<string> > DeBruijn_graph(vector<string> kmers, bool paired = false){
    
    int node_len = 0, index =0;
    string node1, node2, temp_kmer;
    map<string, vector<string> > adjacency_list;
    map<string, vector<string> >::iterator it;
    vector<string> temp; // temp to store the first edge for a node
    
    if(!paired){
        node_len = kmers[0].length() - 1; // length of the string node in the debruijn graph = kmer_len -1
    }
    else{
        index = kmers[0].find("|");
        temp_kmer = kmers[0].substr(0,index);
        node_len = temp_kmer.length() - 1;
    }
    
    for(int i=0;i<kmers.size();i++){ // for each kmer, extracting the two kmer_len - 1 substrings and adding that edge and the edge with the previous kmer into the adj list
        
        if(!paired){
            node1=kmers[i].substr(0,node_len);
            node2=kmers[i].substr(1);
        }
        else{ // if paired reads, the node consists of two kmer_len - 1 strings
            temp_kmer = kmers[i].substr(0,index);
            node1=temp_kmer.substr(0,node_len) + "|";
            node2=temp_kmer.substr(1) + "|";
            
            temp_kmer = kmers[i].substr(index+1);
            node1=node1 + temp_kmer.substr(0,node_len);
            node2=node2 + temp_kmer.substr(1);
        }
        
        //adding the edge between node1 and node2
        it = adjacency_list.find(node1);
        if(it == adjacency_list.end()){ // this node is not in the graph, so inserting it
            temp.clear();
            temp.push_back(node2);
            adjacency_list[node1]=temp;
        }
        else
            adjacency_list[node1].push_back(node2);
    }
    
    return adjacency_list;
}


/*** function to return the adjacency list for the DeBruijn graph for a given text and kmer length (input and output via file streams
 the bool text_given means if we have to read the text and construct substrings (true) or we are already given substrings***/
map<string, vector<string> > DeBruijn_graph(ifstream & input, ofstream & output, bool text_given = true){
    
    vector<string> kmers;
    string line;
    map<string, vector<string> > adjacency_list; // adjacency list of the debruijn graph
    
    //get all the kmers for the text from the input stream
    if(text_given)
        kmers=Composition_substrings(input,output,false,false);
    else{
        while(input >> line && !line.empty())
            kmers.push_back(line);
    }
    
    adjacency_list=DeBruijn_graph(kmers);
    return adjacency_list;
}


/*** function to find the Eulerian cycle from the given adjaceny list of the graph - assuming the graph is Eulerian
***/
vector<string> EulerianCycle(map<string, vector<string> > adjacency_list){
    
    vector<string> cycle;
    map<string, vector<string> >::iterator it;
    map<string,int> edges_left;
    string start_node, prev_node, cur_node = "dummy", last_con_node, temp_node;
    int num_edges = 0;
    
    for(it=adjacency_list.begin();it!=adjacency_list.end();it++){ // counting the total number of edges
        num_edges += it->second.size();
        edges_left[it->first]=it->second.size(); // counting the number of edges for each node
    }
    
    start_node = adjacency_list.begin()->first;
    while(cycle.size() < num_edges){ // while the cycle is not eulerain
        
        if(edges_left[start_node] == 0){ // if you are stuck at the start node
            
            for(int i=cycle.size()-1; i>=0; i--){ // look at the cycle backwards and find the node that has any connections left
                
                if(edges_left[cycle[i]] != 0){
                    last_con_node = cycle[i];
                    break;
                }
            }
            
            start_node = last_con_node; // starting from the new start node
            //updating the cycle so that starts from the last_con_node
            temp_node = cycle.back();
            cycle.pop_back();
            while(temp_node != last_con_node){
                
                cycle.insert(cycle.begin(),temp_node);
                temp_node = cycle.back();
                cycle.pop_back();
            }
            cycle.insert(cycle.begin(),temp_node); // adding the last_con_node to the begining of the cycle
        }
        
        prev_node = start_node;
        cur_node = "dummy";
        while(cur_node != start_node){ // while the path does not come back to the same node
        
            cycle.push_back(prev_node); // add the node to the cycle
            cur_node = adjacency_list[prev_node].back(); // get the connecting node
            adjacency_list[prev_node].pop_back(); // delete the connecting node
            edges_left[prev_node]--; // update the number of edges left
            prev_node = cur_node;
        }
    
    }
    cycle.push_back(start_node); // push the current start node to the end of the list to complete the cycle
    
    return cycle;
}


/*** function to generate a Eulerian path given the graph - assuming that the path exists and that this given graph is unbalanced
 i.e., we have to make it a Eulerian graph by adding an edge ***/
vector<string> EulerianPath(map<string, vector<string> > adj_list){
    
    vector<string>path,temp,temp2;
    string start_node, last_node, node1, node2;
    map<string, int> incoming;
    map<string,int>:: iterator it2;
    map<string, vector<string> >:: iterator it;
    int index = 0;
    
    //this will also add all node names to set
    for(it=adj_list.begin(); it!= adj_list.end(); it++){ // for each node in the graph
        
        // check if the node exists in the incoming list or not, if not insert it
        it2 = incoming.find(it->first);
        if(it2 == incoming.end())
            incoming[it->first] = 0;
        
        temp = it->second;
        for(int i=0; i<temp.size();i++){ // for each outgoing edge, increment the incoming edge for that node
            
            it2 = incoming.find(temp[i]);
            if(it2 == incoming.end())
                incoming[temp[i]] = 0;
            incoming[temp[i]]++;
        }
    }
    
    //check the number of incoming and outgoing edges for each node to detect the start and the end node
    for(it2=incoming.begin(); it2!= incoming.end(); it2++){
        
        it = adj_list.find(it2->first); // add to adj_list if node not found
        if(it == adj_list.end()){
            temp.clear();
            adj_list[it2->first]=temp;
        }
        
        if(it2->second > adj_list[it2->first].size()) // incoming greater than outgoing, last node
            last_node = it2->first;
        else if(it2->second < adj_list[it2->first].size()) // incoming less than outgoing, start node
            start_node = it2->first;
    }
    
    //insert the edge last_node -> start_node
    adj_list[last_node].insert(adj_list[last_node].begin(),start_node);
    //calling the eulerian cycle function
    path = EulerianCycle(adj_list);
    
    //now fix the path so that it starts and ends correctly
    path.pop_back(); // removing the last duplicate node
    if(path[0] != start_node or path.back() != last_node){ // if it does not start and end correctly
        for(index=0; index<path.size()-1; index++){
        
            if(path[index] == last_node and path[index+1] == start_node) // find the place where the extra edge is there
                break;
        }
        
        temp = path;
        path.clear();
        for(int i=index+1;i<temp.size();i++)
            path.push_back(temp[i]);
        for(int i=0; i<=index; i++)
            path.push_back(temp[i]);
    }
    
    return path;
}


/*** function to generate k-universal string for given kmer length, all kmers are binary strings of length k ***/
string kUniversal(int kmer_len){
    
    vector<string> kmers, path;
    string sequence;
    map<string, vector<string> > adj_list;
    
    kmers=GenBinStrings(kmer_len,"",kmers);
    adj_list=DeBruijn_graph(kmers);
    path=EulerianCycle(adj_list);
    for(int i=0; i<kmer_len -1; i++)
        path.pop_back(); // since the string is also cyclic
    sequence=PathToGenome(path);
    
    return sequence;
}


/*** function to return a string that is spelled by the path of paired reads of kmer_len-1 len and gap between them***/
string StringFromGapPattern(vector<string> gappedPairs, int kmer_len, int gap){
 
    string sequence = "", read1_seq, read2_seq, overlap1, overlap2;
    vector<string> reads1, reads2;
    
    for(int i=0; i<gappedPairs.size(); i++){ //separating read1 and read2 from each pair
        
        reads1.push_back(gappedPairs[i].substr(0,kmer_len-1));
        reads2.push_back(gappedPairs[i].substr(kmer_len));
    }
    
    //constructing individual sequences from reads1 and reads2
    read1_seq=PathToGenome(reads1);
    read2_seq=PathToGenome(reads2);
    
    //just for ruddii genome
    //sequence = read1_seq + "|" + read2_seq;
    
    //getting the overlap and checking if it is the same
    overlap1 = read1_seq.substr(kmer_len+gap);
    overlap2 = read2_seq.substr(0,read2_seq.length() - kmer_len - gap);
    
    if(overlap1 == overlap2) // if they are the same, then combine them
        sequence = read1_seq + read2_seq.substr(read2_seq.length() - kmer_len - gap);
    
    return sequence;
}


/*** function to find all maximal non branching paths in the given adjacency list ***/
vector< vector<string> > MaximalNonBranchPaths(map<string, vector<string> > adj_list){
    
    vector<vector<string> > paths;
    vector<string> cur_path, connect_nodes;
    map<string,int> incoming, outgoing;
    map<string, vector<string> >:: iterator it;
    map<string,int>::iterator it2;
    vector<string> OneToOne, Not_OneToOne; // to keep track of one-to-one nodes
    vector<string>::iterator it3;
    string cur_node;
    
    for(it=adj_list.begin();it!=adj_list.end();it++){
        
        //adding all nodes to incoming and outgoing lists
        it2=incoming.find(it->first);
        if(it2==incoming.end())
            incoming[it->first]=0;
        outgoing[it->first]=it->second.size();
        
        for(int i=0; i<it->second.size(); i++){
            
            // increment the total number of incoming edges
            it2=incoming.find(it->second[i]);
            if(it2 == incoming.end())
                incoming[it->second[i]] = 0;
            incoming[it->second[i]]++;
            
            //add the nodes to the outgoing map if not alraedy there
            it2=outgoing.find(it->second[i]);
            if(it2 == outgoing.end())
                outgoing[it->second[i]] = 0;
        }
    }
    
    //talking a tally of which nodes are OneToOne and which ones are not
    for(it2=incoming.begin(); it2!=incoming.end(); it2++){
        
        if(incoming[it2->first] == 1 and outgoing[it2->first] == 1)
            OneToOne.push_back(it2->first);
        else
            Not_OneToOne.push_back(it2->first);
    }
    
    //iterating over all Not_OneToOne nodes as paths will start from here
    for(int i=0;i<Not_OneToOne.size();i++){
        
        if(outgoing[Not_OneToOne[i]] == 0) // no paths if no going edges
            continue;
        
        connect_nodes=adj_list[Not_OneToOne[i]]; // list of nodes connected to this node
        for(int j=0;j<connect_nodes.size();j++){
            
            cur_path.clear();
            cur_path.push_back(Not_OneToOne[i]); // adding the first two nodes of the path
            
            //elongating the paths
            cur_node = connect_nodes[j];
            it3 = find(OneToOne.begin(), OneToOne.end(), cur_node);
            while(it3 != OneToOne.end()){ // while the current node is OneToOne - keep elongating
                
                cur_path.push_back(cur_node);
                OneToOne.erase(it3); // erase this from OneToOne as it has been included in the path
                cur_node=adj_list[cur_node][0]; // get the next connecting node
                it3 = find(OneToOne.begin(), OneToOne.end(), cur_node);
            }
            cur_path.push_back(cur_node); // as this would not be added before
            paths.push_back(cur_path);
        }
    }
    
    //find isolated cycles from the left nodes in OneToOne
    while(!OneToOne.empty()){
        
        cur_path.clear();
        cur_path.push_back(OneToOne[0]); // add the first node
        cur_node=adj_list[OneToOne[0]][0];
        OneToOne.erase(OneToOne.begin()); // erase from OneToOne list
        it3 = find(OneToOne.begin(), OneToOne.end(), cur_node);
        while(it3 != OneToOne.end()){
            
            cur_path.push_back(cur_node);
            OneToOne.erase(it3);
            cur_node=adj_list[cur_node][0];
            it3 = find(OneToOne.begin(), OneToOne.end(), cur_node);
        }
        cur_path.push_back(cur_node);
        paths.push_back(cur_path);
    }
    
    return paths;
}


