//
//  antibioticSequence.cpp
//  Course2_Assembly
//
//  Created by Swati Jain on 5/2/20.
//  Copyright © 2020 Swati Jain. All rights reserved.
//

#include "genomeAssembly.h"

/*** function to find all substrings of given DNA string that encode for the given protein ***/
vector<string> encodeProtein(string DNA, string peptide){
    
    vector<string> encode_strings;
    string cur_string, revComp;
    int sub_len = 0;
    
    sub_len = peptide.length() * 3; // length of the string that will encode the peptide
    
    //checking the substrings of DNA to see if they encode the peptide
    for(int i=0; i< DNA.length() - sub_len + 1; i++){
        
        cur_string = DNA.substr(i,sub_len);
        revComp = reverseCompliment(cur_string);
        if(translate(transcribe(cur_string)) == peptide || translate(transcribe(revComp)) == peptide) // if the cur_string or its reverse compliment encodes for the peptide
            encode_strings.push_back(cur_string);
    }
    
    return encode_strings;
}


/*** function to calculate the linear spectrum of a given peptide ***/
vector<int> LinearSpectrum(string peptide, bool cyclic, bool onlyMasses){
    
    vector<int> spectrum, prefix_mass;
    map<char,int> AAmass;
    int peptide_len = 0, mass =0, index =0, prev_index = 0;
    char aa;
    ifstream file;
    
    //if the given string contains one letter AA codes
    if(!onlyMasses){
        //reading the amino acid masses
        file.open("AminoAcidMass.txt");
        while(file >> aa && !isblank(aa)){
        
            file >> mass;
            AAmass[aa]=mass;
        }
        file.close();
    
        peptide_len = peptide.length();
        prefix_mass.push_back(0);
        for(int i=0;i<peptide_len;i++) // precalculating the prefix mass
            prefix_mass.push_back(prefix_mass[i] + AAmass[peptide[i]]);
    }
    else{ // if the string contains amino acid masses
        
        prefix_mass.push_back(0);
        prev_index = -1;
        index=peptide.find('-');
        while(index != string::npos){ // while the char "-" is found
            
            mass=stoi(peptide.substr(prev_index+1,index-prev_index-1)); // converting the string between two - into an integer
            prefix_mass.push_back(prefix_mass.back() + mass); // adding that mass to the previous cumulative mass
            prev_index=index;
            index=peptide.find('-',prev_index+1);
        }
        mass=stoi(peptide.substr(prev_index+1));
        prefix_mass.push_back(prefix_mass.back() + mass);
        
        peptide_len=prefix_mass.size() -1;
    }
    
    spectrum.push_back(0); // the mass of the empty subpeptide
    spectrum.push_back(prefix_mass[peptide_len]); // mass of the full peptide
    for(int i=1; i<peptide_len; i++){ //length of the subpeptide
    //for(int i=3; i==3; i++){
        if(!cyclic){
            for(int j=0; j<peptide_len - i + 1; j++)
                spectrum.push_back(prefix_mass[j+i] - prefix_mass[j]); // mass of the substring
        }
        else{
            for(int j=0; j<peptide_len; j++){
                
                if(j >= peptide_len - i + 1) // if the substring rolls back to the start of the peptide
                    spectrum.push_back(prefix_mass[peptide_len] - (prefix_mass[j] - prefix_mass[i - peptide_len + j]));
                else
                    spectrum.push_back(prefix_mass[j+i] - prefix_mass[j]); // mass of the substring
            }
        }
    }
    
    sort(spectrum.begin(), spectrum.end());
    return spectrum;
}


/*** function to check if spectrum small_spec is consistent with spectrum full spec, i.e. all elements of small_spec are contained within full_spec ***/
bool SpectrumConsistent(vector<int> small_spec, vector<int> full_spec){
    
    bool flag = true;
    int index = 0; // keep track of small_spec index
    
    for(int i=0;i<full_spec.size();i++){
        
        if (full_spec[i] > small_spec[index]){ // small_spec cur element will not be found in full_spec
            flag=false;
            break;
        }
        
        if(full_spec[i] == small_spec[index]) // current small_spec element found in full_spec, go to next small_spec element
            index++;
        
        if(index == small_spec.size()) // all elements in small_spec are found
            break;
    }
    
    return flag;
}


/*** function to return peptide sequences (in masses) that matches the given spectrum ***/
// this function can be made more efficient by only considering aa in expansion ones whose individual masses are present in the spectrum 
vector<string> CyclopeptideSequencing(vector<int> spectrum){
    
    vector<string> peptides, single_aa, consistent_peptides;
    ifstream file;
    string temp;
    int peptide_len = 0, spec_len = 0;
    
    //read in amino acid masses - both as initial peptides and alone as well for extension
    file.open("UniqueAAMasses.txt");
    while(file >> temp && !temp.empty()){
        peptides.push_back(temp); // initial peptides
        single_aa.push_back(temp); // individual aas used for extension
    }
    file.close();
    
    spec_len = spectrum.size() - 2; // exclude the 0 and the full peptide mass
    if(1+sqrt(4*spec_len+1)/2 > 0)
        peptide_len = int(1+sqrt(4*spec_len+1)/2);
    else
        peptide_len = int(1-sqrt(4*spec_len+1)/2);
    
    for(int i=1;i<peptide_len;i++){ // while the length of the peptide is less than the full length
        
        consistent_peptides.clear();
        for(int j=0;j<peptides.size();j++){
            
            if(SpectrumConsistent(LinearSpectrum(peptides[j],false,true), spectrum)) // if the spectrum of the current peptide is consistent, keep it
                consistent_peptides.push_back(peptides[j]);
        }
        
        //extend the consistent peptides by one amino acid
        peptides.clear();
        for(int j=0; j<consistent_peptides.size();j++){
            for(int k=0;k<single_aa.size();k++){
                peptides.push_back(consistent_peptides[j] + "-" + single_aa[k]);
            }
        }
    }
    
    //Now all peptides are the required length, now check consistency for the full peptide
    consistent_peptides.clear();
    for(int i=0;i<peptides.size();i++){
        if(SpectrumConsistent(LinearSpectrum(peptides[i],true,true), spectrum))
            consistent_peptides.push_back(peptides[i]);
    }
    
    return consistent_peptides;
}


/*** function to return the score of the cyclic peptide with respect to the given spectrum
 this assumes that the two spectrums are sorted ***/
int SpecScore(string peptide, vector<int> spectrum, bool cyclic, bool onlyMasses){
    
    int score = 0, spec_index = 0;
    vector<int> peptide_spectrum;
    
    peptide_spectrum=LinearSpectrum(peptide,cyclic,onlyMasses); // calculate the linear/cyclic spectrum of the given peptide
    
    for(int i = 0; i< spectrum.size();i++){
     
        if(spectrum[i] == peptide_spectrum[spec_index]){ // found the peptide spectrum cur element in the main spectrum
            score++;
            spec_index++;
        }
        else if (spectrum[i] > peptide_spectrum[spec_index]){ // if it is greater, then cur element of peptide spectrum will not be found
            spec_index++; // go to the next element
            i--; // so that the current element can be compared with the new peptide spectrum cur element as well
        }
        
        if(spec_index == peptide_spectrum.size()) // all elements of the peptide specturm have been identified
            break;
    }
    
    return score;
}


/*** function to return the mass of a given peptide ***/
int peptideMass(string peptide, bool onlyMasses){
    
    map<char,int> AAmass;
    int total_mass = 0, mass = 0, index = 0, prev_index = 0;
    char aa;
    ifstream file;
    
    //if the given string contains one letter AA codes
    if(!onlyMasses){
        //reading the amino acid masses
        file.open("AminoAcidMass.txt");
        while(file >> aa && !isblank(aa)){
            
            file >> mass;
            AAmass[aa]=mass;
        }
        file.close();
        
        for(int i=0;i<peptide.length();i++) // calculating the total mass
            total_mass += AAmass[peptide[i]];
    }
    else{ // if the string contains amino acid masses
        
        prev_index = -1;
        index=peptide.find('-');
        while(index != string::npos){ // while the char "-" is found
            
            mass=stoi(peptide.substr(prev_index+1,index-prev_index-1)); // converting the string between two - into an integer
            total_mass+=mass;
            prev_index=index;
            index=peptide.find('-',prev_index+1);
        }
        mass=stoi(peptide.substr(prev_index+1));
        total_mass+=mass;
    }
    
    return total_mass;
}

/*** function to return the peptides that has the highest score with respect to the spectrum using leader board also where at each step you only keep given number of candidates
 the peptide is given and grown in terms of aa masses ***/
vector<string> LeaderBoardCycloPepSeq(vector<int> spectrum, int keepN, int topAA){
    
    string temp, new_peptide;
    vector<string> single_aa, leaderPeptide;
    ifstream file;
    map<int,int> convolution;
    map<int,int>::iterator it;
    int massSpectrum = 0, pepScore = 0, leaderScore = 0, pepMass = 0, count = 0, NScore = 0, pepScore_cyclic = 0, AAcount = 0;
    priority_queue<pair<int,string> > leaderBoard, newLeaderBoard;
    priority_queue<pair<int,int> > countAA;
    pair<int, string> peptide;
    pair<int,int> AA;
    
    // read individual aas used for extension
    if(topAA == 0){ // consider all amino acids
        file.open("UniqueAAMasses.txt");
        while(file >> temp && !temp.empty())
            single_aa.push_back(temp);
        file.close();
        //for(int i=57;i<=200;i++) // expanded list of amino acids
        //    single_aa.push_back(to_string(i));
    }
    else{ // consider all masses that are frequent in the spectrum
        convolution=SpecConvolution(spectrum);
        for(it=convolution.begin();it!=convolution.end();it++)// make a queue for all masses within 57 and 200
            if(it->first >=57 && it->first <=200)
                countAA.push(make_pair(it->second, it->first));
        
        //taking the topAA number of masses with ties
        count = 0;
        AAcount = 0;
        while(!countAA.empty()){
            
            AA=countAA.top();
            countAA.pop();
            count++;
            
            if(count <= topAA || AA.first == AAcount){
                single_aa.push_back(to_string(AA.second));
                AAcount = AA.first;
            }
        }
    }
    
    massSpectrum=spectrum.back(); // the mass of the spectrum is the last element in the spectrum
    leaderBoard.push(make_pair(1, "")); // initialize the leader board with an empty peptide
    // initialize the leader peptide with an empty peptide
    leaderScore = 1;
    
    while(!leaderBoard.empty()){
    //for(int t=0;t<4;t++){
        
        while(!leaderBoard.empty()){ // for each peptide in the leaderBoard, pop the peptide and expand it
        
            peptide=leaderBoard.top();
            leaderBoard.pop();
            for(int j=0; j<single_aa.size(); j++){
                
                // expand the peptide, gets its score and mass
                new_peptide = peptide.second + "-" + single_aa[j];
                pepScore = SpecScore(new_peptide.substr(1), spectrum, false, true); // substring to remove the extra "-"
                pepScore_cyclic = SpecScore(new_peptide.substr(1), spectrum, true, true); // substring to remove the extra "-"
                pepMass=peptideMass(new_peptide.substr(1),true);
                
                if(pepMass <= massSpectrum) // if the mass is less than or equal to the spectrum mass, insert it into the new leader board
                    newLeaderBoard.push(make_pair(pepScore, new_peptide));
                
                if(pepMass == massSpectrum){
                    
                    if(pepScore_cyclic > leaderScore){ // if the mass if equal to the spectrum and score is better than leader peptide, then update leader peptide
                        leaderScore=pepScore_cyclic;
                        leaderPeptide.clear();
                        leaderPeptide.push_back(new_peptide.substr(1)); // substring to remove the extra "-"
                    }
                    else if (pepScore_cyclic == leaderScore)
                        leaderPeptide.push_back(new_peptide.substr(1)); // multiple peptides with the same score
                }
            }
        }
        // leaderBoard is empty, now have to trim the newLeaderBoard by adding top elements to leaderBoard
        count = 0;
        while(!newLeaderBoard.empty()){
                
            peptide=newLeaderBoard.top();
            newLeaderBoard.pop();
            count++;
                
            if(count <= keepN || peptide.first == NScore){
                leaderBoard.push(peptide);
                NScore = peptide.first;
            }
        }
    }
    
    cout << leaderScore << endl;
    return leaderPeptide;
}


/*** function return the convolution of the spectrum, with each element counted as many times as it occurs in the convolution ***/
map<int,int> SpecConvolution(vector<int> spectrum){
    
    int diff = 0;
    map<int,int> count;
    map<int,int>:: iterator it;
    
    for(int i=0;i<spectrum.size();i++){ // take the difference between cur element and all elements before it, and add/increment to the map
        
        for(int j=0; j<i;j++){
            
            diff = spectrum[i] - spectrum[j];
            it=count.find(diff);
            if(it==count.end())
                count[diff]=1;
            else
                count[diff]++;
        }
    }
    
    return count;
}
