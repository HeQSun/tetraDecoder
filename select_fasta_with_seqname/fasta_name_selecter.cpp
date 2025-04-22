#include  <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  <fstream>
#include <iostream>
#include  <sstream>
#include <map>
#include <vector>

#include "split_string.h"

using namespace std;

int main(int argc, char* argv[])
{
    // g++ fasta_name_selecter.cpp -o fasta_name_selecter -O3
    if(argc < 3) {printf("Usage: fasta_name_selecter seq.fa seqname.txt\n"); return false;}
    std::string genome_file = argv[1];
    
    ifstream snfp;
    snfp.open(argv[2]);
    if(!snfp.is_open())
    {
        cout << "   Error: cannot open seq name file: " << argv[2] << endl;
        return false;
    }
    map<string, int> seqname;
    while(snfp.good())
    {
        string line("");
        getline(snfp, line);
        if(line.size()==0) continue;
        
        seqname.insert(std::pair<string, int>(line, 0));
    }
    cout << seqname.size() << " seq names have been collected. " << endl;
    snfp.close();
    
    // prepare output file
    vector<string> split_genomefilename = split_string(genome_file, '/');
    
    std::stringstream ss;
    ss.str("");
    int replace_len = 0;
    size_t pos;
    string ofilename = split_genomefilename[split_genomefilename.size()-1];
    if( ofilename.find(".fasta") != std::string::npos)
    {
        ss << "_selected" << ".fasta";
        replace_len = 6;    
        pos = ofilename.find(".fasta");
    }
    else
    if( ofilename.find(".fas") != std::string::npos)
    {
        ss << "_selected" << ".fas";
        replace_len = 4;
        pos = ofilename.find(".fas");
    }
    else
    if( ofilename.find(".fa") != std::string::npos)
    {
        ss << "_selected" << ".fa";
        replace_len = 3;
        pos = ofilename.find(".fa");
    }
    else
    {
        cout << "   Error: input file is not with .fa, .fas or .fasta suffix. " << endl;
        return false;
    }
    ofilename.replace(pos, pos+replace_len-1, ss.str());
    std::ofstream ofp;
    ofp.open(ofilename.c_str(), std::ios::out);
    if(!ofp.is_open())
    {
        cout << "   Error: cannot open " << ofilename << " to write selected sequences. " << endl;
        return false;
    }
    
    // open input fasta to read sequences
    std::ifstream infp (genome_file.c_str());
    if(!infp.is_open())
    {
        printf("Cannot open file \'%s\' to read sequences in fasta format. Exited.\n", genome_file.c_str());
        exit(1);
    }
    printf("Reading sequence info from file:\t%s...\n", genome_file.c_str());
    
    // traverse input fasta file
    int total_seq_num      = 0;
    int total_seq_selected = 0;
    std::string line("");
    getline(infp, line);
    while(infp.good())
    {
        if(line.find(">") != std::string::npos)
        {
            string seq_name = line;
            total_seq_num ++;
            string seq("");
            int seq_line_id = 0;
            while(infp.good())
            {
                std::string chrseq("");                
                getline(infp, chrseq);
                seq_line_id ++;
                line    = chrseq;
                if(line.find(">")!=std::string::npos) break; // next seq
                if(seq_line_id > 1) seq += '\n';
                seq    += chrseq;
            }
            if(seqname.find(seq_name.substr(1)) != seqname.end())                                                 // seq length
            {
                total_seq_selected ++;
                if(total_seq_selected >= 2)
                ofp << endl;
                ofp << seq_name << endl;
                ofp << seq;
            }
        }
    }
    cout << "   Info: total number of sequences in give fasta: " << total_seq_num      << endl;
    cout << "         total number of sequences selected:      " << total_seq_selected << endl;
    infp.close();
    ofp.close();
}
