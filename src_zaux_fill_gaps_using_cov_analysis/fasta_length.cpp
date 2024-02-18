#include  <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    // g++ fasta_length.cpp -o fasta_length -O3
    if(argc < 2) {printf("Usage: fasta_length filename.fa\n"); exit(1);}
    std::string genome_file = argv[1];
    
    std::ifstream infp (genome_file.c_str());
    if(!infp.is_open())
    {
        printf("Cannot open file \'%s\' to read sequences in fasta format. Exited.\n", genome_file.c_str());
        exit(1);
    }
    printf("Reading sequence info from file:\t%s...\n", genome_file.c_str());
    
    std::string line("");
    getline(infp, line);
    while(infp.good())
    {
        if(line.size()==0) 
        {
           getline(infp, line);
           continue;
        }
        if(line.find(">") != std::string::npos)
        {
            unsigned long chrlen = 0;
            cout << line << "\t"; // seq id
            while(infp.good())
            {
                std::string chrseq("");                
                getline(infp, chrseq);
                line    = chrseq;
                if(line.find(">")!=std::string::npos) break; // next seq 
                chrlen += chrseq.size();
            }
            cout << chrlen << endl;                                                    // seq length
        }
    }
    infp.close();
}
