#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctype.h>
#include "split_string.h"
//
using namespace std;
bool read_chromosome_seq(char* fgenome, char* chr_id,
                         unsigned long max_chr_len, std::string* chr_seq);
//
bool extract_substr_chr(string fasta, string chr, long lmutMargin, long rmutMargin, string* subseq_name, long* subseq_size)
{
    bool extract_from_ref_seq = true; // reading consensus from ref seq instead of consensus_summary.txt
    if(extract_from_ref_seq)
    {
        /* read in chromosome seq */
        std::string chr_seq;
        if(!read_chromosome_seq((char*)fasta.c_str(), (char*)chr.c_str(), 100000000, &chr_seq) )
        {
            printf("   Error: cannot find specified chromosome id. Exited (in function extract_substr_chr(...)).\n");
            printf("   Hint: check if all IDs of chromosomes in all files are consistent.\n");
            return false;
        }
        /* prepare output file*/
        long reg_begin = lmutMargin; // cnt of seq position starting from 1
        long reg_end   = rmutMargin; // cnt of seq position starting from 1
        cout << "   Info: given region to extract sequence: " << reg_begin << "-" << reg_end << endl;
        if(reg_end > chr_seq.size())
        {
            printf("   warning: not enough bases on the right side of the mutation position. Resetting...\n");
            reg_end = chr_seq.size();
        }
        if(reg_begin <= 0)
        {
            printf("   Info: left starting position is <=0. Reset as 1.\n");
            reg_begin = 1;
        }
        char outFile[64];
        //sprintf(outFile, "subseq_%s_%ld_%ld.fasta", (char*)chr.c_str(), reg_begin, reg_end);
        sprintf(outFile, "subseq_tmp.fasta");
        FILE* outfp = fopen(outFile, "w");
        if(outfp == NULL)
        {
            printf("   Error: cannot open output file.\n");
            return false;
        }
        /* output substring to file */
        long vstart = reg_begin - 1; // c/c++ strings start from 0
        long len    = reg_end   - reg_begin + 1;
        fprintf(outfp, ">%s_%ld_%ld_cp\n", (char*)chr.c_str(), reg_begin, reg_end);
        fprintf(outfp, "%s", chr_seq.substr(vstart, len).c_str());
        char tmpname[100];
        sprintf(tmpname,  "%s_%ld_%ld_cp", (char*)chr.c_str(), reg_begin, reg_end);
        // return values
        *subseq_name = (string)tmpname;
        *subseq_size = len;
        cout << "   Info: extracted sub-seq name " << *subseq_name
             << " with size "                      << *subseq_size
             << endl;
        /* close file */
        fclose(outfp);
        cout << "   Info: subsequence extracted into subseq_tmp.fasta. " << endl;
    }
    return true;
}
/* read in a chromosome sequence: max_chr_len might be obtained from chrsizes files */
bool read_chromosome_seq(char* fgenome, char* chr_id, unsigned long max_chr_len, std::string* chr_seq)
{
    (*chr_seq).clear();
    unsigned long chrlen = 0;
    std::ifstream infp (fgenome);
    if(!infp.is_open())
    {
        printf("   Error: cannot open file \'%s\' to read sequences in fasta format. Exited.\n", fgenome);
        return false;
    }
    printf("   Info: reading sequence info from file:\t%s...\n", fgenome);
    while(infp.good())
    {
        std::string line("");
        getline(infp, line);
        if(line.size()==0) continue;
        if(line.substr(1).compare((string)chr_id) == 0)
        {
            cout << "   Info: this chr is " << line.substr(1) << "\t"; // seq id
            while(infp.good())
            {
                std::string chrseq("");
                getline(infp, chrseq);
                line    = chrseq;
                if(line.find(">")!=std::string::npos) break; // next seq
                chrlen += chrseq.size();
                (*chr_seq) += chrseq;
            }
            cout << chrlen << endl;                         // seq length
        }
        if(chrlen>0) break;
    }
    infp.close();
    if(chrlen==0) return false;
    return true;
}
