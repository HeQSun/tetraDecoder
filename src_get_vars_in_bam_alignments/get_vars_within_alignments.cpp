/* Function: check how many variations in an aligment
             to determine whether if there is share between reference and query sequences.

   Written by Hequan Sun
   Email: hequan.sun@xjtu.edu.cn/sunhequan@gmail.com
   20230522
*/
#include   <iostream>
#include    <fstream>
#include    <sstream>
#include        <map>
#include     <vector>
#include     <string>
#include    <iomanip>
#include   <stdlib.h>
#include  <algorithm>
#include <sys/stat.h>
#include   <dirent.h>
#include   <assert.h>
#include     <time.h>  /* clock_t, clock, CLOCKS_PER_SEC */
#include "split_string.h"
struct CIGAR{
   unsigned long qry_len;
   unsigned long num_X;
   unsigned long num_E;
   unsigned long num_D;
   unsigned long num_I;
   unsigned long len_S;
   unsigned long len_H;
};
int decipher_cigar(string       cigar,
                 vector<char>*  operation,
                 vector<int>*   count,
                 CIGAR*         cigar_dec);

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        cout << "\n Function: check variations in alignments from a bam, \n"
             << "           determine if share exists between query and reference sequences.  \n"
             << " Usage...: samtools view -F 3840 query_to_ref.bam | get_vars_within_alignments - out_prefix_str "
             << endl;
        cout << "           note: -F 3840==-F 256 -F 512 -F 1024 -F 2048, excluding all non-primary alignment!"
             << endl << endl;
        return -1;
    }
    clock_t tbeg;
    tbeg = clock();
    string outprefix = (string)argv[2];
    if(outprefix.size()==0) outprefix = "fun_convert";
    //
    bool verbose = true; // to add future option
    // open files under folder
    stringstream ofilename;
    ofilename.str("");
    ofilename << outprefix << "_vars_between_query_and_ref.txt";
    ofstream ofp;
    ofp.open((ofilename.str()).c_str());
    if(!ofp.good())
    {
         cout << "Cannot open new file " << ofilename.str() << " to write aln-var info. Exited.\n";
         return -1;
    }
    //
    cout << "   Info: variations between query reads and reference sequence will be collected in " << endl
         << "      " << ofilename.str() << endl;
    ofp  << "#qry_name"   << "\tres\t"
         << "l_qry"       << "\t"
         << "n_snp"       << "\t"
         << "n_mat"       << "\t"
         << "n_del"       << "\t"
         << "n_ins"       << "\t"
         << "l_socl"      << "\t"
         << "l_hacl"      << endl; //
    //
    unsigned long numraw   = 0;
    unsigned long total_R1 = 0;
    string    exampleNoSEQ = "";
    unsigned long numNoSEQ = 0;
    string     exampleNoCG = "";
    unsigned long  numNoCG = 0;
    string     exampleNoMA = "";
    unsigned long numNoMA  = 0; // number of reads showing multiple alignment: only one kept.
    //
    std::string line;
    while (std::getline(std::cin, line))
    {
        if(line.compare("quit")==0 || line.compare("q")==0 || line.compare("exit")==0) break;
        if(line.size()==0 || line[0]=='#') continue;
        numraw ++;
        if(numraw%100000000 == 0)
        {
            cout << "   info: " << numraw << "th alignment..." << endl;
        }
        //
        vector<string> lineinfo = split_string(line, '\t');
        //
        string pbname    = lineinfo[0]; // read name
        string thisflag  = lineinfo[1]; // flag value
        string thisctg   = lineinfo[2]; // reference contig name; there are special case
        string thispos   = lineinfo[3]; // first matching reference coordinate
        string thiscigar = lineinfo[5]; // cigar string
        string thisseq   = lineinfo[9]; // read sequence
        if(thisseq.compare("*") == 0)
        {
            numNoSEQ ++;
            if(exampleNoSEQ.size()==0)
            {
                exampleNoSEQ = line;
                /*
                cout << endl
                     << "   Warning: there are alignments without clear sequence, secondary skipped, e.g.: "
                     << line
                     << endl;
                */
                cout << "   Info: you will get a total number of such alignments in the end of program. "
                     << endl;
            }
            continue; // caution
        }
        // special case 1: cigar string as star: not aligned
        if(thiscigar.compare("*") == 0)
        {
            if(exampleNoCG.size()==0)
            {
                exampleNoCG = line;
                /*
                cout << endl
                     << "   Warning: there are alignments without explicit CIGAR string, skipped., e.g.: "
                     << line
                     << endl;
                */
                cout << "   Info: you will get a total number of such alignments in the end of program. "
                     << endl;
            }
            numNoCG ++;
            continue;
        }
        // special case 2:
        int hexflag = strtol(thisflag.c_str(), NULL, 0);
        if(hexflag > 255)  // do not skip such lines???????? to discuss!!!
        {
            numNoMA ++;
            if(exampleNoMA.size()==0)
            {
                exampleNoMA = line;
                /*
                cout << endl
                     << "   Warning: there are alignments being secondary/supplementary etc, skipped, e.g.: "
                     << line
                     << endl;
                */
                cout << "   Info: you will get a total number of such alignments in the end of program. "
                     << endl;
            }
            continue;
        }
        // cout << "   Check: read: " << pbname << " aligned to " << thisctg << "\t" << thispos << endl;
        //
        vector<char>  operation;
        vector<int>   count;
        CIGAR         cigar_dec;
        // covered len by alignment
        int           covlen  = decipher_cigar(thiscigar, &operation, &count, &cigar_dec);
        // alignment start and end position on ref seq
        unsigned long firstRefsMatch = strtoul(thispos.c_str(), NULL, 0);
        unsigned long spansecond     = firstRefsMatch+covlen-1; // 1-based
        // cout << "          span: " << firstRefsMatch              << "-"   << spansecond
        //     << " => "             << spansecond-firstRefsMatch+1 << " bp" << endl; // 1-based
        double snp_nonzero = cigar_dec.num_X;
        if(cigar_dec.num_X == 0)
        {
          snp_nonzero = 1.0;
        }
        ofp  << pbname            << "\tres\t"
             << cigar_dec.qry_len << "\t"
             << cigar_dec.num_X   << "\t"
             << cigar_dec.num_E   << "\t"
             << cigar_dec.num_D   << "\t"
             << cigar_dec.num_I   << "\t"
             << cigar_dec.len_S   << "\t"
             << cigar_dec.len_H   << "\t"
             << std::fixed
             << setprecision(4)
             << cigar_dec.num_E*1.0 / cigar_dec.qry_len << "\t"
             << setprecision(0)
             << cigar_dec.qry_len*1.0/snp_nonzero       << endl;
        //
        line.clear();
    }
    // close files
    ofp.close();
    cout << "   Info: time on analyzing streaming bam info: "
         << (float)(clock()-tbeg)/CLOCKS_PER_SEC
         << " seconds.\n"
         << endl;
    return 0;
}
int decipher_cigar(string cigar, vector<char>* operation, vector<int>* count, CIGAR* cigar_dec)
{
    /*
	https://samtools.github.io/hts-specs/SAMv1.pdf:
        //
	Op	BAM	Description						Consumes_query	Consumes_reference
	M	0	alignment match (can be a sequence match or mismatch)	yes		yes
	I	1	insertion to the reference				                    yes		no
	D	2	deletion from the reference				                     no		yes
	N	3	skipped region from the reference		                 	 no		yes
	S	4	soft clipping (clipped sequences present in	SEQ)	    yes		no
	H	5	hard clipping (clipped sequences NOT present in	SEQ)	 no		no
	P	6	padding (silent deletion from padded reference)		     no		no
	=	7	sequence match						                            yes		yes
	X	8	sequence mismatch					                            yes		yes
	NOTE:
        // CIGAR alphabet: [0-9]+[MIDNSHPX=] and *: 23S23M1D8M1I16M1D29M38H
        // Sum of lengths of the MIS=X operations shall equal the length of SEQ in bam
        // H can only be present as the first and/or last operation.
        // S may only have H operations between them and the ends of the CIGAR string
        // or mRNA-to-genome alignment, an N operation represents an intron.
    */
    // check number of variations
    unsigned long num_X  = 0; // get length of mismatches
    unsigned long num_E  = 0; // get length of    matches
    unsigned long num_D  = 0; // get time of deletions  occur, not length of all deletions.
    unsigned long num_I  = 0; // get time of insertions occur, not length of all insertions.
    unsigned long len_I  = 0; // get length of insertions.
    unsigned long len_S  = 0; // get length of soft clipping
    unsigned long len_H  = 0; // get length of hard clipping
    //
    char *cstr = (char*)cigar.c_str(); //
    string numstr("");
    unsigned long covlen = 0;
    for (int i=0; i<cigar.size(); i++)
    {
        if(cstr[i]>='0' && cstr[i]<='9')
        {
            numstr += cstr[i];
        }
        else
        if(cstr[i] == 'M' || cstr[i] == 'I' || cstr[i] == 'D' ||
           cstr[i] == 'N' ||
           cstr[i] == 'S' || cstr[i] == 'H' ||
           cstr[i] == 'P' ||
           cstr[i] == 'X' || cstr[i] == '=')
        {
            unsigned long opt_len = strtoul(numstr.c_str(), NULL, 0);
            (*operation).push_back(cstr[i]);
            (*count).push_back(opt_len);
            //
            if(cstr[i] == 'M' || cstr[i] == 'D' || cstr[i] == 'N' || cstr[i] == 'X' || cstr[i] == '=')
            {
                covlen += opt_len; // positions of ref covered by read
            }
            if(cstr[i] == 'P' || cstr[i] == 'N')
            {
                cout << "   Warning: you have special operation \'" << cstr[i] << "\' in Cigar: " << cigar << endl;
            }
            if(cstr[i] == 'X')
            {
                // mismatch found
                num_X += opt_len;
            }
            if(cstr[i] == '=')
            {
                // match found
                num_E += opt_len;
            }
            if(cstr[i] == 'D')
            {
                // one more deletion found
                num_D += 1; // atoi(numstr.c_str());
            }
            if(cstr[i] == 'I')
            {
                // one more insertion found
                num_I += 1; // atoi(numstr.c_str());
                len_I += opt_len;
            }
            if(cstr[i] == 'S')
            {
                // one more insertion found
                len_S += opt_len;
            }
            if(cstr[i] == 'H')
            {
                // one more insertion found
                len_H += opt_len;
            }
            numstr.clear();
        }
    }
    // get length of query of sequence
    unsigned long qry_len = num_X + num_E + len_I + len_S + len_H;
    (*cigar_dec).num_X    = num_X;
    (*cigar_dec).num_E    = num_E;
    (*cigar_dec).num_D    = num_D;
    (*cigar_dec).num_I    = num_I;
    (*cigar_dec).len_S    = len_S;
    (*cigar_dec).len_H    = len_H;
    (*cigar_dec).qry_len  = qry_len;
    //
    return covlen;
}
