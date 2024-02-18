/*
    Given 
    
        a list of phased snp markers within a genomic region
        a bam file with reads aligned,
        
    separate reads in that region into 2 clusters.
    
    2022-09-23 Started: Hequan Sun, MPIPZ
    Email: sun@mpipz.mpg.de/sunhequan@gmail.com
*/
#include         <map>
#include      <string>
#include      <vector>
#include     <fstream>
#include     <sstream>
#include    <iostream>
#include   <algorithm>
#include    <string.h>
#include    <stdlib.h>
#include    <assert.h>
#include  <sys/stat.h>
#include    <dirent.h>
#include    <unistd.h>
#include "split_string.h"
using namespace std;
//
struct PAMA
{
    string patallele;
    string matallele;
};
double minCOscore  = 0.64; // 
//
bool read_marker(string fmarker,
                 map<string, map<unsigned long, PAMA> >* mkr);
int decipher_cigar(string                                cigar, 
                 vector<char>*                           operation, 
                 vector<int>*                            count);
bool align_read_to_ref(vector<char>                      operation, 
                 vector<int>                             count, 
                 string*                                 orignal_read, 
                 string*                                 updated_read); 
bool genotype_pbread(string*                             updated_read, 
                 map<unsigned long, PAMA>*               this_marker, 
                 unsigned long                           start, 
                 unsigned long                           end, 
                 string*                                 geno,
                 map<unsigned long, PAMA>*               intersectMkr);
//
int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        cout << "\n   Function: extract reads carrying alternative alleles in a given genomic region. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "samtools view sample.bam[or sample.cram] region:ctg_id | bam_extracter snp_list.txt out_prefix_str -"
             << endl << endl;
        return 1;
    }
    double startT= clock();
    //
    string marker_file = (string)argv[1];
    string out_prefix  = (string)argv[2];
    // s1. read markers
    cout << "   Info: reading snp markers from " << marker_file << endl;
    map<string, map<unsigned long, PAMA> > snp_marker; // map<chr, map<position, {mat-allele, pat-allele} > >
    if(!read_marker(marker_file, &snp_marker) )
    {
        cout << "   Error: read snp marker failed. " << endl;
        return false;
    }      
    cout << "   Info: reading snp marker done. " << endl;  
    //
    string out_file = out_prefix + "_extracted_reads.fq";
    ofstream ofp;
    ofp.open(out_file.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open output file " << out_file << endl;
        return 1;
    }
    //
    unsigned long numraw   = 0;
    unsigned long numNoSEQ = 0;
    unsigned long numNoCG  = 0;
    unsigned long numNoMA  = 0;
    string exampleNoSEQ("");
    string exampleNoCG("");
    string exampleNoMA("");
    std::string line;
    while (std::getline(std::cin, line)) 
    {
        if(line.compare("quit")==0 || line.compare("q")==0 || line.compare("exit")==0) break;
        if(line.size()==0) continue;
        if(line[0] == '@') continue;
        //
        numraw ++;
        if(numraw%100000 == 0)
        {
            cout << "   info: " << numraw << "th aligm..." << endl;
        }        
        //        
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size()<11)
        {
            cout << "   Warning: insufficient line info, skipped: " << line << endl;
            continue;
        }
        string pbname    = lineinfo[0]; // read name
        string thisflag  = lineinfo[1]; // flag value
        string thisctg   = lineinfo[2]; // reference contig name; there are special case
        string thispos   = lineinfo[3]; // first matching reference coordinate
        string thiscigar = lineinfo[5]; // cigar string
        string thisseq   = lineinfo[9]; // read sequence  
        string thisqua   = lineinfo[10];
        if(thisseq.compare("*") == 0) 
        {
            numNoSEQ ++;
            if(exampleNoSEQ.size()==0)
            {
                exampleNoSEQ = line;
                cout << endl
                     << "   Warning: there are alignments without clear sequence, secondary skipped, e.g.: " 
                     << line 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;                
            }
            continue; // caution
        }
        // special case 1: cigar string as star: not aligned
        if(thiscigar.compare("*") == 0) 
        {
            numNoCG ++;
            if(exampleNoCG.size()==0)
            {
                exampleNoCG = line;
                cout << endl
                     << "   Warning: there are alignments without explicit CIGAR string, skipped., e.g.: " 
                     << line 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;
            }
            //     
            continue;
        }       
        // special case 2:
        int hexflag = strtol(thisflag.c_str(), NULL, 0);
        if(hexflag > 255)  // do not skip such lines?
        {
            numNoMA ++;
            if(exampleNoMA.size()==0)
            {
                exampleNoMA = line;
                cout << endl
                     << "   Warning: there are alignments being secondary/supplementary etc, skipped, e.g.: " 
                     << line 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;                
            }
            continue;
        }         
        //
        string alignedtogrp = "";
        cout << "   Check: read: " << pbname << " aligned to " << thisctg << "\t" << thispos << endl;
        // If it contains 0x10, SEQ in BAM/CRAM/SAM has been reversed+complemented from read in fastq
        string reversed = "nrc";
        if((hexflag & 0x10) == 0x10) reversed = "rc";
        // covered len by alignment
        vector<char>  operation;
        vector<int>   count;        
        int           covlen  = decipher_cigar(thiscigar, &operation, &count);                   
        // alignment start and end position on ref seq: these would intersect interested markers for current pb read
        unsigned long firstRefsMatch = strtoul(thispos.c_str(), NULL, 0);                       
        unsigned long spansecond     = firstRefsMatch+covlen-1; // 1-based         
        cout << "          span: " << firstRefsMatch              << "-"   << spansecond 
             << " => "             << spansecond-firstRefsMatch+1 << " bp" << endl; // 1-based 
        //
        string updated_read("");
        if(!align_read_to_ref(operation, count, &thisseq, &updated_read))
        {
            cout << "   Warning: unexpected read " << thisseq << endl;
            continue;
        }
        cout << "        : original read: " 
             << thisseq.size()      << " bp; " << endl;
        cout << "        : updated read with \'I\' removed, and \'D:-\' added: " 
             << updated_read.size() << " bp. " << endl;            
        // cout << "        " << updated_read << endl;         
        // prepare all markers along this contig        
        map<unsigned long, PAMA> this_contig_marker;
        this_contig_marker.clear();
        map<string, map<unsigned long, PAMA> >::iterator cpitr = snp_marker.find(thisctg);     
        if(cpitr != snp_marker.end())
        {
            this_contig_marker = snp_marker[thisctg];
        }
        cout << "        : there are " << this_contig_marker.size() << " snp markers for " << thisctg << endl;
        // snp markers: get genotype of the read along intersected snp markers
        string                   thisgeno; // PMu pattern
        map<unsigned long, PAMA> intersectMkr;
        genotype_pbread(&updated_read, 
                        &this_contig_marker, 
                        firstRefsMatch, 
                        spansecond, 
                        &thisgeno,
                        &intersectMkr);
        //cout << "   Check this geno: " << thisgeno << endl;
        size_t pcnt = std::count(thisgeno.begin(), thisgeno.end(), 'P'); // alt
        size_t mcnt = std::count(thisgeno.begin(), thisgeno.end(), 'M'); // ref
        double alt_ratio = pcnt*1.0/(pcnt + mcnt + 0.0000001);         // caution: avoid 0
        if(alt_ratio > 0.65) // mostly alt, at least 2/3
        {
            ofp << "@"      << pbname 
                << "/"      << thisctg 
                << ":"      << firstRefsMatch 
                << "-"      << spansecond 
                << "/"      << pcnt
                << "/"      << pcnt + mcnt
                << endl;
            ofp << thisseq  << endl;
            ofp << "+"      << endl;
            ofp << thisqua  << endl;
        }
        cout << "   Info: ratio of alt allele on this read: " << alt_ratio << endl << endl;
    }
    //
    ofp.close();
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;
}
//
bool genotype_pbread(string*                   updated_read, 
                     map<unsigned long, PAMA>* this_marker, 
                     unsigned long             start, 
                     unsigned long             end, 
                     string*                   geno,
                     map<unsigned long, PAMA>* intersectMkr)
{
    (*geno) = "";
    map<unsigned long, PAMA>::iterator positr;
    map<unsigned long, PAMA>::iterator positr_end;
    positr     = (*this_marker).begin();
    positr_end = (*this_marker).end();
    while(positr != positr_end)
    {
        unsigned long thispos = (*positr).first;
        if(thispos>=start && thispos<=end)
        {
            unsigned long this_dist           = thispos - start;
            unsigned long read_correspondence = 0 + this_dist;
            string        this_read_base      = (*updated_read).substr(read_correspondence, 1);            
            string        this_mat_base       = (*positr).second.matallele; // ref
            string        this_pat_base       = (*positr).second.patallele; // alt
            //
            // cout << "        : read " << this_read_base 
            //      << " ref "           << this_mat_base 
            //      << " alt "           << this_pat_base << endl;
            if(this_read_base.compare(this_mat_base) == 0)
            {
                (*geno) += "M"; // ref
            }else
            if(this_read_base.compare(this_pat_base) == 0)
            {
                (*geno) += "P"; // alt
            }else
            {
                (*geno) += "u";
            }  
            (*intersectMkr).insert(std::pair<unsigned long, PAMA>((*positr).first, (*positr).second));                           
        }else
        if(thispos>end) 
        {
            break;
        }
        else ;
        //
        positr ++;
    }
    //    
    return true;
}
//
bool align_read_to_ref(vector<char> operation, vector<int> count, string* orignal_read, string* updated_read)
{
    /*
        align the reads to reference according to cigar info:
        //          firstRefsMatch|
        //                        x
        //                        |4M  3D  6M  2I    11M    2D 4M  2d 3M 3S                           CIGAR
        // GCTATTTTGCGGAACAACGAATTCCTGGATCCACCA--GAAAATAAAGTTTGCTGCAGGACTTTTGCCAGCCATAAGCTTCGGTCAGGCT REF
        //                        CCTG---CCACCAGAGAAAATGAAGT--GTTGC--GACT                             READ
        //                        ||||DDD||||||II|||||||||||DD|R|||DD||||                             MM/INDELs
        //                        x         |                        |  x
        //          firstReadMatch|                                     |(secondReadMatch)

        The read would become: 
        //                        CCTG---CCACCA  GAAAATGAAGT--GTTGC--GACT                             READ        
                                               GA
        //
	https://samtools.github.io/hts-specs/SAMv1.pdf: 
        //
	Op	BAM	Description						Consumes_query	Consumes_reference	
	M	0	alignment match (can be a sequence match or mismatch)	yes		yes	
	I	1	insertion to the reference				yes		no	
	D	2	deletion from the reference				no		yes	
	N	3	skipped region from the reference			no		yes	
	S	4	soft clipping (clipped sequences present in	SEQ)	yes		no	
	H	5	hard clipping (clipped sequences NOT present in	SEQ)	no		no	
	P	6	padding (silent deletion from padded reference)		no		no	
	=	7	sequence match						yes		yes	
	X	8	sequence mismatch					yes		yes  	                                                       
    */
    (*updated_read) = "";         // delete 'SI', add 'MDX=' on real read; no action with 'H'; 'NP'
    unsigned long next_start = 0; // on real read
    for(int oi=0; oi<operation.size(); oi++)
    {
        if(operation[oi]=='H') 
        {
            // need not change read sequence
            (*updated_read) += "";
            next_start      += 0;
        }else
        if(operation[oi]=='S') 
        {
            // clip the read sequence 
            (*updated_read) += "";
            next_start      += count[oi];
        }else
        if(operation[oi]=='M' || operation[oi]=='X' || operation[oi]=='=') 
        {
            (*updated_read) += (*orignal_read).substr(next_start, count[oi]);
            next_start      += count[oi];
        }else
        if(operation[oi]=='D') 
        {
            // extend the read with "-"
            std::string s(count[oi], '-');            
            (*updated_read) += s;
            next_start      += 0;
        }else
        if(operation[oi]=='I') 
        {
            // delete the read subsequence - so update nothing
            (*updated_read) += "";
            next_start      += count[oi];
        }else
        {
            ; // N and P
        }            
    }    
    return true;
}
// read markers
bool read_marker(string fmarker,
                 map<string, map<unsigned long, PAMA> >* mkr)
{
    /*
      format: contig_id   position .   maternal paternal  optional
              utg000103lc 727428   .   A        T         0.2 RefCall .....
    */
    ifstream ifp;
    ifp.open(fmarker.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open marker file " << fmarker << endl;
        return false;
    }
    map<string, bool>::iterator cpitr;
    map<string, map<unsigned long, PAMA> >::iterator cmitr;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0]  =='#') continue;        
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 5) 
        {
            cout << "   Warning: unexpected line: " << line;
            continue;
        }
        string        thisctg   = lineinfo[0];
        unsigned long thispos   = strtoul(lineinfo[1].c_str(), NULL, 0);
        string        matallele = lineinfo[3]; // M - ref
        string        patallele = lineinfo[4]; // P - alt
        //
        PAMA tpama;
        tpama.matallele = matallele;
        tpama.patallele = patallele;
        //
        cmitr = (*mkr).find(thisctg);
        if(cmitr != (*mkr).end())
        {
            (*cmitr).second.insert(std::pair<unsigned long, PAMA>(thispos, tpama));
        }else
        {
            map<unsigned long, PAMA> posmkr;
            posmkr.insert(std::pair<unsigned long, PAMA>(thispos, tpama));
            (*mkr).insert(std::pair<string, map<unsigned long, PAMA> >(thisctg, posmkr));
        }
    }
    ifp.close();
    return true;
}
//
int decipher_cigar(string cigar, vector<char>* operation, vector<int>* count)
{
    /*
	https://samtools.github.io/hts-specs/SAMv1.pdf: 
        //
	Op	BAM	Description						Consumes_query	Consumes_reference	
	M	0	alignment match (can be a sequence match or mismatch)	yes		yes	
	I	1	insertion to the reference				yes		no	
	D	2	deletion from the reference				no		yes	
	N	3	skipped region from the reference			no		yes	
	S	4	soft clipping (clipped sequences present in	SEQ)	yes		no	
	H	5	hard clipping (clipped sequences NOT present in	SEQ)	no		no	
	P	6	padding (silent deletion from padded reference)		no		no	
	=	7	sequence match						yes		yes	
	X	8	sequence mismatch					yes		yes  	
	NOTE:
        // CIGAR alphabet: [0-9]+[MIDNSHPX=] and *: 23S23M1D8M1I16M1D29M38H 		
        // Sum of lengths of the MIS=X operations shall equal the length of SEQ in bam 
        // H can only be present as the first and/or last operation.
        // S may only have H operations between them and the ends of the CIGAR string
        // or mRNA-to-genome alignment, an N operation represents an intron.  	  
    */
    char *cstr = (char*)cigar.c_str(); // 
    string numstr("");
    int covlen = 0;
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
            (*operation).push_back(cstr[i]);    
            (*count).push_back(atoi(numstr.c_str()));
            //                    
            if(cstr[i] == 'M' || cstr[i] == 'D' || cstr[i] == 'N' || cstr[i] == 'X' || cstr[i] == '=')
            {
                covlen += atoi(numstr.c_str()); // positions of ref covered by read
            }
            numstr.clear();
        }
        if(cstr[i] == 'X' || cstr[i] == '=' || cstr[i] == 'P' || cstr[i] == 'N') 
        {            
            cout << "   Warning: you have special operation \'" << cstr[i] << "\' in Cigar: " << cigar << endl;
        }
    }
    // check
    if(true)
    {
        //cout << endl << "   check: cigar=" << cigar << endl;
        for(int ci = 0; ci < (*operation).size(); ci ++)
        {
            if(ci!=0 && ci!=(*operation).size()-1)
            {
                char tmpc = (*operation)[ci];
                if(tmpc=='S' || tmpc=='H')
                {
                    cout << "   Warning: \'" << tmpc << "\' operation happened in middle of alignment with cigar " 
                         << cigar            << endl;
                }            
            }
            // cout << "   check: operation " << (*operation)[ci] << "\t" << (*count)[ci] << endl;
        }
        // cout << "   check: reference length covered: " << covlen << " bp. " << endl; 
    }
    //
    return covlen;
}
/* 0.QNAME 1.FLAG 2.RNAME 3.POS 4.MAPQ 5.CIGAR 6.RNEXT 7.PNEXT 8.TLEN 9.SEQ 10.QUAL

 0.QNAME:   M05453:196:000000000-CGNKH:1:1101:16586:1168    
 1.FLAG:    163 
 2.REFNAME: refname    
 3.POS:     1   
 4.MAPQ:    42  
 5.CIGAR:   129M    
 6.RNEXT:   =   
 7.PNEXT:   208 
 8.REFLEN:  431 
 9.SEQ:     GGTCAAGGCAAGACGATATAACTGAACTCCGTTGTAGCATTAGAGCTGAAATGTTCTGTGGTTGAATTAATTTGTTTCTGGCAAATAATTAAAGTTGTTGCTGTTGGATTTACGTTGTAGGTATTTGGG   
10.QUAL:    GF9@EFFFCFGGGGG+,CFCCCFFFFFGCGED@GGCCF<ACE96CECG@<,CFF<FF@6EFC8FEGDFGFCFFCC@,CE@EE@F8EF9FFC<@,CEEFF8FFF:B<8AFFGGG,,BEGFGGFEAFFC?7   
 ...
*/
