/*
    Function: given     
        a list of contigs,
        a bam file with pair-end omnic short reads aligned,        
    extract pair-end short reads aligned to the contig markers.
    
    Written by: Hequan Sun
    Address...: Carl-von-linne-weg 10, 50829 Koeln (MPIPZ)
    Email.....: sunhequan@gmail.com/sun@mpipz.mpg.de
    Date......: 2022-10-05
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
bool verbose = false;
//
bool read_ctg(string ctg_file, map<string, unsigned long>* ctg_list);
bool create_sam_header(map<string, unsigned long> ctg_list, string* sam_header);
string get_R12(int thisflag_int);
int decipher_cigar(string cigar, vector<char>* operation, vector<int>* count);
//
int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        // g++ omnic_read_extracter.cpp split_string.cpp -O3 -o omnic_read_extracter
        cout << "\n   Function: extract pair-end short reads aligned to contigs in the given list.";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "samtools view -F 3840 omnic_read_alignment.bam | omnic_read_extracter - target_contig_sizes.txt out_prefix_str"
             << endl            << endl;
        cout << "   Note 1: bam must read name sorted. "             << endl;
        cout << "   Note 2: -F 3840 only check primary alignments. " << endl;
        cout << "   Note 3: output will be in a sam format. "        << endl;
        cout << "   Note 4: only if both read pairs aligned to the given contigs, they will be collected. " << endl
             << endl;
        return 1;
    }
    //
    string checkflag = (string)argv[1];
    if(checkflag.compare("-") != 0)
    {
        cout << "\n   Function: extract pair-end short reads aligned to contigs in the given list. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "samtools view -F 3840 omnic_read_alignment.bam | omnic_read_extracter - target_contig_sizes.txt out_prefix_str"
             << endl            << endl;
        cout << "   Error: \"omnic_read_extracter -\" must be provided. " << endl << endl;
        return 1;
    }
    //
    string ctg_file = (string)argv[2];
    string out_pref = (string)argv[3];
    // s1. read contig and sizes.
    map<string, unsigned long> ctg_list;
    if(!read_ctg(ctg_file, &ctg_list) )
    {
        cout << "   Error: reading ctg file failed. " << endl;
        return 1;
    }
    // s2. prepare sam file and header
    string sam_header = "";
    if(!create_sam_header(ctg_list, &sam_header) )
    {
        return 1;
    }
    string ofilename = out_pref + "_extract.sam";
    fstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    ofp << sam_header;
    // s3. extract reads into sam
    unsigned long lenNoCG  = 0; // no cigar string unmapped - length
    unsigned long numNoCG  = 0; // no cigar string unmapped - number
    string exampleNoCG     = "";
    unsigned long numraw             = 0;
    unsigned long num_pair_collected = 0;   
    unsigned long num_pair_inter_ctg = 0; 
    unsigned long num_pair_intra_ctg = 0;
    unsigned long num_pair_inter_ctg_p = 0;     
    unsigned long num_pair_skip      = 0;         
    std::string read_a; // first  row
    std::string read_b; // second row
    while (std::getline(std::cin, read_a)) 
    {
        if(read_a.size()==0) continue;
        if(read_a[0] == '@') continue;
        numraw += 1;         
        // get second read
        std::getline(std::cin, read_b);
        if(read_b.size()==0)
        {
            cout << "   Warning: you have read_a = " << read_a << " but read_b does not exist. " << endl;
            num_pair_skip ++;
            continue;
        }
        numraw += 1;         
        //
        if(numraw%1000000 == 0)
        {
            cout << "   Info: " << numraw << "th aligm..." << endl;
        }
        //
        vector<string> read_a_info = split_string(read_a, '\t');
        vector<string> read_b_info = split_string(read_b, '\t');
        if(read_a_info.size()<11 || read_b_info.size()<11)
        {
            cout << "   Warning: insufficient read_a info, skipped: " << read_a << endl;
            cout << "   Warning: insufficient read_b info, skipped: " << read_b << endl;
            num_pair_skip ++;
            continue;
        }        
        //
        string readname_a   = read_a_info[0]; // a read name
        string thisflag_a   = read_a_info[1]; // a flag value
        string this_a_cigar = read_a_info[5]; // a cigar string        
        string this_a_seq   = read_a_info[9]; // a read sequence         
        string readname_b   = read_b_info[0]; // b read name
        string thisflag_b   = read_b_info[1]; // b flag value   
        string this_b_cigar = read_b_info[5]; // b cigar string        
        string this_b_seq   = read_b_info[9]; // a read sequence
        int hexflag_a = strtol(thisflag_a.c_str(), NULL, 0); // a flag value 
        int hexflag_b = strtol(thisflag_b.c_str(), NULL, 0); // b flag value    
        // special case 1: cigar string as star: not aligned
        if(this_a_cigar.compare("*")==0 || this_b_cigar.compare("*")==0) 
        {
            //
            if(exampleNoCG.size()==0)
            {
                exampleNoCG = read_a;
                cout << endl
                     << "   Warning: there are paired-alignments without explicit CIGAR string, skipped., e.g.: " 
                     << read_a 
                     << endl
                     << read_b
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;
            }            
            numNoCG += 2; // individual reads; read pair number = numNoCG/2
            lenNoCG  += this_a_seq.size();
            lenNoCG  += this_b_seq.size();
            num_pair_skip ++;
            //
            continue;
        }        
        if(readname_a.compare(readname_b) != 0)
        {
            cout << "   Error: probably unsorted bam provided? Please sort it with read name. " << endl;
            cout << read_a << endl;
            cout << read_b << endl;
            return 1;
        }
        // read a: covered len by alignment
        vector<char>  operation_a;
        vector<int>   count_a;        
        int           covlen_a  = decipher_cigar(this_a_cigar, &operation_a, &count_a); 
        // read b: covered len by alignment
        vector<char>  operation_b;
        vector<int>   count_b;        
        int           covlen_b  = decipher_cigar(this_b_cigar, &operation_b, &count_b);        
        // check if read is R1 or R2
        string read_a_R12 = get_R12( hexflag_a );
        string read_b_R12 = get_R12( hexflag_b );
        if(verbose)
        cout << "   check: read a: " << read_a_R12 << ", and read b: " << read_b_R12 << endl;
        if(!(read_a_R12.compare("R1")==0 && read_b_R12.compare("R2")==0) &&
           !(read_a_R12.compare("R2")==0 && read_b_R12.compare("R1")==0) )
        {
            cout << "   Warning: unclear read pair info: " << endl;
            cout << "      read a = " << read_a << endl;
            cout << "      read b = " << read_b << endl;
            num_pair_skip ++;
            continue;
        }
        string this_a_ctg   = read_a_info[2];  // a reference contig name
        string this_b_ctg   = read_b_info[2];  // b reference contig name
        //
        if( ctg_list.find(this_a_ctg) != ctg_list.end() && ctg_list.find(this_b_ctg) != ctg_list.end() )
        {
            ofp << read_a << endl;
            ofp << read_b << endl;
            num_pair_collected ++;
            if(this_a_ctg.compare(this_b_ctg) == 0)
            {
                // intra-ctg
                size_t this_a_cigar_len = operation_a.size();
                size_t this_b_cigar_len = operation_b.size();                
                vector<char>::iterator aL_opr_itr = operation_a.begin();
                vector<char>::iterator aR_opr_itr = operation_a.end();                
                aR_opr_itr --;
                vector<int>::iterator aL_cnt_itr  = count_a.begin();
                vector<int>::iterator aR_cnt_itr  = count_a.end();                
                aR_cnt_itr --;        
                vector<char>::iterator bL_opr_itr = operation_b.begin();
                vector<char>::iterator bR_opr_itr = operation_b.end();                
                bR_opr_itr --;
                vector<int>::iterator bL_cnt_itr  = count_b.begin();
                vector<int>::iterator bR_cnt_itr  = count_b.end();                
                bR_cnt_itr --;                
                if( ( (*aL_opr_itr)=='S' && (*aL_cnt_itr)>=10 ) ||
                    ( (*aR_opr_itr)=='S' && (*aR_cnt_itr)>=10 ) ||
                    ( (*bL_opr_itr)=='S' && (*bL_cnt_itr)>=10 ) ||
                    ( (*bR_opr_itr)=='S' && (*bR_cnt_itr)>=10 )                                      
                  )
                {
                    // clipping                
                    num_pair_inter_ctg_p ++;
                }
                else
                {
                    // no clipping
                    num_pair_intra_ctg ++;
                }
            }else
            {
                // inter-ctg
                num_pair_inter_ctg ++;
            }            
        }else
        {
            num_pair_skip ++;
        }
    }
    //
    unsigned long num_pair_kept = num_pair_intra_ctg + num_pair_inter_ctg + num_pair_inter_ctg_p;
    cout << "   Info: summary on " << numraw               << " read alignment lines checked, among these, "     << endl
         << "         "            << num_pair_skip        << " read pairs skipped (/aligned to other contigs), "<< endl
         << "         "            << num_pair_kept        << " collected into given group: "            << endl
         << "         "            << num_pair_intra_ctg   << " with R1-R2 aligned within a contig, "    << endl
         << "         "            << num_pair_inter_ctg   << " with R1-R2 aligned to two contigs,"      << endl
         << "         "            << num_pair_inter_ctg_p << " with R1-R2 aligned to two contigs potentially. "<< endl;         
    //
    cout << "   Info: extract reads from bam/sam into linkage done. "     << endl;    
    //
    return 0;    
}
//
bool create_sam_header(map<string, unsigned long> ctg_list, string* sam_header)
{
    if(ctg_list.size()==0)
    {
        cout << "   Error: ctg list is empty." << endl;
        return false;
    }
    // initialize 
    std::stringstream ss;
    ss.str("");
    ss << "@HD\tVN:1.0\tSO:coordinate" << endl;    
    //
    map<string, unsigned long>::iterator citr;
    map<string, unsigned long>::iterator citr_end;
    citr     = ctg_list.begin();
    citr_end = ctg_list.end();
    while(citr != citr_end)
    {
        string this_ctg = (*citr).first;
        unsigned long this_ctg_size = (*citr).second;
        ss << "@SQ"                  << "\t" 
           << "SN:" << this_ctg      << "\t" 
           << "LN:" << this_ctg_size << endl;
        citr ++;
    }
    //
    ss << "@PG"        << "\t"
       << "ID:bowtie2" << "\t"
       << "PN:bowtie2" << "\t"
       << "VN:2.4.4"   << "\t"
       << "CL:\"/opt/share/software/packages/bowtie2-2.4.4/bin/bowtie2-align-s --wrapper basic-0 -p 20"
       << " -x clipped2_A_hifiasm.p_utg.gfa.fa -1 ln_A_BGI_DNB_R1.fq.gz -2 ln_A_BGI_DNB_R2.fq.gz\""
       << endl;
    *sam_header = ss.str();
    //
    return true;
}
//
bool read_ctg(string ctg_file, map<string, unsigned long>* ctg_list)
{
    fstream ifp;
    ifp.open(ctg_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open ctg file " << ctg_file << endl;
        return false;
    }
    //
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 2)
        {
            cout << "   Warning: skipping line with insufficient info: " << line << endl;
            continue;
        }
        string this_ctg_id = lineinfo[0];
        unsigned long this_ctg_size = strtoul(lineinfo[1].c_str(), NULL, 0);
        if( (*ctg_list).find(this_ctg_id) != (*ctg_list).end() )
        {
            if( (*ctg_list)[this_ctg_id] != this_ctg_size)
            {
                cout << "   Warning: different sizes with the same contig id, please check: " << this_ctg_id << endl;
            }else
            {
                cout << "   Warning: repeated contig id: " << this_ctg_id << endl;
            }
        }else
        {
            (*ctg_list).insert(std::pair<string, unsigned long>(this_ctg_id, this_ctg_size));
        }
    }
    //
    return true;
}
//
string get_R12(int thisflag_int)
{
    /* function: check if read is R1 or R2*/
    string read_R12 = "RU";
    if((thisflag_int & 0x40) == 0x40 && (thisflag_int & 0x80) == 0)
    {
        read_R12 = "R1";
    }
    else
    if((thisflag_int & 0x40) == 0    && (thisflag_int & 0x80) == 0x80)
    {
        read_R12 = "R2";
    }
    else
    {
        read_R12 = "RU";
        cout << "   Warning: unclear read info: "    << endl;
        cout << "      flag = "      << thisflag_int << endl;
        cout << "      flag & 0x40 == 0x40 = " << (thisflag_int & 0x40) << endl;
        cout << "      flag & 0x80 == 0    = " << (thisflag_int & 0x80) << endl;            
    } 
    //
    return read_R12;
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
        }else
        if(cstr[i] == 'H') 
        {            
            cout << "   Warning: hard clipping occurred: \'"    << cstr[i] << "\' in Cigar: " << cigar << endl;
        }        
    }
    //
    return covlen;
}
                               
