/*
    Function: given     
        a list of linkage-group wise window markers,
        a bam file with pair-end short reads aligned,        
    separate pair-end short reads into ploidy*4 clusters (as indicated by window marker info).    
    
    Written by: Hequan Sun
    Address...: Carl-von-linne-weg 10, 50829 Koeln (MPIPZ)
    Email.....: sunhequan@gmail.com/sun@mpipz.mpg.de
    Date......: 2022-11-09
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
#include      <cerrno> // Include the errno header file
#include "split_string.h"
#include   "./gzlib/gzstream.h"
using namespace std;
// window marker
struct winMARKER
{
//  string         ctg; // ctd id
//  unsigned long  sta; // window marker start 
    unsigned long  end; // window marker end
    string        type; // hap/dip/trip/tetrap/rep
    vector<string> wlg; // linkage group id 1:4
    string       homlg; // homo-lg id...... 1:12
};
bool verbose = false;
bool complex_readname = true; // this is needed as some tools does not recognize read pairs with complex names
//
bool read_gb_window_marker_grouping(string                 gbmarker_file, 
         map<string, map<unsigned long, winMARKER > >*     gbmarker,
         map<string, map<string, string> >*                ctg2lg,
         map<string, string>*                              homLGs);   
bool compare_lgs(vector<string> this_lgs, vector<string>   last_lgs);         
bool merge_gb_window_marker_grouping(
         map<string, map<unsigned long, winMARKER > >      gbmarker,
         map<string, map<unsigned long, winMARKER > >*     gbmarker_updated);                                                  
int decipher_cigar(string                                  cigar, 
         vector<char>*                                     operation, 
         vector<int>*                                      count);
bool align_read_to_ref(vector<char>                        operation, 
         vector<int>                                       count, 
         string*                                           orignal_read, 
         string*                                           updated_read); 
bool overlap_read_and_win_marker(map<unsigned long, unsigned long> this_win_marker, 
         unsigned long                                     read_sta, 
         unsigned long                                     read_end,
         string                                            this_ctg,
         map<unsigned long, unsigned long>*                this_win_marker_overlapped);
unsigned long overlap_read_and_win_marker_v2(unsigned long this_win_sta, 
         unsigned long                                     this_win_end, 
         unsigned long                                     read_sta, 
         unsigned long                                     read_end,
         string                                            this_ctg);   
string get_R12(int                                         thisflag_int);         
bool get_target_lg_R12(vector<string>                      best_overlap_lg_a, 
         string                                            best_overlap_type_a,
         vector<string>                                    best_overlap_lg_b, 
         string                                            best_overlap_type_b,
         vector<string>*                                   common_lg_togo);
//
int main(int argc, char* argv[])
{
    if(argc < 5)
    {
        // g++ short_read_separator.cpp split_string.cpp -O3 -o short_read_separator
        cout << "\n   Function: separate pair-end short reads using phased contig markers. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "samtools view -F 3840 read_alignment.bam | hic_short_read_separator - win_marker.txt out_prefix_str informative_readname "
             << endl            << endl;
        cout << "   Note 1: bam must read name sorted. " << endl;
        cout << "   Note 2: -F 3840 only check primary alignments. " << endl;
        cout << "   Note 3: informative_readname is 0: give simple original read name; otherwise, given read name with alignment info to markers. "
             << endl;
        return 1;
    }
    double startT= clock();
    //
    string checkflag      = (string)argv[1];
    if(checkflag.compare("-") != 0)
    {
        cout << "\n   Function: separate pair-end short reads using phased contig markers. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "samtools view -F 3840 read_alignment.bam | hic_short_read_separator - win_marker.txt out_prefix_str "
             << endl            << endl;
        cout << "   Error: \"short_read_separator -\" must be provided. " << endl << endl;
        return 1;
    }
    string gbmarker_file  = (string)argv[2];
    string outprefix      = (string)argv[3];
    string complex_readname_singal = (string)argv[4];
    if(complex_readname_singal.compare("0") != 0)
    { 
        // non-0
        cout << "   Info: read names would be more informative with alignment info to markers as asked." << endl;                
    }else
    {
        complex_readname = false;    // 0
        cout << "   Info: read names would be original/raw as asked." << endl;            
    }
    //
    cout << "   Info: reading gamete_binning grouped window marker info.. " << endl;
    map<string, map<unsigned long, winMARKER > > gbmarker;// <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
    map<string, map<string, string> >            ctg2lg;  // <ctg_id, <lg_id=1:4, 1:12> >
    map<string, string>                          homLGs;  // <lg_id=1:4,  homLG_id=1:12>
    if(!read_gb_window_marker_grouping(gbmarker_file, &gbmarker, &ctg2lg,& homLGs))
    {
        return 1;
    }
    cout << "   Info: reading gamete_binning grouped window marker done. " << endl;    
    //
    cout << "   Info: merging gamete_binning grouped window marker info.. " << endl;
    map<string, map<unsigned long, winMARKER > > gbmarker_updated;
    if(!merge_gb_window_marker_grouping(gbmarker, &gbmarker_updated))
    {
        cout << "   Error: failed in merging window markers. " << endl;
        return 1;
    }
    bool stopnow=false;  
    if(stopnow)
    {
        cout << "    Testing: returned. " << endl;
        return 1;
    }    
    cout << "   Info: merging gamete_binning grouped window marker info done. " << endl;
    cout << "   Info: extract read from bam/sam into linkage groups... " << endl;
    // output files
    // create an intermediate folder for collecting details about a CO-molecule
    string tmpfolder = outprefix + "_window_marker_separated_reads";
    DIR* dir = opendir(tmpfolder.c_str());
    if (dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        const int dir_err = mkdir(tmpfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << tmpfolder << endl;
            return 1;
        }
    }
    else;      
    // file and variable preparation
    // file handles on output files: <"xxLG_R1/2.fq.gz", *ofp >    
    map<string, ogzstream*> allofp;         
    // collecting group-specific unique read names: <"yyLG_R1/2.fq.gz", <pbreadname, cnt> >               
    map<string, unsigned long> groupedillureadnames;  // pair-end read count in each linkage group <lg, read1/2_cnt>
    // no real data - just for initializing groupedillureadnames.    
    map<string, int> init_tmp;     
    map<string, string>::iterator grpfileitr;
    map<string, string>::iterator grpfileitr_end;
    grpfileitr     = homLGs.begin(); 
    grpfileitr_end = homLGs.end();
    while(grpfileitr != grpfileitr_end)
    {
        string this_lg    = (*grpfileitr).first; // -1, 1:4
        string this_homLG = (*grpfileitr).second;// 1:12
        if(this_lg.compare("-1") == 0)
        {
             grpfileitr ++;
             continue;
        }
        string this_lg_readfile("");        
        // R1
        this_lg_readfile = "homLG_" + this_homLG + "_LG_" + this_lg + "_reads_R1.fq.gz";
        if(verbose)        
        cout << "   info: " << this_lg_readfile << " to be created. " << endl;
        ogzstream *ofp_R1 = new ogzstream( (tmpfolder + "/" + outprefix + "_" + "" + this_lg_readfile).c_str(), ios::out );
        if(!(*ofp_R1).good())
        {
            cout << "   Error: cannot set up output file for " << tmpfolder + "/" + outprefix + "_" + "" + this_lg_readfile << endl;
            return 1;
        }
        allofp.insert(std::pair<string, ogzstream*>(this_lg_readfile, ofp_R1));
        groupedillureadnames.insert(std::pair<string, unsigned long>(this_lg_readfile, 0));
        // R2
        this_lg_readfile = "homLG_" + this_homLG + "_LG_" + this_lg + "_reads_R2.fq.gz";
        if(verbose)        
        cout << "   info: " << this_lg_readfile << " to be created. " << endl;        
        ogzstream *ofp_R2 = new ogzstream( (tmpfolder + "/" + outprefix + "_" + "" + this_lg_readfile).c_str(), ios::out);
        if(!(*ofp_R2).good())
        {
            cout << "   Error: cannot set up output file for " << tmpfolder + "/" + outprefix + "_" + "" + this_lg_readfile << endl;
            return 1;
        }
        allofp.insert(std::pair<string, ogzstream*>(this_lg_readfile, ofp_R2));
        groupedillureadnames.insert(std::pair<string, unsigned long>(this_lg_readfile, 0));        
        //
        grpfileitr ++;
    } 
    // short reads collected    
    map<string, int> illureadnames;                    // all collected illu read names (unique of 7043689)
    map<string, vector<string> > illureadnames_detail; // all collected illu read names (unique of 7043689) <readname, <grp, align_len> >   
    map<string, int> illureadnames_nolg;     // all collected illu read names without linkage group info (unique of )
    map<string, int> illureadnames_lg;       // all collected illu read names with    linkage group info and with    marker (unique of )
    map<string, int> illureadnames_lg_nomkr; // all collected illu read names with    linkage group info and without marker (unique of )    
    // unmapped pair-end short reads, any of them with cigar as "*"
    string unmappedfile_R1   = tmpfolder + "/" + outprefix + "_" + "unmapped_starcigar_illu_reads_R1.fq.gz"; 
    ogzstream *ofp_R1 = new ogzstream(unmappedfile_R1.c_str(), ios::out);
    allofp.insert(std::pair<string, ogzstream*>("unmapped_starcigar_illu_reads_R1.fq.gz", ofp_R1));
    groupedillureadnames.insert(std::pair<string, unsigned long>("unmapped_starcigar_illu_reads_R1.fq.gz", 0));  
    //
    string unmappedfile_R2   = tmpfolder + "/" + outprefix + "_" + "unmapped_starcigar_illu_reads_R2.fq.gz"; 
    ogzstream *ofp_R2 = new ogzstream(unmappedfile_R2.c_str(), ios::out);
    allofp.insert(std::pair<string, ogzstream*>("unmapped_starcigar_illu_reads_R2.fq.gz", ofp_R2));
    groupedillureadnames.insert(std::pair<string, unsigned long>("unmapped_starcigar_illu_reads_R2.fq.gz", 0));        
    // mapped paired-end short reads, any of them without linkage group info
    string mappednolgfile_R1 = tmpfolder + "/" + outprefix + "_" + "mapped_no_lg_illu_reads_R1.fq.gz"; 
    ogzstream *nolgofp_R1 = new ogzstream(mappednolgfile_R1.c_str(), ios::out);
    allofp.insert(std::pair<string, ogzstream*>("mapped_no_lg_illu_reads_R1.fq.gz", nolgofp_R1));  
    groupedillureadnames.insert(std::pair<string, unsigned long>("mapped_no_lg_illu_reads_R1.fq.gz", 0)); 
    //
    string mappednolgfile_R2 = tmpfolder + "/" + outprefix + "_" + "mapped_no_lg_illu_reads_R2.fq.gz"; 
    ogzstream *nolgofp_R2 = new ogzstream(mappednolgfile_R2.c_str(), ios::out);
    allofp.insert(std::pair<string, ogzstream*>("mapped_no_lg_illu_reads_R2.fq.gz", nolgofp_R2));  
    groupedillureadnames.insert(std::pair<string, unsigned long>("mapped_no_lg_illu_reads_R2.fq.gz", 0));              
    // mapped paired-end short reads, but reads aligned to different linkage groups -- not informative for scaffolding
    string mappedinterlgfile_R1 = tmpfolder + "/" + outprefix + "_" + "mapped_inter_lg_illu_reads_R1.fq.gz"; 
    ogzstream *interlgofp_R1 = new ogzstream(mappedinterlgfile_R1.c_str(), ios::out);
    allofp.insert(std::pair<string, ogzstream*>("mapped_inter_lg_illu_reads_R1.fq.gz", interlgofp_R1));  
    groupedillureadnames.insert(std::pair<string, unsigned long>("mapped_inter_lg_illu_reads_R1.fq.gz", 0)); 
    //
    string mappedinterlgfile_R2 = tmpfolder + "/" + outprefix + "_" + "mapped_inter_lg_illu_reads_R2.fq.gz"; 
    ogzstream *interlgofp_R2 = new ogzstream(mappedinterlgfile_R2.c_str(), ios::out);
    allofp.insert(std::pair<string, ogzstream*>("mapped_inter_lg_illu_reads_R2.fq.gz", interlgofp_R2));  
    groupedillureadnames.insert(std::pair<string, unsigned long>("mapped_inter_lg_illu_reads_R2.fq.gz", 0));                    
    //
    unsigned long lenNoCG  = 0; // no cigar string unmapped - length
    unsigned long numNoCG  = 0; // no cigar string unmapped - number
    string exampleNoCG     = "";
    unsigned long numNoMA  = 0; // number of reads showing multiple alignment: only one kept.
    string exampleNoMA     = "";
    unsigned long numraw   = 0;    
    string exampleNoSEQ    = "";
    unsigned long numNoSEQ = 0;
    unsigned long lengrped = 0; // illu alignments clearly related to linkage groups - length 
    unsigned long grpedNo  = 0; // illu alignments clearly related to linkage groups; note 1 read may have >1 alignments
    unsigned long lennonlg = 0; // illu alignments not related to linkage groups - length       
    unsigned long nonlgNo  = 0; // illu alignments not related to linkage groups - number 
    // m64078_200122_223326/73138406/ccs   16  utg000001l_pilon    458658  60  9555M1D5579M3S  *   0   0   GTGAAAGTTGTTGTCGTCGTTCTAGCTCATCG    
    srand (2); // to randomly separate reads if the marker relates to several linkage groups of 1:4
    std::string read_a; // first  row 
    std::string read_b; // second row
    while (std::getline(std::cin, read_a)) 
    {
        if(read_a.size()==0) continue;
        if(read_a[0] == '@') continue;
        // get second read 
        std::getline(std::cin, read_b);
        if(read_b.size()==0)
        {
            cout << "   Warning: you have read_a = " << read_a << " but read_b does not exist. " << endl;
            continue;
        }
        //
        if(verbose)
        {
            cout << endl;
            cout << "   check: read_a = " << read_a << endl;
            cout << "          read_b = " << read_b << endl;
        }
        //
        numraw += 2; 
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
            continue;
        }        
        //
        string readname_a  = read_a_info[0]; // a read name
        string thisflag_a  = read_a_info[1]; // a flag value
        string readname_b  = read_b_info[0]; // b read name
        string thisflag_b  = read_b_info[1]; // b flag value   
        int hexflag_a = strtol(thisflag_a.c_str(), NULL, 0); // a flag value 
        int hexflag_b = strtol(thisflag_b.c_str(), NULL, 0); // b flag value    
        if(readname_a.compare(readname_b) != 0)
        {
            cout << "   Error: probably unsorted bam provided? Please sort it with read name. " << endl;
            return 1;
        }    
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
            continue;
        }
        /*
		1	0x1	template having multiple segments in sequencing
		2	0x2	each segment properly aligned according to the aligner
		4	0x4	segment unmapped
		8	0x8	next segment in the template unmapped
		16	0x10	SEQ being reverse complemented
		32	0x20	SEQ of the next segment in the template being reverse complemented
		64	0x40	the first segment in the template
		128	0x80	the last segment in the template
		256	0x100	secondary alignment
		512	0x200	not passing filters, such as platform/vendor quality controls
		1024	0x400	PCR or optical duplicate
		2048	0x800	supplementary alignment        
        */
        string this_a_ctg   = read_a_info[2];  // a reference contig name
        string this_a_pos   = read_a_info[3];  // a first matching reference coordinate
        string this_a_cigar = read_a_info[5];  // a cigar string
        string this_a_seq   = read_a_info[9];  // a read sequence 
        string this_a_qual  = read_a_info[10]; // a read sequence quality
        //
        string this_b_ctg   = read_b_info[2];  // b reference contig name
        string this_b_pos   = read_b_info[3];  // b first matching reference coordinate
        string this_b_cigar = read_b_info[5];  // b cigar string
        string this_b_seq   = read_b_info[9];  // b read sequence   
        string this_b_qual  = read_b_info[10]; // b read sequence quality
        // special case 1: cigar string as star: not aligned
        if(this_a_cigar.compare("*")==0 || this_b_cigar.compare("*")==0) 
        {
            // READ SET I
            if(read_a_R12.compare("R1") == 0)
            {
                (*allofp["unmapped_starcigar_illu_reads_R1.fq.gz"])  << "@" << readname_a  << endl;
                (*allofp["unmapped_starcigar_illu_reads_R1.fq.gz"])  <<        this_a_seq  << endl;
                (*allofp["unmapped_starcigar_illu_reads_R1.fq.gz"])  << "+"                << endl;                
                (*allofp["unmapped_starcigar_illu_reads_R1.fq.gz"])  <<        this_a_qual << endl;                
                //
                (*allofp["unmapped_starcigar_illu_reads_R2.fq.gz"])  << "@" << readname_b  << endl;
                (*allofp["unmapped_starcigar_illu_reads_R2.fq.gz"])  <<        this_b_seq  << endl;  
                (*allofp["unmapped_starcigar_illu_reads_R2.fq.gz"])  << "+"                << endl;                
                (*allofp["unmapped_starcigar_illu_reads_R2.fq.gz"])  <<        this_b_qual << endl;                                
            }else
            {
                (*allofp["unmapped_starcigar_illu_reads_R1.fq.gz"])  << "@" << readname_b  << endl;
                (*allofp["unmapped_starcigar_illu_reads_R1.fq.gz"])  <<        this_b_seq  << endl;
                (*allofp["unmapped_starcigar_illu_reads_R1.fq.gz"])  << "+"                << endl;                
                (*allofp["unmapped_starcigar_illu_reads_R1.fq.gz"])  <<        this_b_qual << endl;                  
                //
                (*allofp["unmapped_starcigar_illu_reads_R2.fq.gz"])  << "@" << readname_a  << endl;
                (*allofp["unmapped_starcigar_illu_reads_R2.fq.gz"])  <<        this_a_seq  << endl; 
                (*allofp["unmapped_starcigar_illu_reads_R2.fq.gz"])  << "+"                << endl;                
                (*allofp["unmapped_starcigar_illu_reads_R2.fq.gz"])  <<        this_a_qual << endl;                           
            }
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
            //
            /* 20201024 will use too much mem for hundreds of millions of short reads!
            if(groupedillureadnames["unmapped_starcigar_illu_reads.fq.gz"].find(readname_a) == 
               groupedillureadnames["unmapped_starcigar_illu_reads.fq.gz"].end())
            {
                groupedillureadnames["unmapped_starcigar_illu_reads.fq.gz"].insert(std::pair<string, int>(readname_a, 1));
            }else
            {
                groupedillureadnames["unmapped_starcigar_illu_reads.fq.gz"][readname_a] += 1;
            }
            */
            groupedillureadnames["unmapped_starcigar_illu_reads_R1.fq.gz"] += 1;
            groupedillureadnames["unmapped_starcigar_illu_reads_R2.fq.gz"] += 1;
            //     
            continue;
        }
        // special case 2: skip secondary/
        // note: as "samtools view -F 3840" applied to bam, no non-primary alignment would get here
        if(hexflag_a >= 256 || hexflag_b >= 256)
        {
            numNoMA += 2;
            if(exampleNoMA.size()==0)
            {
                exampleNoMA = read_a;
                cout << endl
                     << "   Warning: there are alignments being secondary/supplementary etc, skipped, e.g.: " 
                     << read_a 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;                
            }
            continue;
        } 
        // collect all informative reads
        /* 20201024 will use too much mem for hundreds of millions of short reads!        
        if(illureadnames.find(readname_a) == illureadnames.end())
        {
            illureadnames.insert(std::pair<string, int>(readname_a, 1));
        }else
        {
            illureadnames[readname_a] += 1;
        }
        */
        string aligned_grp_a = "";
        if(verbose)        
        cout << "   check: read: " << readname_a << " " << read_a_R12 << " aligned to " << this_a_ctg << "\t" << this_a_pos << endl;
        string aligned_grp_b = "";
        if(verbose)        
        cout << "   check: read: " << readname_b << " " << read_b_R12 << " aligned to " << this_b_ctg << "\t" << this_b_pos << endl;                             
        // If it contains 0x10, SEQ in BAM/SAM has been reversed+complemented from read in fastq
        string read_a_reversed = "nrc";
        if((hexflag_a & 0x10) == 0x10) read_a_reversed = "rc";
        string read_b_reversed = "nrc";
        if((hexflag_a & 0x10) == 0x10) read_b_reversed = "rc";         
        // new 20210824
        string this_a_seq_rc("");
        string this_a_qual_rc("");
        if(read_a_reversed.compare("rc")==0)
        {
            this_a_seq_rc = this_a_seq;                    
            /* reverse sequence    */
            std::reverse(this_a_seq_rc.begin(), this_a_seq_rc.end());
            /* lowercase sequence  */
            std::transform(this_a_seq_rc.begin(), this_a_seq_rc.end(),this_a_seq_rc.begin(), ::tolower);
            /* complement sequence */
            std::replace(this_a_seq_rc.begin(), this_a_seq_rc.end(), 'a', 'T');
            std::replace(this_a_seq_rc.begin(), this_a_seq_rc.end(), 'c', 'G'); // caution1: here c:G
            std::replace(this_a_seq_rc.begin(), this_a_seq_rc.end(), 'g', 'C'); // caution2: here g:C
            std::replace(this_a_seq_rc.begin(), this_a_seq_rc.end(), 't', 'A'); // if caution 1&2 is not consistent,
            std::replace(this_a_seq_rc.begin(), this_a_seq_rc.end(), 'u', 'A'); // they will result in c:G+G:C=c:C
            //
            std::transform(this_a_seq_rc.begin(), this_a_seq_rc.end(),this_a_seq_rc.begin(), ::toupper);
            //
            // also reverse quality
            this_a_qual_rc = this_a_qual;
            std::reverse(this_a_qual_rc.begin(), this_a_qual_rc.end());
        }
        string this_b_seq_rc("");
        string this_b_qual_rc("");
        if(read_b_reversed.compare("rc")==0)
        {
            this_b_seq_rc = this_b_seq;                    
            /* reverse sequence    */
            std::reverse(this_b_seq_rc.begin(), this_b_seq_rc.end());
            /* lowercase sequence  */
            std::transform(this_b_seq_rc.begin(), this_b_seq_rc.end(),this_b_seq_rc.begin(), ::tolower);
            /* complement sequence */
            std::replace(this_b_seq_rc.begin(), this_b_seq_rc.end(), 'a', 'T');
            std::replace(this_b_seq_rc.begin(), this_b_seq_rc.end(), 'c', 'G'); // caution1: here c:G
            std::replace(this_b_seq_rc.begin(), this_b_seq_rc.end(), 'g', 'C'); // caution2: here g:C
            std::replace(this_b_seq_rc.begin(), this_b_seq_rc.end(), 't', 'A'); // if caution 1&2 is not consistent,
            std::replace(this_b_seq_rc.begin(), this_b_seq_rc.end(), 'u', 'A'); // they will result in c:G+G:C=c:C
            //
            std::transform(this_b_seq_rc.begin(), this_b_seq_rc.end(),this_b_seq_rc.begin(), ::toupper);
            //
            // also reverse quality
            this_b_qual_rc = this_b_qual;
            std::reverse(this_b_qual_rc.begin(), this_b_qual_rc.end());
        }          
        // new 20210824               
        // read a: covered len by alignment
        vector<char>  operation_a;
        vector<int>   count_a;        
        int           covlen_a  = decipher_cigar(this_a_cigar, &operation_a, &count_a); 
        // read b: covered len by alignment
        vector<char>  operation_b;
        vector<int>   count_b;        
        int           covlen_b  = decipher_cigar(this_b_cigar, &operation_b, &count_b);                   
        // a alignment start and end position on ref seq: these would intersect interested markers for current read
        unsigned long firstRefsMatch_a = strtoul(this_a_pos.c_str(), NULL, 0);                       
        unsigned long spansecond_a     = firstRefsMatch_a+covlen_a-1; // 1-based 
        if(verbose)                
        cout << "          span: " << firstRefsMatch_a              << "-"   << spansecond_a 
             << " => "             << spansecond_a-firstRefsMatch_a+1 << " bp" << endl; // 1-based     
        // b alignment start and end position on ref seq: these would intersect interested markers for current read
        unsigned long firstRefsMatch_b = strtoul(this_b_pos.c_str(), NULL, 0);                       
        unsigned long spansecond_b     = firstRefsMatch_b+covlen_b-1; // 1-based 
        if(verbose)                
        cout << "          span: " << firstRefsMatch_b              << "-"   << spansecond_b 
             << " => "             << spansecond_b-firstRefsMatch_b+1 << " bp" << endl; // 1-based                   
        // a/b read updated with operation_a/b
        string updated_read_a("");
        string updated_read_b("");        
        if(!align_read_to_ref(operation_a, count_a, &this_a_seq, &updated_read_a) ||
           !align_read_to_ref(operation_b, count_b, &this_b_seq, &updated_read_b) )
        {
            cout << "   Warning: unexpected read " << this_a_seq << ", or ";
            cout << "            unexpected read " << this_b_seq << endl;            
            continue;
        }
        /*
        cout << "   check: original read in bam/sam: " 
             << this_a_seq.size()      << " bp; " << endl;
        cout << "   check: updated read with \'I\' removed, and \'D:-\' added: " 
             << updated_read_a.size() << " bp. " << endl;
        cout << "   check: original read in bam/sam: " 
             << this_b_seq.size()      << " bp; " << endl;
        cout << "   check: updated read with \'I\' removed, and \'D:-\' added: " 
             << updated_read_b.size() << " bp. " << endl;
        */
        // map<string, map<unsigned long, winMARKER > > gbmarker_updated; // <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
        // map<string, map<string, string> >            ctg2lg; // <ctg_id, <lg_id=1:4, 1:12> >
        // map<string, string>                          homLGs; // <lg_id=-1,1:4,  homLG_id=1:12>
        assert( gbmarker_updated.find(this_a_ctg) != gbmarker_updated.end() );
        assert( gbmarker_updated.find(this_b_ctg) != gbmarker_updated.end() );  
        // find best matching linkage group for read a      
        map<unsigned long, winMARKER > this_ctg_win_list_a = gbmarker_updated[this_a_ctg]; // can be several large windows
        vector<string> best_overlap_lg_a;
        unsigned long  best_overlap_len_a  = 0;
        string         best_overlap_type_a = "";
        map<unsigned long, winMARKER >::iterator winitr_a;
        map<unsigned long, winMARKER >::iterator winitr_a_end;
        winitr_a     = this_ctg_win_list_a.begin();
        winitr_a_end = this_ctg_win_list_a.end();
        while(winitr_a != winitr_a_end)
        {
            unsigned long this_win_sta = (*winitr_a).first;
            winMARKER tmp_pos_marker   = (*winitr_a).second;
            unsigned long this_win_end = tmp_pos_marker.end;                        
            //
            if(this_win_sta > spansecond_a)
            {
                break;
            }
            if(this_win_end < firstRefsMatch_a)
            {
                winitr_a ++;
                continue;
            }
            //
            unsigned long tmp_overlap_len = overlap_read_and_win_marker_v2(this_win_sta, 
                                                                           this_win_end, 
                                                                           firstRefsMatch_a, 
                                                                           spansecond_a,
                                                                           this_a_ctg);   
            if(tmp_overlap_len > best_overlap_len_a)
            {
                best_overlap_len_a  = tmp_overlap_len;
                best_overlap_lg_a   = tmp_pos_marker.wlg; // vector<lg>
                best_overlap_type_a = tmp_pos_marker.type;
                if(verbose)                
                cout << "   check: best matching window marker updated as " << this_win_sta 
                     << "-"                                                 << this_win_end
                     << " for read a. "                                     << endl; 
            }
            //
            winitr_a ++;
        }
        // find best matching linkage group for read b      
        map<unsigned long, winMARKER > this_ctg_win_list_b = gbmarker_updated[this_b_ctg]; // can be several large windows
        vector<string> best_overlap_lg_b;
        unsigned long  best_overlap_len_b  = 0;
        string         best_overlap_type_b = "";
        map<unsigned long, winMARKER >::iterator winitr_b;
        map<unsigned long, winMARKER >::iterator winitr_b_end;
        winitr_b     = this_ctg_win_list_b.begin();
        winitr_b_end = this_ctg_win_list_b.end();
        while(winitr_b != winitr_b_end)
        {
            unsigned long this_win_sta = (*winitr_b).first;
            winMARKER tmp_pos_marker   = (*winitr_b).second;
            unsigned long this_win_end = tmp_pos_marker.end;                        
            //
            if(this_win_sta > spansecond_b)
            {
                break;
            }
            if(this_win_end < firstRefsMatch_b)
            {
                winitr_b ++;
                continue;
            }
            //
            unsigned long tmp_overlap_len = overlap_read_and_win_marker_v2(this_win_sta, 
                                                                           this_win_end, 
                                                                           firstRefsMatch_b, 
                                                                           spansecond_b,
                                                                           this_b_ctg);   
            if(tmp_overlap_len > best_overlap_len_b)
            {
                best_overlap_len_b  = tmp_overlap_len;
                best_overlap_lg_b   = tmp_pos_marker.wlg;
                best_overlap_type_b = tmp_pos_marker.type;
                if(verbose)                
                cout << "   check: best matching window marker updated as " << this_win_sta 
                     << "-"                                                 << this_win_end
                     << " for read b. "                                     << endl; 
            }
            //
            winitr_b ++;
        }  
        // non-signal reads
        if(best_overlap_type_a.compare("rep")==0 || best_overlap_type_b.compare("rep")==0)
        {
            string this_lg_readfile_a = "mapped_no_lg_illu_reads_" + read_a_R12 + ".fq.gz";
            if(complex_readname)
            {
                (*allofp[this_lg_readfile_a]) << "@"   
                                   << readname_a 
                                   << " "
                                   << read_a_reversed                                   
                                   << " + lg="
                                   << "-1/48"
                                   << " homLG="
                                   << "-1/12"
                                   << " rep.span="         
                                   << this_a_ctg 
                                   << ":"   
                                   << firstRefsMatch_a   
                                   << "-"          
                                   << spansecond_a
                                   << endl;   // read name  
            }else
            {
                (*allofp[this_lg_readfile_a]) << "@"   
                                   << readname_a 
                                   << " "
                                   << read_a_reversed
                                   << endl;   // read name                  
            }
            (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;   // read seq
            (*allofp[this_lg_readfile_a]) << "+" << endl;
            (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;                                   
            aligned_grp_a = this_lg_readfile_a;              
            //                                   
            string this_lg_readfile_b = "mapped_no_lg_illu_reads_" + read_b_R12 + ".fq.gz";
            if(complex_readname)
            {
                (*allofp[this_lg_readfile_b]) << "@"   
                                   << readname_b 
                                   << " "
                                   << read_b_reversed                                   
                                   << " + lg="
                                   << "-1/48"
                                   << " homLG="
                                   << "-1/12"
                                   << " rep.span="         
                                   << this_b_ctg 
                                   << ":"   
                                   << firstRefsMatch_b   
                                   << "-"          
                                   << spansecond_b
                                   << endl;   // read name  
            }else
            {
                (*allofp[this_lg_readfile_b]) << "@"   
                                   << readname_b 
                                   << " "
                                   << read_b_reversed                                   
                                   << endl;   // read name              
            }
            (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
            (*allofp[this_lg_readfile_b]) << "+" << endl;
            (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl;                                            
            aligned_grp_b = this_lg_readfile_b;              
            //                                   
            nonlgNo  += 2;           
            lennonlg += this_a_seq.size();  
            lennonlg += this_b_seq.size();              
            // without lg
            /* 20201024 will use too much mem for hundreds of millions of short reads!            
            if(illureadnames_nolg.find(readname_a) == illureadnames_nolg.end())
            {
                illureadnames_nolg.insert(std::pair<string,int>(readname_a, 1));
            } 
            */
            //
            /* 20201024 will use too much mem for hundreds of millions of short reads!            
            if(groupedillureadnames[this_lg_readfile].find(readname_a) == 
               groupedillureadnames[this_lg_readfile].end())
            {
                groupedillureadnames[this_lg_readfile].insert(std::pair<string, int>(readname_a, 1));
            }else
            {
                groupedillureadnames[this_lg_readfile][readname_a] += 1;
            }
            */
            groupedillureadnames[this_lg_readfile_a] += 1;
            groupedillureadnames[this_lg_readfile_b] += 1;            
            //                                                           
            continue; // caution: next read from bam                       
        }
        // linkage groups of signal read a 
        string this_wlg_a("");
        for(int jj = 0; jj < best_overlap_lg_a.size(); jj ++)
        {
            if(jj > 0)
            this_wlg_a += ",";
            this_wlg_a += best_overlap_lg_a[jj];
        }
        // linkage groups of signal read b
        string this_wlg_b("");
        for(int jj = 0; jj < best_overlap_lg_b.size(); jj ++)
        {
            if(jj > 0)
            this_wlg_b += ",";
            this_wlg_b += best_overlap_lg_b[jj];
        }         
        // find common lgs between two reads         
        vector<string> common_lg_togo;
        if(!get_target_lg_R12(best_overlap_lg_a,
                              best_overlap_type_a,
                              best_overlap_lg_b, 
                              best_overlap_type_b,
                              &common_lg_togo) )
        {
            return 1;
        }
        // no common linkage group found between the two reads; inter-lg read pairs - will not be used in scaffolding
        if(common_lg_togo.size() == 0)
        {
           if(verbose)        
            cout << "   check: read a would go to lg:" << this_wlg_a  
                 << " while read b would go to lg:"    << this_wlg_b << endl;
            if(verbose)                 
            cout << "          no common lg to go."    << endl;
            string this_lg_readfile_a = "mapped_inter_lg_illu_reads_" + read_a_R12 + ".fq.gz"; 
            assert(allofp.find(this_lg_readfile_a) != allofp.end() );           
            if(complex_readname)
            {            
                (*allofp[this_lg_readfile_a]) << "@"   
                                   << readname_a 
                                   << " "
                                   << read_a_reversed                                   
                                   << " + lg="
                                   << this_wlg_a
                                   << "/48"
                                   << " homLG="
                                   << "unchecked/12"
                                   << " rep.span="         
                                   << this_a_ctg 
                                   << ":"   
                                   << firstRefsMatch_a   
                                   << "-"          
                                   << spansecond_a
                                   << endl;   // read name  
            }else
            {
                (*allofp[this_lg_readfile_a]) << "@"   
                                   << readname_a 
                                   << " "
                                   << read_a_reversed                                   
                                   << endl;   // read name             
            }                      
            (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;   // read seq
            (*allofp[this_lg_readfile_a]) << "+" << endl;
            (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;                                            
            aligned_grp_a = this_lg_readfile_a; 
            //
            string this_lg_readfile_b = "mapped_inter_lg_illu_reads_" + read_b_R12 + ".fq.gz"; 
            assert(allofp.find(this_lg_readfile_b) != allofp.end() );           
            if(complex_readname)
            {             
                (*allofp[this_lg_readfile_b]) << "@"   
                                   << readname_b 
                                   << " "
                                   << read_b_reversed                                   
                                   << " + lg="
                                   << this_wlg_b                                   
                                   << "/48"
                                   << " homLG="
                                   << "unchecked/12"
                                   << " rep.span="         
                                   << this_b_ctg 
                                   << ":"   
                                   << firstRefsMatch_b   
                                   << "-"          
                                   << spansecond_b
                                   << endl;   // read name  
            }else
            {
                (*allofp[this_lg_readfile_b]) << "@"   
                                   << readname_b 
                                   << " "
                                   << read_b_reversed                                   
                                   << endl;   // read name                  
            }
            (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
            (*allofp[this_lg_readfile_b]) << "+" << endl;
            (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl;                                   
            aligned_grp_b = this_lg_readfile_b; 
            //
            groupedillureadnames[this_lg_readfile_a] += 1;
            groupedillureadnames[this_lg_readfile_b] += 1;                         
            //
            continue;                                   
        }
        // linkage groups of signal read a and b
        string this_wlg_ab("");
        for(int jj = 0; jj < common_lg_togo.size(); jj ++)
        {
            if(jj > 0)
            this_wlg_ab += ",";
            this_wlg_ab += common_lg_togo[jj];
        } 
        if(verbose)                
        cout << "   check: read a would go to lg:" << this_wlg_a  
             << " while read b would go to lg:"    << this_wlg_b  << endl;
        if(verbose)             
        cout << "          common lg to go "       << this_wlg_ab << endl;        
        double randra = rand()%1000000000*1.0 / 1000000000.0;
        if(common_lg_togo.size() == 1)
        {
            string tmp_goto_lg = common_lg_togo[0];
            assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
            string this_homLG = homLGs[tmp_goto_lg];
            /*
            cout << "   check: read "  
                 << readname_a    << " " << read_a_R12 << " would go to LG " 
                 << this_wlg_a    << ", and selected "
                 << tmp_goto_lg << "/48 in " 
                 << this_homLG  << "/12" 
                 << endl; 
            */           
            string this_lg_readfile_a = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_a_R12 +".fq.gz";
            string this_lg_readfile_b = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_b_R12 +".fq.gz";            
            if(tmp_goto_lg.compare("-1") == 0)
            {
                this_lg_readfile_a = "mapped_no_lg_illu_reads_" + read_a_R12 + ".fq.gz"; 
                this_lg_readfile_b = "mapped_no_lg_illu_reads_" + read_b_R12 + ".fq.gz";                 
                nonlgNo  += 2;           
                lennonlg += this_a_seq.size();            
                lennonlg += this_b_seq.size();                            
                // without lg
                // 20201024 will use too much mem for hundreds of millions of short reads!     
                /*           
                if(illureadnames_nolg.find(readname_a) == illureadnames_nolg.end())
                {
                    illureadnames_nolg.insert(std::pair<string,int>(readname_a, 1));
                } 
                */
             }
             // read a
            if(complex_readname)
            {                 
                 (*allofp[this_lg_readfile_a]) << "@"   
                                    << readname_a 
                                    << " "
                                    << read_a_reversed                                    
                                    << " + lg="
                                    << tmp_goto_lg << "/48"
                                    << " homLG="
                                    << this_homLG  << "/12"
                                    << " span="         
                                    << this_a_ctg 
                                    << ":"   
                                    << firstRefsMatch_a   
                                    << "-"          
                                    << spansecond_a
                                    << endl;   // read name  
             }else
             {
                 (*allofp[this_lg_readfile_a]) << "@"   
                                    << readname_a 
                                    << " "
                                    << read_a_reversed                                               
                                    << endl;   // read name               
             }
            (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;   // read seq
            (*allofp[this_lg_readfile_a]) << "+" << endl;
            (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;   
             // read b            
            if(complex_readname)
            {             
                 (*allofp[this_lg_readfile_b]) << "@"   
                                    << readname_b 
                                    << " "
                                    << read_b_reversed                                               
                                    << " + lg="
                                    << tmp_goto_lg << "/48"
                                    << " homLG="
                                    << this_homLG  << "/12"
                                    << " span="         
                                    << this_b_ctg 
                                    << ":"   
                                    << firstRefsMatch_b  
                                    << "-"          
                                    << spansecond_b
                                    << endl;   // read name  
             }else
             {
                 (*allofp[this_lg_readfile_b]) << "@"   
                                    << readname_b 
                                    << " "
                                    << read_b_reversed                                               
                                    << endl;   // read name               
             }
            (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
            (*allofp[this_lg_readfile_b]) << "+" << endl;
            (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl;                                   
            //
            /* 20201024 will use too much mem for hundreds of millions of short reads!            
            if(groupedillureadnames[this_lg_readfile].find(readname_a) == 
               groupedillureadnames[this_lg_readfile].end())
            {
                groupedillureadnames[this_lg_readfile].insert(std::pair<string, int>(readname_a, 1));
            }else
            {
                groupedillureadnames[this_lg_readfile][readname_a] += 1;
            } 
            */
            groupedillureadnames[this_lg_readfile_a] += 1;
            groupedillureadnames[this_lg_readfile_b] += 1;                 
            //
            aligned_grp_a = this_lg_readfile_a;
            aligned_grp_b = this_lg_readfile_b;
        }else
        if(common_lg_togo.size() == 2)
        {
            if(randra <= 0.5)
            {
                string tmp_goto_lg = common_lg_togo[0];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                /*
                cout << "   check: read "  
                     << readname_a    << " " << read_a_R12 << " would go to LG " 
                     << this_wlg_a    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in [0, 0.5]"
                     << endl;         
                */   
                string this_lg_readfile_a = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_a_R12 +".fq.gz";
                string this_lg_readfile_b = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_b_R12 +".fq.gz";            
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile_a = "mapped_no_lg_illu_reads_" + read_a_R12 + ".fq.gz"; 
                    this_lg_readfile_b = "mapped_no_lg_illu_reads_" + read_b_R12 + ".fq.gz";                 
                    nonlgNo  += 2;           
                    lennonlg += this_a_seq.size();            
                    lennonlg += this_b_seq.size();                            
                    // without lg
                    // 20201024 will use too much mem for hundreds of millions of short reads!     
                    /*           
                    if(illureadnames_nolg.find(readname_a) == illureadnames_nolg.end())
                    {
                        illureadnames_nolg.insert(std::pair<string,int>(readname_a, 1));
                    } 
                    */
                 }  
                 // read a
                if(complex_readname)
                {                    
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_a_ctg 
                                        << ":"   
                                        << firstRefsMatch_a   
                                        << "-"          
                                        << spansecond_a
                                        << " "
                                        << randra                                         
                                        << endl;   // read name 
                 }else
                 {
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << endl;   // read name                  
                 }                        
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_a]) << "+" << endl;
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;   
                 // read b   
                if(complex_readname)
                {                            
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_b_ctg 
                                        << ":"   
                                        << firstRefsMatch_b  
                                        << "-"          
                                        << spansecond_b
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << endl;   // read name                   
                 }
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_b]) << "+" << endl;
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl;    
                //
                /* 20201024 will use too much mem for hundreds of millions of short reads!                
                if(groupedillureadnames[this_lg_readfile].find(readname_a) == 
                   groupedillureadnames[this_lg_readfile].end())
                {
                    groupedillureadnames[this_lg_readfile].insert(std::pair<string, int>(readname_a, 1));
                }else
                {
                    groupedillureadnames[this_lg_readfile][readname_a] += 1;
                }
                */ 
                groupedillureadnames[this_lg_readfile_a] += 1;
                groupedillureadnames[this_lg_readfile_b] += 1;                     
                //
                aligned_grp_a = this_lg_readfile_a;                                                      
                aligned_grp_b = this_lg_readfile_b;                                                                      
            }else
            {
                string tmp_goto_lg = common_lg_togo[1];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                /*
                cout << "   check: read "  
                     << readname_a    << " " << read_a_R12 << " would go to LG " 
                     << this_wlg_a    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (0.5, 1]"                     
                     << endl;            
                */
                string this_lg_readfile_a = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_a_R12 +".fq.gz";
                string this_lg_readfile_b = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_b_R12 +".fq.gz";            
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile_a = "mapped_no_lg_illu_reads_" + read_a_R12 + ".fq.gz"; 
                    this_lg_readfile_b = "mapped_no_lg_illu_reads_" + read_b_R12 + ".fq.gz";                 
                    nonlgNo  += 2;           
                    lennonlg += this_a_seq.size();            
                    lennonlg += this_b_seq.size();                            
                    // without lg
                    // 20201024 will use too much mem for hundreds of millions of short reads!     
                    /*           
                    if(illureadnames_nolg.find(readname_a) == illureadnames_nolg.end())
                    {
                        illureadnames_nolg.insert(std::pair<string,int>(readname_a, 1));
                    } 
                    */
                 }  
                 // read a
                if(complex_readname)
                {                   
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_a_ctg 
                                        << ":"   
                                        << firstRefsMatch_a   
                                        << "-"          
                                        << spansecond_a
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << endl;   // read name                   
                 }                     
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_a]) << "+" << endl;
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;  
                 // read b    
                if(complex_readname)
                {                          
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_b_ctg 
                                        << ":"   
                                        << firstRefsMatch_b  
                                        << "-"          
                                        << spansecond_b
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << endl;   // read name                  
                 }
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_b]) << "+" << endl;
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl;  
                //
                /* 20201024 will use too much mem for hundreds of millions of short reads!                
                if(groupedillureadnames[this_lg_readfile].find(readname_a) == 
                   groupedillureadnames[this_lg_readfile].end())
                {
                    groupedillureadnames[this_lg_readfile].insert(std::pair<string, int>(readname_a, 1));
                }else
                {
                    groupedillureadnames[this_lg_readfile][readname_a] += 1;
                }
                */ 
                groupedillureadnames[this_lg_readfile_a] += 1;
                groupedillureadnames[this_lg_readfile_b] += 1;                                     
                //
                aligned_grp_a = this_lg_readfile_a;                                                      
                aligned_grp_b = this_lg_readfile_b;                                                   
            }                       
        }else
        if(common_lg_togo.size() == 3)
        {
            if(randra <= 1.0/3)
            {
                string tmp_goto_lg = common_lg_togo[0];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                /*
                cout << "   check: read "  
                     << readname_a    << " " << read_a_R12 << " would go to LG " 
                     << this_wlg_a    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in [0, 1/3]"                     
                     << endl;            
                */
                string this_lg_readfile_a = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_a_R12 +".fq.gz";
                string this_lg_readfile_b = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_b_R12 +".fq.gz";            
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile_a = "mapped_no_lg_illu_reads_" + read_a_R12 + ".fq.gz"; 
                    this_lg_readfile_b = "mapped_no_lg_illu_reads_" + read_b_R12 + ".fq.gz";                 
                    nonlgNo  += 2;           
                    lennonlg += this_a_seq.size();            
                    lennonlg += this_b_seq.size();                            
                    // without lg
                    // 20201024 will use too much mem for hundreds of millions of short reads!     
                    /*           
                    if(illureadnames_nolg.find(readname_a) == illureadnames_nolg.end())
                    {
                        illureadnames_nolg.insert(std::pair<string,int>(readname_a, 1));
                    } 
                    */
                 }  
                 // read a
                if(complex_readname)
                {                  
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_a_ctg 
                                        << ":"   
                                        << firstRefsMatch_a   
                                        << "-"          
                                        << spansecond_a
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << endl;   // read name                   
                 }                       
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_a]) << "+" << endl;
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;  
                 // read b     
                if(complex_readname)
                {                        
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_b_ctg 
                                        << ":"   
                                        << firstRefsMatch_b  
                                        << "-"          
                                        << spansecond_b
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << endl;   // read name                   
                 }                       
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_b]) << "+" << endl;
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl; 
                //
                /* 20201024 will use too much mem for hundreds of millions of short reads!                
                if(groupedillureadnames[this_lg_readfile].find(readname_a) == 
                   groupedillureadnames[this_lg_readfile].end())
                {
                    groupedillureadnames[this_lg_readfile].insert(std::pair<string, int>(readname_a, 1));
                }else
                {
                    groupedillureadnames[this_lg_readfile][readname_a] += 1;
                }
                */ 
                groupedillureadnames[this_lg_readfile_a] += 1;
                groupedillureadnames[this_lg_readfile_b] += 1;                                     
                //
                aligned_grp_a = this_lg_readfile_a;                                                      
                aligned_grp_b = this_lg_readfile_b;                                                   
            }else
            if(1.0/3<randra && randra<=2.0/3)            
            {
                string tmp_goto_lg = common_lg_togo[1];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                /*
                cout << "   check: read "  
                     << readname_a    << " " << read_a_R12 << " would go to LG " 
                     << this_wlg_a    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (1/3, 2/3]"                     
                     << endl;        
                */    
                string this_lg_readfile_a = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_a_R12 +".fq.gz";
                string this_lg_readfile_b = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_b_R12 +".fq.gz";            
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile_a = "mapped_no_lg_illu_reads_" + read_a_R12 + ".fq.gz"; 
                    this_lg_readfile_b = "mapped_no_lg_illu_reads_" + read_b_R12 + ".fq.gz";                 
                    nonlgNo  += 2;           
                    lennonlg += this_a_seq.size();            
                    lennonlg += this_b_seq.size();                            
                    // without lg
                    // 20201024 will use too much mem for hundreds of millions of short reads!     
                    /*           
                    if(illureadnames_nolg.find(readname_a) == illureadnames_nolg.end())
                    {
                        illureadnames_nolg.insert(std::pair<string,int>(readname_a, 1));
                    } 
                    */
                 }  
                 // read a
                if(complex_readname)
                {                   
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_a_ctg 
                                        << ":"   
                                        << firstRefsMatch_a   
                                        << "-"          
                                        << spansecond_a
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << endl;   // read name                  
                 }
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_a]) << "+" << endl;
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;  
                 // read b   
                if(complex_readname)
                {                             
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_b_ctg 
                                        << ":"   
                                        << firstRefsMatch_b  
                                        << "-"          
                                        << spansecond_b
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << endl;   // read name                  
                 }
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_b]) << "+" << endl;
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl;  
                //
                /* 20201024 will use too much mem for hundreds of millions of short reads!                
                if(groupedillureadnames[this_lg_readfile].find(readname_a) == 
                   groupedillureadnames[this_lg_readfile].end())
                {
                    groupedillureadnames[this_lg_readfile].insert(std::pair<string, int>(readname_a, 1));
                }else
                {
                    groupedillureadnames[this_lg_readfile][readname_a] += 1;
                }
                */
                groupedillureadnames[this_lg_readfile_a] += 1;
                groupedillureadnames[this_lg_readfile_b] += 1;                                      
                //
                aligned_grp_a = this_lg_readfile_a;                                                      
                aligned_grp_b = this_lg_readfile_b;                                                      
            }else
            {
                string tmp_goto_lg = common_lg_togo[2];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                /*
                cout << "   check: read "  
                     << readname_a    << " " << read_a_R12 << " would go to LG " 
                     << this_wlg_a    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (2/3, 1.0]"                     
                     << endl;
                */
                string this_lg_readfile_a = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_a_R12 +".fq.gz";
                string this_lg_readfile_b = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_b_R12 +".fq.gz";            
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile_a = "mapped_no_lg_illu_reads_" + read_a_R12 + ".fq.gz"; 
                    this_lg_readfile_b = "mapped_no_lg_illu_reads_" + read_b_R12 + ".fq.gz";                 
                    nonlgNo  += 2;           
                    lennonlg += this_a_seq.size();            
                    lennonlg += this_b_seq.size();                            
                    // without lg
                    // 20201024 will use too much mem for hundreds of millions of short reads!     
                    /*           
                    if(illureadnames_nolg.find(readname_a) == illureadnames_nolg.end())
                    {
                        illureadnames_nolg.insert(std::pair<string,int>(readname_a, 1));
                    } 
                    */
                 }  
                 // read a
                if(complex_readname)
                {                   
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a
                                        << " "
                                        << read_a_reversed                                                    
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_a_ctg 
                                        << ":"   
                                        << firstRefsMatch_a   
                                        << "-"          
                                        << spansecond_a
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << endl;   // read name                  
                 }
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_a]) << "+" << endl;
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;  
                 // read b  
                if(complex_readname)
                {                            
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_b_ctg 
                                        << ":"   
                                        << firstRefsMatch_b  
                                        << "-"          
                                        << spansecond_b
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << endl;   // read name                   
                 }
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_b]) << "+" << endl;
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl; 
                //
                /* 20201024 will use too much mem for hundreds of millions of short reads!                
                if(groupedillureadnames[this_lg_readfile].find(readname_a) == 
                   groupedillureadnames[this_lg_readfile].end())
                {
                    groupedillureadnames[this_lg_readfile].insert(std::pair<string, int>(readname_a, 1));
                }else
                {
                    groupedillureadnames[this_lg_readfile][readname_a] += 1;
                }
                */
                groupedillureadnames[this_lg_readfile_a] += 1;
                groupedillureadnames[this_lg_readfile_b] += 1;                                      
                //
                aligned_grp_a = this_lg_readfile_a;                                                      
                aligned_grp_b = this_lg_readfile_b;                                                       
            }                       
        }else
        if(common_lg_togo.size() == 4)
        {
            if(randra <= 1.0/4)
            {
                string tmp_goto_lg = common_lg_togo[0];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                /*
                cout << "   check: read "  
                     << readname_a    << " " << read_a_R12 << " would go to LG " 
                     << this_wlg_a    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in [0, 1/4]"                     
                     << endl;         
                */   
                string this_lg_readfile_a = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_a_R12 +".fq.gz";
                string this_lg_readfile_b = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_b_R12 +".fq.gz";            
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile_a = "mapped_no_lg_illu_reads_" + read_a_R12 + ".fq.gz"; 
                    this_lg_readfile_b = "mapped_no_lg_illu_reads_" + read_b_R12 + ".fq.gz";                 
                    nonlgNo  += 2;           
                    lennonlg += this_a_seq.size();            
                    lennonlg += this_b_seq.size();                            
                    // without lg
                    // 20201024 will use too much mem for hundreds of millions of short reads!     
                    /*           
                    if(illureadnames_nolg.find(readname_a) == illureadnames_nolg.end())
                    {
                        illureadnames_nolg.insert(std::pair<string,int>(readname_a, 1));
                    } 
                    */
                 }  
                 // read a
                if(complex_readname)
                {                  
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_a_ctg 
                                        << ":"   
                                        << firstRefsMatch_a   
                                        << "-"          
                                        << spansecond_a
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << endl;   // read name                   
                 }                       
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_a]) << "+" << endl;
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;  
                 // read b  
                 if(complex_readname)
                 {                             
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_b_ctg 
                                        << ":"   
                                        << firstRefsMatch_b  
                                        << "-"          
                                        << spansecond_b
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                  }else
                  {
                      (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                   
                                        << endl;   // read name                   
                  }                       
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_b]) << "+" << endl;
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl;  
                //
                /* 20201024 will use too much mem for hundreds of millions of short reads!                
                if(groupedillureadnames[this_lg_readfile].find(readname_a) == 
                   groupedillureadnames[this_lg_readfile].end())
                {
                    groupedillureadnames[this_lg_readfile].insert(std::pair<string, int>(readname_a, 1));
                }else
                {
                    groupedillureadnames[this_lg_readfile][readname_a] += 1;
                }
                */
                groupedillureadnames[this_lg_readfile_a] += 1;
                groupedillureadnames[this_lg_readfile_b] += 1;                                      
                //
                aligned_grp_a = this_lg_readfile_a;                                                      
                aligned_grp_b = this_lg_readfile_b;                                                   
            }else
            if(1.0/4<randra && randra<=2.0/4)            
            {
                string tmp_goto_lg = common_lg_togo[1];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                /*
                cout << "   check: read "  
                     << readname_a    << " " << read_a_R12 << " would go to LG " 
                     << this_wlg_a    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (1/4, 2/4]"                     
                     << endl;            
                */
                string this_lg_readfile_a = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_a_R12 +".fq.gz";
                string this_lg_readfile_b = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_b_R12 +".fq.gz";            
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile_a = "mapped_no_lg_illu_reads_" + read_a_R12 + ".fq.gz"; 
                    this_lg_readfile_b = "mapped_no_lg_illu_reads_" + read_b_R12 + ".fq.gz";                 
                    nonlgNo  += 2;           
                    lennonlg += this_a_seq.size();            
                    lennonlg += this_b_seq.size();                            
                    // without lg
                    // 20201024 will use too much mem for hundreds of millions of short reads!     
                    /*           
                    if(illureadnames_nolg.find(readname_a) == illureadnames_nolg.end())
                    {
                        illureadnames_nolg.insert(std::pair<string,int>(readname_a, 1));
                    } 
                    */
                 }  
                 // read a
                 if(complex_readname)
                 {                      
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                   
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_a_ctg 
                                        << ":"   
                                        << firstRefsMatch_a   
                                        << "-"          
                                        << spansecond_a
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a     
                                        << " "
                                        << read_a_reversed                                                                                                        
                                        << endl;   // read name                   
                 }
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_a]) << "+" << endl;
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;  
                 // read b    
                 if(complex_readname)
                 {                         
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                                      
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_b_ctg 
                                        << ":"   
                                        << firstRefsMatch_b  
                                        << "-"          
                                        << spansecond_b
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                                      
                                        << endl;   // read name                   
                 }                       
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_b]) << "+" << endl;
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl; 
                //
                /* 20201024 will use too much mem for hundreds of millions of short reads!                
                if(groupedillureadnames[this_lg_readfile].find(readname_a) == 
                   groupedillureadnames[this_lg_readfile].end())
                {
                    groupedillureadnames[this_lg_readfile].insert(std::pair<string, int>(readname_a, 1));
                }else
                {
                    groupedillureadnames[this_lg_readfile][readname_a] += 1;
                }
                */ 
                groupedillureadnames[this_lg_readfile_a] += 1;
                groupedillureadnames[this_lg_readfile_b] += 1;                                     
                //
                aligned_grp_a = this_lg_readfile_a;                                                      
                aligned_grp_b = this_lg_readfile_b;                                                      
            }else
            if(2.0/4<randra && randra<=3.0/4)            
            {
                string tmp_goto_lg = common_lg_togo[2];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                /*
                cout << "   check: read "  
                     << readname_a    << " " << read_a_R12 << " would go to LG " 
                     << this_wlg_a    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (2/4, 3/4]"                     
                     << endl;     
                */       
                string this_lg_readfile_a = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_a_R12 +".fq.gz";
                string this_lg_readfile_b = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_b_R12 +".fq.gz";            
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile_a = "mapped_no_lg_illu_reads_" + read_a_R12 + ".fq.gz"; 
                    this_lg_readfile_b = "mapped_no_lg_illu_reads_" + read_b_R12 + ".fq.gz";                 
                    nonlgNo  += 2;           
                    lennonlg += this_a_seq.size();            
                    lennonlg += this_b_seq.size();                            
                    // without lg
                    // 20201024 will use too much mem for hundreds of millions of short reads!     
                    /*           
                    if(illureadnames_nolg.find(readname_a) == illureadnames_nolg.end())
                    {
                        illureadnames_nolg.insert(std::pair<string,int>(readname_a, 1));
                    } 
                    */
                 }  
                 // read a
                 if(complex_readname)
                 {                   
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                                      
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_a_ctg 
                                        << ":"   
                                        << firstRefsMatch_a   
                                        << "-"          
                                        << spansecond_a
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                                      
                                        << endl;   // read name                   
                 }
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_a]) << "+" << endl;
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;  
                 // read b     
                 if(complex_readname)
                 {                        
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                                      
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_b_ctg 
                                        << ":"   
                                        << firstRefsMatch_b  
                                        << "-"          
                                        << spansecond_b
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                                      
                                        << endl;   // read name                   
                 }
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_b]) << "+" << endl;
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl; 
                //
                /* 20201024 will use too much mem for hundreds of millions of short reads!                
                if(groupedillureadnames[this_lg_readfile].find(readname_a) == 
                   groupedillureadnames[this_lg_readfile].end())
                {
                    groupedillureadnames[this_lg_readfile].insert(std::pair<string, int>(readname_a, 1));
                }else
                {
                    groupedillureadnames[this_lg_readfile][readname_a] += 1;
                }
                */ 
                groupedillureadnames[this_lg_readfile_a] += 1;
                groupedillureadnames[this_lg_readfile_b] += 1;                                     
                //
                aligned_grp_a = this_lg_readfile_a;                                                      
                aligned_grp_b = this_lg_readfile_b;                                                      
            }else
            {
                string tmp_goto_lg = common_lg_togo[3];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                /*
                cout << "   check: read "  
                     << readname_a    << " " << read_a_R12 << " would go to LG " 
                     << this_wlg_a    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (3/4, 1.0]"                     
                     << endl;       
                */     
                string this_lg_readfile_a = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_a_R12 +".fq.gz";
                string this_lg_readfile_b = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads_" + read_b_R12 +".fq.gz";            
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile_a = "mapped_no_lg_illu_reads_" + read_a_R12 + ".fq.gz"; 
                    this_lg_readfile_b = "mapped_no_lg_illu_reads_" + read_b_R12 + ".fq.gz";                 
                    nonlgNo  += 2;           
                    lennonlg += this_a_seq.size();            
                    lennonlg += this_b_seq.size();                            
                    // without lg
                    // 20201024 will use too much mem for hundreds of millions of short reads!     
                    /*           
                    if(illureadnames_nolg.find(readname_a) == illureadnames_nolg.end())
                    {
                        illureadnames_nolg.insert(std::pair<string,int>(readname_a, 1));
                    } 
                    */
                 }  
                 // read a
                 if(complex_readname)
                 {                   
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                                      
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_a_ctg 
                                        << ":"   
                                        << firstRefsMatch_a   
                                        << "-"          
                                        << spansecond_a
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_a]) << "@"   
                                        << readname_a 
                                        << " "
                                        << read_a_reversed                                                                      
                                        << endl;   // read name                  
                 }                       
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_seq_rc:this_a_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_a]) << "+" << endl;
                (*allofp[this_lg_readfile_a]) << (this_a_seq_rc.size()>0?this_a_qual_rc:this_a_qual) << endl;  
                 // read b 
                 if(complex_readname)
                 {                              
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                                      
                                        << " + lg="
                                        << tmp_goto_lg << "/48"
                                        << " homLG="
                                        << this_homLG  << "/12"
                                        << " span="         
                                        << this_b_ctg 
                                        << ":"   
                                        << firstRefsMatch_b  
                                        << "-"          
                                        << spansecond_b
                                        << " "
                                        << randra                                         
                                        << endl;   // read name  
                 }else
                 {
                     (*allofp[this_lg_readfile_b]) << "@"   
                                        << readname_b 
                                        << " "
                                        << read_b_reversed                                                                      
                                        << endl;   // read name                   
                 }
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_seq_rc:this_b_seq)   << endl;// read seq
                (*allofp[this_lg_readfile_b]) << "+" << endl;
                (*allofp[this_lg_readfile_b]) << (this_b_seq_rc.size()>0?this_b_qual_rc:this_b_qual) << endl; 
                //
                /* 20201024 will use too much mem for hundreds of millions of short reads!                
                if(groupedillureadnames[this_lg_readfile].find(readname_a) == 
                   groupedillureadnames[this_lg_readfile].end())
                {
                    groupedillureadnames[this_lg_readfile].insert(std::pair<string, int>(readname_a, 1));
                }else
                {
                    groupedillureadnames[this_lg_readfile][readname_a] += 1;
                }
                */ 
                groupedillureadnames[this_lg_readfile_a] += 1;
                groupedillureadnames[this_lg_readfile_b] += 1;                                     
                //
                aligned_grp_a = this_lg_readfile_a;                                                      
                aligned_grp_b = this_lg_readfile_b;                                                     
            } 
        } 
        //
        /* 20201024 will use too much mem for hundreds of millions of short reads!        
        // collect informative reads with details
        std::stringstream readss;
        readss.str("");
        //if(aligned_grp_a.size() == 0) cout << "   Warning: this grp name is null: check " << readname_a << endl;
        readss << aligned_grp_a << ":" << spansecond_a-firstRefsMatch_a+1;      
        if(illureadnames_detail.find(readname_a) == illureadnames_detail.end())
        {
            vector<string> alignvec;
            alignvec.push_back(readss.str());              
            illureadnames_detail.insert(std::pair<string, vector<string> >(readname_a, alignvec));
        }else
        {
            illureadnames_detail[readname_a].push_back(readss.str());
        }                                              
        */
        // end of parsing bam/sam
    } 
    //
    cout << endl;
    if(numNoCG > 0)
    {
        cout << "   Warning: there are a1="
             << numNoCG 
             << " alignments, totaling v1=" 
             << lenNoCG/1000000000.0 
             << " Gb "
             << " without explicit CIAGR info -- collected in unmapped file."   
             << endl;     
    }
    if(numNoMA > 0)
    {
        cout << "   Warning: there are a2="
             << numNoMA 
             << " alignments being secondary/supplementary alignment, skipped. "
             << endl;
    }
    if(numNoSEQ > 0)
    {
        cout << "   Warning: there are a3="
             << numNoSEQ 
             << " alignments without seq field - secondary/supplementary alignment"
             << " or not passing filters, skipped." 
             << endl;
    }
    cout << "   Info: in total all=" 
         << numraw                  
         << " aligment lines collected (<-header line not counted; raw alignment including secondary etc - skipped in read sep). " 
         << endl;
    // no lg
    cout << "   Info: number of illu alignment seqs WITHOUT linkage Info: "   
         << nonlgNo       
         << ", totaling v2=" 
         << lennonlg/1000000000.0      
         << " Gb" 
         << endl;
    // read statistics 
    // map<string, map<string, int> > groupedillureadnames;    
    cout << endl 
         << "   Info: distribution of short reads (alignment) in linkage groups: "
         << endl << endl;
    unsigned long total_reads = 0;
    map<string, unsigned long>::iterator gpbitr;
    map<string, unsigned long>::iterator gpbitr_end;
    gpbitr     = groupedillureadnames.begin();
    gpbitr_end = groupedillureadnames.end();
    while(gpbitr != gpbitr_end)
    {
        cout << "\t" << (*gpbitr).first 
             << "\t" << (*gpbitr).second 
             << endl;
        total_reads += (*gpbitr).second;
        gpbitr ++;
    }       
    cout << "   Info: " << total_reads << " R1+R2 reads refer to unique read counts if samtools view -F 3840 is applied. " << endl;                   
    // close output files
    map<string, ogzstream*>::iterator ofileitr;
    map<string, ogzstream*>::iterator ofileitr_end;
    ofileitr     = allofp.begin();
    ofileitr_end = allofp.end();
    while(ofileitr != ofileitr_end)
    {
        ogzstream *ofp = (*ofileitr).second;
        (*ofp).close();
        //
        ofileitr ++;
    }    
    cout << "   Info: extract reads from bam/sam into linkage done. "     << endl;
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;    
}
//
bool get_target_lg_R12(vector<string>  best_overlap_lg_a, 
                       string          best_overlap_type_a,
                       vector<string>  best_overlap_lg_b, 
                       string          best_overlap_type_b,
                       vector<string>* common_lg_togo)
{
    /* function: find common linkage group x/48 the paired reads should go to */
    (*common_lg_togo).clear();
    for(int i = 0; i < best_overlap_lg_a.size(); i ++)
    {
        vector<string>::iterator tmp_itr = std::find(best_overlap_lg_b.begin(), 
                                                     best_overlap_lg_b.end(),
                                                     best_overlap_lg_a[i]);
        if(tmp_itr != best_overlap_lg_b.end() )
        {
            (*common_lg_togo).push_back( best_overlap_lg_a[i] );
        }
    }
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
unsigned long overlap_read_and_win_marker_v2(unsigned long    this_win_sta, 
                           unsigned long                      this_win_end, 
                           unsigned long                      read_sta, 
                           unsigned long                      read_end,
                           string                             this_ctg)
{
    /*
        Input: 
           this_win_marker sta and end
           read_sta        current long read alignment starting position along the reference 
           read_end        current long read alignment ending   position along the reference            
        Output:            length of current window marker overlapping with the current short read.    
    */
    // score is the overlapping length between window marker and read
    unsigned long win_sta = this_win_sta;
    unsigned long win_end = this_win_end;        
    unsigned long olspan = 0; // overlapping span
    if(win_sta<=read_sta && read_sta<=win_end)
    {
        /*
                  sta---------win-makrer--------end         window marker span
                        sta-----read-----end                short read case 1: right_min = read_end
                        sta-----read------------------end   short read case 2: right_min = win_end
        */
        unsigned long right_min = read_end<win_end?read_end:win_end;
        olspan = right_min - read_sta + 1;
    }else
    if(win_sta<=read_end && read_end<=win_end)
    {
        /*
                  sta---------win-marker--------end         window marker span
            sta--------------read-----end                   short read case 2: right_min = win_end
        */              
        unsigned long left_max = read_sta>win_sta?read_sta:win_sta;
        olspan = read_end - left_max + 1;
    }else
    if(win_sta>=read_sta && read_end>=win_end)
    {
        /*
                  sta---------win-marker--------end         window marker span
            sta--------------read---------------------end   short read case 2: right_min = win_end
        */    
        olspan = win_end - win_sta + 1;         
    }
    // do we need a cutoff for considering overlapping as sufficient? current as least 10-20 bp.
    if(olspan >= 20 || (olspan>=10 && read_end - read_sta + 1 <= 80))
    {
        if(verbose)    
        cout << "          : window-maker " 
             << this_ctg                << ":" 
             << win_sta                 << "-"
             << win_end                 << " overlap short read spanning "
             << read_sta                << "-" 
             << read_end                << " by " 
             << olspan                  << " bp. "
             << endl;
    }else
    {
        if(0 && verbose)
        cout << "          : window-marker " 
             << this_ctg                << ":" 
             << win_sta                 << "-" 
             << win_end                 << " does not overlap short read spaning "                               
             << read_sta                << "-" 
             << read_end                << "." 
             << endl; 
    }
    //
    return olspan;
}
//
bool overlap_read_and_win_marker(map<unsigned long, unsigned long> this_win_marker, 
                           unsigned long                           read_sta, 
                           unsigned long                           read_end,
                           string                                  this_ctg,
                           map<unsigned long, unsigned long>*      this_win_marker_overlapped)
{
    /*
        Input: 
           this_win_marker lists of window markers along a contig
           read_sta        current long read alignment starting position along the reference 
           read_end        current long read alignment ending   position along the reference            
        Output:            window markers overlapping with the current short read.    
    */
    map<unsigned long, unsigned long>::iterator witr;
    map<unsigned long, unsigned long>::iterator witr_end;
    witr     = this_win_marker.begin();
    witr_end = this_win_marker.end();
    // score is the overlapping length between window marker and read
    // if both genotypes are collected, we choose the one with larger overlapping region.
    map<string, unsigned long> genotype_max_overlap; 
    while(witr != witr_end)
    {
        unsigned long win_sta = (*witr).first;
        unsigned long win_end = (*witr).second;        
        unsigned long olspan = 0; // overlapping span
        if(win_sta<=read_sta && read_sta<=win_end)
        {
            /*
                      sta---------win-makrer--------end         window marker span
                            sta-----read-----end                short read case 1: right_min = read_end
                            sta-----read------------------end   short read case 2: right_min = win_end
            */
            unsigned long right_min = read_end<win_end?read_end:win_end;
            olspan = right_min - read_sta + 1;
        }else
        if(win_sta<=read_end && read_end<=win_end)
        {
            /*
                      sta---------win-marker--------end         window marker span
                sta--------------read-----end                   short read case 2: right_min = win_end
            */              
            unsigned long left_max = read_sta>win_sta?read_sta:win_sta;
            olspan = read_end - left_max + 1;
        }else
        if(win_sta>=read_sta && read_end>=win_end)
        {
            /*
                      sta---------win-marker--------end         window marker span
                sta--------------read---------------------end   short read case 2: right_min = win_end
            */    
            olspan = win_end - win_sta + 1;         
        }
        // do we need a cutoff for considering overlapping as sufficient? current as least 500 bp.
        if(olspan >= 20 || (olspan>=10 && read_end - read_sta + 1 <= 80))
        {
            (*this_win_marker_overlapped).insert(std::pair<unsigned long, unsigned long>((*witr).first, (*witr).second));
            if(verbose)            
            cout << "          : window-maker " 
                 << this_ctg                << ":" 
                 << win_sta                 << "-"
                 << win_end                 << " overlap illu read spanning "
                 << read_sta                << "-" 
                 << read_end                << " by " 
                 << olspan                  << " bp. "
                 << endl;
        }else
        {
            if(0 && verbose)
            cout << "          : window-marker " 
                 << this_ctg                << ":" 
                 << win_sta                 << "-" 
                 << win_end                 << " does not overlap illu read spaning "                               
                 << read_sta                << "-" 
                 << read_end                << "." 
                 << endl; 
        }
        witr ++;
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
        if(operation[oi]=='D') // TODO N can be added here! 20210822
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
//
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
//
bool merge_gb_window_marker_grouping(
                 map<string, map<unsigned long, winMARKER > >  gbmarker,
                 map<string, map<unsigned long, winMARKER > >* gbmarker_updated)
{
    /*
       function: merge smallers window markers into large ones, if they are contigous; same type and same <LG-list>
       
	utg000005l_pilon	2260001	2310000	dip	34	11
	utg000005l_pilon	2310001	2360000	dip	1	11
	utg000005l_pilon	2310001	2360000	dip	34	11
	utg000005l_pilon	2360001	2400000	dip	1	11
	utg000005l_pilon	2360001	2400000	dip	34	11
	utg000005l_pilon	2400001	2410000	trip	1	11
	utg000005l_pilon	2400001	2410000	trip	33	11
	utg000005l_pilon	2400001	2410000	trip	44	11
	utg000005l_pilon	2410001	2460000	dip	1	11
	utg000005l_pilon	2410001	2460000	dip	34	11
	utg000005l_pilon	2460001	2510000	dip	1	11
	utg000005l_pilon	2460001	2510000	dip	34	11
	
	=> 3 larger markers
	
	utg000005l_pilon	2260001	2400000	dip	1,34	11
	utg000005l_pilon	2400001	2410000	trip	1,33,34	11
	utg000005l_pilon	2410001	2510000	dip	1,34
       
       return..: gbmarker_updated: <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<1:4>, 1:12 } > >	
    */
    (*gbmarker_updated).clear();
    //
    map<string, map<unsigned long, winMARKER > >::iterator mcitr;
    map<string, map<unsigned long, winMARKER > >::iterator mcitr_end;
    mcitr     = gbmarker.begin();
    mcitr_end = gbmarker.end();   
    while(mcitr != mcitr_end)
    {
        string            last_ctg  = "";
        unsigned long last_win_sta  = 0;
        unsigned long last_win_end  = 0;
        string            last_type = "";
        vector<string>    last_lgs;
        last_lgs.clear();
        string          last_homLG  = "";     
        //
        string this_ctg = (*mcitr).first;
        map<unsigned long, winMARKER> tmp_win_list = (*mcitr).second;
        map<unsigned long, winMARKER>::iterator mpitr;
        map<unsigned long, winMARKER>::iterator mpitr_end;
        mpitr     = tmp_win_list.begin();
        mpitr_end = tmp_win_list.end();
        while(mpitr != mpitr_end)
        {
            unsigned long this_win_sta = (*mpitr).first;
            winMARKER      tmp_marker  = (*mpitr).second;
            unsigned long this_win_end = tmp_marker.end;
            string           this_type = tmp_marker.type;
            vector<string>    this_lgs = tmp_marker.wlg;
            string          this_homLG = tmp_marker.homlg;
            //
            if(last_ctg.size() == 0)
            {
                // initialize
                last_ctg      = this_ctg;
                last_win_sta  = this_win_sta;
                last_win_end  = this_win_end;
                last_type     = this_type;
                last_lgs      = this_lgs;
                last_homLG    = this_homLG;             
            }else
            {
                // continous coordinate at same ctg; same type; same lgs to go
                if(this_win_sta == last_win_end + 1 && 
                   this_type.compare(last_type)==0  && 
                   compare_lgs(this_lgs, last_lgs) )
                {
                    // merge with last window
                    last_win_end = this_win_end;
                }else
                {
                    // collect current merged windows
                    if( (*gbmarker_updated).find(last_ctg) == (*gbmarker_updated).end() )
                    {
                        winMARKER tmp_pos_marker;
                        tmp_pos_marker.end   = last_win_end;
                        tmp_pos_marker.type  = last_type;
                        tmp_pos_marker.wlg   = last_lgs;
                        tmp_pos_marker.homlg = last_homLG;
                        map<unsigned long, winMARKER> merged_win_list;
                        merged_win_list.insert(std::pair<unsigned long, winMARKER>(last_win_sta, tmp_pos_marker));
                        (*gbmarker_updated).insert(std::pair<string, map<unsigned long, winMARKER> >
                                                   (last_ctg, merged_win_list) );                   
                    }else
                    {
                        winMARKER tmp_pos_marker;
                        tmp_pos_marker.end   = last_win_end;
                        tmp_pos_marker.type  = last_type;
                        tmp_pos_marker.wlg   = last_lgs;
                        tmp_pos_marker.homlg = last_homLG;                  
                        assert( (*gbmarker_updated)[last_ctg].find(last_win_sta) == (*gbmarker_updated)[last_ctg].end() );
                        (*gbmarker_updated)[last_ctg].insert(std::pair<unsigned long, winMARKER >
                                                             (last_win_sta, tmp_pos_marker) );
                    }
                    // update next
                    last_ctg      = this_ctg;
                    last_win_sta  = this_win_sta;
                    last_win_end  = this_win_end;
                    last_type     = this_type;
                    last_lgs      = this_lgs;
                    last_homLG    = this_homLG;                          
                }                
            }                            
            // next window marker
            mpitr ++;
        }
        // collect last (merged) window(s) of current ctg
        if((*gbmarker_updated).find(last_ctg) == (*gbmarker_updated).end())
        {
            winMARKER tmp_pos_marker;
            tmp_pos_marker.end   = last_win_end;
            tmp_pos_marker.type  = last_type;
            tmp_pos_marker.wlg   = last_lgs;
            tmp_pos_marker.homlg = last_homLG;
            map<unsigned long, winMARKER> merged_win_list;
            merged_win_list.insert(std::pair<unsigned long, winMARKER>(last_win_sta, tmp_pos_marker));
            (*gbmarker_updated).insert(std::pair<string, map<unsigned long, winMARKER> >
                                       (last_ctg, merged_win_list) );
        }else
        {
            winMARKER tmp_pos_marker;
            tmp_pos_marker.end   = last_win_end;
            tmp_pos_marker.type  = last_type;
            tmp_pos_marker.wlg   = last_lgs;
            tmp_pos_marker.homlg = last_homLG;
            assert( (*gbmarker_updated)[last_ctg].find(last_win_sta) == (*gbmarker_updated)[last_ctg].end() );
            (*gbmarker_updated)[last_ctg].insert(std::pair<unsigned long, winMARKER >
                                                 (last_win_sta, tmp_pos_marker) );
        }        
        // next ctg
        mcitr ++;
    }
    //
    bool check_this = true;
    if(check_this && verbose)
    {
        cout << "   check: merged ctg win marker distribution in LGs: " << endl;
        map<string, map<unsigned long, winMARKER > >::iterator mcitr;
        map<string, map<unsigned long, winMARKER > >::iterator mcitr_end;
        mcitr     = (*gbmarker_updated).begin();
        mcitr_end = (*gbmarker_updated).end();
        while(mcitr != mcitr_end)
        {
            string this_ctg = (*mcitr).first;
            map<unsigned long, winMARKER > tmp_win_list = (*mcitr).second;
            map<unsigned long, winMARKER >::iterator mpitr;
            map<unsigned long, winMARKER >::iterator mpitr_end;
            mpitr     = tmp_win_list.begin();
            mpitr_end = tmp_win_list.end();
            while(mpitr != mpitr_end)
            {
                unsigned long win_sta = (*mpitr).first;
                winMARKER tmp_marker  = (*mpitr).second;
                cout << "        : " << this_ctg 
                     << "\t"         <<  win_sta 
                     << "\t"         << tmp_marker.end 
                     << "\t"         << tmp_marker.type 
                     << "\t"         << tmp_marker.homlg
                     << "\t";
                vector<string> this_LGs = tmp_marker.wlg;
                vector<string>::iterator lgitr;
                vector<string>::iterator lgitr_end;
                lgitr     = this_LGs.begin();
                lgitr_end = this_LGs.end();
                while(lgitr != lgitr_end)
                {
                    if(lgitr != this_LGs.begin())
                    {
                        cout << "," << *lgitr;
                    }else
                    {
                        cout <<        *lgitr;
                    }
                    lgitr ++;
                }
                cout << endl;                
                //
                mpitr ++;
            }
            //
            mcitr ++;
        }
    }    
    //
    return true;
}   
//
bool compare_lgs(vector<string> this_lgs, vector<string> last_lgs)
{
    if(this_lgs.size() != last_lgs.size())
    {
        return false;
    }
    for(int i=0; i<this_lgs.size(); i++)
    {
        vector<string>::iterator lgitr = std::find(last_lgs.begin(), last_lgs.end(), this_lgs[i]);
        if(lgitr == last_lgs.end() )
        {
            return false;
        }
    }
    return true;
}              
//
bool read_gb_window_marker_grouping(string gbmarker_file, 
         map<string, map<unsigned long, winMARKER > >* gbmarker,
         map<string, map<string, string> >*            ctg2lg,
         map<string, string>*                          homLGs)
{
    /*
       function: read gamete_binning grouped window markers       
       return..: gbmarker: <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
                 ctg2lg..: <ctg_id, <lg_id=1:4, 1:12> >
                 homLGs..: <lg_id=1:4,  homLG_id=1:12>
       //          
	struct winMARKER
	{
	    unsigned long  end; // window marker end
	    string        type; // hap/dip/trip/tetrap/rep
	    vector<string> wlg; // linkage group id 
	};
    */
    ifstream ifp;
    ifp.open(gbmarker_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << gbmarker_file << endl;
        return false;
    }
    while(ifp.good() )
    {
        string line("");
        getline(ifp, line);
        if(line.size() == 0 || line[0]=='#') continue;  
        // utg000380l_pilon    510001  560000  dip 1   11  utg001654l_pilon    0.958721    1   0   linked  case_2.2      
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 6)
        {
            cout << "   Warning: skipped insufficient info line of " << line << endl;
            continue;
        }  
        string this_ctg       = lineinfo[0];
        unsigned long win_sta = strtoul(lineinfo[1].c_str(), NULL, 0);
        unsigned long win_end = strtoul(lineinfo[2].c_str(), NULL, 0);
        string this_type      = lineinfo[3];
        string this_lg        = lineinfo[4]; // 1:4
        string this_homLG     = lineinfo[5]; // 1:12
        //
        if( (*gbmarker).find(this_ctg) == (*gbmarker).end() )
        {
            winMARKER tmp_marker;
            tmp_marker.end   = win_end;
            tmp_marker.type  = this_type;
            tmp_marker.wlg.push_back(this_lg);
            tmp_marker.homlg = this_homLG;
            map<unsigned long, winMARKER > tmp_win_list;
            tmp_win_list.insert(std::pair<unsigned long, winMARKER>(win_sta, tmp_marker));
            (*gbmarker).insert(std::pair<string, map<unsigned long, winMARKER > >(this_ctg, tmp_win_list));
        }else
        {
            if( (*gbmarker)[this_ctg].find( win_sta ) == (*gbmarker)[this_ctg].end() )
            {
                winMARKER tmp_marker;
                tmp_marker.end   = win_end;
                tmp_marker.type  = this_type;
                tmp_marker.wlg.push_back(this_lg);
                tmp_marker.homlg = this_homLG;                
                (*gbmarker)[this_ctg].insert(std::pair<unsigned long, winMARKER>(win_sta, tmp_marker));
            }else
            {
                assert(this_type.compare( (*gbmarker)[this_ctg][win_sta].type ) == 0 );
                (*gbmarker)[this_ctg][win_sta].wlg.push_back(this_lg);
            }
        }
        // update ctg2lg
        if( (*ctg2lg).find(this_ctg) == (*ctg2lg).end() )
        {
            // lg list 
            map<string, string> lgtmp;
            lgtmp.insert(std::pair<string, string>(this_lg, this_homLG));
            // ctg:: lg list
            (*ctg2lg).insert(std::pair<string, map<string, string> >(this_ctg, lgtmp));
        }else
        {
            if((*ctg2lg)[this_ctg].find(this_lg) == (*ctg2lg)[this_ctg].end() )
            {
                (*ctg2lg)[this_ctg].insert(std::pair<string, string>(this_lg, this_homLG));
            }
        }        
        //
        if( (*homLGs).find(this_lg) == (*homLGs).end() )
        {
            (*homLGs).insert(std::pair<string, string>(this_lg, this_homLG) );
        }else
        {
            assert( this_homLG.compare( (*homLGs)[this_lg] ) == 0);
        }
        //
    }
    bool check_this = true;
    if(check_this && verbose)
    {
        cout << "   check: ctg win marker distribution in LGs: " << endl;
        map<string, map<unsigned long, winMARKER > >::iterator mcitr;
        map<string, map<unsigned long, winMARKER > >::iterator mcitr_end;
        mcitr     = (*gbmarker).begin();
        mcitr_end = (*gbmarker).end();
        while(mcitr != mcitr_end)
        {
            string this_ctg = (*mcitr).first;
            map<unsigned long, winMARKER > tmp_win_list = (*mcitr).second;
            map<unsigned long, winMARKER >::iterator mpitr;
            map<unsigned long, winMARKER >::iterator mpitr_end;
            mpitr     = tmp_win_list.begin();
            mpitr_end = tmp_win_list.end();
            while(mpitr != mpitr_end)
            {
                unsigned long win_sta = (*mpitr).first;
                winMARKER tmp_marker  = (*mpitr).second;
                cout << "        : " << this_ctg 
                     << "\t"         <<  win_sta 
                     << "\t"         << tmp_marker.end 
                     << "\t"         << tmp_marker.type 
                     << "\t"         << tmp_marker.homlg
                     << "\t";
                vector<string> this_LGs = tmp_marker.wlg;
                vector<string>::iterator lgitr;
                vector<string>::iterator lgitr_end;
                lgitr     = this_LGs.begin();
                lgitr_end = this_LGs.end();
                while(lgitr != lgitr_end)
                {
                    if(lgitr != this_LGs.begin())
                    {
                        cout << "," << *lgitr;
                    }else
                    {
                        cout <<        *lgitr;
                    }
                    lgitr ++;
                }
                cout << endl;                
                //
                mpitr ++;
            }
            //
            mcitr ++;
        }     
        //
        cout << "   check: LG distribution in CTGs: " << endl;        
        map<string, map<string, string> >::iterator citr;
        map<string, map<string, string> >::iterator citr_end;
        citr     = (*ctg2lg).begin();
        citr_end = (*ctg2lg).end();
        while(citr != citr_end)
        {
            cout << "   check: CTG=" << (*citr).first << " in LGs: ";
            map<string, string>::iterator lgitr2;
            map<string, string>::iterator lgitr2_end;
            lgitr2     = ((*citr).second).begin();
            lgitr2_end = ((*citr).second).end();
            while(lgitr2 != lgitr2_end)
            {
                if(lgitr2 != ((*citr).second).begin() )
                {
                    cout << "," << (*lgitr2).first;
                }else
                {
                    cout <<        (*lgitr2).first;                    
                }
                lgitr2 ++;
            }
            cout << endl;
            //
            citr ++;
        }
        //
        cout << "   check:  homologous LGs distribution 1-48 versus 1-12: " << endl;
        map<string, string>::iterator hlgitr;
        map<string, string>::iterator hlgitr_end;
        hlgitr     = (*homLGs).begin();
        hlgitr_end = (*homLGs).end();
        while(hlgitr != hlgitr_end)
        {
            cout << "        : lg = " << (*hlgitr).first << " in " << (*hlgitr).second << endl;
            hlgitr ++;
        }
    }
    //
    ifp.close();
    //
    return true;     
}                                         




                               
