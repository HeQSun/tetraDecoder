/*
    Given     
        a list of linkage-group wise window markers,
        a bam/cram file with long reads aligned,        
    separate pacbio reads into 2*n clusters (as indicated by window marker info).    
    
    Note: there is 500 bp cutoff to consider a window marker overlapping a long read.
    Note: some contigs in bam might not be in marker file, due to missed alignment to ref genome.
          they need to be recovered with hic contact after chrs' assembly
    
    Written by: Hequan Sun
    Address...: Carl-von-linne-weg 10, 50829 Koeln (MPIPZ)
    Email.....: sunhequan@gmail.com/sun@mpipz.mpg.de
    Date......: 2022-12-30
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
// window marker
struct winMARKER
{
//  string         ctg; // ctd id
//  unsigned long  sta; // window marker start 
    unsigned long  end; // window marker end
    string        type; // hap/dip/trip/tetrap/rep
    vector<string> wlg; // linkage group id 1:48
    string       homlg; // homo-lg id...... 1:12
};
double minCOscore  = 0.64; // 
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
//
int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        // g++ long_read_separator.cpp split_string.cpp -O3 -o long_read_separator
        cout << "\n   Function: separate long reads using genotypes according to phased contig-snp/del markers. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "samtools view read_alignment.bam[cram] | long_read_separator_v2 - win_marker.txt out_prefix_str "
             << endl            << endl;
        return 1;
    }
    double startT= clock();
    //
    string checkflag      = (string)argv[1];
    if(checkflag.compare("-") != 0)
    {
        cout << "\n   Function: separate long reads using genotypes according to phased contig-snp/del markers. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "samtools view read_alignment.bam[cram] | long_read_separator_v2 - win_marker.txt out_prefix_str "
             << endl            << endl;
        cout << "   Error: \"pacbio_genotyper -\" must be provided. " << endl << endl;
        return 1;
    }
    string gbmarker_file  = (string)argv[2];
    string outprefix      = (string)argv[3];
    //
    cout << "   Info: reading gamete_binning grouped window marker info.. " << endl;
    map<string, map<unsigned long, winMARKER > > gbmarker;// <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
    map<string, map<string, string> >            ctg2lg;  // <ctg_id, <lg_id=1:48, 1:12> >
    map<string, string>                          homLGs;  // <lg_id=1:48,  homLG_id=1:12>
    if(!read_gb_window_marker_grouping(gbmarker_file, &gbmarker, &ctg2lg,& homLGs))
    {
        return 1;
    }
    cout << "   Info: reading gamete_binning grouped window marker done. "  << endl;    
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
    // file handles on output files: <infilename+"_PPP/MMM_pbreads.fa", *ofp >    
    map<string, ofstream*> allofp;         
    // collecting group-specific unique read names: <infilename+"_PPP/MMM_pbreads.fa", <pbreadname, cnt> >               
    map<string, map<string, int> > groupedpbreadnames;  
    // no real data - just for initializing groupedpbreadnames.    
    map<string, int> init_tmp;     
    map<string, string>::iterator grpfileitr;
    map<string, string>::iterator grpfileitr_end;
    grpfileitr     = homLGs.begin(); 
    grpfileitr_end = homLGs.end();
    while(grpfileitr != grpfileitr_end)
    {
        string this_lg    = (*grpfileitr).first; // -1, 1:48
        string this_homLG = (*grpfileitr).second;// 1:12
        if(this_lg.compare("-1") == 0)
        {
             grpfileitr ++;
             continue;
        }
        //
        string this_lg_readfile = "homLG_" + this_homLG + "_LG_" + this_lg + "_reads.fa";
        ofstream *ofp = new ofstream( (tmpfolder + "/" + this_lg_readfile).c_str(), ios::out);
        if(!(*ofp).good())
        {
            cout << "   Error: cannot set up output file for " << tmpfolder + "/" + this_lg_readfile << endl;
            return 1;
        }
        allofp.insert(std::pair<string, ofstream*>(this_lg_readfile, ofp));
        groupedpbreadnames.insert(std::pair<string, map<string, int> >(this_lg_readfile, init_tmp));
        //
        grpfileitr ++;
    } 
    // pacbio reads collected    
    map<string, int> pbreadnames;                    // all collected pb read names (unique of 7043689)
    map<string, vector<string> > pbreadnames_detail; // all collected pb read names (unique of 7043689) <readname, <grp, align_len> >   
    map<string, int> pbreadnames_nolg;     // all collected pb read names without linkage group info (unique of )
    map<string, int> pbreadnames_lg;       // all collected pb read names with    linkage group info and with    marker (unique of )
    map<string, int> pbreadnames_lg_nomkr; // all collected pb read names with    linkage group info and without marker (unique of )    
    // unmapped pacbio reads with cigar as "*"
    string unmappedfile   = tmpfolder + "/unmapped_starcigar_pb_reads.fa"; 
    ofstream *ofp         = new ofstream(unmappedfile.c_str(), ios::out);
    allofp.insert(std::pair<string, ofstream*>("unmapped_starcigar_pb_reads.fa", ofp));
    groupedpbreadnames.insert(std::pair<string, map<string, int> >("unmapped_starcigar_pb_reads.fa", init_tmp));    
    // mapped pacbio reads but no linkage group info
    string mappednolgfile = tmpfolder + "/mapped_no_lg_pb_reads.fa"; 
    ofstream *nolgofp     = new ofstream(mappednolgfile.c_str(), ios::out);
    allofp.insert(std::pair<string, ofstream*>("mapped_no_lg_pb_reads.fa", nolgofp));  
    groupedpbreadnames.insert(std::pair<string, map<string, int> >("mapped_no_lg_pb_reads.fa", init_tmp));              
    //
    unsigned long lenNoCG  = 0; // no cigar string unmapped - length
    unsigned long numNoCG  = 0; // no cigar string unmapped - number
    string exampleNoCG     = "";
    unsigned long numNoMA  = 0; // number of reads showing multiple alignment: only one kept.
    string exampleNoMA     = "";
    unsigned long numraw   = 0;    
    string exampleNoSEQ    = "";
    unsigned long numNoSEQ = 0;
    unsigned long lengrped = 0; // pb alignments clearly related to linkage groups - length 
    unsigned long grpedNo  = 0; // pb alignments clearly related to linkage groups; note 1 read may have >1 alignments
    unsigned long lennonlg = 0; // pb alignments not related to linkage groups - length       
    unsigned long nonlgNo  = 0; // pb alignments not related to linkage groups - number 
    unsigned long ctg_miss = 0; // contigs in bam/cram but not in marker list (which are not aligned to the ref genome)
    // m64078_200122_223326/73138406/ccs   16  utg000001l_pilon    458658  60  9555M1D5579M3S  *   0   0   GTGAAAGTTGTTGTCGTCGTTCTAGCTCATCG    
    srand (2); // to randomly separate reads if the marker relates to several linkage groups of 1:48
    std::string line;    
    while (std::getline(std::cin, line)) 
    {
        if(line.size()==0) continue;
        if(line[0] == '@') continue;
        //
        numraw ++;
        if(numraw%1000000 == 0)
        {
            cout << "   Info: " << numraw << "th aligm..." << endl;
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
        string thisctg   = lineinfo[2]; // reference contig name
        string thispos   = lineinfo[3]; // first matching reference coordinate
        string thiscigar = lineinfo[5]; // cigar string
        string thisseq   = lineinfo[9]; // read sequence      
        // special case 1: cigar string as star: not aligned
        if(thiscigar.compare("*") == 0) 
        {
            // READ SET I
            (*allofp["unmapped_starcigar_pb_reads.fa"])  << ">" << pbname  << endl;
            (*allofp["unmapped_starcigar_pb_reads.fa"])  <<        thisseq << endl;
            //
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
            numNoCG ++;
            lenNoCG  += thisseq.size();
            //
            if(groupedpbreadnames["unmapped_starcigar_pb_reads.fa"].find(pbname) == 
               groupedpbreadnames["unmapped_starcigar_pb_reads.fa"].end())
            {
                groupedpbreadnames["unmapped_starcigar_pb_reads.fa"].insert(std::pair<string, int>(pbname, 1));
            }else
            {
                groupedpbreadnames["unmapped_starcigar_pb_reads.fa"][pbname] += 1;
            }
            //     
            continue;
        }    
        // special case 2: skip secondary/
        int hexflag = strtol(thisflag.c_str(), NULL, 0);
        if(hexflag >= 256)
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
        // collect all informative reads
        if(pbreadnames.find(pbname) == pbreadnames.end())
        {
            pbreadnames.insert(std::pair<string, int>(pbname, 1));
        }else
        {
            pbreadnames[pbname] += 1;
        }
        cout << endl;
        string alignedtogrp = "";
        cout << "   Check: read: " << pbname << " aligned to " << thisctg << "\t" << thispos << endl;                     
        // If it contains 0x10, SEQ in BAM/SAM has been reversed+complemented from read in fastq
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
        cout << "   Check: original read in bam/sam: " 
             << thisseq.size()      << " bp; " << endl;
        cout << "   Check: updated read with \'I\' removed, and \'D:-\' added: " 
             << updated_read.size() << " bp. " << endl;            
        // map<string, map<unsigned long, winMARKER > > gbmarker_updated; // <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
        // map<string, map<string, string> >            ctg2lg; // <ctg_id, <lg_id=1:48, 1:12> >
        // map<string, string>                          homLGs; // <lg_id=1:48,  homLG_id=1:12>           
        if( gbmarker_updated.find(thisctg) == gbmarker_updated.end() )
        {
            cout << "   warning: " << thisctg << " not in marker list. " << endl;
            string this_lg_readfile = "mapped_no_lg_pb_reads.fa";
            (*allofp[this_lg_readfile]) << ">"   
                                   << pbname 
                                   << " + lg="
                                   << "-1/48"
                                   << " homLG="
                                   << "-1/12"
                                   << " rep.span="         
                                   << thisctg 
                                   << ":"   
                                   << firstRefsMatch   
                                   << "-"          
                                   << spansecond
                                   << endl;   // read name  
            (*allofp[this_lg_readfile]) << thisseq         
                                   << endl;   // read seq
            nonlgNo  ++;           
            lennonlg += thisseq.size();  
            // without lg
            if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
            {
                pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
            } 
            //
            if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
               groupedpbreadnames[this_lg_readfile].end())
            {
                groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
            }else
            {
                groupedpbreadnames[this_lg_readfile][pbname] += 1;
            }
            //
            alignedtogrp = this_lg_readfile;              
            //
            ctg_miss ++;
            continue;
        }
        map<unsigned long, winMARKER > this_ctg_win_list = gbmarker_updated[thisctg]; // can be several large windows
        vector<string> best_overlap_lg;
        unsigned long  best_overlap_len  = 0;
        string         best_overlap_type = "";
        map<unsigned long, winMARKER >::iterator winitr;
        map<unsigned long, winMARKER >::iterator winitr_end;
        winitr     = this_ctg_win_list.begin();
        winitr_end = this_ctg_win_list.end();
        while(winitr != winitr_end)
        {
            unsigned long this_win_sta = (*winitr).first;
            winMARKER tmp_pos_marker   = (*winitr).second;
            unsigned long this_win_end = tmp_pos_marker.end;                        
            //
            if(this_win_sta > spansecond)
            {
                break;
            }
            if(this_win_end < firstRefsMatch)
            {
                winitr ++;
                continue;
            }
            //
            unsigned long tmp_overlap_len = overlap_read_and_win_marker_v2(this_win_sta, 
                                                                           this_win_end, 
                                                                           firstRefsMatch, 
                                                                           spansecond,
                                                                           thisctg);   
            if(tmp_overlap_len > best_overlap_len)
            {
                best_overlap_len  = tmp_overlap_len;
                best_overlap_lg   = tmp_pos_marker.wlg;
                best_overlap_type = tmp_pos_marker.type;
                cout << "   Check: best matching window marker updated as " << this_win_sta 
                     << "-"                                                 << this_win_end
                     << endl; 
            }
            //
            winitr ++;
        }
        // non-signal reads
        if(best_overlap_type.compare("rep") == 0)
        {
            string this_lg_readfile = "mapped_no_lg_pb_reads.fa";
            (*allofp[this_lg_readfile]) << ">"   
                                   << pbname 
                                   << " + lg="
                                   << "-1/48"
                                   << " homLG="
                                   << "-1/12"
                                   << " rep.span="         
                                   << thisctg 
                                   << ":"   
                                   << firstRefsMatch   
                                   << "-"          
                                   << spansecond
                                   << endl;   // read name  
            (*allofp[this_lg_readfile]) << thisseq         
                                   << endl;   // read seq
            nonlgNo  ++;           
            lennonlg += thisseq.size();  
            // without lg
            if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
            {
                pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
            } 
            //
            if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
               groupedpbreadnames[this_lg_readfile].end())
            {
                groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
            }else
            {
                groupedpbreadnames[this_lg_readfile][pbname] += 1;
            }
            //
            alignedtogrp = this_lg_readfile;              
            //                                                           
            continue; // caution: next read from bam                       
        }
        // signal reads
        string this_wlg("");
        for(int jj = 0; jj < best_overlap_lg.size(); jj ++)
        {
            if(jj > 0)
            this_wlg += ",";
            this_wlg += best_overlap_lg[jj];
        }
        double randra = rand()%1000000000*1.0 / 1000000000.0;
        if(best_overlap_lg.size() == 1)
        {
            string tmp_goto_lg = best_overlap_lg[0];
            assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
            string this_homLG = homLGs[tmp_goto_lg];
            cout << "   Check: read "  
                 << pbname      << " would go to LG " 
                 << this_wlg    << ", and selected "
                 << tmp_goto_lg << "/48 in " 
                 << this_homLG  << "/12" 
                 << endl;            
            string this_lg_readfile = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads.fa";
            if(tmp_goto_lg.compare("-1") == 0)
            {
                this_lg_readfile = "mapped_no_lg_pb_reads.fa"; 
                nonlgNo  ++;           
                lennonlg += thisseq.size();            
                // without lg
                if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
                {
                    pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
                }                         
            }  
            (*allofp[this_lg_readfile]) << ">"   
                                   << pbname 
                                   << " + lg="
                                   << tmp_goto_lg << "/48"
                                   << " homLG="
                                   << this_homLG  << "/12"
                                   << " span="         
                                   << thisctg 
                                   << ":"   
                                   << firstRefsMatch   
                                   << "-"          
                                   << spansecond
                                   << endl;   // read name  
            (*allofp[this_lg_readfile]) << thisseq         
                                   << endl;   // read seq   
            //
            if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
               groupedpbreadnames[this_lg_readfile].end())
            {
                groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
            }else
            {
                groupedpbreadnames[this_lg_readfile][pbname] += 1;
            } 
            //
            alignedtogrp = this_lg_readfile;                                                       
        }else
        if(best_overlap_lg.size() == 2)
        {
            if(randra <= 0.5)
            {
                string tmp_goto_lg = best_overlap_lg[0];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                cout << "   Check: read "  
                     << pbname      << " would go to LG " 
                     << this_wlg    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in [0, 0.5]"
                     << endl;            
                string this_lg_readfile = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads.fa";
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile = "mapped_no_lg_pb_reads.fa"; 
                    nonlgNo  ++;           
                    lennonlg += thisseq.size();      
                    // without lg
                    if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
                    {
                        pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
                    }                                      
                }  
                (*allofp[this_lg_readfile]) << ">"   
                                       << pbname 
                                       << " + lg="
                                       << tmp_goto_lg << "/48"
                                       << " homLG="
                                       << this_homLG  << "/12"
                                       << " span="         
                                       << thisctg 
                                       << ":"   
                                       << firstRefsMatch   
                                       << "-"          
                                       << spansecond
                                       << " "
                                       << randra
                                       << endl;   // read name  
                (*allofp[this_lg_readfile]) << thisseq         
                                       << endl;   // read seq 
                //
                if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
                   groupedpbreadnames[this_lg_readfile].end())
                {
                    groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
                }else
                {
                    groupedpbreadnames[this_lg_readfile][pbname] += 1;
                } 
                //
                alignedtogrp = this_lg_readfile;                                                      
            }else
            {
                string tmp_goto_lg = best_overlap_lg[1];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                cout << "   Check: read "  
                     << pbname      << " would go to LG " 
                     << this_wlg    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (0.5, 1]"                     
                     << endl;            
                string this_lg_readfile = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads.fa";
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile = "mapped_no_lg_pb_reads.fa"; 
                    nonlgNo  ++;           
                    lennonlg += thisseq.size();        
                    // without lg
                    if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
                    {
                        pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
                    }                                    
                }  
                (*allofp[this_lg_readfile]) << ">"   
                                       << pbname 
                                       << " + lg="
                                       << tmp_goto_lg << "/48"
                                       << " homLG="
                                       << this_homLG  << "/12"
                                       << " span="         
                                       << thisctg 
                                       << ":"   
                                       << firstRefsMatch   
                                       << "-"          
                                       << spansecond
                                       << " "
                                       << randra                                       
                                       << endl;   // read name  
                (*allofp[this_lg_readfile]) << thisseq         
                                       << endl;   // read seq 
                //
                if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
                   groupedpbreadnames[this_lg_readfile].end())
                {
                    groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
                }else
                {
                    groupedpbreadnames[this_lg_readfile][pbname] += 1;
                }
                //
                alignedtogrp = this_lg_readfile;                                                       
            }                       
        }else
        if(best_overlap_lg.size() == 3)
        {
            if(randra <= 1.0/3)
            {
                string tmp_goto_lg = best_overlap_lg[0];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                cout << "   Check: read "  
                     << pbname      << " would go to LG " 
                     << this_wlg    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in [0, 1/3]"                     
                     << endl;            
                string this_lg_readfile = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads.fa";
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile = "mapped_no_lg_pb_reads.fa"; 
                    nonlgNo  ++;           
                    lennonlg += thisseq.size();   
                    // without lg
                    if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
                    {
                        pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
                    }                                         
                }  
                (*allofp[this_lg_readfile]) << ">"   
                                       << pbname 
                                       << " + lg="
                                       << tmp_goto_lg << "/48"
                                       << " homLG="
                                       << this_homLG  << "/12"
                                       << " span="         
                                       << thisctg 
                                       << ":"   
                                       << firstRefsMatch   
                                       << "-"          
                                       << spansecond
                                       << " "
                                       << randra
                                       << endl;   // read name  
                (*allofp[this_lg_readfile]) << thisseq         
                                       << endl;   // read seq 
                //
                if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
                   groupedpbreadnames[this_lg_readfile].end())
                {
                    groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
                }else
                {
                    groupedpbreadnames[this_lg_readfile][pbname] += 1;
                }  
                //
                alignedtogrp = this_lg_readfile;                                                     
            }else
            if(1.0/3<randra && randra<=2.0/3)            
            {
                string tmp_goto_lg = best_overlap_lg[1];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                cout << "   Check: read "  
                     << pbname      << " would go to LG " 
                     << this_wlg    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (1/3, 2/3]"                     
                     << endl;            
                string this_lg_readfile = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads.fa";
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile = "mapped_no_lg_pb_reads.fa"; 
                    nonlgNo  ++;           
                    lennonlg += thisseq.size();  
                    // without lg
                    if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
                    {
                        pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
                    }                                          
                }  
                (*allofp[this_lg_readfile]) << ">"   
                                       << pbname 
                                       << " + lg="
                                       << tmp_goto_lg << "/48"
                                       << " homLG="
                                       << this_homLG  << "/12"
                                       << " span="         
                                       << thisctg 
                                       << ":"   
                                       << firstRefsMatch   
                                       << "-"          
                                       << spansecond
                                       << " "
                                       << randra                                       
                                       << endl;   // read name  
                (*allofp[this_lg_readfile]) << thisseq         
                                       << endl;   // read seq 
                //
                if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
                   groupedpbreadnames[this_lg_readfile].end())
                {
                    groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
                }else
                {
                    groupedpbreadnames[this_lg_readfile][pbname] += 1;
                }
                //
                alignedtogrp = this_lg_readfile;                                                       
            }else
            {
                string tmp_goto_lg = best_overlap_lg[2];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                cout << "   Check: read "  
                     << pbname      << " would go to LG " 
                     << this_wlg    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (2/3, 1.0]"                     
                     << endl;            
                string this_lg_readfile = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads.fa";
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile = "mapped_no_lg_pb_reads.fa"; 
                    nonlgNo  ++;           
                    lennonlg += thisseq.size();  
                    // without lg
                    if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
                    {
                        pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
                    }                                          
                }  
                (*allofp[this_lg_readfile]) << ">"   
                                       << pbname 
                                       << " + lg="
                                       << tmp_goto_lg << "/48"
                                       << " homLG="
                                       << this_homLG  << "/12"
                                       << " span="         
                                       << thisctg 
                                       << ":"   
                                       << firstRefsMatch   
                                       << "-"          
                                       << spansecond
                                       << " "
                                       << randra                                       
                                       << endl;   // read name  
                (*allofp[this_lg_readfile]) << thisseq         
                                       << endl;   // read seq 
                //
                if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
                   groupedpbreadnames[this_lg_readfile].end())
                {
                    groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
                }else
                {
                    groupedpbreadnames[this_lg_readfile][pbname] += 1;
                } 
                //
                alignedtogrp = this_lg_readfile;                                                      
            }                       
        }else
        if(best_overlap_lg.size() == 4)
        {
            if(randra <= 1.0/4)
            {
                string tmp_goto_lg = best_overlap_lg[0];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                cout << "   Check: read "  
                     << pbname      << " would go to LG " 
                     << this_wlg    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in [0, 1/4]"                     
                     << endl;            
                string this_lg_readfile = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads.fa";
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile = "mapped_no_lg_pb_reads.fa"; 
                    nonlgNo  ++;           
                    lennonlg += thisseq.size();      
                    // without lg
                    if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
                    {
                        pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
                    }                                      
                }  
                (*allofp[this_lg_readfile]) << ">"   
                                       << pbname 
                                       << " + lg="
                                       << tmp_goto_lg << "/48"
                                       << " homLG="
                                       << this_homLG  << "/12"
                                       << " span="         
                                       << thisctg 
                                       << ":"   
                                       << firstRefsMatch   
                                       << "-"          
                                       << spansecond
                                       << " "
                                       << randra
                                       << endl;   // read name  
                (*allofp[this_lg_readfile]) << thisseq         
                                       << endl;   // read seq 
                //
                if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
                   groupedpbreadnames[this_lg_readfile].end())
                {
                    groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
                }else
                {
                    groupedpbreadnames[this_lg_readfile][pbname] += 1;
                } 
                //
                alignedtogrp = this_lg_readfile;                                                      
            }else
            if(1.0/4<randra && randra<=2.0/4)            
            {
                string tmp_goto_lg = best_overlap_lg[1];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                cout << "   Check: read "  
                     << pbname      << " would go to LG " 
                     << this_wlg    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (1/4, 2/4]"                     
                     << endl;            
                string this_lg_readfile = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads.fa";
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile = "mapped_no_lg_pb_reads.fa"; 
                    nonlgNo  ++;           
                    lennonlg += thisseq.size();   
                    // without lg
                    if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
                    {
                        pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
                    }                     
                }  
                (*allofp[this_lg_readfile]) << ">"   
                                       << pbname 
                                       << " + lg="
                                       << tmp_goto_lg << "/48"
                                       << " homLG="
                                       << this_homLG  << "/12"
                                       << " span="         
                                       << thisctg 
                                       << ":"   
                                       << firstRefsMatch   
                                       << "-"          
                                       << spansecond
                                       << " "
                                       << randra                                       
                                       << endl;   // read name  
                (*allofp[this_lg_readfile]) << thisseq         
                                       << endl;   // read seq 
                //
                if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
                   groupedpbreadnames[this_lg_readfile].end())
                {
                    groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
                }else
                {
                    groupedpbreadnames[this_lg_readfile][pbname] += 1;
                } 
                //
                alignedtogrp = this_lg_readfile;                                                      
            }else
            if(2.0/4<randra && randra<=3.0/4)            
            {
                string tmp_goto_lg = best_overlap_lg[2];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                cout << "   Check: read "  
                     << pbname      << " would go to LG " 
                     << this_wlg    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (2/4, 3/4]"                     
                     << endl;            
                string this_lg_readfile = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads.fa";
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile = "mapped_no_lg_pb_reads.fa"; 
                    nonlgNo  ++;           
                    lennonlg += thisseq.size(); 
                    // without lg
                    if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
                    {
                        pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
                    }                                           
                }  
                (*allofp[this_lg_readfile]) << ">"   
                                       << pbname 
                                       << " + lg="
                                       << tmp_goto_lg << "/48"
                                       << " homLG="
                                       << this_homLG  << "/12"
                                       << " span="         
                                       << thisctg 
                                       << ":"   
                                       << firstRefsMatch   
                                       << "-"          
                                       << spansecond
                                       << " "
                                       << randra                                       
                                       << endl;   // read name  
                (*allofp[this_lg_readfile]) << thisseq         
                                       << endl;   // read seq 
                //
                if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
                   groupedpbreadnames[this_lg_readfile].end())
                {
                    groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
                }else
                {
                    groupedpbreadnames[this_lg_readfile][pbname] += 1;
                } 
                //
                alignedtogrp = this_lg_readfile;                                                      
            }else
            {
                string tmp_goto_lg = best_overlap_lg[3];
                assert(homLGs.find(tmp_goto_lg) != homLGs.end() );
                string this_homLG = homLGs[tmp_goto_lg];
                cout << "   Check: read "  
                     << pbname      << " would go to LG " 
                     << this_wlg    << ", and selected "                     
                     << tmp_goto_lg << "/48 in " 
                     << this_homLG  << "/12 " 
                     << randra      << " in (3/4, 1.0]"                     
                     << endl;            
                string this_lg_readfile = "homLG_" + this_homLG + "_LG_" + tmp_goto_lg + "_reads.fa";
                if(tmp_goto_lg.compare("-1") == 0)
                {
                    this_lg_readfile = "mapped_no_lg_pb_reads.fa"; 
                    nonlgNo  ++;           
                    lennonlg += thisseq.size(); 
                    // without lg
                    if(pbreadnames_nolg.find(pbname) == pbreadnames_nolg.end())
                    {
                        pbreadnames_nolg.insert(std::pair<string,int>(pbname, 1));
                    }                                           
                }  
                (*allofp[this_lg_readfile]) << ">"   
                                       << pbname 
                                       << " + lg="
                                       << tmp_goto_lg << "/48"
                                       << " homLG="
                                       << this_homLG  << "/12"
                                       << " span="         
                                       << thisctg 
                                       << ":"   
                                       << firstRefsMatch   
                                       << "-"          
                                       << spansecond
                                       << " "
                                       << randra                                       
                                       << endl;   // read name  
                (*allofp[this_lg_readfile]) << thisseq         
                                       << endl;   // read seq 
                //
                if(groupedpbreadnames[this_lg_readfile].find(pbname) == 
                   groupedpbreadnames[this_lg_readfile].end())
                {
                    groupedpbreadnames[this_lg_readfile].insert(std::pair<string, int>(pbname, 1));
                }else
                {
                    groupedpbreadnames[this_lg_readfile][pbname] += 1;
                } 
                //
                alignedtogrp = this_lg_readfile;                                                      
            } 
        } 
        // collect informative reads with details
        std::stringstream readss;
        readss.str("");
        //if(alignedtogrp.size() == 0) cout << "   Warning: this grp name is null: check " << pbname << endl;
        readss << alignedtogrp << ":" << spansecond-firstRefsMatch+1;      
        if(pbreadnames_detail.find(pbname) == pbreadnames_detail.end())
        {
            vector<string> alignvec;
            alignvec.push_back(readss.str());              
            pbreadnames_detail.insert(std::pair<string, vector<string> >(pbname, alignvec));
        }else
        {
            pbreadnames_detail[pbname].push_back(readss.str());
        }                                              
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
    cout << "   Info: in total a4=" 
         << pbreadnames.size() 
         << " reads from all=" 
         << numraw                  
         << " aligment lines collected (<-header line not counted; raw alignment including secondary etc - skipped in read sep). " 
         << endl;
    // no lg
    cout << "   Info: number of pb alignment seqs WITHOUT linkage Info: "   
         << nonlgNo       
         << ", totaling v2=" 
         << lennonlg/1000000000.0      
         << " Gb" 
         << endl;
    if(ctg_miss > 0)
    {
        cout << "       : including " << ctg_miss << " reads missed due to contigs not in marker list. " << endl;
    }         
    cout << "         (u0=" 
         << pbreadnames_nolg.size() 
         << " unique reads in this cluster of no lg or not grouped. )"    
         << endl; 
    // read statistics 
    // map<string, map<string, int> > groupedpbreadnames;    
    cout << endl 
         << "   Info: distribution of pacbio reads in linkage groups: "
         << endl << endl;
    map<string, map<string, int> >::iterator gpbitr;
    map<string, map<string, int> >::iterator gpbitr_end;
    gpbitr     = groupedpbreadnames.begin();
    gpbitr_end = groupedpbreadnames.end();
    while(gpbitr != gpbitr_end)
    {
        cout << "\t" << (*gpbitr).first 
             << "\t" << (*gpbitr).second.size() 
             << "\t";
        map<string, int>::iterator nameitr;
        map<string, int>::iterator nameitr_end;
        nameitr     = (*gpbitr).second.begin();
        nameitr_end = (*gpbitr).second.end(); 
        unsigned long redundcount = 0; // one read may have >1 alignments
        while(nameitr != nameitr_end)
        {
            redundcount += (*nameitr).second;
            nameitr ++;
        }
        cout << redundcount << endl;
        gpbitr ++;
    }                          
    // close output files
    map<string, ofstream*>::iterator ofileitr;
    map<string, ofstream*>::iterator ofileitr_end;
    ofileitr     = allofp.begin();
    ofileitr_end = allofp.end();
    while(ofileitr != ofileitr_end)
    {
        ofstream *ofp = (*ofileitr).second;
        (*ofp).close();
        //
        ofileitr ++;
    }    
    cout << "   Info: extract reads from bam/sam into linkage done. "     << endl;
    // for check purpose
    string ofilenamecheck = tmpfolder + "/zcheck_read_group_distribution.txt";
    ofstream checkofp;
    checkofp.open(ofilenamecheck.c_str(), ios::out);
    if(!checkofp.good())
    {
        cout << "   Error: cannot open file to write read-group info file " << ofilenamecheck << endl;
        return 1;
    }
    checkofp << "#pbreadname\tgroup_assigned:score[;...]\tnum_assignment\tbest_assigned_group" << endl;
    map<string, vector<string> >::iterator rgitr;
    map<string, vector<string> >::iterator rgitr_end;
    rgitr     = pbreadnames_detail.begin();
    rgitr_end = pbreadnames_detail.end();
    while(rgitr != rgitr_end)
    {
        unsigned long best_align_score = 0;
        string        best_align_group = "";        
        //cout << "   checking " << (*rgitr).first << endl;        
        checkofp << (*rgitr).first << "\t";
        vector<string> alignvec = (*rgitr).second;
        vector<string>::iterator alitr;
        vector<string>::iterator alitr_end;
        alitr     = alignvec.begin();
        alitr_end = alignvec.end();
        while(alitr != alitr_end)
        {
            checkofp << *alitr;
            //cout << "   checking " << *alitr << endl;
            vector<string> alinfo = split_string(*alitr, ':');
            unsigned long  alisco = strtoul(alinfo[1].c_str(), NULL, 0);
            if(alisco > best_align_score)
            {
                best_align_score = alisco;
                best_align_group = alinfo[0];
            }
            alitr ++;
            if(alitr != alitr_end)
            {
                checkofp << ";";
            }
        }
        checkofp << "\t" << alignvec.size() << "\t" << best_align_group << endl;
        rgitr ++;
    }
    checkofp.close();    
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;    
}
//
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
        Output:            length of current window marker overlapping with the current pacbio read.    
    */
    // score is the overlapping length between window marker and read
    unsigned long win_sta = this_win_sta;
    unsigned long win_end = this_win_end;        
    unsigned long olspan = 0; // overlapping span
    if(win_sta<=read_sta && read_sta<=win_end)
    {
        /*
                  sta---------win-makrer--------end         window marker span
                        sta-----read-----end                pacbio read case 1: right_min = read_end
                        sta-----read------------------end   pacbio read case 2: right_min = win_end
        */
        unsigned long right_min = read_end<win_end?read_end:win_end;
        olspan = right_min - read_sta + 1;
    }else
    if(win_sta<=read_end && read_end<=win_end)
    {
        /*
                  sta---------win-marker--------end         window marker span
            sta--------------read-----end                   pacbio read case 2: right_min = win_end
        */              
        unsigned long left_max = read_sta>win_sta?read_sta:win_sta;
        olspan = read_end - left_max + 1;
    }else
    if(win_sta>=read_sta && read_end>=win_end)
    {
        /*
                  sta---------win-marker--------end         window marker span
            sta--------------read---------------------end   pacbio read case 2: right_min = win_end
        */    
        olspan = win_end - win_sta + 1;         
    }
    // do we need a cutoff for considering overlapping as sufficient? current as least 500 bp.
    if(olspan >= 500 || (olspan>=250 && read_end - read_sta + 1 <= 1000))
    {
        cout << "          : window-maker " 
             << this_ctg                << ":" 
             << win_sta                 << "-"
             << win_end                 << " overlap pb read spanning "
             << read_sta                << "-" 
             << read_end                << " by " 
             << olspan                  << " bp. "
             << endl;
    }else
    {
        if(0)
        cout << "          : window-marker " 
             << this_ctg                << ":" 
             << win_sta                 << "-" 
             << win_end                 << " does not overlap pb read spaning "                               
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
        Output:            window markers overlapping with the current pacbio read.    
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
                            sta-----read-----end                pacbio read case 1: right_min = read_end
                            sta-----read------------------end   pacbio read case 2: right_min = win_end
            */
            unsigned long right_min = read_end<win_end?read_end:win_end;
            olspan = right_min - read_sta + 1;
        }else
        if(win_sta<=read_end && read_end<=win_end)
        {
            /*
                      sta---------win-marker--------end         window marker span
                sta--------------read-----end                   pacbio read case 2: right_min = win_end
            */              
            unsigned long left_max = read_sta>win_sta?read_sta:win_sta;
            olspan = read_end - left_max + 1;
        }else
        if(win_sta>=read_sta && read_end>=win_end)
        {
            /*
                      sta---------win-marker--------end         window marker span
                sta--------------read---------------------end   pacbio read case 2: right_min = win_end
            */    
            olspan = win_end - win_sta + 1;         
        }
        // do we need a cutoff for considering overlapping as sufficient? current as least 500 bp.
        if(olspan >= 500 || (olspan>=250 && read_end - read_sta + 1 <= 1000))
        {
            (*this_win_marker_overlapped).insert(std::pair<unsigned long, unsigned long>((*witr).first, (*witr).second));
            cout << "          : window-maker " 
                 << this_ctg                << ":" 
                 << win_sta                 << "-"
                 << win_end                 << " overlap pb read spanning "
                 << read_sta                << "-" 
                 << read_end                << " by " 
                 << olspan                  << " bp. "
                 << endl;
        }else
        {
            if(1)
            cout << "          : window-marker " 
                 << this_ctg                << ":" 
                 << win_sta                 << "-" 
                 << win_end                 << " does not overlap pb read spaning "                               
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
       
       return..: gbmarker_updated: <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<1:48>, 1:12 } > >	
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
    if(check_this)
    {
        cout << "   CHECK: merged ctg win marker distribution in LGs: " << endl;
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
                 ctg2lg..: <ctg_id, <lg_id=1:48, 1:12> >
                 homLGs..: <lg_id=1:48,  homLG_id=1:12>
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
        string this_lg        = lineinfo[4]; // 1:48
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
    if(check_this)
    {
        cout << "   CHECK: ctg win marker distribution in LGs: " << endl;
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
        cout << "   CHECK: LG distribution in CTGs: " << endl;        
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
        cout << "   CHECK:  homologous LGs distribution 1-48 versus 1-12: " << endl;
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




                               
