/* this function  
       given a cnv_winsize1000_step1000_hq.txt and a fasta file, split ends of contigs based coverage analysis.
         from coverage patten 50.40.50.30.40.100.100.100.100.100.40.50.30.40.50.50.40x => 3 segments from the contig:
           50.40.50.30.40x
           100.100.100.100.100x
           40.50.30.40.50.50.40x
                     
   Written by Hequan Sun, MPIPZ
   Email: sunhequan@gmail.com/sun@mpipz.mpg.de
   Date.: 20220706     
*/
#include        <stdio.h>
#include       <string.h>
#include       <stdlib.h>
#include        <fstream>
#include       <iostream>
#include        <iomanip>
#include        <sstream>
#include         <vector>
#include            <map>
#include         <time.h>
#include       <assert.h>
#include "split_string.h"
//
bool find_cliping_site(string covfile, int min_cov, map<string, map<unsigned long, unsigned long> >* clip_regions);
bool read_and_split_contigs(string fafile, map<string, map<unsigned long, unsigned long> > clip_regions);
//
int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        printf("\nFunction: split contigs into [1, sta], (sta, end), [end, ctg_size] according to coverage analysis.\n");
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl << endl;
        printf("Usage: ctg_end_cov_splitter assembly.fasta cnv_winsize1000_step1000_hq.txt min_cov\n\n");
        cout << "       Note, cnv_winsize1000_step1000_hq.txt is output of CNV_HQ_v3."    << endl;
        cout << "             split windows from both ends of a contig until window depth > min_cov." << endl;
        cout << "        : after splitting, re-align reads and purge haplotigs using existing tools."
             << endl        
             << endl
             << endl;
        exit(1);
    }
    // get files
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    printf("\nSplitting ends of contigs started on %s\n", asctime(tinfo));    
    // step 0. read inputs
    string fafile  = (string)argv[1]; // file name of fasta file 
    string covfile = (string)argv[2]; // file name of coverage analysis
    int min_cov    = atoi(argv[3]);   // cutoff on selecting windows
    //
    // step 1. read coverage info and define splitting sites along contigs:
    //         format: ctg win-sta win-en 1 normalized_cnt raw_avg_cnt ctg-size
    map<string, map<unsigned long, unsigned long> > clip_regions;
    if(!find_cliping_site(covfile, min_cov, &clip_regions))
    {
        cout << "   Error: analyzing coverage depth of contigs failed. " << endl;
        return -1;
    }
    // check split status
    if(1)
    {
        map<string, map<unsigned long, unsigned long> >::iterator citr;
        map<string, map<unsigned long, unsigned long> >::iterator citr_end;
        citr     = clip_regions.begin();
        citr_end = clip_regions.end();
        while(citr != citr_end)
        {
            map<unsigned long, unsigned long> regions = (*citr).second;
            map<unsigned long, unsigned long>::iterator ritr;
            map<unsigned long, unsigned long>::iterator ritr_end;
            ritr = regions.begin();
            ritr_end = regions.end();
            while(ritr != ritr_end)
            {
                cout << (*citr).first << "\t" << (*ritr).first << "\t" << (*ritr).second << endl;
                //
                ritr ++;
            }
            //
            citr ++;
        }
    }
    // step 2. extract sequences from fasta 
    if(!read_and_split_contigs(fafile, clip_regions))
    {
        cout << "   Error: split the fasta file failed. " << endl;
        return -1;
    }
    //
    time_t endftime;
    struct tm* endtinfo;
    time(&endftime);
    endtinfo = localtime(&endftime);
    printf("\nSplitting ends of contigs successfully finished on %s\n", asctime(endtinfo));
    return 0;
}
//
bool read_and_split_contigs(string fafile, map<string, map<unsigned long, unsigned long> > clip_regions)
{
    ifstream ifp;
    ifp.open(fafile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open fasta file " << fafile << endl;
        return false;
    }
    //
    size_t psuffix = fafile.find(".fa");
    if(psuffix==std::string::npos)
    {
        cout << "   Error: please ensure fasta file with .fa or .fasta suffix." << endl;
        return false;
    }
    string ofilename      = fafile.substr(0, psuffix) + "_clipped.fa";
    string ofilename_size = fafile.substr(0, psuffix) + "_clipped.ctgsizes";       
    //
    ofstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open file for output: " << ofilename << endl;
        return false;
    }
    ofstream ofp_size;
    ofp_size.open(ofilename_size.c_str(), ios::out);
    if(!ofp_size.good())
    {
        cout << "   Error: cannot open file for output: " << ofilename_size << endl;
        return false;
    }    
    // traverse input fasta file
    std::string line("");
    getline(ifp, line);
    while(ifp.good())
    {
        if(line.find(">") != std::string::npos)
        {
            string seq_name = line.substr(1);
            string seq("");
            while(ifp.good())
            {
                std::string chrseq("");                
                getline(ifp, chrseq);
                if(chrseq.size()==0) continue;
                line    = chrseq;
                if(line.find(">")!=std::string::npos) break; // next seq 
                seq    += chrseq;
            }
            if(clip_regions.find(seq_name) != clip_regions.end())
            {
                map<unsigned long, unsigned long> regions = clip_regions[seq_name];
                map<unsigned long, unsigned long>::iterator ritr;
                map<unsigned long, unsigned long>::iterator ritr_end;
                ritr     = regions.begin();
                ritr_end = regions.end();
                int si   = 0;
                while(ritr != ritr_end)
                {
                    // get region of subseq
                    unsigned long sta = (*ritr).first  - 1; // 0-based
                    unsigned long end = (*ritr).second - 1; // 0-based
                    unsigned long len = end - sta + 1;
                    // get sub seq                    
                    string this_seq = seq.substr(sta, len);                    
                    //
                    si ++;
                    ofp_size << ""  << seq_name << "_" << si << "\t" << this_seq.size() << "\t" << sta+1 << "\t" << end+1 << endl; // window=1-based now
                    ofp      << ">" << seq_name << "_" << si << endl; // window=1-based now
                    //
                    int line_bp = 80; // 80 bp per line
                    int win_num = len / line_bp;
                    for(int j = 0; j < win_num; j ++)
                    {
                        ofp << this_seq.substr(j*line_bp, line_bp) << endl;
                    }
                    if(len > win_num * line_bp)
                    {
                        ofp << this_seq.substr(win_num*line_bp, len - win_num*line_bp) << endl;
                    }
                    //
                    //
                    ritr ++;
                } 
            }
        }
    }
    ifp.close();
    ofp.close();
    ofp_size.close();
    //
    return true;
}
//
bool find_cliping_site(string covfile, int min_cov, map<string, map<unsigned long, unsigned long> >* clip_regions)
{
    ifstream ifp;
    ifp.open(covfile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Warning: cannot open file " << covfile << endl;
        return false;
    }
    string            last_ctg = "";
    unsigned long last_ctg_siz = 0;
    map<string, double> ctg_cov;
    vector<string> ctg_win_order;
    // for each ctg, the {sta, end} should make the full contig clip_regions: <ctg, <sta, end> >
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 7)
        {
            cout << "   Warning: skipping insufficient line info at " << line << endl;
            continue;
        }
        //
        string this_ctg            = lineinfo[0];
        string this_win_sta_str    = lineinfo[1];
        string this_win_end_str    = lineinfo[2];
        string this_key            = this_ctg + "#" + this_win_sta_str + "#" + this_win_end_str;
        unsigned long this_win_sta = strtoul(this_win_sta_str.c_str(), NULL, 0);
        unsigned long this_win_end = strtoul(this_win_end_str.c_str(), NULL, 0);
        double        this_win_cov = atof(lineinfo[5].c_str());
        unsigned long this_win_siz = strtoul(lineinfo[6].c_str(), NULL, 0);
        //
        if(last_ctg.size() != 0 && last_ctg.compare(this_ctg) != 0)
        {
            // analyze last_ctg
            // cout << "   check: ctg " << last_ctg << " with " << ctg_cov.size() << " windows. " << endl;
            // case 1. left end of the contig
            vector<string>::iterator witr;
            vector<string>::iterator witr_end;
            witr     = ctg_win_order.begin();
            witr_end = ctg_win_order.end();
            bool win_low_cov = true;
            unsigned long clip_left_sta = 1;
            unsigned long clip_left_end = 0;
            while(witr != witr_end)
            {
                double tmp_win_cov = ctg_cov[(*witr)];
                // cout << "   check: " << *witr << ": " << tmp_win_cov << endl;
                if(win_low_cov && tmp_win_cov < min_cov)
                {
                    vector<string> keyinfo = split_string(*witr, '#');
                    clip_left_end = strtoul(keyinfo[2].c_str(), NULL, 0);
                }else
                {
                    win_low_cov = false; // region showing good coverage happens
                    //break;
                }
                witr ++;
            }
            // cout << "        : ctg " << last_ctg << " needs to split [" << clip_left_sta << ", " << clip_left_end << "] from left end." << endl;
            // case 2. right end of the contig
            vector<string>::reverse_iterator r_witr     = ctg_win_order.rbegin();
            vector<string>::reverse_iterator r_witr_end = ctg_win_order.rend();
            win_low_cov = true;
            unsigned long clip_right_sta = last_ctg_siz;
            unsigned long clip_right_end = last_ctg_siz;
            while(r_witr != r_witr_end)
            {
                double tmp_win_cov = ctg_cov[*r_witr];
                if(win_low_cov && tmp_win_cov < min_cov)
                {
                    vector<string> keyinfo = split_string(*r_witr, '#');                    
                    clip_right_sta = strtoul(keyinfo[1].c_str(), NULL, 0);
                }else
                {
                    win_low_cov = false; // region showing good coverage happens
                    //break;
                }
                r_witr ++;
            }
            // cout << "        : ctg " << last_ctg << " needs to split [" << clip_right_sta << ", " << clip_right_end << "] from right end." << endl; 
            // check if left and right clippings have overlaps
            if(clip_left_end>1 && clip_right_sta<last_ctg_siz && clip_right_sta>1)
            {
                // clipping from both ends
                map<unsigned long, unsigned long> tmp;
                tmp.insert(std::pair<unsigned long, unsigned long>(clip_left_sta, clip_left_end));
                (*clip_regions).insert(std::pair<string, map<unsigned long, unsigned long> >(last_ctg, tmp) );
                (*clip_regions)[last_ctg].insert(std::pair<unsigned long, unsigned long>(clip_left_end+1, clip_right_sta-1));
                (*clip_regions)[last_ctg].insert(std::pair<unsigned long, unsigned long>(clip_right_sta, clip_right_end));
            }else
            if(clip_left_end>1 && clip_left_end<last_ctg_siz && clip_right_sta==last_ctg_siz)
            {
                // clipping from left end only 
                map<unsigned long, unsigned long> tmp;
                tmp.insert(std::pair<unsigned long, unsigned long>(clip_left_sta, clip_left_end));
                (*clip_regions).insert(std::pair<string, map<unsigned long, unsigned long> >(last_ctg, tmp) );
                (*clip_regions)[last_ctg].insert(std::pair<unsigned long, unsigned long>(clip_left_end+1, clip_right_end));
            }else
            if(clip_left_end==0 && clip_right_sta<last_ctg_siz && clip_right_sta>1)
            {
                // clipping from right end only
                map<unsigned long, unsigned long> tmp;
                tmp.insert(std::pair<unsigned long, unsigned long>(clip_right_sta, clip_right_end));
                (*clip_regions).insert(std::pair<string, map<unsigned long, unsigned long> >(last_ctg, tmp) );
                (*clip_regions)[last_ctg].insert(std::pair<unsigned long, unsigned long>(clip_left_sta, clip_right_sta-1));
            }else
            {
                // no clipping: either whole contig low coverage (to purge later) or good coverage (to keep)
                map<unsigned long, unsigned long> tmp;
                tmp.insert(std::pair<unsigned long, unsigned long>(clip_left_sta, clip_right_end));
                (*clip_regions).insert(std::pair<string, map<unsigned long, unsigned long> >(last_ctg, tmp) );
            }                       
            //
            // update
            ctg_cov.clear();
            ctg_win_order.clear();
        }
        //
        ctg_cov.insert(std::pair<string, double>(this_key, this_win_cov));
        ctg_win_order.push_back(this_key); // order of windows
        last_ctg_siz = this_win_siz;
        last_ctg     = this_ctg;
    }
    // check last contig 
    // case 1. left end of the contig
    vector<string>::iterator witr;
    vector<string>::iterator witr_end;
    witr     = ctg_win_order.begin();
    witr_end = ctg_win_order.end();
    bool win_low_cov = true;
    unsigned long clip_left_sta = 1;
    unsigned long clip_left_end = 0;
    while(witr != witr_end)
    {
        double tmp_win_cov = ctg_cov[(*witr)];
        // cout << "   check: " << *witr << ": " << tmp_win_cov << endl;
        if(win_low_cov && tmp_win_cov < min_cov)
        {
            vector<string> keyinfo = split_string(*witr, '#');
            clip_left_end = strtoul(keyinfo[2].c_str(), NULL, 0);
        }else
        {
            win_low_cov = false; // region showing good coverage happens
            //break;
        }
        witr ++;
    }
    // cout << "        : ctg " << last_ctg << " needs to split [" << clip_left_sta << ", " << clip_left_end << "] from left end." << endl;
    // case 2. right end of the contig
    vector<string>::reverse_iterator r_witr     = ctg_win_order.rbegin();
    vector<string>::reverse_iterator r_witr_end = ctg_win_order.rend();
    win_low_cov = true;
    unsigned long clip_right_sta = last_ctg_siz;
    unsigned long clip_right_end = last_ctg_siz;
    while(r_witr != r_witr_end)
    {
        double tmp_win_cov = ctg_cov[*r_witr];
        if(win_low_cov && tmp_win_cov < min_cov)
        {
            vector<string> keyinfo = split_string(*r_witr, '#');                    
            clip_right_sta = strtoul(keyinfo[1].c_str(), NULL, 0);
        }else
        {
            win_low_cov = false; // region showing good coverage happens
            //break;
        }
        r_witr ++;
    }
    // cout << "        : ctg " << last_ctg << " needs to split [" << clip_right_sta << ", " << clip_right_end << "] from right end." << endl; 
    // check if left and right clippings have overlaps
    if(clip_left_end>1 && clip_right_sta<last_ctg_siz && clip_right_sta>1)
    {
        // clipping from both ends
        map<unsigned long, unsigned long> tmp;
        tmp.insert(std::pair<unsigned long, unsigned long>(clip_left_sta, clip_left_end));
        (*clip_regions).insert(std::pair<string, map<unsigned long, unsigned long> >(last_ctg, tmp) );
        (*clip_regions)[last_ctg].insert(std::pair<unsigned long, unsigned long>(clip_left_end+1, clip_right_sta-1));
        (*clip_regions)[last_ctg].insert(std::pair<unsigned long, unsigned long>(clip_right_sta, clip_right_end));
    }else
    if(clip_left_end>1 && clip_left_end<last_ctg_siz && clip_right_sta==last_ctg_siz)
    {
        // clipping from left end only 
        map<unsigned long, unsigned long> tmp;
        tmp.insert(std::pair<unsigned long, unsigned long>(clip_left_sta, clip_left_end));
        (*clip_regions).insert(std::pair<string, map<unsigned long, unsigned long> >(last_ctg, tmp) );
        (*clip_regions)[last_ctg].insert(std::pair<unsigned long, unsigned long>(clip_left_end+1, clip_right_end));
    }else
    if(clip_left_end==0 && clip_right_sta<last_ctg_siz && clip_right_sta>1)
    {
        // clipping from right end only
        map<unsigned long, unsigned long> tmp;
        tmp.insert(std::pair<unsigned long, unsigned long>(clip_right_sta, clip_right_end));
        (*clip_regions).insert(std::pair<string, map<unsigned long, unsigned long> >(last_ctg, tmp) );
        (*clip_regions)[last_ctg].insert(std::pair<unsigned long, unsigned long>(clip_left_sta, clip_right_sta-1));
    }else
    {
        // no clipping: either whole contig low coverage (to purge later) or good coverage (to keep)
        map<unsigned long, unsigned long> tmp;
        tmp.insert(std::pair<unsigned long, unsigned long>(clip_left_sta, clip_right_end));
        (*clip_regions).insert(std::pair<string, map<unsigned long, unsigned long> >(last_ctg, tmp) );
    }  
    //
    ifp.close();
    return true;
}










