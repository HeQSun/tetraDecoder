/*  this function take hic_binning grouped markers and a fasta file,
       extract/cut fasta sequences into groups for parental kmer validation.
    #
    #
    Written by: Hequan Sun
    Address...: Carl-von-linne-weg 10, 50829 Koeln (MPIPZ)
    Email.....: sunhequan@gmail.com/sun@mpipz.mpg.de
    Date......: 2022-11-04    
*/
#include        <map>
#include     <string>
#include     <vector>
#include    <fstream>
#include    <sstream>
#include   <iostream>
#include  <algorithm>
#include   <string.h>
#include   <stdlib.h>
#include   <assert.h>
#include     <math.h>
#include <sys/stat.h>
#include   <dirent.h>
#include   <iomanip>      // std::setprecision
#include "split_string.h"

bool read_hb_window_marker_grouping(string hbmarker_file, 
         map<string, map<string, map<unsigned long, unsigned long> > >* hbmarker,
         map<string, map<string, int> >*     ctg2lg);
bool merge_hb_window_marker_grouping(map<string, map<string, map<unsigned long, unsigned long> > >  hbmarker,
                                     map<string, map<string, map<unsigned long, unsigned long> > >* hbmarker_updated);  
bool extract_fasta( string fasta_file, 
                    map<string, map<string, map<unsigned long, unsigned long> > > hbmarker_updated,
                    map<string, map<string, int> > ctg2lg);                                            
//
using namespace std;
//
int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        cout << "\nFunction: extract fasta according gamete_binning grouped window markers. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "Usage: extract_fasta_with_hbMarker_grouping fasta gamete_binning_window_marker_grouping.txt"
             << endl << endl;
        return 1;
    }
    double startT= clock();
    //
    string fasta_file    = (string)argv[1];
    string hbmarker_file = (string)argv[2];
    //
    cout << "   info: reading gamete_binning grouped window marker info.. " << endl;    
    map<string, map<string, map<unsigned long, unsigned long> > > hbmarker; // <lg_id,  <ctg_id, <win_sta, win_end> > >
    map<string, map<string, int> > ctg2lg;                                  // <ctg_id, <lg_id, 1> >
    if( !read_hb_window_marker_grouping(hbmarker_file, &hbmarker, &ctg2lg))
    {
        cout << "   Error: cannot read gamete_binning grouped window marker info. " << endl;
        return 1;
    }
    cout << "   info: reading gamete_binning grouped window marker done. " << endl;    
    //
    cout << "   info: merging gamete_binning grouped window marker info.. " << endl;
    map<string, map<string, map<unsigned long, unsigned long> > > hbmarker_updated;
    if(!merge_hb_window_marker_grouping(hbmarker, &hbmarker_updated) )
    {
        cout << "   Error: cannot merge gamete_binning grouped window marker info. " << endl;
        return 1;        
    }
    cout << "   info: merging gamete_binning grouped window marker info done. " << endl;    
    //
    cout << "   info: extract ctg sequences from fasta into linkage groups... " << endl;
    if(!(extract_fasta( fasta_file, hbmarker_updated, ctg2lg)))
    {
        cout << "   Error: failed in extract ctg sequences from fasta to linkage groups. " << endl;
    }
    cout << "   info: extract ctg sequences from fasta into linkage done. "     << endl;    
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    //
    return 0;
}
//
bool extract_fasta( string fasta_file, 
                    map<string, map<string, map<unsigned long, unsigned long> > > hbmarker_updated,
                    map<string, map<string, int> > ctg2lg)
{
    /*
        function: extract ctg-window based sequences into linkage groups
          hbmarker; <lg_id,  <ctg_id, <win_sta, win_end> > >
          ctg2lg..; <ctg_id, <lg_id, 1> >                
    */
    // initialize 48+1 files to collect reads, where there is not-assigned group
    map<string, ofstream*> allofp;
    map<string, map<string, map<unsigned long, unsigned long> > >::iterator lgitr;
    map<string, map<string, map<unsigned long, unsigned long> > >::iterator lgitr_end;
    lgitr     = hbmarker_updated.begin();
    lgitr_end = hbmarker_updated.end();
    while(lgitr != lgitr_end)
    {
        string this_lg = (*lgitr).first; // note that '-1' can occur here, referring to not-assigned group!
        std::stringstream ss;
        ss.str("");
        ss << "ctg_window_sequence_for_lg_" << this_lg << ".fa";
        ofstream *ofp = new ofstream( (ss.str()).c_str(), ios::out);
        if(!(*ofp).good())
        {
            cout << "   Error: cannot set up output file for " << ss.str() << endl;
            return 1;
        }
        allofp.insert(std::pair<string, ofstream*>(ss.str(), ofp));
        //
        lgitr ++;
    }
    // collect reads from fasta
    ifstream ifp;
    ifp.open(fasta_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << fasta_file << endl;
        return false;
    }
    std::string line("");
    getline(ifp, line);    
    while(ifp.good())
    {
        if(line.size()==0) 
        {
            getline(ifp, line);
            continue;
        }
        if(line.find(">") != std::string::npos)
        {
            string this_ctg = line.substr(1);
            string seq("");
            int seq_line_id = 0;
            while(ifp.good())
            {
                std::string chrseq("");                
                getline(ifp, chrseq);
                if(chrseq.size()==0) continue;
                seq_line_id ++;
                line    = chrseq;
                if(line.find(">")!=std::string::npos) break; // next seq
                // if(seq_line_id > 1) seq += '\n'; // get one liner sequence 
                seq    += chrseq;
            }
            // extract regions of this seq to different linkage groups 
            cout << "   check-size: contig " << this_ctg << " " << seq.size() << endl;
            // find lg(s)
            if( ctg2lg.find(this_ctg) == ctg2lg.end() )
            {
                cout << "   warning: skipping " << this_ctg << " not in marker list. " << endl;
            }
            map<string, int> related_lgs = ctg2lg[this_ctg];
            map<string, int>::iterator olgitr;
            map<string, int>::iterator olgitr_end;
            olgitr     = related_lgs.begin();
            olgitr_end = related_lgs.end();
            while(olgitr != olgitr_end)
            {
                string this_lg = (*olgitr).first;
                string olgfilename = "ctg_window_sequence_for_lg_" + this_lg + ".fa";                
                // 
                assert( hbmarker_updated.find(this_lg) != hbmarker_updated.end() );
                map<string, map<unsigned long, unsigned long> > ctg_win_sta_end = hbmarker_updated[ this_lg ];
                //
                assert( ctg_win_sta_end.find(this_ctg) != ctg_win_sta_end.end() );
                map<unsigned long, unsigned long> win_sta_end = ctg_win_sta_end[this_ctg];
                map<unsigned long, unsigned long>::iterator witr;
                map<unsigned long, unsigned long>::iterator witr_end;
                witr     = win_sta_end.begin();
                witr_end = win_sta_end.end();
                while(witr != witr_end)
                {
                    unsigned long this_sta = (*witr).first;                    // 1-based system
                    unsigned long this_end = (*witr).second;
                    unsigned long this_len = this_end - this_sta + 1;
                    string extracted_seq   = seq.substr(this_sta-1, this_len); // 0-based system
                    //
                    assert( allofp.find(olgfilename) != allofp.end() );
                    (*allofp[olgfilename]) << ">" << this_ctg << "_" << this_sta << "_" << this_end << endl;
                    //
                    unsigned long line_num = ceil(this_len * 1.0 / 80);        // 80 bases per line 
                    cout << "   check: " << this_ctg 
                         << "_"          << this_sta 
                         << "_"          << this_end 
                         << " with "     << line_num 
                         << " lines. "   << endl;
                    for(int li = 0; li < line_num; li ++)
                    {
                        unsigned long line_seq_sta = li     * 80;    // 0-based
                        unsigned long line_seq_end = (li+1) * 80 - 1;// 0-based
                        if(line_seq_end > extracted_seq.size()-1)
                        {
                            cout << "   check: line end changed from " << line_seq_end << " to " << extracted_seq.size()-1 << endl;
                            line_seq_end = extracted_seq.size()-1;
                        }
                        unsigned long line_seq_len = line_seq_end - line_seq_sta + 1;
                        (*allofp[olgfilename]) << extracted_seq.substr(line_seq_sta, line_seq_len) << "\n";                        
                    }
                    //
                    witr ++;
                }
                //
                // 
                olgitr ++;
            }
            // 
        }
        else
        {
            getline(ifp, line);
        }        
    }
    ifp.close();
    //
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
    //
    return true;
}
//
bool merge_hb_window_marker_grouping(map<string, map<string, map<unsigned long, unsigned long> > >  hbmarker,
                                     map<string, map<string, map<unsigned long, unsigned long> > >* hbmarker_updated)
{
    /*
       function: merge smallers windows into large ones, if they are contigous
       e.g.,
          : LG=1, CTG=utg000001l_pilon, WIN=1-50000
	  : LG=1, CTG=utg000001l_pilon, WIN=50001-100000
	  : LG=1, CTG=utg000001l_pilon, WIN=100001-150000	  
	  : LG=1, CTG=utg000001l_pilon, WIN=300001-350000
	  : LG=1, CTG=utg000001l_pilon, WIN=350001-400000
	  will be merged as two larger windows:
          : LG=1, CTG=utg000001l_pilon, WIN=1-150000
	  : LG=1, CTG=utg000001l_pilon, WIN=300001-400000	  
    */
    bool merge_verbose = false;
    //
    map<string, map<string, map<unsigned long, unsigned long> > >::iterator lgitr;
    map<string, map<string, map<unsigned long, unsigned long> > >::iterator lgitr_end;
    lgitr     = hbmarker.begin();
    lgitr_end = hbmarker.end();
    while(lgitr != lgitr_end)
    {
        string this_lg = (*lgitr).first;
        map<string, map<unsigned long, unsigned long> > ctg_win_sta_end = (*lgitr).second; 
        //
        map<string, map<unsigned long, unsigned long> >::iterator citr;
        map<string, map<unsigned long, unsigned long> >::iterator citr_end;
        citr     = ctg_win_sta_end.begin();
        citr_end = ctg_win_sta_end.end();
        while(citr != citr_end)
        {
            string this_ctg = (*citr).first;
            map<unsigned long, unsigned long> win_sta_end = (*citr).second;
            map<unsigned long, unsigned long>::iterator witr;
            map<unsigned long, unsigned long>::iterator witr_end;
            witr     = win_sta_end.begin();
            witr_end = win_sta_end.end();
            unsigned long new_sta = 0;
            unsigned long new_end = 0;
            unsigned long last_sta = 0;
            unsigned long last_end = 0;            
            while(witr != witr_end)
            {
                unsigned long this_sta = (*witr).first;
                unsigned long this_end = (*witr).second;
                if(merge_verbose)
                cout << "   check: LG=" << this_lg 
                     << ", CTG="        << this_ctg 
                     << ", WIN="        << this_sta 
                     << "-"             << this_end 
                     << endl;
                if(last_sta == 0)
                {
                    last_sta = this_sta;
                    last_end = this_end;
                }else
                if(this_sta == last_end + 1)
                {
                    // connect 
                    last_end = this_end;
                }else
                {
                    // collect current merged widows
                    if((*hbmarker_updated).find(this_lg) == (*hbmarker_updated).end())
                    {
                        // window marker 
                        map<unsigned long, unsigned long> win_sta_end;
                        win_sta_end.insert(std::pair<unsigned long, unsigned long>(last_sta, last_end));
                        // ctg:: window marker 
                        map<string, map<unsigned long, unsigned long> > ctg_win_sta_end;
                        ctg_win_sta_end.insert(std::pair<string, map<unsigned long, unsigned long> >(this_ctg, win_sta_end));
                        // lg:: ctg:: window marker 
                        (*hbmarker_updated).insert(std::pair<string, map<string, map<unsigned long, unsigned long> > >(this_lg, ctg_win_sta_end) );
                    }else
                    {
                        if( (*hbmarker_updated)[this_lg].find(this_ctg) == (*hbmarker_updated)[this_lg].end() )
                        {
                            map<unsigned long, unsigned long> win_sta_end;
                            win_sta_end.insert(std::pair<unsigned long, unsigned long>(last_sta, last_end));      
                            (*hbmarker_updated)[this_lg].insert(std::pair<string, map<unsigned long, unsigned long> >(this_ctg, win_sta_end));                          
                        }else
                        {
                            (*hbmarker_updated)[this_lg][this_ctg].insert( std::pair<unsigned long, unsigned long>(last_sta, last_end) );
                        }
                    }
                    if(merge_verbose)
                    cout << "        : inter new WIN=" << last_sta
                         << "-"                        << last_end
                         << endl;
                    // initialize new large window 
                    last_sta = this_sta;
                    last_end = this_end;                    
                }
                // next window marker 
                witr ++;
            }
            // last window on current contig 
            // collect current merged widows
            if((*hbmarker_updated).find(this_lg) == (*hbmarker_updated).end())
            {
                // window marker 
                map<unsigned long, unsigned long> win_sta_end;
                win_sta_end.insert(std::pair<unsigned long, unsigned long>(last_sta, last_end));
                // ctg:: window marker 
                map<string, map<unsigned long, unsigned long> > ctg_win_sta_end;
                ctg_win_sta_end.insert(std::pair<string, map<unsigned long, unsigned long> >(this_ctg, win_sta_end));
                // lg:: ctg:: window marker 
                (*hbmarker_updated).insert(std::pair<string, map<string, map<unsigned long, unsigned long> > >(this_lg, ctg_win_sta_end) );
            }else
            {
                if( (*hbmarker_updated)[this_lg].find(this_ctg) == (*hbmarker_updated)[this_lg].end() )
                {
                    map<unsigned long, unsigned long> win_sta_end;
                    win_sta_end.insert(std::pair<unsigned long, unsigned long>(last_sta, last_end));      
                    (*hbmarker_updated)[this_lg].insert(std::pair<string, map<unsigned long, unsigned long> >(this_ctg, win_sta_end));                          
                }else
                {
                    (*hbmarker_updated)[this_lg][this_ctg].insert( std::pair<unsigned long, unsigned long>(last_sta, last_end) );
                }
            }
            if(merge_verbose)
            cout << "        : lastt new WIN=" << last_sta
                 << "-"                        << last_end
                 << endl;            
            // next ctg 
            citr ++;
        }
        // next lg
        lgitr ++;
    } 
    //
    bool check_this = true;
    if(check_this)
    {
        cout << "   CHECK: ctg win marker distribution in LGs, after merging into larger windows: " << endl;
        map<string, map<string, map<unsigned long, unsigned long> > >::iterator lgitr;
        map<string, map<string, map<unsigned long, unsigned long> > >::iterator lgitr_end;
        lgitr     = (*hbmarker_updated).begin();
        lgitr_end = (*hbmarker_updated).end();
        while(lgitr != lgitr_end)
        {
            string this_lg = (*lgitr).first;
            map<string, map<unsigned long, unsigned long> > ctg_win_sta_end = (*lgitr).second; 
            //
            map<string, map<unsigned long, unsigned long> >::iterator citr;
            map<string, map<unsigned long, unsigned long> >::iterator citr_end;
            citr     = ctg_win_sta_end.begin();
            citr_end = ctg_win_sta_end.end();
            while(citr != citr_end)
            {
                string this_ctg = (*citr).first;
                map<unsigned long, unsigned long> win_sta_end = (*citr).second;
                map<unsigned long, unsigned long>::iterator witr;
                map<unsigned long, unsigned long>::iterator witr_end;
                witr     = win_sta_end.begin();
                witr_end = win_sta_end.end();
                while(witr != witr_end)
                {
                    cout << "   check-merged: LG=" << this_lg 
                         << ", CTG="               << this_ctg 
                         << ", WIN="               << (*witr).first 
                         << "-"                    << (*witr).second 
                         << endl;
                    // next window marker 
                    witr ++;
                }
                // next ctg 
                citr ++;
            }
            // next lg
            lgitr ++;
        }
    }       
    
    //
    return true;
}                                    
//
bool read_hb_window_marker_grouping(string              hbmarker_file, 
                                    map<string, map<string, map<unsigned long, unsigned long> > >* hbmarker,
                                    map<string, map<string, int> > * ctg2lg)
{
    /*
       function: read gamete_binning grouped window markers       
       return..: hbmarker: <lg_id, <ctg_id, <win_sta, win_end> > >
                 ctg2lg..: <ctg_id, <lg_id, 1> >
    */
    ifstream ifp;
    ifp.open(hbmarker_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << hbmarker_file << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size() == 0 || line[0]=='#') continue;  
        // utg000380l_pilon    510001  560000  dip 1   11  utg001654l_pilon    0.958721    1   0   linked  case_2.2      
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 5)
        {
            cout << "   Warning: skipped insufficient info line of " << line << endl;
            continue;
        }
        string this_lg  = lineinfo[4];        
        string this_ctg = lineinfo[0];
        unsigned long win_sta = strtoul(lineinfo[1].c_str(), NULL, 0);
        unsigned long win_end = strtoul(lineinfo[2].c_str(), NULL, 0);
        // update hbmarker
        if((*hbmarker).find(this_lg) == (*hbmarker).end())
        {
            // window marker 
            map<unsigned long, unsigned long> win_sta_end;
            win_sta_end.insert(std::pair<unsigned long, unsigned long>(win_sta, win_end));
            // ctg:: window marker 
            map<string, map<unsigned long, unsigned long> > ctg_win_sta_end;
            ctg_win_sta_end.insert(std::pair<string, map<unsigned long, unsigned long> >(this_ctg, win_sta_end));
            // lg:: ctg:: window marker 
            (*hbmarker).insert(std::pair<string, map<string, map<unsigned long, unsigned long> > >(this_lg, ctg_win_sta_end) );
        }else
        {
            if( (*hbmarker)[this_lg].find(this_ctg) == (*hbmarker)[this_lg].end() )
            {
                map<unsigned long, unsigned long> win_sta_end;
                win_sta_end.insert(std::pair<unsigned long, unsigned long>(win_sta, win_end));      
                (*hbmarker)[this_lg].insert(std::pair<string, map<unsigned long, unsigned long> >(this_ctg, win_sta_end));                          
            }else
            {
                (*hbmarker)[this_lg][this_ctg].insert( std::pair<unsigned long, unsigned long>(win_sta, win_end) );
            }
        }
        // update ctg2lg
        if( (*ctg2lg).find(this_ctg) == (*ctg2lg).end() )
        {
            // lg list 
            map<string, int> lgtmp;
            lgtmp.insert(std::pair<string, int>(this_lg, 1));
            // ctg:: lg list
            (*ctg2lg).insert(std::pair<string, map<string, int> >(this_ctg, lgtmp));
        }else
        {
            if((*ctg2lg)[this_ctg].find(this_lg) == (*ctg2lg)[this_ctg].end() )
            {
                (*ctg2lg)[this_ctg].insert(std::pair<string, int>(this_lg, 1));
            }
        }
    }
    ifp.close();
    //
    bool check_this = true;
    if(check_this)
    {
        cout << "   CHECK: ctg win marker distribution in LGs: " << endl;
        map<string, map<string, map<unsigned long, unsigned long> > >::iterator lgitr;
        map<string, map<string, map<unsigned long, unsigned long> > >::iterator lgitr_end;
        lgitr     = (*hbmarker).begin();
        lgitr_end = (*hbmarker).end();
        while(lgitr != lgitr_end)
        {
            string this_lg = (*lgitr).first;
            map<string, map<unsigned long, unsigned long> > ctg_win_sta_end = (*lgitr).second; 
            //
            map<string, map<unsigned long, unsigned long> >::iterator citr;
            map<string, map<unsigned long, unsigned long> >::iterator citr_end;
            citr     = ctg_win_sta_end.begin();
            citr_end = ctg_win_sta_end.end();
            while(citr != citr_end)
            {
                string this_ctg = (*citr).first;
                map<unsigned long, unsigned long> win_sta_end = (*citr).second;
                map<unsigned long, unsigned long>::iterator witr;
                map<unsigned long, unsigned long>::iterator witr_end;
                witr     = win_sta_end.begin();
                witr_end = win_sta_end.end();
                while(witr != witr_end)
                {
                    cout << "   check: LG=" << this_lg 
                         << ", CTG="        << this_ctg 
                         << ", WIN="        << (*witr).first 
                         << "-"             << (*witr).second 
                         << endl;
                    // next window marker 
                    witr ++;
                }
                // next ctg 
                citr ++;
            }
            // next lg
            lgitr ++;
        }
        //
        cout << "   CHECK: LG distribution in CTGs: " << endl;        
        map<string, map<string, int> >::iterator citr;
        map<string, map<string, int> >::iterator citr_end;
        citr     = (*ctg2lg).begin();
        citr_end = (*ctg2lg).end();
        while(citr != citr_end)
        {
            cout << "   check: CTG=" << (*citr).first << " in LGs: ";
            map<string, int>::iterator lgitr2;
            map<string, int>::iterator lgitr2_end;
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
    }
    //
    return true;
}



























