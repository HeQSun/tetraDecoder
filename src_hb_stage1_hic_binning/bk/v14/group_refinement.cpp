/* 
   Refine hic-based haplotyping: 
      
      input:  s8_grouping_window_markers.txt
      output: s8_grouping_window_markers_refined.txt
      
   Written by Hequan Sun, MPIPZ (Germany), 2022.........................................................................
   Email: sun@mpipz.mpg.de/sunhequan@gmail.com..........................................................................
*/
#include      <algorithm>
#include            <map>
#include         <string>
#include         <vector>
#include        <fstream>
#include        <sstream>
#include       <iostream>
#include        <iomanip>  // std::setprecision
#include     <sys/stat.h>
#include       <string.h>
#include       <stdlib.h>
#include       <assert.h>
#include         <math.h>
#include       <dirent.h>
#include "split_string.h"
using namespace std;
//
// window marker
struct winMARKER
{
//  string         ctg; // ctd id
//  unsigned long  sta; // window marker start 
    unsigned long  end; // window marker end
    string        type; // hap/dip/trip/tetrap/rep
    vector<string> wlg; // linkage group id 1:4
    vector<string> ctg; // linking contig 
    vector<string> raw; // original line    
    vector<double> hic; // linking hic contact
    string       homlg; // homo-lg id...... 1:12
};
vector<string> all_lgs; // 1:4
//
bool read_hb_window_marker_grouping(string             hbmarker_file, 
         vector<string>*                               ctg_order,
         map<string, map<unsigned long, winMARKER > >* hbmarker,
         map<string, map<string, string> >*            ctg2lg,
         map<string, string>*                          homLGs); 
bool refine_hb_window_marker_intra_ctg(vector<string>  ctg_order,
         map<string, map<unsigned long, winMARKER > >  hbmarker,
         map<string, map<string, string> >             ctg2lg,
         string                                        out_folder);      
bool refine_hb_window_marker_inter_grp(vector<string>  ctg_order,
         map<string, map<unsigned long, winMARKER > >  hbmarker,
         map<string, map<string, string> >             ctg2lg,
         string                                        out_folder);            
//
bool group_refinement(string hbmarker_file, string out_folder)
{
    double startT= clock();
    //
    cout << "   Info: reading hic_binning grouped window marker info.. " << endl;
    vector<string>                               ctg_order;// <ctg_id> keep input order of contigs
    map<string, map<unsigned long, winMARKER > > hbmarker; // <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
    map<string, map<string, string> >            ctg2lg;   // <ctg_id, <lg_id=1:4, 1:12> >
    map<string, string>                          homLGs;   // <lg_id=1:4,  homLG_id=1:12>
    if(!read_hb_window_marker_grouping(hbmarker_file, &ctg_order, &hbmarker, &ctg2lg,& homLGs))
    {
        return false;
    }
    cout << "   Info: reading hic_binning grouped window marker done. "  << endl;      
    cout << endl;
    cout << "   Info: refining hic_binning grouped window marker info.." << endl;
    if(!refine_hb_window_marker_intra_ctg(ctg_order, hbmarker, ctg2lg, out_folder))
    {
        return false;
    }
    cout << "   Info: refining hic_binning grouped window marker done."  << endl;    
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    //
    return true;
}
bool refine_hb_window_marker_inter_grp(vector<string>  ctg_order,
         map<string, map<unsigned long, winMARKER > >  hbmarker,
         map<string, map<string, string> >             ctg2lg,
         string                                        out_folder)
{
    return true;
}
bool refine_hb_window_marker_intra_ctg(vector<string> ctg_order,
         map<string, map<unsigned long, winMARKER > > hbmarker,
         map<string, map<string, string> >            ctg2lg,
         string                                       out_folder)
{
    /*
       function: refine grouping of window markers within each contig 
    */
    if(ctg_order.size()==0 || hbmarker.size()==0 || ctg2lg.size()==0)
    {
        cout << "   Error: incorrect marker grouping data. " << endl;
        return false;
    }
    string out_file = out_folder + "/s8_grouping_window_markers_refined.txt\0";
    ofstream ofp;
    ofp.open(out_file.c_str(), ios::out);
    if(!ofp.good())
    {
        cout   << "   Error: cannot open file " << out_file << endl;
        return false;
    }        
    //
    vector<string>::iterator citr;
    vector<string>::iterator citr_end;
    citr     = ctg_order.begin();
    citr_end = ctg_order.end();
    while(citr != citr_end)
    {
        string this_ctg = (*citr);
        if(hbmarker.find(this_ctg) == hbmarker.end())
        {
            cout << "   Error: marker grouping info not matched. " << endl;
            return false;
        }
        //
        map<string, map<string, double> > type_lg_hic;// <type, <lg, hic> >// how each lg is linked under each type 
        map<string, unsigned long>        type_along_ctg; // <type, total_size>
        map<string, double> link_lg_hic;    // <linked_lg, hic_contact>    // how much hic  the ctg linked to this lg
        map<string, int>    link_lg_cnt;    // <linked_lg, num_win_marker> // how many time the ctg linked to this lg
        //
        map<unsigned long, winMARKER > tmp_win_list = hbmarker[this_ctg];
        map<unsigned long, winMARKER >::iterator mpitr;
        map<unsigned long, winMARKER >::iterator mpitr_end;
        mpitr     = tmp_win_list.begin();
        mpitr_end = tmp_win_list.end();
        while(mpitr != mpitr_end)
        {
            unsigned long win_sta = (*mpitr).first;
            winMARKER tmp_marker  = (*mpitr).second;
            cout << "       : " << this_ctg 
                 << "\t"        << win_sta 
                 << "\t"        << tmp_marker.end 
                 << "\t"        << tmp_marker.type 
                 << "\t"        << tmp_marker.homlg
                 << "\t";
            unsigned long size_win  = tmp_marker.end - win_sta + 1;
            vector<string> this_ctg = tmp_marker.ctg; // linking contig 
            vector<double> this_hic = tmp_marker.hic; // linking hic contact                     
            vector<string> this_LGs = tmp_marker.wlg; // linking group id 
            assert(this_ctg.size()==this_hic.size() && this_ctg.size()==this_LGs.size());                
            vector<string>::iterator lgitr;
            vector<string>::iterator lgitr_end;
            lgitr     = this_LGs.begin();
            lgitr_end = this_LGs.end();
            int i = 0;
            while(lgitr != lgitr_end)
            {
                if(lgitr != this_LGs.begin())
                {
                    cout << ";{" << this_LGs[i] << "," << this_ctg[i] << "," << this_hic[i] << "}";
                }else
                {
                    cout <<  "{" << this_LGs[i] << "," << this_ctg[i] << "," << this_hic[i] << "}";
                }
                // how much hic the ctg linked to this lg and how many time the ctg linked to this lg
                if(link_lg_hic.find( this_LGs[i] ) == link_lg_hic.end() )
                {
                   link_lg_hic.insert(std::pair<string, double>( this_LGs[i], this_hic[i] ) );
                   link_lg_cnt.insert(std::pair<string, int>( this_LGs[i], 1 ) );
                }else
                {
                   link_lg_hic[ this_LGs[i] ] += this_hic[i];
                   link_lg_cnt[ this_LGs[i] ] += 1;                    
                }
                // how each lg is linked with current ctg under each type 
                if( type_lg_hic.find(tmp_marker.type) == type_lg_hic.end() )
                {
                    map<string, double> tmp_lg_hic;
                    tmp_lg_hic.insert(std::pair<string, double>( this_LGs[i], this_hic[i]  ));
                    type_lg_hic.insert(std::pair<string, map<string, double> >(tmp_marker.type, tmp_lg_hic));
                }else
                {
                    if( type_lg_hic[tmp_marker.type].find( this_LGs[i] ) ==  type_lg_hic[tmp_marker.type].end() )
                    {
                        type_lg_hic[tmp_marker.type].insert(std::pair<string, double>( this_LGs[i], this_hic[i]  ));
                    }else
                    {
                        type_lg_hic[tmp_marker.type][ this_LGs[i] ] += this_hic[i];
                    }
                }
                //
                lgitr ++;
                i ++;
            }
            cout << endl;   
            // marker type along contig 
            if(type_along_ctg.find(tmp_marker.type) == type_along_ctg.end())
            {
                type_along_ctg.insert(std::pair<string, unsigned long>(tmp_marker.type, size_win));
            }else
            {
                type_along_ctg[tmp_marker.type] += size_win;
            }
            //
            mpitr ++;
        } 
        // analyze type 
        cout << "       : " << this_ctg << " with " << type_along_ctg.size() << " types of markers: " << endl; 
        // analyze total hic-score/cnt 
        map<string, unsigned long>::iterator titr; // <type, 1>
        map<string, unsigned long>::iterator titr_end;
        titr     = type_along_ctg.begin();
        titr_end = type_along_ctg.end();
        while(titr != titr_end)
        {
            string this_mkr_type = (*titr).first;
            cout << "         | " << (*titr).first << ", " << (*titr).second << " bp: " << endl;
            assert( type_lg_hic.find( this_mkr_type ) != type_lg_hic.end() );
            map<string, double> tmp_lg_hic = type_lg_hic[ this_mkr_type ];
            map<string, double>::iterator tlgitr; // lg under this type 
            map<string, double>::iterator tlgitr_end;
            tlgitr     = tmp_lg_hic.begin();
            tlgitr_end = tmp_lg_hic.end();
            while(tlgitr != tlgitr_end)
            {
                cout << "            - LG-"      << (*tlgitr).first 
                     << " with total hic-score " << (*tlgitr).second 
                     << endl;
                //
                tlgitr ++;
            }
            //
            titr ++;
        }
        // sort the candidate linkage groups according to the hic-score 
        vector<string> final_lgs; // <1st_best_lg,     2nd_best_lg,     ...>
        vector<double> final_hic; // <1st_best_lg_hic, 2nd_best_lg_hic, ...>
        map<string, double>::iterator hitr;
        map<string, double>::iterator hitr_end;
        hitr     = link_lg_hic.begin();
        hitr_end = link_lg_hic.end();
        while(hitr != hitr_end)
        {
            string linked_lg = (*hitr).first;
            assert( link_lg_cnt.find(linked_lg) != link_lg_cnt.end() );
            cout << "         | total hic-link score with LG-"<< linked_lg << ": "<< link_lg_hic[ linked_lg ] << endl;
            cout << "         | total hic-link cnt   with LG-"<< linked_lg << ": "<< link_lg_cnt[ linked_lg ] << endl;
            //
            if(final_lgs.size() == 0)
            {
                final_lgs.push_back(              linked_lg   );
                final_hic.push_back( link_lg_hic[ linked_lg ] );    
            }else
            {
                bool tmp_inserted = false;
                vector<string>::iterator tmp_lg_itr  = final_lgs.begin();
                vector<double>::iterator tmp_hic_itr = final_hic.begin();
                for(int ii = 0; ii < final_lgs.size(); ii ++)
                {
                    if( final_hic[ii] <= link_lg_hic[ linked_lg ] )
                    {
                        final_lgs.insert(tmp_lg_itr, linked_lg);
                        final_hic.insert(tmp_hic_itr, link_lg_hic[ linked_lg ]);
                        tmp_inserted = true;
                        break;
                    }
                    tmp_lg_itr ++;
                    tmp_hic_itr ++;                    
                }
                if(!tmp_inserted)
                {
                    final_lgs.insert(tmp_lg_itr, linked_lg);
                    final_hic.insert(tmp_hic_itr, link_lg_hic[ linked_lg ]);                    
                }
            }
            //
            hitr ++;
        }
        // print to screen
        for(int ii = 0; ii < final_lgs.size(); ii ++)
        {
            cout << "         * final ordering: LG-" << final_lgs[ii] << " with score " << final_hic[ii] << endl;
        }        
        cout << endl;
        // output updated grouping to file: utg000235l_pilon	590001	750000	trip	1	CHR1	utg000148l_pilon	11.3636	assign1
        mpitr     = tmp_win_list.begin();
        mpitr_end = tmp_win_list.end();
        while(mpitr != mpitr_end)
        {
            unsigned long win_sta = (*mpitr).first;
            winMARKER tmp_marker  = (*mpitr).second;
            int out_num = -1;
            if(tmp_marker.type.compare("hap") == 0)
            {
                out_num = 1;
            }else            
            if(tmp_marker.type.compare("dip") == 0)
            {
                out_num = 2;
            }else
            if(tmp_marker.type.compare("trip") == 0)
            {
                out_num = 3;
            }else    
            if(tmp_marker.type.compare("tetrap") == 0)
            {
                out_num = 4;
                final_lgs = all_lgs;
            }else
            if(tmp_marker.type.compare("rep") == 0)
            {
                ofp  <<  this_ctg 
                     << "\t"        << win_sta 
                     << "\t"        << tmp_marker.end 
                     << "\t"        << tmp_marker.type 
                     << "\t"        << "-1"
                     << "\t"        << tmp_marker.homlg
                     << "\traw\t"   << tmp_marker.raw[ 0 ]
                     << endl;  
            }            
            for (int oi = 0; oi < out_num; oi ++)
            {
                if( final_lgs.size() >= out_num)
                {
                    ofp  << this_ctg 
                         << "\t"        << win_sta 
                         << "\t"        << tmp_marker.end 
                         << "\t"        << tmp_marker.type 
                         << "\t"        << final_lgs[ oi ]
                         << "\t"        << tmp_marker.homlg
                         << "\traw\t"   << tmp_marker.raw[ oi ]
                         << endl;   
                }else
                {
                    ofp  << this_ctg 
                         << "\t"        << win_sta 
                         << "\t"        << tmp_marker.end 
                         << "\t"        << tmp_marker.type 
                         << "\t"        << "-1"
                         << "\t"        << tmp_marker.homlg
                         << "\traw\t"   << tmp_marker.raw[ oi ]
                         << endl;                      
                }
            }
            //
            mpitr ++;
        }
        //
        citr ++;
    }
    //
    ofp.close();
    //
    return true;
}
//
bool read_hb_window_marker_grouping(string hbmarker_file, 
         vector<string>*                               ctg_order,
         map<string, map<unsigned long, winMARKER > >* hbmarker,
         map<string, map<string, string> >*            ctg2lg,
         map<string, string>*                          homLGs)
{
    /*
       function: read hic_binning grouped window markers       
       return..: hbmarker.: <ctg_id, <win_sta, { win-end, 
                                                type=hap/dip/...,
                                                wlg=<1/2/..>, 
                                                ctg=<ctg1/ctg2/..>, 
                                                val=<val1/val2/..>, 
                                                homlg=1/../12 
                                               } > >
                 ctg2lg...: <ctg_id, <lg_id=1:4, 1:12> >
                 homLGs...: <lg_id=1:4,  homLG_id=1:12>
                 ctg_order: <ctg_id>, keeping order of contigs in input
       //
	struct winMARKER
	{
	    unsigned long  end; // window marker end
	    string        type; // hap/dip/trip/tetrap/rep
	    vector<string> wlg; // linkage group id 
            vector<string> ctg; // linking contig 
            vector<double> hic; // linking hic contact
            string       homlg; // homo-lg id...... 1:12	    
	};
    */
    unsigned long win_num = 0;
    ifstream ifp;
    ifp.open(hbmarker_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << hbmarker_file << endl;
        return false;
    }
    while(ifp.good() )
    {
        string line("");
        getline(ifp, line);
        if(line.size() == 0 || line[0]=='#') continue;  
        /*
	   utg000041l_pilon          1    10000  hap 4   CHR1    utg000129l_pilon    9.21569 assign1
	   utg000041l_pilon      10001   510000  dip 4   CHR1    utg000129l_pilon    63.3    assign1
	   utg000041l_pilon      10001   510000  dip 1   CHR1    utg001514l_pilon    41.3959 assign2
	   utg000041l_pilon     510001  1120000  dip 4   CHR1    utg000129l_pilon    104.865 assign1
	   utg000041l_pilon     510001  1120000  dip 1   CHR1    utg001514l_pilon    37.7689 assign2
	   utg000041l_pilon    1120001  1139408  hap 1   CHR1    utg000560l_pilon    10.2618 assign1       
        */
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 8)
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
        string link_ctg       = lineinfo[6];
        double link_hic       = atof(lineinfo[7].c_str());
        string raw_line_sub   = "";
        for(int ri = 4; ri < lineinfo.size(); ri ++)
        {
            if(ri > 4)
            raw_line_sub += "\t";
            raw_line_sub += lineinfo[4];
        }
        win_num ++;
        //
        if(std::find(all_lgs.begin(), all_lgs.end(), this_lg) == all_lgs.end() )
        {
            all_lgs.push_back(this_lg);
        }
        //
        if(std::find( (*ctg_order).begin(), (*ctg_order).end(), this_ctg) == (*ctg_order).end() )
        {
            (*ctg_order).push_back(this_ctg);
        }
        //
        if( (*hbmarker).find(this_ctg) == (*hbmarker).end() )
        {
            winMARKER tmp_marker;
            tmp_marker.end   = win_end;
            tmp_marker.type  = this_type;
            tmp_marker.wlg.push_back(this_lg);
            tmp_marker.ctg.push_back(link_ctg);
            tmp_marker.hic.push_back(link_hic);
            tmp_marker.raw.push_back(line);            
            tmp_marker.homlg = this_homLG;
            map<unsigned long, winMARKER > tmp_win_list;
            tmp_win_list.insert(std::pair<unsigned long, winMARKER>(win_sta, tmp_marker));
            (*hbmarker).insert(std::pair<string, map<unsigned long, winMARKER > >(this_ctg, tmp_win_list));
        }else
        {
            if( (*hbmarker)[this_ctg].find( win_sta ) == (*hbmarker)[this_ctg].end() )
            {
                winMARKER tmp_marker;
                tmp_marker.end   = win_end;
                tmp_marker.type  = this_type;
                tmp_marker.wlg.push_back(this_lg);
                tmp_marker.ctg.push_back(link_ctg);
                tmp_marker.hic.push_back(link_hic);                
                tmp_marker.raw.push_back(line);
                tmp_marker.homlg = this_homLG;    
                (*hbmarker)[this_ctg].insert(std::pair<unsigned long, winMARKER>(win_sta, tmp_marker));
            }else
            {
                assert(this_type.compare( (*hbmarker)[this_ctg][win_sta].type ) == 0 );
                (*hbmarker)[this_ctg][win_sta].wlg.push_back(this_lg);
                (*hbmarker)[this_ctg][win_sta].ctg.push_back(link_ctg);
                (*hbmarker)[this_ctg][win_sta].hic.push_back(link_hic);
                (*hbmarker)[this_ctg][win_sta].raw.push_back(line);
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
    //
    cout << "   Info: " << (*ctg_order).size() << " ctgs checked, with " << win_num << " window markers. " << endl;
    //
    bool check_this = false;
    if(check_this)
    {
        cout << "   CHECK: ctg win marker distribution in LGs: " << endl;
        map<string, map<unsigned long, winMARKER > >::iterator mcitr;
        map<string, map<unsigned long, winMARKER > >::iterator mcitr_end;
        mcitr     = (*hbmarker).begin();
        mcitr_end = (*hbmarker).end();
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
                     << "\t"         << win_sta 
                     << "\t"         << tmp_marker.end 
                     << "\t"         << tmp_marker.type 
                     << "\t"         << tmp_marker.homlg
                     << "\t";
                vector<string> this_ctg = tmp_marker.ctg; // linking contig 
                vector<double> this_hic = tmp_marker.hic; // linking hic contact                     
                vector<string> this_LGs = tmp_marker.wlg; // linking group id 
                assert(this_ctg.size()==this_hic.size() && this_ctg.size()==this_LGs.size());                
                vector<string>::iterator lgitr;
                vector<string>::iterator lgitr_end;
                lgitr     = this_LGs.begin();
                lgitr_end = this_LGs.end();
                int i = 0;
                while(lgitr != lgitr_end)
                {
                    if(lgitr != this_LGs.begin())
                    {
                        cout << ";{" << *lgitr << "," << this_ctg[i] << "," << this_hic[i] << "}";
                    }else
                    {
                        cout <<  "{" << *lgitr << "," << this_ctg[i] << "," << this_hic[i] << "}";
                    }
                    lgitr ++;
                    i ++;
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
            cout << "        : CTG=" << (*citr).first << " in LGs: ";
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
        cout << "   CHECK:  homologous LGs distribution 1-4 versus 1-12: " << endl;
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
