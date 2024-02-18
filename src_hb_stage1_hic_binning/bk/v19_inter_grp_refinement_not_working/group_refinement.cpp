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
struct winMARKER2
{
//  string         ctg; // ctd id
//  unsigned long  sta; // window marker start 
    unsigned long  end; // window marker end
    string        type; // hap/dip/trip/tetrap/rep
    string          lg; // currently assigned lg
    vector<string> raw; // original line                    
    map<string, double> lg_hic_sum;    // <1:4, sum_hic_contact>
    map<string, double> lg_hic_sum_cnt;// <1:4, sum_hic_contact_cnt>        
};
struct WHIC 
{
  //unsigned long sta1; // marker-window sta of ctg1
    unsigned long end1; // marker-window end of ctg1
    string        typ1; // type of window of ctg1: hap/dip/...
  //unsigned long sta2; // marker-window sta of ctg2
    unsigned long end2; // marker-window end of ctg2 
    string        typ2; // type of window of ctg2: hap/dip/...
    unsigned long hcnt; // hic-link-cnt between {ctg1:sta1-end1} and {ctg2:sta2-end2}    
};
struct NONALLELIE
{
 // unsigned long sta;    // sta of window marker
    unsigned long end;    // end of window marker
    map<int, double> grp; // <non_allelic_grp_id, Hi-C_score>
};
vector<string> all_lgs; // 1:4
string homLG = "";
//
bool read_hb_window_marker_grouping(string             hbmarker_file, 
         vector<string>*                               ctg_order,
         map<string, map<unsigned long, winMARKER > >* hbmarker,
         map<string, map<string, string> >*            ctg2lg,
         map<string, string>*                          homLGs);
bool refine_hb_window_marker_intra_ctg(vector<string>  ctg_order,
         map<string, map<unsigned long, winMARKER > >  hbmarker,
         map<string, map<string, string> >             ctg2lg,
         map<string, unsigned long>                    target_contig_size,
         string                                        out_label,
         string                                        out_folder);      
bool refine_hb_window_marker_inter_grp(vector<string>  ctg_order,
         map<string, map<unsigned long, winMARKER > >  hbmarker,
         map<string, map<string, string> >             ctg2lg,
         string                                        hic_matrix_file,
         map<string, map<unsigned long, NONALLELIE> >  ctg_win_clear_allelic_group,
         map<string, map<int, map<string, unsigned long> > > ctg_allelic_group,
         map<string, unsigned long>                    allelic_map_seq_size,
         map<string, unsigned long>                    target_contig_size,         
         string                                        out_folder);  
bool add_vet_to_group(string                           vet, 
         int                                           hic_grpid_old,
         int                                           hic_grpid_new,
         map<int, map<string, int> >*                  hic_group_2ctg,
         map<string, int>*                             ctg2_hic_group);
bool further_extend_hic_binning_lg(string              new_hbmarker_file, 
         map<string, unsigned long>                    target_contig_size,
         map<string, double>                           ctg_hap_ratio,
         map<string, int>*                             ctg2_hic_group,
         map<int, map<string, int> >*                  hic_group_2ctg);  
bool define_group_with_allelic(map<string, int>        ctg2_hic_group,
         map<int, map<string, int> >                   hic_group_2ctg,
         map<string, string>                           allelic_map, 
         map<string, string>                           allelic_map_seq,
         map<string, unsigned long>                    target_contig_size,
         map<string, map<int, map<string, unsigned long> > >* ctg_allelic_group);              
//
bool group_refinement(string                           hbmarker_file, 
         string                                        out_folder, 
         map<string, map<unsigned long, map<unsigned long, WHIC> > > hic_matrix,
         map<string, unsigned long>                    target_contig_size,
         map<string, double>                           ctg_hap_ratio,
         map<string, map<unsigned long, NONALLELIE> >  ctg_win_clear_allelic_group,
         map<string, int>                              ctg2_hic_group,
         map<int, map<string, int> >                   hic_group_2ctg,
         map<string, string>                           allelic_map, 
         map<string, string>                           allelic_map_seq,
         map<string, unsigned long>                    allelic_map_seq_size)
{
    double startT= clock();
    //
    vector<string>                               ctg_order;// <ctg_id> keep input order of contigs
    map<string, map<unsigned long, winMARKER > > hbmarker; // <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
    map<string, map<string, string> >            ctg2lg;   // <ctg_id, <lg_id=1:4, 1:12> >
    map<string, string>                          homLGs;   // <lg_id=1:4,  homLG_id=1:12>    
    //
    for(int intra_i = 0; intra_i < 1; intra_i ++)
    {
        if(intra_i > 0)
        {
            hbmarker_file = out_folder + "/s8_grouping_window_markers_refined_1st.txt\0"; 
        }
        // s1: 1st intra-ctg refinement
        cout << "   Info: reading hic_binning grouped window marker info for 1st intra-ctg refinement.. " << endl;    
        ctg_order.clear();// <ctg_id> keep input order of contigs
        hbmarker.clear(); // <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
        ctg2lg.clear();   // <ctg_id, <lg_id=1:4, 1:12> >
        homLGs.clear();   // <lg_id=1:4,  homLG_id=1:12>        
        if(!read_hb_window_marker_grouping(hbmarker_file, &ctg_order, &hbmarker, &ctg2lg,& homLGs))
        {
            // note, hbmarker_file = tmpfolder_s6 + "/s8_grouping_window_markers.txt\0"
            return false;
        }
        cout << "   Info: reading hic_binning grouped window marker info for 1st intra-ctg refinement done. " << endl;      
        cout << endl;
        cout << "   Info: 1st-refining hic_binning grouped window marker info.." << endl;
        if(!refine_hb_window_marker_intra_ctg(ctg_order, hbmarker, ctg2lg, target_contig_size, "1st", out_folder))
        {
            return false;
        }
        cout << "   Info: 1st-refining hic_binning grouped window marker done."           << endl;
    }
    //
    string hic_matrix_file = out_folder + "/zadd_hic_matrix.txt\0";     
    string new_hbmarker_file = out_folder + "/s8_grouping_window_markers_refined_1st.txt\0";   
    cout << "   Info: 1st-refining marker collected in " << new_hbmarker_file << endl;
    cout << "   Info: Hi-C matrix file: "                << hic_matrix_file   << endl;
    cout << endl; 
    cout << "   Info: further extending the backbone groups with non-pure haplotigs... " << endl;
    if(!further_extend_hic_binning_lg(new_hbmarker_file,
                                      target_contig_size,
                                      ctg_hap_ratio,
                                      &ctg2_hic_group,
                                      &hic_group_2ctg) )
    {
        cout << "   Error: failed in extending the haplotig groups. " << endl;
        return false;
    }
    cout << "   Info: further extending the backbone groups with non-pure haplotigs done. " << endl;    
    cout << "   Info: re-building allelic groups for contigs..."                            << endl;
    /* Rebuild allelic groups
       <ctg, <allelic_grp, <allelic_ctg, allelic_ctg_size > > >
       ctg has a lower chance to be assigned to allelic_grp, where the prob can be defined by number of allelic ctgs.
    */
    map<string, map<int, map<string, unsigned long> > > ctg_allelic_group; 
    if(!define_group_with_allelic(ctg2_hic_group, 
                                  hic_group_2ctg, 
                                  allelic_map, 
                                  allelic_map_seq, 
                                  target_contig_size,
                                  &ctg_allelic_group))
    {
        // using the info below: 
        // map<string, string> allelic_map;       // allelic relationship between contigs <ctg1-ctg2, "0"> - gmap/gene based infer!
        // map<string, string> allelic_map_seq;   // allelic relationship between contigs <ctg1-ctg2, "0"> - ref-align based infer!
        // ctg2_hic_group.clear(); // <ctg_id, clusterid_old>
        // hic_group_2ctg.clear(); // <clusterid_old, <ctg_id, clusterid_new> >  
        cout << "   Error: failed in defining grouping based on allelic info. " << endl;
        return false; 
    }     
    cout << "   Info: re-building allelic groups for contigs done."             << endl;                           
    // s2: inter-grp refinement  
    cout << "   Info: reading hic_binning grouped window marker info for 1st inter-grp refinement.. " << endl;
    ctg_order.clear();// <ctg_id> keep input order of contigs
    hbmarker.clear(); // <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
    ctg2lg.clear();   // <ctg_id, <lg_id=1:4, 1:12> >
    homLGs.clear();   // <lg_id=1:4,  homLG_id=1:12>
    if(!read_hb_window_marker_grouping(new_hbmarker_file, &ctg_order, &hbmarker, &ctg2lg, &homLGs))
    {
        return false;
    }
    cout << "   Info: reading hic_binning grouped window marker info for 1st inter-grp refinement done. " << endl;
    cout << endl;
    if(!refine_hb_window_marker_inter_grp(ctg_order, 
                                          hbmarker, 
                                          ctg2lg, 
                                          hic_matrix_file, 
                                          ctg_win_clear_allelic_group,
                                          ctg_allelic_group,
                                          allelic_map_seq_size,
                                          target_contig_size,                                          
                                          out_folder) )
    {
        return false;
    }            
    // s3: 2nd intra-ctg refinement
    new_hbmarker_file = out_folder + "/s8_grouping_window_markers_tmp2_overall_hic.txt\0";
    cout << "   Info: reading hic_binning grouped window marker info.. " << endl;
    ctg_order.clear();// <ctg_id> keep input order of contigs
    hbmarker.clear(); // <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
    ctg2lg.clear();   // <ctg_id, <lg_id=1:4, 1:12> >
    homLGs.clear();   // <lg_id=1:4,  homLG_id=1:12>
    if(!read_hb_window_marker_grouping(new_hbmarker_file, &ctg_order, &hbmarker, &ctg2lg, &homLGs))
    {
        return false;
    }
    cout << "   Info: reading hic_binning grouped window marker done. "  << endl;      
    cout << endl;  
    cout << "   Info: refining hic_binning grouped window marker info.." << endl;
    if(!refine_hb_window_marker_intra_ctg(ctg_order, hbmarker, ctg2lg, target_contig_size, "2nd", out_folder))
    {
        return false;
    }
    cout << "   Info: refining hic_binning grouped window marker done."  << endl;    
    
    
    
    
    //////////////////////////////////////
    for(int ri = 0; ri < 0; ri ++)
    {
        new_hbmarker_file = out_folder + "/s8_grouping_window_markers_refined_2nd.txt\0";   
        cout << endl;
        // inter-grp refinement  
        cout << "   Info: reading hic_binning grouped marker for " << ri + 3 << "th inter-grp refinement.. " << endl;
        ctg_order.clear();// <ctg_id> keep input order of contigs
        hbmarker.clear(); // <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
        ctg2lg.clear();   // <ctg_id, <lg_id=1:4, 1:12> >
        homLGs.clear();   // <lg_id=1:4,  homLG_id=1:12>
        if(!read_hb_window_marker_grouping(new_hbmarker_file, &ctg_order, &hbmarker, &ctg2lg, &homLGs))
        {
            return false;
        }
        cout << "   Info: reading hic_binning grouped marker for " << ri + 3 << "th inter-grp refinement done. " << endl;
        cout << endl;
        if(!refine_hb_window_marker_inter_grp(ctg_order, 
                                              hbmarker, 
                                              ctg2lg, 
                                              hic_matrix_file, 
                                              ctg_win_clear_allelic_group,
                                              ctg_allelic_group,
                                              allelic_map_seq_size,
                                              target_contig_size,                                              
                                              out_folder))
        {
            return false;
        }          
        // intra-ctg refinement
        new_hbmarker_file = out_folder + "/s8_grouping_window_markers_tmp2_overall_hic.txt\0";
        cout << "   Info: reading hic_binning grouped window marker info.. " << endl;
        ctg_order.clear();// <ctg_id> keep input order of contigs
        hbmarker.clear(); // <ctg_id, <win_sta, {win-end, type=hap/dip/..., wlg=<lgs>, 1:12 } > >
        ctg2lg.clear();   // <ctg_id, <lg_id=1:4, 1:12> >
        homLGs.clear();   // <lg_id=1:4,  homLG_id=1:12>
        if(!read_hb_window_marker_grouping(new_hbmarker_file, &ctg_order, &hbmarker, &ctg2lg, &homLGs))
        {
            return false;
        }
        cout << "   Info: reading hic_binning grouped window marker done. "  << endl;      
        cout << endl;  
        cout << "   Info: refining hic_binning grouped window marker info.." << endl;
        if(!refine_hb_window_marker_intra_ctg(ctg_order, hbmarker, ctg2lg, target_contig_size, "2nd", out_folder))
        {
            return false;
        }
        cout << "   Info: refining hic_binning grouped window marker done."  << endl;        
        //////////////////////////////////////
    }   
       
       
        
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    //
    return true;
}
bool refine_hb_window_marker_inter_grp(vector<string> ctg_order,
         map<string, map<unsigned long, winMARKER > > hbmarker,
         map<string, map<string, string> >            ctg2lg,
         string                                       hic_matrix_file,
         map<string, map<unsigned long, NONALLELIE> > ctg_win_clear_allelic_group,
         map<string, map<int, map<string, unsigned long> > > ctg_allelic_group,
         map<string, unsigned long>                   allelic_map_seq_size,
         map<string, unsigned long>                   target_contig_size,
         string                                       out_folder)
{
    /* 
       function: refine grouping of window markers among groups 
                 ctg_allelic_group................<ctg, <allelic_grp_id, <allelic_ctg, allelic_ctg_size > > >
       return..: s8_grouping_window_markers_tmp[2]_overall_hic.txt
       //
       struct winMARKER
       {
	   unsigned long  end; // window marker end
	   string        type; // hap/dip/trip/tetrap/rep
	   vector<string> wlg; // linkage group id 
           vector<string> ctg; // linking contig 
           vector<string> raw; // original line                
           vector<double> hic; // linking hic contact
           string       homlg; // homo-lg id...... 1:12	    
       };
       struct winMARKER2
       {
	   unsigned long  end; // window marker end
	   string        type; // hap/dip/trip/tetrap/rep
           string          lg; // currently assigned lg
           vector<string> raw; // original line       
	   map<string, double> lg_hic_sum;// <1:4, sum_hic_contact>
           map<string, double> lg_hic_sum_cnt;// <1:4, sum_hic_contact_cnt>        
       };
    */    
    if(ctg_order.size()==0 || hbmarker.size()==0 || ctg2lg.size()==0)
    {
        cout << "   Error: incorrect marker grouping data. " << endl;
        return false;
    } 
    // total amount of hic to each linkage group
    map<string, map<unsigned long, winMARKER2 > > hbmarker_grp_hic;
    //
    ifstream ifp;
    ifp.open(hic_matrix_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open hic matrix file " << hic_matrix_file << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size() == 0 || line[0]=='#') continue;
        /*
           e.g., utg000007l_pilon    1   500000  hap utg000114l_pilon    4500001 5000000 hap 0.1
        */
        vector<string> lineinfo = split_string(line, '\t');
        string            ctg1 = lineinfo[0];
        unsigned long win_sta1 = strtoul(lineinfo[1].c_str(), NULL, 0);
        unsigned long win_end1 = strtoul(lineinfo[2].c_str(), NULL, 0);        
        string           type1 = lineinfo[3];
        string            ctg2 = lineinfo[4];
        unsigned long win_sta2 = strtoul(lineinfo[5].c_str(), NULL, 0);
        unsigned long win_end2 = strtoul(lineinfo[6].c_str(), NULL, 0);        
        string           type2 = lineinfo[7];        
        double     hic_contact = atof(lineinfo[8].c_str());
        if(type1.compare("hap")==0 && type2.compare("hap")==0) // testing! 20221206
        //if(1)
        {
            //assert(ctg2lg.find(ctg1) != ctg2lg.end());
            //assert(ctg2lg.find(ctg2) != ctg2lg.end());            
            //string ctg1_lg = (*(ctg2lg[ctg1].begin())).first;
            //string ctg2_lg = (*(ctg2lg[ctg2].begin())).first;
            // current haplotype group
            assert( hbmarker.find(ctg1) != hbmarker.end() );
            assert( hbmarker[ctg1].find(win_sta1) != hbmarker[ctg1].end() );            
            assert( hbmarker.find(ctg2) != hbmarker.end() );
            assert( hbmarker[ctg2].find(win_sta2) != hbmarker[ctg2].end() );            
            string ctg1_lg = hbmarker[ctg1][win_sta1].wlg[0];        
            string ctg2_lg = hbmarker[ctg2][win_sta2].wlg[0];
            // catch ctg1
            if(hbmarker_grp_hic.find(ctg1) == hbmarker_grp_hic.end())
            {
                winMARKER2 tmp_marker2;
                tmp_marker2.end  = win_end1;
                tmp_marker2.type = type1;
                tmp_marker2.lg   = ctg1_lg;
                tmp_marker2.raw.push_back(line);
                tmp_marker2.lg_hic_sum.insert(std::pair<string, double>(ctg2_lg, hic_contact));
                tmp_marker2.lg_hic_sum_cnt.insert(std::pair<string, double>(ctg2_lg, 1.0));                        
                map<unsigned long, winMARKER2 > tmp_win2_list;
                tmp_win2_list.insert(std::pair<unsigned long, winMARKER2>(win_sta1, tmp_marker2));
                hbmarker_grp_hic.insert(std::pair<string, map<unsigned long, winMARKER2> >(ctg1, tmp_win2_list));
            }else
            {
                if(hbmarker_grp_hic[ctg1].find(win_sta1) == hbmarker_grp_hic[ctg1].end() )
                {
                    winMARKER2 tmp_marker2;
                    tmp_marker2.end  = win_end1;
                    tmp_marker2.type = type1;
                    tmp_marker2.lg   = ctg1_lg;   
                    tmp_marker2.raw.push_back(line);                                     
                    tmp_marker2.lg_hic_sum.insert(std::pair<string, double>(ctg2_lg, hic_contact));                    
                    tmp_marker2.lg_hic_sum_cnt.insert(std::pair<string, double>(ctg2_lg, 1.0));                                            
                    hbmarker_grp_hic[ctg1].insert(std::pair<unsigned long, winMARKER2>(win_sta1, tmp_marker2));
                }else
                {
                    if(hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum.find(ctg2_lg) == hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum.end() )
                    {
                        hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum.insert(std::pair<string, double>(ctg2_lg, hic_contact));
                        hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum_cnt.insert(std::pair<string, double>(ctg2_lg, 1.0));                                                
                        hbmarker_grp_hic[ctg1][win_sta1].raw.push_back(line);
                    }else
                    {
                        hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum[ctg2_lg] += hic_contact;
                        hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum_cnt[ctg2_lg] += 1.0;
                        hbmarker_grp_hic[ctg1][win_sta1].raw.push_back(line);                        
                    }
                }
            }
            // catch ctg2
            if(hbmarker_grp_hic.find(ctg2) == hbmarker_grp_hic.end())
            {
                winMARKER2 tmp_marker2;
                tmp_marker2.end  = win_end2;
                tmp_marker2.type = type2;
                tmp_marker2.lg   = ctg2_lg;  
                tmp_marker2.raw.push_back(line);                              
                tmp_marker2.lg_hic_sum.insert(std::pair<string, double>(ctg1_lg, hic_contact));
                tmp_marker2.lg_hic_sum_cnt.insert(std::pair<string, double>(ctg1_lg, 1.0));                
                map<unsigned long, winMARKER2 > tmp_win2_list;
                tmp_win2_list.insert(std::pair<unsigned long, winMARKER2>(win_sta2, tmp_marker2));
                hbmarker_grp_hic.insert(std::pair<string, map<unsigned long, winMARKER2> >(ctg2, tmp_win2_list));
            }else
            {
                if(hbmarker_grp_hic[ctg2].find(win_sta2) == hbmarker_grp_hic[ctg2].end() )
                {
                    winMARKER2 tmp_marker2;
                    tmp_marker2.end  = win_end2;
                    tmp_marker2.type = type2;
                    tmp_marker2.lg   = ctg2_lg;  
                    tmp_marker2.raw.push_back(line);                                      
                    tmp_marker2.lg_hic_sum.insert(std::pair<string, double>(ctg1_lg, hic_contact));                    
                    tmp_marker2.lg_hic_sum_cnt.insert(std::pair<string, double>(ctg1_lg, 1.0));                    
                    hbmarker_grp_hic[ctg2].insert(std::pair<unsigned long, winMARKER2>(win_sta2, tmp_marker2));
                }else
                {
                    if(hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum.find(ctg1_lg) == hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum.end() )
                    {
                        hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum.insert(std::pair<string, double>(ctg1_lg, hic_contact));
                        hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum_cnt.insert(std::pair<string, double>(ctg1_lg, 1.0));                        
                        hbmarker_grp_hic[ctg2][win_sta2].raw.push_back(line);                        
                    }else
                    {
                        hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum[ctg1_lg] += hic_contact;
                        hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum_cnt[ctg1_lg] += 1.0;                        
                        hbmarker_grp_hic[ctg2][win_sta2].raw.push_back(line);                        
                    }
                }
            }            
        }else
        if(type1.compare("hap")==0)
        {
            //assert(ctg2lg.find(ctg1) != ctg2lg.end());
            //assert(ctg2lg.find(ctg2) != ctg2lg.end());            
            //string ctg1_lg = (*(ctg2lg[ctg1].begin())).first;
            //string ctg2_lg = (*(ctg2lg[ctg2].begin())).first;
            // current haplotype group
            assert( hbmarker.find(ctg1) != hbmarker.end() );
            assert( hbmarker[ctg1].find(win_sta1) != hbmarker[ctg1].end() );            
            assert( hbmarker.find(ctg2) != hbmarker.end() );
            assert( hbmarker[ctg2].find(win_sta2) != hbmarker[ctg2].end() );            
            string ctg1_lg = hbmarker[ctg1][win_sta1].wlg[0];        
            string ctg2_lg = hbmarker[ctg2][win_sta2].wlg[0]; // .wlg can be with many lgs            
            // catch ctg2
            if(hbmarker_grp_hic.find(ctg2) == hbmarker_grp_hic.end())
            {
                winMARKER2 tmp_marker2;
                tmp_marker2.end  = win_end2;
                tmp_marker2.type = type2;
                tmp_marker2.lg   = ctg2_lg;   
                tmp_marker2.raw.push_back(line);                             
                tmp_marker2.lg_hic_sum.insert(std::pair<string, double>(ctg1_lg, hic_contact));
                tmp_marker2.lg_hic_sum_cnt.insert(std::pair<string, double>(ctg1_lg, 1.0));                
                map<unsigned long, winMARKER2 > tmp_win2_list;
                tmp_win2_list.insert(std::pair<unsigned long, winMARKER2>(win_sta2, tmp_marker2));
                hbmarker_grp_hic.insert(std::pair<string, map<unsigned long, winMARKER2> >(ctg2, tmp_win2_list));
            }else
            {
                if(hbmarker_grp_hic[ctg2].find(win_sta2) == hbmarker_grp_hic[ctg2].end() )
                {
                    winMARKER2 tmp_marker2;
                    tmp_marker2.end  = win_end2;
                    tmp_marker2.type = type2;
                    tmp_marker2.lg   = ctg2_lg;     
                    tmp_marker2.raw.push_back(line);                                   
                    tmp_marker2.lg_hic_sum.insert(std::pair<string, double>(ctg1_lg, hic_contact));                    
                    tmp_marker2.lg_hic_sum_cnt.insert(std::pair<string, double>(ctg1_lg, 1.0));                    
                    hbmarker_grp_hic[ctg2].insert(std::pair<unsigned long, winMARKER2>(win_sta2, tmp_marker2));
                }else
                {
                    if(hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum.find(ctg1_lg) == hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum.end() )
                    {
                        hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum.insert(std::pair<string, double>(ctg1_lg, hic_contact));
                        hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum_cnt.insert(std::pair<string, double>(ctg1_lg, 1.0));                        
                        hbmarker_grp_hic[ctg2][win_sta2].raw.push_back(line);                        
                    }else
                    {
                        hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum[ctg1_lg] += hic_contact;
                        hbmarker_grp_hic[ctg2][win_sta2].lg_hic_sum_cnt[ctg1_lg] += 1.0;                        
                        hbmarker_grp_hic[ctg2][win_sta2].raw.push_back(line);                        
                    }
                }
            }    
        }else
        if(type2.compare("hap")==0)
        {
            //assert(ctg2lg.find(ctg1) != ctg2lg.end());
            //assert(ctg2lg.find(ctg2) != ctg2lg.end());            
            //string ctg1_lg = (*(ctg2lg[ctg1].begin())).first;
            //string ctg2_lg = (*(ctg2lg[ctg2].begin())).first;
            // current haplotype group
            assert( hbmarker.find(ctg1) != hbmarker.end() );
            assert( hbmarker[ctg1].find(win_sta1) != hbmarker[ctg1].end() );            
            assert( hbmarker.find(ctg2) != hbmarker.end() );
            assert( hbmarker[ctg2].find(win_sta2) != hbmarker[ctg2].end() );            
            string ctg1_lg = hbmarker[ctg1][win_sta1].wlg[0];        
            string ctg2_lg = hbmarker[ctg2][win_sta2].wlg[0];               
            // catch ctg1
            if(hbmarker_grp_hic.find(ctg1) == hbmarker_grp_hic.end())
            {
                winMARKER2 tmp_marker2;
                tmp_marker2.end  = win_end1;
                tmp_marker2.type = type1;
                tmp_marker2.lg   = ctg1_lg;                     
                tmp_marker2.raw.push_back(line);                                                   
                tmp_marker2.lg_hic_sum.insert(std::pair<string, double>(ctg2_lg, hic_contact));
                tmp_marker2.lg_hic_sum_cnt.insert(std::pair<string, double>(ctg2_lg, 1.0));                
                map<unsigned long, winMARKER2 > tmp_win2_list;
                tmp_win2_list.insert(std::pair<unsigned long, winMARKER2>(win_sta1, tmp_marker2));
                hbmarker_grp_hic.insert(std::pair<string, map<unsigned long, winMARKER2> >(ctg1, tmp_win2_list));
            }else
            {
                if(hbmarker_grp_hic[ctg1].find(win_sta1) == hbmarker_grp_hic[ctg1].end() )
                {
                    winMARKER2 tmp_marker2;
                    tmp_marker2.end  = win_end1;
                    tmp_marker2.type = type1;
                    tmp_marker2.lg   = ctg1_lg;
                    tmp_marker2.raw.push_back(line);                                                       
                    tmp_marker2.lg_hic_sum.insert(std::pair<string, double>(ctg2_lg, hic_contact));
                    tmp_marker2.lg_hic_sum_cnt.insert(std::pair<string, double>(ctg2_lg, 1.0));                    
                    hbmarker_grp_hic[ctg1].insert(std::pair<unsigned long, winMARKER2>(win_sta1, tmp_marker2));
                }else
                {
                    if(hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum.find(ctg2_lg) == hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum.end() )
                    {
                        hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum.insert(std::pair<string, double>(ctg2_lg, hic_contact));
                        hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum_cnt.insert(std::pair<string, double>(ctg2_lg, 1.0));                        
                        hbmarker_grp_hic[ctg1][win_sta1].raw.push_back(line);                        
                    }else
                    {
                        hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum[ctg2_lg] += hic_contact;
                        hbmarker_grp_hic[ctg1][win_sta1].lg_hic_sum_cnt[ctg2_lg] += 1.0;                        
                        hbmarker_grp_hic[ctg1][win_sta1].raw.push_back(line);                        
                    }
                }
            }
        }else ;
    }
    ifp.close();
    // output 1
    string out_file = out_folder + "/s8_grouping_window_markers_tmp1_overall_hic.txt\0";
    ofstream ofp;
    ofp.open(out_file.c_str(), ios::out);
    if(!ofp.good())
    {
        cout   << "   Error: cannot open file " << out_file << endl;
        return false;
    }
    map<string, map<unsigned long, winMARKER2 > >::iterator sum_itr;
    map<string, map<unsigned long, winMARKER2 > >::iterator sum_itr_end;
    sum_itr     = hbmarker_grp_hic.begin();
    sum_itr_end = hbmarker_grp_hic.end();
    while(sum_itr != sum_itr_end)
    {
        string this_ctg = (*sum_itr).first;
        // define allelic groups for this contig
        cout << "   checkig: phasing windows of " << this_ctg << endl;
        /* check allelic contigs in groups .... */
        map<int, map<string, unsigned long> > allelic_grps;
        if(ctg_allelic_group.find(this_ctg) != ctg_allelic_group.end() )
        {
            allelic_grps = ctg_allelic_group[this_ctg];
        }        
        // sort the linkage groups according to the number of allelic groups: grps with more ctgs with less priority
        vector<int>    allelic_final_lgs;      // <1st_allelic_lg,     2nd_allelic_lg,     ...>
        vector<double> allelic_final_ctg_cnt;  // <1st_allelic_lg_hic, 2nd_allelic_lg_hic, ...>
        vector<unsigned long> allelic_final_ctg_size; // <1st_allelic_lg_size,2nd_allelic_lg_size, ...>
        vector<unsigned long> allelic_final_ctg_size_overlap; // <1st_allelic_lg_overlap_size,2nd_allelic_lg_overlap_size, ...>
        vector<double> allelic_final_ctg_size_overlap_ratio; // <1st_allelic_lg_overlap_ratio,2nd_allelic_lg_overlap_rato, ...>        
        map<int, map<string, unsigned long> >::iterator alitr; // <allelic_grp_id, <allelic_ctg, allelic_ctg_size > >
        map<int, map<string, unsigned long> >::iterator alitr_end;
        alitr     = allelic_grps.begin();
        alitr_end = allelic_grps.end();
        while(alitr != alitr_end)
        {
            int linked_lg = (*alitr).first;
            int allelic_ctg_num = (*alitr).second.size();
            unsigned long allelic_ctg_size = 0;
            unsigned long allelic_ctg_size_overlap = 0;            
            double  allelic_ctg_size_overlap_ratio = 0;
            map<string, unsigned long>::iterator sitr;
            map<string, unsigned long>::iterator sitr_end;
            sitr     = (*alitr).second.begin();
            sitr_end = (*alitr).second.end();
            while(sitr != sitr_end)
            {
                allelic_ctg_size += (*sitr).second;
                string allelic_ctg_id = (*sitr).first;
                string allelic_key1 = this_ctg       + "-" + allelic_ctg_id;
                string allelic_key2 = allelic_ctg_id + "-" + this_ctg;
                if( allelic_map_seq_size.find( allelic_key1 ) != allelic_map_seq_size.end() )
                {
                    allelic_ctg_size_overlap += allelic_map_seq_size[ allelic_key1 ];
                }else
                if( allelic_map_seq_size.find( allelic_key2 ) != allelic_map_seq_size.end() )
                {
                    allelic_ctg_size_overlap += allelic_map_seq_size[ allelic_key2 ];
                }else
                {
                    allelic_ctg_size_overlap += 0;
                }
                //
                sitr ++;
            }
            // allelic_ctg_size_overlap_ratio = allelic_ctg_size_overlap * 1.0 / (allelic_ctg_size + 1); // target_contig_size[this_ctg]
            // allelic_ctg_size_overlap_ratio = allelic_ctg_size_overlap * 2.0 / (allelic_ctg_size + target_contig_size[this_ctg]);
            // how much of this_ctg itself is allelic to any given group!
            allelic_ctg_size_overlap_ratio = allelic_ctg_size_overlap * 1.0 / (0 + target_contig_size[this_ctg]);            
            cout << "   | num of allelic ctgs in LG-"<< linked_lg 
                 << ", cnt: "                        << allelic_ctg_num 
                 << ", size: "                       << allelic_ctg_size 
                 << ", overlap size: "               << allelic_ctg_size_overlap 
                 << ", overlap ratio: "              << allelic_ctg_size_overlap_ratio << endl;            
            if(allelic_ctg_size_overlap_ratio < 0.1) // TODO: set up options on this threshold
            {
                // if overlap too short, remove allelic relationship
                allelic_ctg_num = 0;
                allelic_ctg_size_overlap = 0;
            }
            //
            if(allelic_final_lgs.size() == 0 && allelic_ctg_num > 0)
            {
                allelic_final_lgs.push_back(           linked_lg );
                allelic_final_ctg_cnt.push_back( allelic_ctg_num );
                allelic_final_ctg_size.push_back( allelic_ctg_size );
                allelic_final_ctg_size_overlap.push_back( allelic_ctg_size_overlap );
                allelic_final_ctg_size_overlap_ratio.push_back(allelic_ctg_size_overlap_ratio);
            }else
            if(allelic_ctg_num > 0)            
            {
                bool tmp_inserted = false;
                vector<int>::iterator                tmp_lg_itr = allelic_final_lgs.begin();
                vector<double>::iterator            tmp_hic_itr = allelic_final_ctg_cnt.begin();
                vector<unsigned long>::iterator    tmp_size_itr = allelic_final_ctg_size.begin();
                vector<unsigned long>::iterator tmp_size_ol_itr = allelic_final_ctg_size_overlap.begin();                
                vector<double>::iterator      tmp_size_ol_r_itr = allelic_final_ctg_size_overlap_ratio.begin();                                
                for(int ii = 0; ii < allelic_final_lgs.size(); ii ++)
                {
                    if( allelic_final_ctg_cnt[ii] == allelic_ctg_num )
                    {
                        if(allelic_final_ctg_size[ii] <= allelic_ctg_size)
                        {
                            allelic_final_lgs.insert(tmp_lg_itr, linked_lg);
                            allelic_final_ctg_cnt.insert(tmp_hic_itr, allelic_ctg_num);
                            allelic_final_ctg_size.insert(tmp_size_itr, allelic_ctg_size);
                            allelic_final_ctg_size_overlap.insert(tmp_size_ol_itr, allelic_ctg_size_overlap);
                            allelic_final_ctg_size_overlap_ratio.insert(tmp_size_ol_r_itr, allelic_ctg_size_overlap_ratio);                            
                            tmp_inserted = true;
                            break;
                        }
                    }else                
                    if( allelic_final_ctg_cnt[ii] <  allelic_ctg_num )
                    {
                        allelic_final_lgs.insert(tmp_lg_itr, linked_lg);
                        allelic_final_ctg_cnt.insert(tmp_hic_itr, allelic_ctg_num);
                        allelic_final_ctg_size.insert(tmp_size_itr, allelic_ctg_size);
                        allelic_final_ctg_size_overlap.insert(tmp_size_ol_itr, allelic_ctg_size_overlap);                        
                        allelic_final_ctg_size_overlap_ratio.insert(tmp_size_ol_r_itr, allelic_ctg_size_overlap_ratio);                                                    
                        tmp_inserted = true;
                        break;
                    }
                    tmp_lg_itr ++;
                    tmp_hic_itr ++;
                    tmp_size_itr ++;
                    tmp_size_ol_itr ++;
                    tmp_size_ol_r_itr ++;
                }
                if(!tmp_inserted)
                {
                    allelic_final_lgs.insert(tmp_lg_itr, linked_lg);
                    allelic_final_ctg_cnt.insert(tmp_hic_itr, allelic_ctg_num);
                    allelic_final_ctg_size.insert(tmp_size_itr, allelic_ctg_size);   
                    allelic_final_ctg_size_overlap.insert(tmp_size_ol_itr, allelic_ctg_size_overlap);                                     
                    allelic_final_ctg_size_overlap_ratio.insert(tmp_size_ol_r_itr, allelic_ctg_size_overlap_ratio);
                }
            }else
            {
                ; // no action if allelic_ctg_num is 0
            }
            //
            alitr ++;
        }        
        // print to screen: larger allelic_final_ctg_cnt[ii] => lower prob to join this_ctg into the group
        for(int ii = 0; ii < allelic_final_lgs.size(); ii ++)
        {
            cout << "   * allelic: LG-"  << allelic_final_lgs[ii] 
                 << " with #ctgs = "     << allelic_final_ctg_cnt[ii]
                 << " with bp = "        << allelic_final_ctg_size[ii]
                 << " with ovelap bp = " << allelic_final_ctg_size_overlap[ii]
                 << " with overlap r = " << allelic_final_ctg_size_overlap_ratio[ii]
                 << endl;
        } 
        if(allelic_final_lgs.size() > 0)
        {       
            cout << endl;
        }else
        {
            cout << "   * no allelic LG found. " << endl;
        }
        //
        map<unsigned long, winMARKER2> tmp_win2_list = (*sum_itr).second;
        map<unsigned long, winMARKER2>::iterator witr;
        map<unsigned long, winMARKER2>::iterator witr_end;
        witr     = tmp_win2_list.begin();
        witr_end = tmp_win2_list.end();
        while(witr != witr_end)
        {
            unsigned long this_win_sta = (*witr).first;
            winMARKER2 tmp_marker2 = (*witr).second;
            map<string, double> lg_hic_sum     = tmp_marker2.lg_hic_sum;
            map<string, double> lg_hic_sum_cnt = tmp_marker2.lg_hic_sum_cnt;            
            cout << "   inter-grp: "
                 << this_ctg        << "\t"
                 << this_win_sta    << "\t"
                 << tmp_marker2.end << "\t"
                 << tmp_marker2.type<< ": "
                 << endl;
            // clear allelic groups 
            map<int, double> skipping_allelic_grps; // <grp_id, hic_score>
            // 1st allelic cases
            if(ctg_win_clear_allelic_group.find(this_ctg) != ctg_win_clear_allelic_group.end() )
            {
                if(ctg_win_clear_allelic_group[this_ctg].find( this_win_sta) != ctg_win_clear_allelic_group[this_ctg].end())
                {
                    skipping_allelic_grps = ctg_win_clear_allelic_group[this_ctg][this_win_sta].grp;
                }
            }
            // 2nd allelic cases
            for(int ii = 0; ii < allelic_final_lgs.size(); ii ++)
            {
                if(allelic_final_ctg_size_overlap[ii] > 50000 && allelic_final_ctg_size_overlap_ratio[ii] >= 0.1)
                {
                    if(skipping_allelic_grps.find( allelic_final_lgs[ii]  )==skipping_allelic_grps.end() )
                    {
                        skipping_allelic_grps.insert(std::pair<int, double>( allelic_final_lgs[ii], 0) );
                    }
                }
            }            
            // sort the candidate linkage groups according to the hic-score 
            vector<string> final_lgs; // <1st_best_lg,     2nd_best_lg,     ...>
            vector<double> final_hic; // <1st_best_lg_hic, 2nd_best_lg_hic, ...>
            map<string, double>::iterator hitr_cnt;
            map<string, double>::iterator hitr;
            map<string, double>::iterator hitr_end;
            hitr_cnt = lg_hic_sum_cnt.begin();
            hitr     = lg_hic_sum.begin();
            hitr_end = lg_hic_sum.end();
            while(hitr != hitr_end)
            {
                string linked_lg = (*hitr).first;
                assert( lg_hic_sum_cnt.find(linked_lg) != lg_hic_sum_cnt.end() );
                if(lg_hic_sum_cnt[ linked_lg ] == 0) lg_hic_sum_cnt[ linked_lg ] = 1;
                double avg_lg_hic = lg_hic_sum[ linked_lg ] * 1.0 / lg_hic_sum_cnt[ linked_lg ];
                cout << "      | total hic-link score with LG-" << linked_lg 
                     <<                                    ": " << lg_hic_sum[ linked_lg ] 
                     <<                               ", avg: " << avg_lg_hic
                     << endl;
                //
                if(skipping_allelic_grps.find( atoi(linked_lg.c_str()) ) != skipping_allelic_grps.end() )
                {
                    cout << "         | skipping clear allelic grp LG-" << linked_lg << endl;
                }else
                if(linked_lg.compare("-1") != 0)
                {
                    if(final_lgs.size() == 0)
                    {
                        final_lgs.push_back( linked_lg  );
                        final_hic.push_back( avg_lg_hic );
                    }else
                    {
                        bool tmp_inserted = false;
                        vector<string>::iterator tmp_lg_itr  = final_lgs.begin();
                        vector<double>::iterator tmp_hic_itr = final_hic.begin();
                        for(int ii = 0; ii < final_lgs.size(); ii ++)
                        {
                            if( final_hic[ii] <= avg_lg_hic )
                            {
                                final_lgs.insert(tmp_lg_itr, linked_lg);
                                final_hic.insert(tmp_hic_itr, avg_lg_hic);
                                tmp_inserted = true;
                                break;
                            }
                            tmp_lg_itr ++;
                            tmp_hic_itr ++;                    
                        }
                        if(!tmp_inserted)
                        {
                            final_lgs.insert(tmp_lg_itr, linked_lg);
                            final_hic.insert(tmp_hic_itr, avg_lg_hic);                    
                        }
                    }
                }else ;
                //
                hitr ++;
                hitr_cnt ++;
            }
            //
            hbmarker_grp_hic[this_ctg][this_win_sta].lg_hic_sum.clear();
            hbmarker_grp_hic[this_ctg][this_win_sta].lg_hic_sum_cnt.clear();
            // print to screen            
            for(int ii = 0; ii < final_lgs.size(); ii ++)
            {
                cout << "      * final-2 ordering: LG-" << final_lgs[ii] << " with score " << final_hic[ii] << endl;
            }        
            cout << endl;            
            // output             
            int out_num = -1;
            if(tmp_marker2.type.compare("hap") == 0)
            {
                out_num = 1;
            }else            
            if(tmp_marker2.type.compare("dip") == 0)
            {
                out_num = 2;
            }else
            if(tmp_marker2.type.compare("trip") == 0)
            {
                out_num = 3;
            }else    
            if(tmp_marker2.type.compare("tetrap") == 0)
            {
                out_num = 4;
            }else
            if(tmp_marker2.type.compare("rep") == 0)
            {
                ofp  <<  this_ctg 
                     << "\t"        << this_win_sta 
                     << "\t"        << tmp_marker2.end 
                     << "\t"        << tmp_marker2.type 
                     << "\t"        << "-1"
                     << "\t"        << homLG
                   //<< "\traw\t"   << tmp_marker.raw[ 0 ]
                   //<< "\t"        << raw_line_part
                     << endl;  
                hbmarker_grp_hic[this_ctg][this_win_sta].lg_hic_sum.insert(std::pair<string, double>("-1", 0));
            }            
            for (int oi = 0; oi < out_num; oi ++)
            {
                if(tmp_marker2.type.compare("tetrap") == 0)
                {
                    for(int ti = 0; ti<4; ti ++)
                    {
                        ofp  << this_ctg 
                             << "\t"        << this_win_sta 
                             << "\t"        << tmp_marker2.end 
                             << "\t"        << tmp_marker2.type 
                             << "\t"        << all_lgs[ ti ] // ca
                           //<< "\t"        << homLG
                           //<< "\t"        << tmp_marker2.raw[ 0 ]
                             << endl;
                     }
                     break;
                }else          
                if( tmp_marker2.type.compare("hap") == 0 ) 
                {   
                    string this_hap_grp = "-1";
                    if(tmp_marker2.type.find("assignbb") != std::string::npos)              
                    {
                         this_hap_grp = tmp_marker2.lg;
                    }else
                    if(final_lgs.size() > 0)
                    {
                         this_hap_grp = final_lgs[0];
                    }else ;
                    ofp  << this_ctg 
                         << "\t"        << this_win_sta 
                         << "\t"        << tmp_marker2.end 
                         << "\t"        << tmp_marker2.type 
                         << "\t"        << this_hap_grp
                         << "\t"        << homLG
                       //<< "\traw\t"   << tmp_marker2.raw[ oi ]
                       //<< "\t"        << raw_line_part
                         << endl;   
                    hbmarker_grp_hic[this_ctg][this_win_sta].lg_hic_sum.insert(std::pair<string, double>(this_hap_grp, 111 ));                         
                }else            
                if( final_lgs.size() >= out_num) // to be replaced wih "oi < final_lgs.size()"? 
                {    
                    ofp  << this_ctg 
                         << "\t"        << this_win_sta 
                         << "\t"        << tmp_marker2.end 
                         << "\t"        << tmp_marker2.type 
                         << "\t"        << final_lgs[ oi ]
                         << "\t"        << homLG
                       //<< "\traw\t"   << tmp_marker2.raw[ oi ]
                       //<< "\t"        << raw_line_part
                         << endl;   
                    hbmarker_grp_hic[this_ctg][this_win_sta].lg_hic_sum.insert(std::pair<string, double>(final_lgs[ oi ], final_hic[ oi ] ));                         
                }else
                {
                    ofp  << this_ctg 
                         << "\t"        << this_win_sta 
                         << "\t"        << tmp_marker2.end 
                         << "\t"        << tmp_marker2.type 
                         << "\t"        << "-1"
                         << "\t"        << homLG
                       //<< "\traw\t"   << tmp_marker2.raw[ oi ]
                       //<< "\t"        << raw_line_part                         
                         << endl;  
                    hbmarker_grp_hic[this_ctg][this_win_sta].lg_hic_sum.insert(std::pair<string, double>("-1", 0 ));
                }
            }                                                
            //
            witr ++;
        }
        //
        sum_itr ++;
    }
    ofp.close();
    // output 2
    string out_file2 = out_folder + "/s8_grouping_window_markers_tmp2_overall_hic.txt\0";
    ofstream ofp2;
    ofp2.open(out_file2.c_str(), ios::out);
    if(!ofp2.good())
    {
        cout   << "   Error: cannot open file " << out_file2 << endl;
        return false;
    }           
    map<string, map<unsigned long, winMARKER > >::iterator preitr;
    map<string, map<unsigned long, winMARKER > >::iterator preitr_end;
    preitr     = hbmarker.begin();
    preitr_end = hbmarker.end();
    while(preitr != preitr_end)
    {
        string this_ctg = (*preitr).first;
        map<unsigned long, winMARKER > tmp_win_list = (*preitr).second;
        map<unsigned long, winMARKER >::iterator mpitr;
        map<unsigned long, winMARKER >::iterator mpitr_end;
        mpitr     = tmp_win_list.begin();
        mpitr_end = tmp_win_list.end();
        while(mpitr != mpitr_end)
        {
            unsigned long win_sta = (*mpitr).first;
            winMARKER tmp_marker  = (*mpitr).second;
            for(int ri = 0; ri < tmp_marker.raw.size(); ri ++)
            {
                /*
                    utg000007l_pilon	1000001	1420000	hap	4	CHR1	utg000129l_pilon	35.3261	assign1_c
                    utg000007l_pilon	1420001	1440000	dip	4	CHR1	utg000129l_pilon	10	assign1_c
                    utg000007l_pilon	1420001	1440000	dip	2	CHR1	utg001439l_pilon	4.77637	assign2_c             
                */
                vector<string> rawinfo = split_string(tmp_marker.raw[ri], '\t');
                string raw_ctg = rawinfo[0];
                unsigned long raw_win_sta   = strtoul(rawinfo[1].c_str(), NULL, 0);
                unsigned long raw_win_end   = strtoul(rawinfo[2].c_str(), NULL, 0);                
                string        raw_win_type  = rawinfo[3];
                string        raw_curr_lg   = rawinfo[4];
                string        raw_homLG     = rawinfo[5];
                string        raw_line_part = rawinfo[6] + "\t" + rawinfo[7] + "\t" + rawinfo[8]+""; // "2" turned off
                // check if this exists in hbmarker_grp_hic
                if(hbmarker_grp_hic.find(raw_ctg) != hbmarker_grp_hic.end() )
                {
                    if(hbmarker_grp_hic[raw_ctg].find(raw_win_sta) != hbmarker_grp_hic[raw_ctg].end() )
                    {
                        map<string, double> lg_hic_sum = hbmarker_grp_hic[raw_ctg][raw_win_sta].lg_hic_sum;
                        map<string, double>::iterator uplgitr;
                        map<string, double>::iterator uplgitr_end;
                        uplgitr     = lg_hic_sum.begin();
                        uplgitr_end = lg_hic_sum.end();
                        int rri     = 0;
                        while(rri < ri)
                        {
                            rri ++;
                            uplgitr ++;
                        }
                        if(rri < lg_hic_sum.size() )
                        {
                            ofp2 << raw_ctg          << "\t"
                                 << raw_win_sta      << "\t"
                                 << raw_win_end      << "\t"
                                 << raw_win_type     << "\t"
                                 << (*uplgitr).first << "\t"
                                 << raw_homLG        << "\t"
                                 << raw_line_part    << endl;    
                        }else
                        {
                            ofp2 << tmp_marker.raw[ri] << "" << endl;
                        }                    
                                         
                    }else
                    {
                        ofp2 << tmp_marker.raw[ri] << "" << endl;                                            
                    }
                }else
                {
                    ofp2 << tmp_marker.raw[ri] << "" << endl;                    
                }
            }
            mpitr ++;
        }                
        //
        preitr ++;
    }
    ofp2.close();
    //
    return true;
}
//
bool further_extend_hic_binning_lg(string               new_hbmarker_file, 
                           map<string, unsigned long>   target_contig_size,
                           map<string, double>          ctg_hap_ratio,
                           map<string, int>*            ctg2_hic_group,
                           map<int, map<string, int> >* hic_group_2ctg)
{
    /*  function: read cluster of haplotigs from a dot file
            ctg2_hic_group........: <ctg_id, clusterid_old>
            hic_group_2ctg........: <clusterid_old, <ctg_id, clusterid_new> >                 
    */
    // s0. find map of old and new group id 
    map<int, int> new2old_grp;
    map<int, map<string, int> >::iterator gitr;
    map<int, map<string, int> >::iterator gitr_end;
    gitr     = (*hic_group_2ctg).begin();
    gitr_end = (*hic_group_2ctg).end();
    while(gitr != gitr_end)
    {
        int hic_grpid_old = (*gitr).first;
        map<string, int> tmp_ctg2grp = (*gitr).second;
        map<string, int>::iterator tgitr = tmp_ctg2grp.begin();
        int hic_grpid_new = (*tgitr).second;
        cout << "   checkl: old group id " << hic_grpid_old  << " => new group id " << hic_grpid_new << endl;
        new2old_grp.insert( std::pair<int, int>(hic_grpid_new, hic_grpid_old) );
        //
        gitr ++;
    }
    // s1. read lg info
    ifstream ifp;
    ifp.open(new_hbmarker_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << new_hbmarker_file << endl;
        return false;
    }
    map<string, map<string, int> > ctg_hap_win_lg_cnt; //<ctg, <lg, win_cnt> >
    while(ifp.good() )
    {
        string line("");
        getline(ifp, line);
        if(line.size() == 0 || line[0]=='#') continue;  
        /*
	    utg000003l_pilon	1	500000	hap	2	CHR8	utg000003l_pilon	111	assignbb
	    utg000003l_pilon	500001	1000000	hap	2	CHR8	utg000003l_pilon	111	assignbb
	    utg000003l_pilon	1000001	1500000	hap	2	CHR8	utg000003l_pilon	111	assignbb
	    utg000003l_pilon	1500001	2000000	hap	2	CHR8	utg000003l_pilon	111	assignbb
	    utg000003l_pilon	2000001	2500000	hap	2	CHR8	utg000003l_pilon	111	assignbb
	    utg000003l_pilon	2500001	2847836	hap	2	CHR8	utg000003l_pilon	111	assignbb
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
        //
        assert(ctg_hap_ratio.find(this_ctg) != ctg_hap_ratio.end() );
        if(ctg_hap_ratio[this_ctg]>0.9 && this_type.compare("hap")==0)
        {
            if(ctg_hap_win_lg_cnt.find(this_ctg) == ctg_hap_win_lg_cnt.end() )
            {
                map<string, int> tmp_lg_cnt;
                tmp_lg_cnt.insert(std::pair<string, int>(this_lg, 1));
                ctg_hap_win_lg_cnt.insert(std::pair<string, map<string, int> >(this_ctg, tmp_lg_cnt));
            }else
            {
                if( ctg_hap_win_lg_cnt[this_ctg].find(this_lg) == ctg_hap_win_lg_cnt[this_ctg].end() )
                {
                    ctg_hap_win_lg_cnt[this_ctg].insert(std::pair<string, int>(this_lg, 1));
                }else
                {
                    ctg_hap_win_lg_cnt[this_ctg][this_lg] += 1;
                }
            }
        }
    }   
    ifp.close();  
    // s2. analyze lg info 
    map<string, map<string, int> >::iterator chitr;
    map<string, map<string, int> >::iterator chitr_end;
    chitr     = ctg_hap_win_lg_cnt.begin();
    chitr_end = ctg_hap_win_lg_cnt.end();
    while(chitr != chitr_end)
    {
        string this_ctg  = (*chitr).first;        
        cout << "   checkl: " << this_ctg << ": ";  
        //
        int    max_cnt = 0;
        string max_lg  = "-1";
        int    hic_grpid_new = -1;
        int    hic_grpid_old = -1;
        map<string, int> tmp_lg_cnt = (*chitr).second;
        map<string, int>::iterator lgitr;
        map<string, int>::iterator lgitr_end;
        lgitr    = tmp_lg_cnt.begin();
        lgitr_end = tmp_lg_cnt.end();
        while(lgitr != lgitr_end)
        {
            int this_cnt = (*lgitr).second;
            cout << (*lgitr).first << " with cnt " << this_cnt; 
            if(this_cnt > max_cnt)
            {
                max_cnt = this_cnt;
                max_lg  = (*lgitr).first;
                hic_grpid_new = atoi ( (*lgitr).first.c_str() );
                assert( new2old_grp.find(hic_grpid_new) != new2old_grp.end() );
                hic_grpid_old = new2old_grp[ hic_grpid_new ];                
            }
            lgitr ++;
        }
        cout << ", => lg " << max_lg << ", win-cnt " << max_cnt << endl;
        //
        if((*ctg2_hic_group).find(this_ctg) == (*ctg2_hic_group).end() )
        {
            if(!add_vet_to_group(this_ctg, 
                                 hic_grpid_old, 
                                 hic_grpid_new, 
                                 hic_group_2ctg, 
                                 ctg2_hic_group) )
            {
                return false;
            } 
        }       
        //
        chitr ++;
    }
    // s3. if output 
    bool checkl = true;
    if(checkl)
    {
        gitr     = (*hic_group_2ctg).begin();
        gitr_end = (*hic_group_2ctg).end();
        while(gitr != gitr_end)
        {
            int            hic_grpid_old = (*gitr).first;
            map<string, int> tmp_ctg2grp = (*gitr).second;
            map<string, int>::iterator tgitr     = tmp_ctg2grp.begin();
            map<string, int>::iterator tgitr_end = tmp_ctg2grp.end();
            while(tgitr != tgitr_end)
            {
                string    this_ctg = (*tgitr).first ;
                int hic_grpid_new = (*tgitr).second;
                cout << "   checkl: " << this_ctg
                     << "\told="      << hic_grpid_old
                     << "\tnew="      << hic_grpid_new
                     << "\thapr="     << ctg_hap_ratio[ this_ctg ]
                     << endl;
                //
                tgitr ++;
            }
            //
            gitr ++;
        }    
    }    
    //  
    return true;
}
bool refine_hb_window_marker_intra_ctg(vector<string> ctg_order,
         map<string, map<unsigned long, winMARKER > > hbmarker,
         map<string, map<string, string> >            ctg2lg,
         map<string, unsigned long>                   target_contig_size,
         string                                       out_label,
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
    string out_file = out_folder + "/s8_grouping_window_markers_refined_" + out_label + ".txt\0";
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
        map<string, string> link_lg_best_ctg;    // ctg in this lg showing best link 
        map<string, double> link_lg_best_ctg_hic;// ctg in this lg showing best link of value
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
            unsigned long size_win   = tmp_marker.end - win_sta + 1;
            vector<string> this_ctgs = tmp_marker.ctg; // linking contig 
            vector<double> this_hic  = tmp_marker.hic; // linking hic contact                     
            vector<string> this_LGs  = tmp_marker.wlg; // linking group id 
            assert(this_ctgs.size()==this_hic.size() && this_ctgs.size()==this_LGs.size());                
            vector<string>::iterator lgitr;
            vector<string>::iterator lgitr_end;
            lgitr     = this_LGs.begin();
            lgitr_end = this_LGs.end();
            int i = 0;
            while(lgitr != lgitr_end)
            {
                if(target_contig_size.find( this_ctgs[i] ) != target_contig_size.end() && 
                   target_contig_size[ this_ctgs[i] ] > 0) // only use long contig to phase other windows
                {
                    if(lgitr != this_LGs.begin())
                    {
                        cout << ";{" << this_LGs[i] << "," << this_ctgs[i] << "," << this_hic[i]*size_win << "}";
                    }else
                    {
                        cout <<  "{" << this_LGs[i] << "," << this_ctgs[i] << "," << this_hic[i]*size_win << "}";
                    }
                    // how much hic the ctg linked to this lg and how many time the ctg linked to this lg
                    if(link_lg_hic.find( this_LGs[i] ) == link_lg_hic.end() )
                    {
                       link_lg_hic.insert(std::pair<string, double>( this_LGs[i], this_hic[i]*size_win ) );
                       link_lg_cnt.insert(   std::pair<string, int>( this_LGs[i], 1 ) );
                       link_lg_best_ctg.insert(    std::pair<string, string>(this_LGs[i], this_ctgs[i]));
                       link_lg_best_ctg_hic.insert(std::pair<string, double>(this_LGs[i], this_hic[i]));
                    }else
                    {
                       link_lg_hic[ this_LGs[i] ] += this_hic[i]*size_win;
                       link_lg_cnt[ this_LGs[i] ] += 1; 
                       // update the best linked ctg  
                       if(this_hic[i] > link_lg_best_ctg_hic[ this_LGs[i] ] )
                       {
                           link_lg_best_ctg_hic[ this_LGs[i] ] = this_hic[i];
                           link_lg_best_ctg[     this_LGs[i] ] = this_ctgs[i];
                       }
                    }
                    // how each lg is linked with current ctg under each type 
                    if( type_lg_hic.find(tmp_marker.type) == type_lg_hic.end() )
                    {
                        map<string, double> tmp_lg_hic;
                        tmp_lg_hic.insert(std::pair<string, double>( this_LGs[i], this_hic[i]*size_win  ));
                        type_lg_hic.insert(std::pair<string, map<string, double> >(tmp_marker.type, tmp_lg_hic));
                    }else
                    {
                        if( type_lg_hic[tmp_marker.type].find( this_LGs[i] ) ==  type_lg_hic[tmp_marker.type].end() )
                        {
                            type_lg_hic[tmp_marker.type].insert(std::pair<string, double>( this_LGs[i], this_hic[i]*size_win  ));
                        }else
                        {
                            type_lg_hic[tmp_marker.type][ this_LGs[i] ] += this_hic[i] * size_win;
                        }
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
        cout << "       : " << out_label << endl;
        cout << "       : " << this_ctg  << " with " << type_along_ctg.size() << " types of markers: " << endl; 
        // analyze total hic-score/cnt 
        map<string, unsigned long>::iterator titr; // <type, 1>
        map<string, unsigned long>::iterator titr_end;
        titr     = type_along_ctg.begin();
        titr_end = type_along_ctg.end();
        while(titr != titr_end)
        {
            string this_mkr_type = (*titr).first;
            if(type_lg_hic.find( this_mkr_type ) != type_lg_hic.end() )          
            {
                cout << "         | " << (*titr).first << ", " << (*titr).second << " bp: " << endl;
                //assert( type_lg_hic.find( this_mkr_type ) != type_lg_hic.end() );
                map<string, double> tmp_lg_hic = type_lg_hic[ this_mkr_type ];
                map<string, double>::iterator tlgitr; // lg under this type 
                map<string, double>::iterator tlgitr_end;
                tlgitr     = tmp_lg_hic.begin();
                tlgitr_end = tmp_lg_hic.end();
                while(tlgitr != tlgitr_end)
                {
                    cout << "            - LG-"      << (*tlgitr).first 
                         << " with total hic-score " << (*tlgitr).second * 1.0
                         << endl;
                    //
                    tlgitr ++;
                }
            }else
            {
                cout << "         | " << (*titr).first << ", " << (*titr).second << " bp: no linked long ctg." << endl;                
            }
            //
            titr ++;
        }
        // sort the candidate linkage groups according to the hic-score 
        vector<string> final_lgs;   // <1st_best_lg,     2nd_best_lg,     ...>
        vector<double> final_hic;   // <1st_best_lg_hic, 2nd_best_lg_hic, ...>
        vector<string> best_ctg;    // <1st_best_ctg,    2nd_best_ctg,    ...>
        vector<double> best_ctg_hic;// <1st_best_ctg_hic,2nd_best_ctg_hic,...>        
        map<string, double>::iterator hitr;
        map<string, double>::iterator hitr_end;
        hitr     = link_lg_hic.begin();
        hitr_end = link_lg_hic.end();
        while(hitr != hitr_end)
        {
            string linked_lg = (*hitr).first;
            assert( link_lg_cnt.find(linked_lg) != link_lg_cnt.end() );
            cout << "         | total   hic-link score with LG-"<< linked_lg << ": "<< link_lg_hic[ linked_lg ] << endl;
            cout << "         | total   hic-link cnt   with LG-"<< linked_lg << ": "<< link_lg_cnt[ linked_lg ] << endl;
            if(link_lg_cnt[ linked_lg ] == 0) link_lg_cnt[ linked_lg ] = 1;
            double avg_to_grp_hic   = link_lg_hic[ linked_lg ]*1.0 / link_lg_cnt[ linked_lg ]; 
            double total_to_grp_hic = link_lg_hic[ linked_lg ]*1.0;
            cout << "         | average hic-link score with LG-" << linked_lg 
                 << ": "                                         << avg_to_grp_hic
                 << endl;            
            //
            if(final_lgs.size() == 0)
            {
                final_lgs.push_back(    linked_lg   );
                final_hic.push_back( total_to_grp_hic );    
                //
                best_ctg.push_back( link_lg_best_ctg[ linked_lg ] );
                best_ctg_hic.push_back( link_lg_best_ctg_hic[ linked_lg ] );
            }else
            {
                bool tmp_inserted = false;
                vector<string>::iterator tmp_lg_itr      = final_lgs.begin();    // lg
                vector<double>::iterator tmp_lg_hic_itr  = final_hic.begin();    // lg-hic
                vector<string>::iterator tmp_ctg_itr     = best_ctg.begin();     // lg-ctg
                vector<double>::iterator tmp_ctg_hic_itr = best_ctg_hic.begin(); // lg-ctg-hic
                for(int ii = 0; ii < final_lgs.size(); ii ++)
                {
                    if( final_hic[ii] <= total_to_grp_hic )
                    {
                        final_lgs.insert(tmp_lg_itr, linked_lg);
                        final_hic.insert(tmp_lg_hic_itr, total_to_grp_hic);
                        //
                        best_ctg.insert(tmp_ctg_itr, link_lg_best_ctg[ linked_lg ] );
                        best_ctg_hic.insert(tmp_ctg_hic_itr, link_lg_best_ctg_hic[ linked_lg ] );             
                        //
                        tmp_inserted = true;                        
                        break;
                    }
                    tmp_lg_itr ++;
                    tmp_lg_hic_itr ++;    
                    tmp_ctg_itr ++;  
                    tmp_ctg_hic_itr ++;              
                }
                if(!tmp_inserted)
                {
                    final_lgs.insert(tmp_lg_itr, linked_lg);
                    final_hic.insert(tmp_lg_hic_itr, total_to_grp_hic);
                    //
                    best_ctg.insert(tmp_ctg_itr, link_lg_best_ctg[ linked_lg ] );
                    best_ctg_hic.insert(tmp_ctg_hic_itr, link_lg_best_ctg_hic[ linked_lg ] );                                  
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
        // output updated grouping to file: utg000235l_pilon 590001 750000 trip 1 CHR1 utg000148l_pilon 11.3636 assign1
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
            }else
            if(tmp_marker.type.compare("rep") == 0)
            {
                vector<string> rawinfo = split_string(tmp_marker.raw[ 0 ], '\t');
                string raw_line_part  = rawinfo[6] + "\t" + rawinfo[7] + "\t" + rawinfo[8]; // caution: this is raw link!                
                ofp  <<  this_ctg 
                     << "\t"        << win_sta 
                     << "\t"        << tmp_marker.end 
                     << "\t"        << tmp_marker.type 
                     << "\t"        << "-1"
                     << "\t"        << tmp_marker.homlg
                   //<< "\traw\t"   << tmp_marker.raw[ 0 ]
                     << "\t"        << raw_line_part
                     << endl;
                out_num = -1;
            }
            for (int oi = 0; oi < out_num; oi ++)
            {
                /*
                cout << "        checkz: output refinement for " << tmp_marker.raw[ oi ] 
                     << ", final_lgs.size()="               << final_lgs.size()
                     << ", best_ctg.size()="                << best_ctg.size()
                     << ", best_ctg_hic.size()="            << best_ctg_hic.size()
                     << endl;
                */
                vector<string> rawinfo = split_string(tmp_marker.raw[ oi ], '\t');
                string raw_line_part  = rawinfo[6] + "\t" + rawinfo[7] + "\t" + rawinfo[8];            
                if(tmp_marker.type.compare("tetrap") == 0)
                {
                    for(int ti = 0; ti<4; ti ++)
                    {
                        ofp  << this_ctg 
                             << "\t"        << win_sta 
                             << "\t"        << tmp_marker.end 
                             << "\t"        << tmp_marker.type 
                             << "\t"        << all_lgs[ ti ] // ca
                             << "\t"        << rawinfo[5]
                             << "\t"        << raw_line_part
                             << endl;
                     }
                     break;             
                }else
                if( final_lgs.size() >= out_num) // to be replaced wih "oi < final_lgs.size()"? 
                {   
                    ofp  << this_ctg 
                         << "\t"        << win_sta 
                         << "\t"        << tmp_marker.end 
                         << "\t"        << tmp_marker.type 
                         << "\t"        << final_lgs[ oi ]
                         << "\t"        << tmp_marker.homlg
                         << "\t"        << best_ctg[ oi ]
                         << "\t"        << best_ctg_hic[ oi ]
                         << "\t"        << rawinfo[8]
                       //<< "\t"        << raw_line_part
                         << endl;
                }else
                {
                    ofp  << this_ctg 
                         << "\t"        << win_sta 
                         << "\t"        << tmp_marker.end 
                         << "\t"        << tmp_marker.type 
                         << "\t"        << "-1"
                         << "\t"        << tmp_marker.homlg
                       //<< "\traw\t"   << tmp_marker.raw[ oi ]
                         << "\t"        << raw_line_part                         
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
                                                hic=<val1/val2/..>, 
                                                raw=<line1/line2/...>,
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
            vector<string> raw; // raw line
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
        if(this_lg.compare("-1") != 0)
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
        cout << "   checkr1: ctg win marker distribution in LGs: " << endl;
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
                vector<string> this_ctgs = tmp_marker.ctg; // linking contig 
                vector<double> this_hic  = tmp_marker.hic; // linking hic contact                     
                vector<string> this_LGs  = tmp_marker.wlg; // linking group id 
                assert(this_ctgs.size()==this_hic.size() && this_ctgs.size()==this_LGs.size());                
                vector<string>::iterator lgitr;
                vector<string>::iterator lgitr_end;
                lgitr     = this_LGs.begin();
                lgitr_end = this_LGs.end();
                int i = 0;
                while(lgitr != lgitr_end)
                {
                    if(lgitr != this_LGs.begin())
                    {
                        cout << ";{" << *lgitr << "," << this_ctgs[i] << "," << this_hic[i] << "}";
                    }else
                    {
                        cout <<  "{" << *lgitr << "," << this_ctgs[i] << "," << this_hic[i] << "}";
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
        cout << "   checkr2: LG distribution in CTGs: " << endl;        
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
        cout << "   checkr3:  homologous LGs distribution 1-4 versus 1-12: " << endl;
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
