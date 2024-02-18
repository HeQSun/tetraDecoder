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
bool read_group_wise_hb_window_marker(string           hbmarker_file, 
         map<string, map<string, map<unsigned long, unsigned long> > >* lg_ctg_win_markers);  
bool calculate_inter_group_hic(map<string, map<unsigned long, map<unsigned long, WHIC> > > hic_matrix,
         map<string, map<string, map<unsigned long, unsigned long> > > lg_ctg_win_markers,
         double                                        normalized_reads_scale,
         string                                        out_folder);           
bool refine_hb_window_marker_intra_ctg(vector<string>  ctg_order,
         map<string, map<unsigned long, winMARKER > >  hbmarker,
         map<string, map<string, string> >             ctg2lg,
         map<string, unsigned long>                    target_contig_size,
         string                                        out_label,
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
         map<string, unsigned long>                    allelic_map_seq_size,
         double                                        normalized_reads_scale)
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
    hbmarker_file = out_folder + "/s8_grouping_window_markers_refined_1st.txt\0"; 
    map<string, map<string, map<unsigned long, unsigned long> > > lg_ctg_win_markers;
    if(!read_group_wise_hb_window_marker(hbmarker_file, &lg_ctg_win_markers))
    {
        cout << "   Error: failed in collected group-wise markers. " << endl;
        return false;
    }
    if(!calculate_inter_group_hic(hic_matrix, lg_ctg_win_markers, normalized_reads_scale, out_folder))
    {
        /*
            map<string, map<unsigned long, map<unsigned long, WHIC> > >   hic_matrix
            map<string, map<string, map<unsigned long, unsigned long> > > lg_ctg_win_markers
        */
        cout << "   Error: failed in calculating Hi-C within/between linkage groups. " << endl;
        return false;
    }
    //
    /* TODO: refinement within groups not use yet!
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
    // Rebuild allelic groups
    //   <ctg, <allelic_grp, <allelic_ctg, allelic_ctg_size > > >
    //   ctg has a lower chance to be assigned to allelic_grp, where the prob can be defined by number of allelic ctgs.
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
    */      
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    //
    return true;
}
//
bool calculate_inter_group_hic(map<string, map<unsigned long, map<unsigned long, WHIC> > > hic_matrix,
         map<string, map<string, map<unsigned long, unsigned long> > > lg_ctg_win_markers,
         double normalized_reads_scale,
         string out_folder)
{
    /*  Function: calculate total Hi-C contact within each linkage group and between pairs of linkage group
                  lg_ctg_win_markers.: <lg_id, <ctg_id, <win_sta, win_end> > >
                  hic_matrix.........: < ctg1-ctg2, <sta1, <sta2, {end1,typ1,end2,typ2,hic-cnt} > > >                  
    */
    if(lg_ctg_win_markers.find("-1") != lg_ctg_win_markers.end() )
    {
        lg_ctg_win_markers.erase("-1");
    }
    string hic_sum_file = out_folder + "/s8_grouping_window_markers_refined_1st_group_hic_stats.txt\0"; 
    ofstream ofp;
    ofp.open(hic_sum_file.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: failed in opening file for writing Hi-C stats: " << hic_sum_file << endl;
        return false;
    }
    int N_cluster = 4;
    if(lg_ctg_win_markers.size() < N_cluster) // 
    {
        cout << "   warning: the number of clusters is only " << lg_ctg_win_markers.size() 
             << ", smaller than expected value of "           << N_cluster
             << ". you should increase initial Hi-C score to get more clusters." << endl;
        for(int gi = 1; gi <= N_cluster; gi ++)
        {
            for(int gj = gi; gj <= N_cluster; gj ++)
            {
                ofp << "group\t" << gi << "\t" << gj << "\t0\t" << "avg\t0" << endl; 
            }
        }
        ofp.close();
        return false;
    }
    // result variable
    map<string, map<string, double> > grp_contact; // <grp1, <grp2, sum_grp12_hic > >
    //
    map<string, map<string, map<unsigned long, unsigned long> > >::iterator lgitr_1;
    map<string, map<string, map<unsigned long, unsigned long> > >::iterator lgitr_end;
    map<string, map<string, map<unsigned long, unsigned long> > >::iterator lgitr_2;
    lgitr_1   = lg_ctg_win_markers.begin();
    lgitr_end = lg_ctg_win_markers.end();
    while(lgitr_1 != lgitr_end)
    {
        string lg_1 = (*lgitr_1).first;
        map<string, map<unsigned long, unsigned long> > tmp_ctg_win_1 = (*lgitr_1).second;
        // initialize result - a
        map<string, double> sum_grp12_hic;
        grp_contact.insert(std::pair<string, map<string, double> >(lg_1, sum_grp12_hic));        
        //
        lgitr_2 = lgitr_1;        
        while(lgitr_2 != lgitr_end)
        {
            string lg_2 = (*lgitr_2).first;
            map<string, map<unsigned long, unsigned long> > tmp_ctg_win_2 = (*lgitr_2).second;
            // initialize result - b
            grp_contact[lg_1].insert(std::pair<string, double>(lg_2, 0));
            int individual_cnts = 0;            
            //
            map<string, map<unsigned long, unsigned long> >::iterator citr_1;
            map<string, map<unsigned long, unsigned long> >::iterator citr_1_end;
            citr_1     = tmp_ctg_win_1.begin();
            citr_1_end = tmp_ctg_win_1.end();
            while(citr_1 != citr_1_end)
            {
                string ctg_1 = (*citr_1).first;
                map<unsigned long, unsigned long> this_wins_1 = (*citr_1).second;
                if(tmp_ctg_win_2.find(ctg_1) != tmp_ctg_win_2.end() &&
                   lg_2.compare(lg_1) != 0
                  )
                {
                    // skipping contigs in both linkage groups
                    citr_1 ++;
                    continue;
                }
                //
                map<string, map<unsigned long, unsigned long> >::iterator citr_2;
                map<string, map<unsigned long, unsigned long> >::iterator citr_2_end;
                citr_2     = tmp_ctg_win_2.begin();
                citr_2_end = tmp_ctg_win_2.end();                
                while(citr_2 != citr_2_end)
                {
                    string ctg_2 = (*citr_2).first;
                    map<unsigned long, unsigned long> this_wins_2 = (*citr_2).second;
                    //
                    if(tmp_ctg_win_1.find(ctg_2) != tmp_ctg_win_1.end() &&
                       lg_2.compare(lg_1) != 0 )
                    {
                        // skipping contigs in both linkage groups                    
                        citr_2 ++;
                        continue;
                    }
                    //
                    string key1 = ctg_1 + "-" + ctg_2;
                    string key2 = ctg_2 + "-" + ctg_1;          
                    map<unsigned long, map<unsigned long, WHIC> > ctg_pair_hic;
                    string key_found = key1;                        
                    if( hic_matrix.find(key1) != hic_matrix.end() )
                    {
                        ctg_pair_hic = hic_matrix[key1];
                    }else
                    if( hic_matrix.find(key2) != hic_matrix.end() )
                    {
                        ctg_pair_hic = hic_matrix[key2];    
                        key_found    = key2;                    
                    }else ;
                    //
                    map<unsigned long, map<unsigned long, WHIC> >::iterator hitr;
                    map<unsigned long, map<unsigned long, WHIC> >::iterator hitr_end;
                    hitr     = ctg_pair_hic.begin();
                    hitr_end = ctg_pair_hic.end();
                    while(hitr != hitr_end)
                    {
                        unsigned long sta1 = (*hitr).first;
                        map<unsigned long, WHIC> tmp_hic_val = (*hitr).second;
                        //
                        map<unsigned long, WHIC>::iterator whitr;
                        map<unsigned long, WHIC>::iterator whitr_end;
                        whitr     = tmp_hic_val.begin();
                        whitr_end = tmp_hic_val.end();
                        while(whitr != whitr_end)
                        {   
                            unsigned long sta2 = (*whitr).first;                                                         
                            WHIC tmp_whic = (*whitr).second;
                            unsigned long ctg1_win_size = tmp_whic.end1 - sta1 + 1;
                            unsigned long ctg2_win_size = tmp_whic.end2 - sta2 + 1;
                            double normalized_reads = tmp_whic.hcnt*1.0 / (ctg1_win_size+ctg2_win_size) * normalized_reads_scale;                                
                            grp_contact[lg_1][lg_2] += normalized_reads;
                            individual_cnts ++;                            
                            //
                            whitr ++;
                        }
                        //
                        hitr ++;
                    }                    
                    //
                    citr_2 ++;
                }
                //
                citr_1 ++;
            }
            if(individual_cnts == 0) individual_cnts = 1;
            cout << "   checkgh: hic of group " << lg_1 
                 <<                         "-" << lg_2 
                 <<                        ": " << grp_contact[lg_1][lg_2] 
                 <<               ", avg hic: " << grp_contact[lg_1][lg_2] / individual_cnts
                 << endl;
            ofp  << std::fixed
                 <<      "group\t" << lg_1 
                 <<           "\t" << lg_2 
                 <<           "\t" << grp_contact[lg_1][lg_2] 
                 <<      "\tavg\t" << grp_contact[lg_1][lg_2] / individual_cnts
                 << endl;
            //
            lgitr_2 ++;
        }
        //
        lgitr_1 ++;
    }
    //
    ofp.close();
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
    /*  Function: read cluster of haplotigs from a dot file
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
       Function: refine grouping of window markers within each contig 
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
                        string this_lgi = all_lgs[0];
                        if(all_lgs.size() > ti)
                        {
                            this_lgi = all_lgs[ti];
                        }
                        ofp  << this_ctg 
                             << "\t"        << win_sta 
                             << "\t"        << tmp_marker.end 
                             << "\t"        << tmp_marker.type 
                             << "\t"        << this_lgi // ca
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
       Function: read hic_binning grouped window markers       
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
    cout << "       : " << (*ctg_order).size() << " ctgs checked, with " << win_num << " window markers. " << endl;
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
bool read_group_wise_hb_window_marker(string hbmarker_file, 
         map<string, map<string, map<unsigned long, unsigned long> > >* lg_ctg_win_markers)
{
    /*
       Function: read hic_binning grouped window markers       
       return..: lg_ctg_win_markers.: <lg_id, <ctg_id, <win_sta, win_end> > >
    */
    map<string, map<string, map<unsigned long, unsigned long> > > tmp_lg_ctg_win;
    map<string, vector<string> > tmp_lg_ctg_order; // <lg_id, <ctg_id> >
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
        //
        if(tmp_lg_ctg_win.find(this_lg) == tmp_lg_ctg_win.end() )
        {
            map<unsigned long, unsigned long> tmp_win;
            tmp_win.insert(std::pair<unsigned long, unsigned long>(win_sta, win_end));
            map<string, map<unsigned long, unsigned long> > tmp_ctg_win;
            tmp_ctg_win.insert(std::pair<string, map<unsigned long, unsigned long> >(this_ctg, tmp_win));
            tmp_lg_ctg_win.insert(std::pair<string, map<string, map<unsigned long, unsigned long> > >(this_lg, 
                                                                                                      tmp_ctg_win));            
            vector<string> tmp_ctgs;
            tmp_ctgs.push_back(this_ctg);
            tmp_lg_ctg_order.insert(std::pair<string, vector<string> >(this_lg, tmp_ctgs));
        }else
        {
            if(tmp_lg_ctg_win[this_lg].find(this_ctg) == tmp_lg_ctg_win[this_lg].end() )
            {
                map<unsigned long, unsigned long> tmp_win;
                tmp_win.insert(std::pair<unsigned long, unsigned long>(win_sta, win_end));
                map<string, map<unsigned long, unsigned long> > tmp_ctg_win;   
                tmp_lg_ctg_win[this_lg].insert(std::pair<string, map<unsigned long, unsigned long> >(this_ctg, tmp_win));                             
                //
                tmp_lg_ctg_order[this_lg].push_back(this_ctg);
            }else
            {
                if(tmp_lg_ctg_win[this_lg][this_ctg].find(win_sta) == tmp_lg_ctg_win[this_lg][this_ctg].end() )
                {
                    tmp_lg_ctg_win[this_lg][this_ctg].insert(std::pair<unsigned long, unsigned long>(win_sta, win_end));
                }
            }
        }
    }
    ifp.close();   
    // output: note, the number of wins here might be smaller than given file as "-1" cases collapse.
    *lg_ctg_win_markers = tmp_lg_ctg_win;
    //
    bool check_this = true;
    if(check_this)
    {
        map<string, vector<string> >::iterator lgitr;
        map<string, vector<string> >::iterator lgitr_end;
        lgitr     = tmp_lg_ctg_order.begin();
        lgitr_end = tmp_lg_ctg_order.end();
        while(lgitr != lgitr_end)
        {
            string this_lg = (*lgitr).first;            
            vector<string> this_group_ctgs = (*lgitr).second;
            //
            assert((*lg_ctg_win_markers).find(this_lg) != (*lg_ctg_win_markers).end() );
            map<string, map<unsigned long, unsigned long> > tmp_ctg_win = (*lg_ctg_win_markers)[this_lg];
            //
            vector<string>::iterator citr;
            vector<string>::iterator citr_end;
            citr     = this_group_ctgs.begin();
            citr_end = this_group_ctgs.end();
            while(citr != citr_end)
            {
                string this_ctg = *citr;
                assert(tmp_ctg_win.find(this_ctg) != tmp_ctg_win.end() );
                map<unsigned long, unsigned long> tmp_win = tmp_ctg_win[this_ctg];
                map<unsigned long, unsigned long>::iterator witr;
                map<unsigned long, unsigned long>::iterator witr_end;
                witr     = tmp_win.begin();
                witr_end = tmp_win.end();
                while(witr != witr_end)
                {
                    cout << "   checklc: " << this_ctg
                         << "\t"           << (*witr).first
                         << "\t"           << (*witr).second 
                         << "\t"           << this_lg
                         << endl;
                    witr ++;
                }
                //
                citr ++;
            }
            //
            lgitr ++;
        }
    }
    //
    return true;
}
