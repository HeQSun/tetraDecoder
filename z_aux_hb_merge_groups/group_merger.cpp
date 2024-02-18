/* 
   Merge many clusters into four only: 
      
      input:  s6_haplotig_hic_contact_matrix_subclusters_raw.dot and zadd_hic_matrix.txt
      output: s6_haplotig_hic_contact_matrix_subclusters_merged.dot
      
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
struct NODE 
{
    unsigned long         end; // window end
    string               type; // hap/dip/trip/tetrap/rep
    unsigned long    win_size; // size of this marker
    unsigned long    ctg_size; // size of the contig having this marker
    map<int, double> grp_link_val; // the best link between this marker and each hap-extended-group
    map<int, string> grp_link_ctg; // the best ctg linking this marker and each hap-extended-group    
};
struct WHIC 
{
  //unsigned long sta1; // marker-window sta of ctg1
    unsigned long end1; // marker-window end of ctg1
    string        typ1; // type of window of ctg1: hap/dip/...
  //unsigned long sta2; // marker-window sta of ctg2
    unsigned long end2; // marker-window end of ctg2 
    string        typ2; // type of window of ctg2: hap/dip/...
    double        hcnt; // hic-link-cnt between {ctg1:sta1-end1} and {ctg2:sta2-end2}    
};
map<string, map<unsigned long, NODE> > win_marker;      // all window markers along given ctgs!!!!!!!!!!!!!!!!!!!!!!!!!!
map<string, map<unsigned long, map<unsigned long, WHIC> > > hic_matrix;//<ctg1-ctg2,<sta1,<sta2,{end1,end2,hic-cnt}> > >
//
bool read_hb_cluster(string                    hb_group_file, 
                  map<string, int>*            ctg2_hb_group,
                  map<int, map<string, int> >* hb_group_2ctg);
bool add_vet_to_group(string                   vet, 
                  int                          hb_grpid_old,
                  int                          hb_grpid_new,
                  map<int, map<string, int> >* hb_group_2ctg,
                  map<string, int>*            ctg2_hb_group);
bool read_win_marker_ctg(string                            win_marker_file,
                  map<string, unsigned long>               target_contig_size,
                  map<string, map<unsigned long, NODE> >*  win_marker);                  
bool read_hic_matrix(string hic_matrix_file);                  
bool write_hic_matrix(string tmpfolder);                  
bool read_target_ctg(string ctg_file, map<string, unsigned long>* target_contig_size);
bool merge_clusters_v2(map<int, map<string, int> >      hb_group_2ctg);
//
int main(int argc, char* argv[])
{
    if(argc < 5)
    {
        cout << "\nFunction: merge many clusters into four only.";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl
             << "\nUsage: hb_group_refiner s6_haplotig_hic_contact_matrix_subclusters_raw.dot zadd_hic_matrix.txt target_ctg_file win_marker_file"
             << endl;
        cout << endl;
        return 1;
    }
    double startT= clock();
    //    
    string hb_group_file   = (string)argv[1];
    string hic_matrix_file = (string)argv[2];
    string ctg_file        = (string)argv[3]; 
    string win_marker_file = (string)argv[4]; 
    // s1. read raw clusters
    map<string, int> ctg2_hb_group;
    map<int, map<string, int> > hb_group_2ctg;
    if(!read_hb_cluster(hb_group_file, &ctg2_hb_group, &hb_group_2ctg))
    {
        cout << "   Error: failed in reading raw clusters: " << endl;
        return 1;
    }
    // s2. read window marker info along target contigs
    cout << "   Info: step 2 - read info on window markers... " << endl;
    map<string, unsigned long> target_contig_size; // <target_ctg, size>
    if(!read_target_ctg(ctg_file, &target_contig_size))
    {
        cout << "   Error: failed in collecting target contig and size info." << endl;
        return 1;
    }    
    if(!read_win_marker_ctg(win_marker_file, target_contig_size, &win_marker) )
    {
        cout << "   Error: failed in reading coverage-defined window markers. " << endl;
        return 1;
    }    
    // s3. read hic contacts
    if(!read_hic_matrix(hic_matrix_file))
    {
        cout << "   Error: failed in reading hic contact matrix from file " << hic_matrix_file << endl;
        return 1;
    }
    // s4. find hic contact between clusters 
    if(!merge_clusters_v2(hb_group_2ctg))
    {
        cout << "   Error: failed in finding hic contact between raw clusters. " << endl;
        return 1;
    }
    
    
    
    
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    //
    return 0;
}
bool merge_clusters_v2(map<int, map<string, int> > hb_group_2ctg)
{
    /* Function: check the hic links between raw clusters. 
               hb_group_2ctg........: <old_group_id, <ctg_id, new_group_id> >
    */
    while(hb_group_2ctg.size() > 4)
    {
        cout << "   Info: ---- number of current clusters " << hb_group_2ctg.size()
             << " > expected cluster number "               << 4 
             << endl;
        int    to_merge_grp1_id = -1;
        int    to_merge_grp2_id = -1;
        double to_merge_hic_max = 0;
        //
        map<int, map<int, double> > grp_contact; // <grp1, <grp2, hic_val> >
        //
        map<int, map<string, int> >::iterator gitr1;
        map<int, map<string, int> >::iterator gitr_end;
        gitr1     = hb_group_2ctg.begin();
        gitr_end = hb_group_2ctg.end();
        while(gitr1 != gitr_end)
        {
            // get cluster 1
            int group1_id = (*gitr1).first;
            // cout << "   check: old group1 id = " << group1_id << endl;
            map<string, int> group1 = (*gitr1).second;
            // initialize result - a
            map<int, double> sum_grp12_hic;
            grp_contact.insert(std::pair<int, map<int, double> >(group1_id, sum_grp12_hic));
            //
            map<int, map<string, int> >::iterator gitr2 = gitr1;
            gitr2 ++;        
            while(gitr2 != gitr_end)
            {
                // get cluster 2
                int group2_id = (*gitr2).first;
                // cout << "        : old group2 id = " << group2_id << endl;
                map<string, int> group2 = (*gitr2).second;
                // initialize result - b
                grp_contact[group1_id].insert(std::pair<int, double>(group2_id, 0));
                int individual_cnts = 0;
                // check conctact between ctg1's and ctg2's
                map<string, int>::iterator citr1;
                map<string, int>::iterator citr1_end;
                citr1     = group1.begin();
                citr1_end = group1.end();
                while(citr1 != citr1_end)
                {
                    string ctg1 = (*citr1).first;
                    //
                    map<string, int>::iterator citr2;
                    map<string, int>::iterator citr2_end;
                    citr2     = group2.begin();
                    citr2_end = group2.end();
                    while(citr2 != citr2_end)
                    {
                        string ctg2 = (*citr2).first;
                        //
                        string key1 = ctg1 + "-" + ctg2;
                        string key2 = ctg2 + "-" + ctg1;          
                        map<unsigned long, map<unsigned long, WHIC> > ctg_pair_hic;
                        if( hic_matrix.find(key1) != hic_matrix.end() )
                        {
                            ctg_pair_hic = hic_matrix[key1];
                        }else
                        if( hic_matrix.find(key2) != hic_matrix.end() )
                        {
                            ctg_pair_hic = hic_matrix[key2];                        
                        }else ;
                        //
                        map<unsigned long, map<unsigned long, WHIC> >::iterator hitr;
                        map<unsigned long, map<unsigned long, WHIC> >::iterator hitr_end;
                        hitr     = ctg_pair_hic.begin();
                        hitr_end = ctg_pair_hic.end();
                        while(hitr != hitr_end)
                        {
                            map<unsigned long, WHIC> tmp_hic_val = (*hitr).second;
                            //
                            map<unsigned long, WHIC>::iterator whitr;
                            map<unsigned long, WHIC>::iterator whitr_end;
                            whitr     = tmp_hic_val.begin();
                            whitr_end = tmp_hic_val.end();
                            while(whitr != whitr_end)
                            {
                                grp_contact[group1_id][group2_id] += (*whitr).second.hcnt;
                                individual_cnts ++;
                                //
                                whitr ++;
                            }
                            //
                            hitr ++;
                        }
                        //
                        citr2 ++;
                    }
                    //
                    citr1 ++;
                } 
                //
                if(individual_cnts == 0) individual_cnts ++;
                cout << "       : hic contact between group-" << group1_id 
                     <<                         " and group-" << group2_id 
                     <<                                  ": " << grp_contact[group1_id][group2_id] 
                     <<                             ", avg: " << grp_contact[group1_id][group2_id] / individual_cnts
                     << endl;
                if(grp_contact[group1_id][group2_id] / individual_cnts > to_merge_hic_max)
                {
                    to_merge_hic_max = grp_contact[group1_id][group2_id] / individual_cnts;
                    to_merge_grp1_id = group1_id;
                    to_merge_grp2_id = group2_id;
                }
                //
                gitr2 ++;
            }                    
            //
            gitr1 ++;
        }
        // merge 
        if(to_merge_grp1_id != -1 && to_merge_grp1_id != -1)
        {
            map<string, int> group2 = hb_group_2ctg[to_merge_grp2_id];
            hb_group_2ctg[to_merge_grp1_id].insert(group2.begin(), group2.end() ); 
            hb_group_2ctg.erase( to_merge_grp2_id );
            cout << "   Info: group-"              << to_merge_grp2_id 
                 << " has been merged into group-" << to_merge_grp1_id 
                 << ", and thus removed list."     << endl;
        }
    }
    //
    if(1)
    {
        map<int, map<string, int> >::iterator gitr;
        map<int, map<string, int> >::iterator gitr_end;
        gitr     = hb_group_2ctg.begin();
        gitr_end = hb_group_2ctg.end();
        while(gitr != gitr_end)
        {
            int this_old_group_id = (*gitr).first;
            cout << "   merged: old group id = " << this_old_group_id << " <= new group id = ";
            map<string, int> this_group = (*gitr).second;
            map<string, int>::iterator citr;
            map<string, int>::iterator citr_end;
            citr     = this_group.begin();
            citr_end = this_group.end();
            bool first_ctg = true;
            while(citr != citr_end)
            {
                string this_ctg = (*citr).first;
                int    this_new_group_id = (*citr).second;
                if(first_ctg)
                {
                    cout << this_new_group_id << endl;
                    first_ctg = false;
                }
                cout << "        : " << this_ctg << endl;
                //
                citr ++;
            }
            //
            gitr ++;
        }
    }
    //
    return true;
}
//
bool read_hic_matrix(string hic_matrix_file)
{
    /* Function: read hic links between ctg1 window and ctg2 window
         collect info to global hic_matrix:
           map<string, map<unsigned long, map<unsigned long, WHIC> > > hic_matrix; < ctg1-ctg2, <sta1, <sta2, {end1,typ1,end2,typ2,hic-cnt} > > >
             WHIC
   	       unsigned long end1; // marker-window end of ctg1
    	       string        typ1; // type of window of ctg1: hap/dip/...
	       unsigned long end2; // marker-window end of ctg2 
	       string        typ2; // type of window of ctg2: hap/dip/...
	       double        hcnt; // hic-link-cnt between {ctg1:sta1-end1} and {ctg2:sta2-end2}	       
    */
    ifstream ifp;
    ifp.open(hic_matrix_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: failed in reading hic matrix file " << hic_matrix_file << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 9) 
        {
            cout << "   warning: insufficient line info " << line << endl;
            continue;
        }
        /*
            utg000011l_pilon	4000001	4500000	hap	utg000073l_pilon	1	10000	hap	0.196078        
        */
        string            ctg1 = lineinfo[0];
        unsigned long win1_sta = strtoul(lineinfo[1].c_str(), NULL, 0);
        unsigned long win1_end = strtoul(lineinfo[2].c_str(), NULL, 0);        
        string        win1_typ = lineinfo[3];
        string            ctg2 = lineinfo[4];
        unsigned long win2_sta = strtoul(lineinfo[5].c_str(), NULL, 0);
        unsigned long win2_end = strtoul(lineinfo[6].c_str(), NULL, 0);        
        string        win2_typ = lineinfo[7];      
        double            hcnt = atof(lineinfo[8].c_str());
        //
        if(0)
        cout << "   checkkk: " 
             << ctg1      << "\t" 
             << win1_sta  << "\t" 
             << win1_end  << "\t" 
             << win1_typ  << "\t"
             << ctg2      << "\t"
             << win2_sta  << "\t" 
             << win2_end  << "\t"
             << win2_typ  << endl;
        // map<ctg1-ctg2, map<sta1, map<sta2, WHIC> > >
        string key1 = ctg1 + "-" + ctg2;
        string key2 = ctg2 + "-" + ctg1;
        if( hic_matrix.find(key1)==hic_matrix.end() && hic_matrix.find(key2)==hic_matrix.end() )
        {
            WHIC tmp_whic;
            tmp_whic.end1 = win1_end;
            tmp_whic.typ1 = win1_typ;
            tmp_whic.end2 = win2_end;
            tmp_whic.typ2 = win2_typ;
            tmp_whic.hcnt = hcnt;
            map<unsigned long, WHIC> tmp2_pos;
            tmp2_pos.insert(std::pair<unsigned long, WHIC>(win2_sta, tmp_whic) );
            map<unsigned long, map<unsigned long, WHIC> > tmp1_pos;
            tmp1_pos.insert(std::pair<unsigned long, map<unsigned long, WHIC> >(win1_sta, tmp2_pos));
            hic_matrix.insert(std::pair<string, map<unsigned long, map<unsigned long, WHIC> > >(key1, tmp1_pos));
        }else
        if( hic_matrix.find(key1)!=hic_matrix.end() )
        {
            if( hic_matrix[key1].find(win1_sta) == hic_matrix[key1].end() )
            {
                WHIC tmp_whic;
                tmp_whic.end1 = win1_end;
                tmp_whic.typ1 = win1_typ;
                tmp_whic.end2 = win2_end;
                tmp_whic.typ2 = win2_typ;
                tmp_whic.hcnt = hcnt;  
                map<unsigned long, WHIC> tmp2_pos;
                tmp2_pos.insert(std::pair<unsigned long, WHIC>(win2_sta, tmp_whic) );
                hic_matrix[key1].insert(std::pair<unsigned long, map<unsigned long, WHIC> >(win1_sta, tmp2_pos));
            }else
            {
                if( hic_matrix[key1][win1_sta].find(win2_sta) == hic_matrix[key1][win1_sta].end() )
                {
                    WHIC tmp_whic;
                    tmp_whic.end1 = win1_end;
                    tmp_whic.typ1 = win1_typ;
                    tmp_whic.end2 = win2_end;
                    tmp_whic.typ2 = win2_typ;
                    tmp_whic.hcnt = hcnt;                  
                    hic_matrix[key1][win1_sta].insert(std::pair<unsigned long, WHIC>(win2_sta, tmp_whic) );                
                }else
                {
                    cout << "   waring: line repeated: " << line << endl;
                }
            }
        }else
        if( hic_matrix.find(key2)!=hic_matrix.end() )
        {
            /* swap ctg2-ctg1 as ctg1-ctg2 */
            string tmp_ctg = ctg1;
            ctg2           = ctg1;
            ctg1           = tmp_ctg;
            unsigned long tmp_sta = win1_sta;
            win1_sta              = win2_sta;
            win2_sta              = tmp_sta;
            //
            if( hic_matrix[key2].find(win1_sta) == hic_matrix[key2].end() )
            {
                WHIC tmp_whic;
                tmp_whic.end1 = win2_end; // swapped
                tmp_whic.typ1 = win2_typ; // swapped
                tmp_whic.end2 = win1_end; // swapped
                tmp_whic.typ2 = win1_typ; // swapped
                tmp_whic.hcnt = hcnt;  
                map<unsigned long, WHIC> tmp2_pos;
                tmp2_pos.insert(std::pair<unsigned long, WHIC>(win2_sta, tmp_whic) );
                hic_matrix[key2].insert(std::pair<unsigned long, map<unsigned long, WHIC> >(win1_sta, tmp2_pos));
            }else
            {
                if( hic_matrix[key2][win1_sta].find(win2_sta) == hic_matrix[key2][win1_sta].end() )
                {
                    WHIC tmp_whic;
                    tmp_whic.end1 = win2_end; // swapped
                    tmp_whic.typ1 = win2_typ; // swapped
                    tmp_whic.end2 = win1_end; // swapped 
                    tmp_whic.typ2 = win1_typ; // swapped
                    tmp_whic.hcnt = hcnt;                  
                    hic_matrix[key2][win1_sta].insert(std::pair<unsigned long, WHIC>(win2_sta, tmp_whic) );                
                }else
                {
                    cout << "   waring: line repeated: " << line << endl;
                }
            }
        } else;
    }
    //
    ifp.close();
    //
    if(0)
    {
        string tmpfolder = "./";
        write_hic_matrix(tmpfolder);
    }
    //
    return true;
}
//
//
bool write_hic_matrix(string tmpfolder)
{
    /* Function: output hic contact among window markers */
    string out_file = "zadd_hic_matrix_to_check.txt\0"; // 
    ofstream ofp;
    ofp.open(out_file.c_str(), ios::out);
    if(!ofp.good())
    {
        cout   << "   Error: cannot open file " << out_file << endl;
        return false;
    }
    //
    map<string, map<unsigned long, map<unsigned long, WHIC> > >::iterator kitr;
    map<string, map<unsigned long, map<unsigned long, WHIC> > >::iterator kitr_end;
    kitr     = hic_matrix.begin();
    kitr_end = hic_matrix.end();
    while(kitr != kitr_end)
    {
        string this_key = (*kitr).first;
        vector<string> keyinfo = split_string(this_key, '-');
        string ctg1 = keyinfo[0];
        string ctg2 = keyinfo[1];
        map<unsigned long, map<unsigned long, WHIC> > tmp1_pos = (*kitr).second;
        //
        map<unsigned long, map<unsigned long, WHIC> >::iterator p1_itr;
        map<unsigned long, map<unsigned long, WHIC> >::iterator p1_itr_end;
        p1_itr     = tmp1_pos.begin();
        p1_itr_end = tmp1_pos.end();
        while(p1_itr != p1_itr_end)
        {
            unsigned long sta1                = (*p1_itr).first;
            map<unsigned long, WHIC> tmp2_pos = (*p1_itr).second;
            //
            map<unsigned long, WHIC>::iterator p2_itr;
            map<unsigned long, WHIC>::iterator p2_itr_end;
            p2_itr     = tmp2_pos.begin();
            p2_itr_end = tmp2_pos.end();
            while(p2_itr != p2_itr_end)
            {
                unsigned long sta2 = (*p2_itr).first;
                WHIC tmp_whic      = (*p2_itr).second;
                unsigned long ctg1_win_size = tmp_whic.end1 - sta1 + 1;
                unsigned long ctg2_win_size = tmp_whic.end2 - sta2 + 1;
                double normalized_reads = tmp_whic.hcnt;                
                ofp << ctg1             << "\t" 
                    << sta1             << "\t" 
                    << tmp_whic.end1    << "\t" 
                    << tmp_whic.typ1    << "\t"
                    << ctg2             << "\t"
                    << sta2             << "\t" 
                    << tmp_whic.end2    << "\t"
                    << tmp_whic.typ2    << "\t"
                    << normalized_reads << endl;                  
                //
                p2_itr ++;
            }
            //
            p1_itr ++;
        }
        //
        kitr ++;
    }
    cout << "   Info: hic matrix wrote out to file " << out_file << endl;
    //
    ofp.close();
    return true;
}
//
bool read_hb_cluster(string                    hb_group_file, 
                  map<string, int>*            ctg2_hb_group,
                  map<int, map<string, int> >* hb_group_2ctg)
{
    /*
        ctg2_hb_group........: <ctg_id, old_group_id>
        hb_group_2ctg........: <old_group_id, <ctg_id, new_group_id> >
    */
    ifstream ifp;
    ifp.open(hb_group_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open gb group file " << hb_group_file << endl;
        return false;
    }
    int hb_grpid_new = -1;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        // "\tutg000001l_pilon -- utg000070l_pilon [color=gold, penwidth=1, arrowsize=1, label=0.969478]; /* cluster 0 */"
        // cout << "   check: " << line << endl;
        if(line.find(" -- ") != std::string::npos)
        {
            vector<string> lineinfo2 = split_string(line, '\t');
            vector<string> lineinfo  = split_string(lineinfo2[0], ' ');            
            string vet1  = lineinfo[0];
            string vet2  = lineinfo[2];            
            int hb_grpid_old = atoi(lineinfo[9].c_str()); // old_group_id   
            //cout << "   check: vet1=" << vet1 << " -- vet2=" << vet2 << ", gb group = " << hb_grpid_old << endl;
            if(!add_vet_to_group(vet1, 
                                 hb_grpid_old, 
                                 hb_grpid_new, 
                                 hb_group_2ctg, 
                                 ctg2_hb_group) )
            {
                return false;
            }
            if(!add_vet_to_group(vet2, 
                                 hb_grpid_old, 
                                 hb_grpid_new, 
                                 hb_group_2ctg, 
                                 ctg2_hb_group) )
            {
                return false;
            }            
        }        
    }
    ifp.close();
    //
    if(1)
    {
        map<int, map<string, int> >::iterator gitr;
        map<int, map<string, int> >::iterator gitr_end;
        gitr     = (*hb_group_2ctg).begin();
        gitr_end = (*hb_group_2ctg).end();
        while(gitr != gitr_end)
        {
            int this_old_group_id = (*gitr).first;
            cout << "   check: old group id = " << this_old_group_id << " <= new group id = ";
            map<string, int> this_group = (*gitr).second;
            map<string, int>::iterator citr;
            map<string, int>::iterator citr_end;
            citr     = this_group.begin();
            citr_end = this_group.end();
            bool first_ctg = true;
            while(citr != citr_end)
            {
                string this_ctg = (*citr).first;
                int    this_new_group_id = (*citr).second;
                if(first_ctg)
                {
                    cout << this_new_group_id << endl;
                    first_ctg = false;
                }
                cout << "        : " << this_ctg << endl;
                //
                citr ++;
            }
            //
            gitr ++;
        }
    }
    //
    return true;
}
//
bool add_vet_to_group(string                       vet, 
                      int                          hb_grpid_old,
                      int                          hb_grpid_new,
                      map<int, map<string, int> >* hb_group_2ctg,
                      map<string, int>*            ctg2_hb_group)
{
    /*
        hb_grpid_old.........: old group id: any number
        hb_grpid_new.........: new group id: 1-4
        hb_group_2ctg........: <old_group_id, <ctg_id, new_group_id> >
    */    
    // update map from linkage group to list of contigs 
    if( (*hb_group_2ctg).find(hb_grpid_old) == (*hb_group_2ctg).end() )
    {
        map<string, int> tmpgroup;
        tmpgroup.insert(std::pair<string, int>(vet, hb_grpid_new));
        (*hb_group_2ctg).insert(std::pair<int, map<string, int> >(hb_grpid_old, tmpgroup));
    }else
    {
        if( (*hb_group_2ctg)[hb_grpid_old].find(vet) == (*hb_group_2ctg)[hb_grpid_old].end() )
        {        
            (*hb_group_2ctg)[hb_grpid_old].insert(std::pair<string, int>(vet, hb_grpid_new));
        }else
        {
            ; // redundant 
        }
    }
    // update map from a contig to its linkage group
    if((*ctg2_hb_group).find(vet) == (*ctg2_hb_group).end())
    {
        (*ctg2_hb_group).insert(std::pair<string, int>(vet, hb_grpid_old));
    }else
    {
        // should never happen
        if( (*ctg2_hb_group)[vet] != hb_grpid_old )
        {
            cout << "   Warning: contig found in group " << (*ctg2_hb_group)[vet] << " and " << hb_grpid_old << endl;
            cout << "            if it is non-haplotig markers, it is fine; otherwise, something wrong!"     << endl;
        }
    }
    //      
    return true;
}
//
bool read_win_marker_ctg(string win_marker_file, 
                         map<string, unsigned long> target_contig_size, 
                         map<string, map<unsigned long, NODE> >* win_marker)
{
    /*
        target_contig_size..... <target_ctg_id, size>
        win_marker............. <ctg, <win_sta, {win_end, hap/dip/trip/tetrap/rep, size_ctg, size_window, ...}> >
        NODE = {
		    unsigned long             end; // window end
		    string                   type; // hap/dip/trip/tetrap/rep
		    unsigned long        win_size; // size of this marker
		    unsigned long        ctg_size; // size of the contig having this marker
   	        }
   	Note, here only contigs in given list of target_contig_size will be collected
    */
    ifstream ifp;
    ifp.open(win_marker_file.c_str(), ios::in);
    if(!ifp.good())
    {
        return false;
    }
    int num_win       = 0; // given group of contigs
    int num_win_total = 0; // full assembly
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        /*
            utg000027lc	1	50000	174.013	156	1551343	hap
	    utg000027lc	50001	100000	174.013	156	1551343	hap
	    utg000027lc	100001	150000	174.013	156	1551343	hap
        */
        num_win_total ++;
        vector<string> lineinfo = split_string(line, '\t');
        string this_ctg         = lineinfo[0];
        if(target_contig_size.find(this_ctg) != target_contig_size.end())
        {
            num_win ++;
            // only collect info in given contig list.
            unsigned long this_win_sta  = strtoul(lineinfo[1].c_str(), NULL, 0);
            unsigned long this_win_end  = strtoul(lineinfo[2].c_str(), NULL, 0);
            unsigned long this_ctg_size = strtoul(lineinfo[5].c_str(), NULL, 0);
            string this_type            = lineinfo[6];
            NODE tmp_node;
            tmp_node.end      = this_win_end;
            tmp_node.type     = this_type;
            tmp_node.win_size = this_win_end - this_win_sta + 1;
            tmp_node.ctg_size = this_ctg_size;
            tmp_node.grp_link_val.clear();
            tmp_node.grp_link_ctg.clear();
            if( (*win_marker).find(this_ctg) == (*win_marker).end() )
            {
                map<unsigned long, NODE> tmp_win;
                tmp_win.insert(std::pair<unsigned long, NODE>(this_win_sta, tmp_node) );
                (*win_marker).insert(std::pair<string, map<unsigned long, NODE> >(this_ctg, tmp_win));
            }else
            {
                (*win_marker)[this_ctg].insert(std::pair<unsigned long, NODE>(this_win_sta, tmp_node));
            }
        }
    }
    ifp.close();
    //
    cout << "   Info: " << target_contig_size.size() << " target ctgs given, "          << endl;
    cout << "       : " << num_win_total             << " windows from full assembly, " << endl;
    cout << "       : " << (*win_marker).size()      << " target ctgs collected with "  << num_win 
         << " windows (note, windows not belonging to given group were skipped). "      << endl;
    //
    return true;
}
//
//
bool read_target_ctg(string ctg_file, map<string, unsigned long>* target_contig_size)
{
    /* read target contigs and their size */
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
        if( (*target_contig_size).find(this_ctg_id) != (*target_contig_size).end() )
        {
            if( (*target_contig_size)[this_ctg_id] != this_ctg_size)
            {
                cout << "   Warning: different sizes with the same contig id, please check: " << this_ctg_id << endl;
            }else
            {
                cout << "   Warning: repeated contig id: " << this_ctg_id << endl;
            }
        }else
        {
            (*target_contig_size).insert(std::pair<string, unsigned long>(this_ctg_id, this_ctg_size));
        }
    }
    //
    return true;
}
