/*
    Given 
    
        two dot files with contigs as vertex, and their hic contact as edges,
        
    compare the two files, and check how each group in file 1 overlaps groups in file 2, regarding contigs.
    
    2022-10-28 Started: Hequan Sun, MPIPZ
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
bool get_dot_contig(string dotfile, map<string, map<string, int> >* group_ctg, map<string, string>* group_ids);
bool compare_group(map<string, map<string, int> > ref_group_ctg, 
                   map<string, string>            ref_group_ids,
                   map<string, map<string, int> > qry_group_ctg, 
                   map<string, string>            qry_group_ids);
//
int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        cout << "\n   Function: given two dot files, with contigs as vertex and their hic contact as edges, " << endl
             << "               check how group in query.dot overlaps groups in ref.dot, regarding contigs."  << endl;
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "dot_comparator ref.dot query.dot"
             << endl << endl;
        return 1;
    }
    double startT= clock();
    // read ref group:contigs
    string ref_dotfile = (string)argv[1];
    map<string, map<string, int> > ref_group_ctg;
    map<string, string> ref_group_ids;
    if(!get_dot_contig(ref_dotfile, &ref_group_ctg, &ref_group_ids))
    {
        cout << "   Error: failed in getting ref group contigs. " << endl;
        return 1;
    }
    // read query group:contigs
    string qry_dotfile = (string)argv[2];
    map<string, map<string, int> > qry_group_ctg;
    map<string, string> qry_group_ids;
    if(!get_dot_contig(qry_dotfile, &qry_group_ctg, &qry_group_ids))
    {
        cout << "   Error: failed in getting query group contigs. " << endl;
        return 1;
    }    
    // compare contigs in query group with reference group 
    if(!compare_group(ref_group_ctg, ref_group_ids,
                      qry_group_ctg, qry_group_ids))
    {
        cout << "   Error: failed in compare contigs of query and reference groups. " << endl;
        return 1;
    }
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;
}
//
bool compare_group(map<string, map<string, int> > ref_group_ctg, 
                   map<string, string>            ref_group_ids,
                   map<string, map<string, int> > qry_group_ctg, 
                   map<string, string>            qry_group_ids)
{
    if(ref_group_ctg.size()==0 || qry_group_ctg.size()==0)
    {
        cout << "   Error: empty group info. " << endl;
        return false;
    }
    map<string, map<string, int> >::iterator qry_grp_itr;
    map<string, map<string, int> >::iterator qry_grp_itr_end;
    qry_grp_itr     = qry_group_ctg.begin();
    qry_grp_itr_end = qry_group_ctg.end(); 
    while(qry_grp_itr != qry_grp_itr_end)
    {
        // current query group info 
        string this_qry_group_id        = (*qry_grp_itr).first;
        assert(qry_group_ids.find(this_qry_group_id) != qry_group_ids.end() );
        string this_qry_group_id_given  = qry_group_ids[ this_qry_group_id ];
        map<string, int> this_qry_group = (*qry_grp_itr).second;
        cout << "   Info: this_qry_group " << this_qry_group_id 
             << "="                        << this_qry_group_id_given 
             << " with "                   << this_qry_group.size() 
             << " contigs: "               << endl; 
        //
        map<string, map<string, int> >::iterator ref_grp_itr;
        map<string, map<string, int> >::iterator ref_grp_itr_end;
        ref_grp_itr     = ref_group_ctg.begin();
        ref_grp_itr_end = ref_group_ctg.end();
        while(ref_grp_itr != ref_grp_itr_end)
        {
            // current reference group info
            string this_ref_group_id        = (*ref_grp_itr).first;
            assert(ref_group_ids.find(this_ref_group_id) != ref_group_ids.end() );
            string this_ref_group_id_given  = ref_group_ids[ this_ref_group_id ];
            map<string, int> this_ref_group = (*ref_grp_itr).second;
            // compare contigs of query group with those in reference group
            int    qry_overlap_ref = 0;
            string qry_overlap_ref_ctgs = "";            
            map<string, int>::iterator qry_ctg_itr;
            map<string, int>::iterator qry_ctg_itr_end;
            qry_ctg_itr     = this_qry_group.begin();
            qry_ctg_itr_end = this_qry_group.end();
            while(qry_ctg_itr != qry_ctg_itr_end)
            {
                string this_qry_ctg = (*qry_ctg_itr).first;
                if( this_ref_group.find(this_qry_ctg) != this_ref_group.end() )
                {
                    qry_overlap_ref ++;
                    qry_overlap_ref_ctgs = qry_overlap_ref_ctgs + "\n              " + this_qry_ctg;
                }
                //
                qry_ctg_itr ++;
            }
            if(qry_overlap_ref > 0)
            cout << "         : overlapping ref group " << this_ref_group_id 
                 << "="                                 << this_ref_group_id_given
                 << " by "                              << qry_overlap_ref
                 << " contigs: "                        << qry_overlap_ref_ctgs
                 << endl;
            //
            ref_grp_itr ++;
        }
        //
        qry_grp_itr ++;
    } 
    //
    return true;
}                   
//
bool get_dot_contig(string dotfile, map<string, map<string, int> >* group_ctg, map<string, string>* group_ids)
{
    ifstream ifp;
    ifp.open(dotfile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: failed in openning dot file " << dotfile << endl;
        return false;
    }
    string group_id;
    string ctg1_id;
    string ctg2_id;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        if(line.find(" -- ")==std::string::npos && line.find("label=\"G")==std::string::npos) continue;        
        vector<string> lineinfo = split_string(line, '\t');  
        if(lineinfo[0].find(" -- ") != std::string::npos)
        {
	    /*
      	    \tutg000055l_pilon -- utg000343l_pilon [color=gold, fontcolor=gold, penwidth=1, label=487.1]; \/* cluster 2 *\/
	    */
	    vector<string> ctginfo = split_string(lineinfo[0], ' ');
	    group_id = ctginfo[8] + " " + ctginfo[9];
	    ctg1_id  = ctginfo[0];
	    ctg2_id  = ctginfo[2];
	    if((*group_ctg).find(group_id) == (*group_ctg).end() )
	    {
	        map<string, int> tmp_ctgs;
	        tmp_ctgs.insert(std::pair<string, int>(ctg1_id, 1));
	        tmp_ctgs.insert(std::pair<string, int>(ctg2_id, 1));
	        (*group_ctg).insert(std::pair<string, map<string, int> >(group_id, tmp_ctgs));
	    }else
	    {
	        if((*group_ctg)[group_id].find(ctg1_id) == (*group_ctg)[group_id].end())
	        {
	            (*group_ctg)[group_id].insert(std::pair<string, int>(ctg1_id, 1));
	        }
	        if((*group_ctg)[group_id].find(ctg2_id) == (*group_ctg)[group_id].end())
	        {
	            (*group_ctg)[group_id].insert(std::pair<string, int>(ctg2_id, 1));
	        }	        
	    }
	    
        }else
        if(lineinfo[0].find("label=\"G") != std::string::npos)
        {
   	    /* "\tlabel="G_2" ;"  */  
   	    vector<string> groupidinfo = split_string(lineinfo[0], '=');
   	    string given_grp_id = groupidinfo[1].substr(0, groupidinfo[1].size()-1);
   	    (*group_ids).insert(std::pair<string, string>(group_id, given_grp_id));
   	    cout << "   Info: group "  << given_grp_id 
   	         << " from "           <<  group_id 
   	         << " collected with " << (*group_ctg)[group_id].size() 
   	         << " contigs. "       << endl;
        }else ;
    }
    ifp.close();
    return true;
}























