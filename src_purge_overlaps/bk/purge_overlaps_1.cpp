/* this tool  

       given a .gfa file, is aimed to 
         purge overlaps between overlapping unitigs.
         
   Written by Hequan Sun, MPIPZ
   Email: sunhequan@gmail.com/sun@mpipz.mpg.de
   Date.: 20220803
   
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
struct NODE
{
    string seq;
    unsigned long len;  
    map<string, unsigned long> out; // <out_ctg_id+direction, overlap_length>  
    map<string, unsigned long> in;  // <in_ctg_id+direction, overlap_length> 
    unsigned long pleft;            // left  clipping position given '+unitig'
    unsigned long pright;           // right clipping position given '+unitig'
    bool rm_ctg;                    // check window-based coverage, if all coverage<cutoff, remove this contig and related edges/overlapping other contigs
};
//
unsigned long final_clipped_off_size = 0;
map<string, unsigned long> clip_all;        // <"ctgX_id+-/lr", overlap length>
map<string, unsigned long> equivalent_clip; // <"ctgX_id+-/lr", overlap length> and <"ctgY_id+-/lr", overlap length>: clip one and keep the other!
//
bool read_gfa(string gfa_file, map<string, NODE>* ctg_info); // <ctg_id, NODE>
bool read_cov(string cov_file, 
              double min_hap_cov, 
              map<string, map<unsigned long, double> >* ctg_cov, 
              map<string, bool>* ctg_low_cov);
bool prune_gfa(map<string, bool> ctg_low_cov, map<string, NODE>* ctg_info);
//
bool find_links(map<string, NODE> ctg_info);
bool find_cliping_site(map<string, unsigned long> this_ctg_relevant_overlaps);
bool remove_equivalent_clip(map<string, NODE> ctg_info, map<string, unsigned long>* clip_all);
//
int main(int argc, char* argv[])
{
    if(argc < 5)
    {
        printf("\nFunction: coverage-based .gfa sequence cleaning.\n");
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl << endl;
        printf("Usage: purge_overlaps unitigs.gfa cnv_winsize1000_step1000_hq.txt hap_cov min_hap_ratio\n\n");
        cout << "       Note, cnv_winsize1000_step1000_hq.txt is output of CNV_HQ_v3."    
             << endl
             << endl;
        exit(1);
    }
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    printf("\nUnitig end cleaning started on %s\n", asctime(tinfo));    
    // s0. get inputs
    string gfa_file    = argv[1];
    string cov_file    = argv[2];
    double hap_cov     = atof(argv[3]); // Expected cov. of haplotigs; 1/2,1/3,1/4,... of hap_cov meaning overlap of 2,3,4,...unitigs
    double min_hap_cov = hap_cov*atof(argv[4]);  // if all windows with coverage less than this value, whole contigs removed.
    cout << "   Info: contigs with coverage < " <<  min_hap_cov << " for all windows will be purged directly. " << endl;
    // s1. get initial graph
    cout << "   Info: reading graph from gfa file..." << endl;
    map<string, NODE> ctg_info; // <ctg_id, NODE>
    if(!read_gfa(gfa_file, &ctg_info))
    {
        cout << "   Error: reading gfa file failed. Check " << gfa_file << endl;
        return false;
    }
    cout << "   Info: reading graph done." << endl;    
    
    
    /* test */
    find_links(ctg_info);
    // return true;
    
    // s2. get coverage of contigs
    cout << "   Info: reading coverage info..." << endl;
    map<string, map<unsigned long, double> > ctg_cov;
    map<string, bool> ctg_low_cov;
    if(!read_cov(cov_file, min_hap_cov, &ctg_cov, &ctg_low_cov))
    {
        cout << "   Error: read coverage file failed. Check " << cov_file << endl;
        return false;
    }
    cout << "   Info: reading coverage done." << endl;     
    // s4. prune graph
    cout << "   Info: contigs in map (before prunning) " << ctg_info.size() << endl;    
    cout << "   Info: pruning initial graph with coverage info..." << endl;    
    if(!prune_gfa(ctg_low_cov, &ctg_info))
    {
        cout << "   Error: pruning assembly graph failed. " << endl;
        return false;
    }
    cout << "   Info: pruning done." << endl;        
    // s3. get potential clipping sites...
    cout << "   Info: get potential clipping sites along contigs..." << endl;
    unsigned long total_ctg_size   = 0;
    unsigned long total_ctg_cnt    = 0;    
    unsigned long clipped_off_size = 0;
    if(1)
    {        
        map<string, NODE>::iterator citr;
        map<string, NODE>::iterator citr_end;
        citr     = ctg_info.begin();
        citr_end = ctg_info.end();
        while(citr != citr_end)
        {
            string this_ctg = (*citr).first;
            NODE tmp_node   = (*citr).second;            
            cout << "      check-ovl: " << this_ctg 
                 << "\t"                << tmp_node.len 
                 << " bp, potential regions to cut: ";
            total_ctg_size += tmp_node.len;
            total_ctg_cnt  ++;
            //
            if(ctg_info[this_ctg].rm_ctg == true)
            {
                // full contig will be removed.
                cout << "\tfull[1-" << tmp_node.len << "]" << endl;
                clipped_off_size += tmp_node.len;
                //
                citr ++;
                continue;
            }
            //
            bool in_overlap = false;
            if(tmp_node.in.size()>0)
            {     
                cout << "\tleft(in)[1-"  << tmp_node.pleft  << "]";
                in_overlap = true;
            }else
            {
                cout << "\tleft(in)[0-"  << 0               << "]";
            }
            bool out_overlap = false;            
            if(tmp_node.out.size()>0)
            {
                cout << "\tright(out):[" << tmp_node.pright << "-" << tmp_node.len << "]";
                out_overlap = true;
            }else
            {
                cout << "\tright(out):[" << tmp_node.len    << "-" << tmp_node.len << "]";             
            }
            cout << endl;            
            if(!out_overlap && !in_overlap)
            {
                cout << "      no clipping on singleton"  << endl;
                citr ++;
                continue;                
            }
            //
            if(tmp_node.in.size()>=2 )
            {
                cout << "       no clipping on left-end"  << endl;
            }else
            if(tmp_node.in.size()==1 )
            {
                cout << "       yes clipping on left-end" << endl;
                clipped_off_size += tmp_node.pleft;
            }else
            {
                cout << "       no clipping on left-end as no overlapping" << endl;
            }
            if(tmp_node.out.size()>=2 )
            {
                cout << "       no clipping on right-end"  << endl;
            }else
            if(tmp_node.out.size()==1 )
            {
                cout << "       yes clipping on right-end" << endl;
                clipped_off_size += (tmp_node.len - tmp_node.pright + 1);                
            }else
            {
                cout << "       no clipping on right-end as no overlapping" << endl;
            }
            //
            map<string, unsigned long>::iterator oitr;
            map<string, unsigned long>::iterator oitr_end;
            oitr     = tmp_node.out.begin();
            oitr_end = tmp_node.out.end();
            while(oitr != oitr_end)
            {
                cout << "       right-end (out) overlapping " << (*oitr).first << " by " << (*oitr).second << endl;
                oitr ++;
            }
            oitr     = tmp_node.in.begin();
            oitr_end = tmp_node.in.end();
            while(oitr != oitr_end)
            {
                cout << "       left-end (in) overlapping " << (*oitr).first << " by " << (*oitr).second << endl;
                oitr ++;
            }         
            //
            citr ++;
        }
    }
    //
    final_clipped_off_size = 0;
    map<string, unsigned long>::iterator clitr     = clip_all.begin();
    map<string, unsigned long>::iterator clitr_end = clip_all.end();
    while(clitr != clitr_end)
    {
        cout << "   Final: clip " << (*clitr).first << " " << (*clitr).second << " bp." << endl;
        final_clipped_off_size += (*clitr).second;
        //
        clitr ++;
    } 
    //
    cout << "   Info: contigs in map (after prunning) "   << ctg_info.size() << endl;        
    cout << "   Info: total contig number........: "      << total_ctg_cnt                           << endl;
    cout << "       : total contig size observed.: "      << total_ctg_size                          << endl;
    cout << "       : total regional size clipped1: "     << clipped_off_size                        << endl;
    cout << "       : total informative size kept1: "     << total_ctg_size - clipped_off_size       << endl;    
    cout << "       : total regional size clipped2: "     << final_clipped_off_size                  << endl;    
    cout << "       : total informative size kept2: "     << total_ctg_size - final_clipped_off_size << endl;    
    cout << "   Info: get potential clipping sites done." << endl; 
    //
    time_t endftime;
    struct tm* endtinfo;
    time(&endftime);
    endtinfo = localtime(&endftime);
    printf("\nUnitig end cleaning successfully finished on %s\n", asctime(endtinfo));
    return 0;
}
// 
bool remove_equivalent_clip(map<string, NODE> ctg_info, map<string, unsigned long>* clip_all)
{
    /*
        given: L   utg000161l  +   utg001792l  -   21838M  L1:i:130099,
               clip utg001792l+r 21838 bp is equivalent to clip utg000161l+r 21838 bp
               we keep the long contig not changed!
        ctg_info.............<ctg_id, NODE>
        clip_all.............<ctg_id + "+/-" + "l/r", overlap_length>
    */ 
    
    map<string, NODE>::iterator citr;
    map<string, NODE>::iterator citr_end;
    citr     = ctg_info.begin();
    citr_end = ctg_info.end();
    while(citr != citr_end)
    {
        string this_ctg = (*citr).first;
        NODE  this_node = (*citr).second;
        
        //
        citr ++;
    }
    
    map<string, unsigned long>::iterator clitr;
    map<string, unsigned long>::iterator clitr_end;
    clitr     = (*clip_all).begin();
    clitr_end = (*clip_all).end();
    while(clitr != clitr_end)
    { 
        string key1 = (*clitr).first; // ctg_id + "+/-" + "l/r"        
        string dir1 = key1.substr( key1.size()-2, 1);        
        string end1 = key1.substr( key1.size()-1, 1);
        //
        string key2 = "";
        if(dir1.compare("+")==0 && end1.compare("l")==0 )
        {
            key2 = key1.substr(0, key1.size()-2) + "-r";
        }else
        if(dir1.compare("+")==0 && end1.compare("r")==0 )
        {
            key2 = key1.substr(0, key1.size()-2) + "-l";
        }else
        if(dir1.compare("-")==0 && end1.compare("l")==0 )
        {
            key2 = key1.substr(0, key1.size()-2) + "+r";
        }else        
        if(dir1.compare("-")==0 && end1.compare("r")==0 )
        {
            key2 = key1.substr(0, key1.size()-2) + "+l";
        }else;        
        //
        clitr ++;
    }
    //
    return true;
}
//
bool find_links(map<string, NODE> ctg_info)
{
    // given a node, find linked-nodes and overlaps    
    map<string, NODE>::iterator citr;
    map<string, NODE>::iterator citr_end;
    citr     = ctg_info.begin();
    citr_end = ctg_info.end();
    while(citr != citr_end)
    {
        string this_ctg = (*citr).first; // the contig to be removed
        cout << endl;
        cout << "check-cls: " << this_ctg << " as focus +" << endl;
        // this collects all related overlaps on the "+" side of this_ctg
        map<string, unsigned long> this_ctg_relevant_overlaps; // < ctg1_id + ">" + ctg1_dir + ">" + ctg2_id + ">" + ctg2_dir, length >
        // find contigs as children of this_ctg, i.e. "this_ctg+" -> "child_ctg+/-"
        map<string, int> connected_ctg; // <"ctg_id+/-", 1>
        //
        NODE tmp_node   = (*citr).second;
        map<string, unsigned long>::iterator oitr;
        map<string, unsigned long>::iterator oitr_end;
        oitr     = tmp_node.out.begin();
        oitr_end = tmp_node.out.end();
        while(oitr != oitr_end)
        {
            string ovl_ctg        = (*oitr).first; // "ctg_id+/-"
            unsigned long ovl_len = (*oitr).second;
            cout << "         : " << this_ctg << "+ >>>>>>>>>> " << ovl_ctg << " by " << ovl_len << " bp" << endl;
            // collect overlap info
            string this_overlap_key = this_ctg                            + ">" + 
                                                                      "+" + ">" + 
                                      ovl_ctg.substr(0, ovl_ctg.size()-1) + ">" + 
                                      ovl_ctg.substr(ovl_ctg.size()-1, 1);
            this_ctg_relevant_overlaps.insert(std::pair<string, unsigned long>(this_overlap_key, ovl_len));
            //
            if(connected_ctg.find(ovl_ctg) == connected_ctg.end())
            {
                connected_ctg.insert(std::pair<string, int>(ovl_ctg, 1));
            }
            //
            oitr ++;
        }
        //
        map<string, int> connected_ctg_extended;
        map<string, int>::iterator ccitr;
        map<string, int>::iterator ccitr_end;
        ccitr     = connected_ctg.begin();
        ccitr_end = connected_ctg.end();
        while(ccitr != ccitr_end)
        {
            string ovl_ctg     = (*ccitr).first; // "ctg_id+/-"
            string ovl_ctg_dir = ovl_ctg.substr(ovl_ctg.size()-1, 1); // "+" or "-"
            //
            if(ovl_ctg_dir.compare("+")==0)
            {
                // ovl_ctg as "+" => no change on direction of other contigs overlapping ovl_ctg
                map<string, NODE>::iterator citr_tmp2 = ctg_info.find(ovl_ctg.substr( 0, ovl_ctg.size()-1) );
                NODE tmp2_node   = (*citr_tmp2).second;
                //
                map<string, unsigned long>::iterator iitr;
                map<string, unsigned long>::iterator iitr_end;
                iitr     = tmp2_node.in.begin();
                iitr_end = tmp2_node.in.end();
                while(iitr != iitr_end)
                {
                    string ovl_ctg2        = (*iitr).first; // "ctg_id+/-"
                    unsigned long ovl_len2 = (*iitr).second;
                    //cout << "         : " << ovl_ctg  << " overlapping " << ovl_ctg2 << " by " << ovl_len2 << " bp" << endl;
                    cout << "         : " << ovl_ctg2 << " >>>>>>>>>> "  << ovl_ctg  << " by " << ovl_len2 << " bp" << endl;                    
                    // collect overlap info                    
                    string this_overlap_key1 = ovl_ctg.substr(0, ovl_ctg.size()-1)   + ">" + 
                                               ovl_ctg.substr(ovl_ctg.size()-1, 1)   + ">" + 
                                               ovl_ctg2.substr(0, ovl_ctg2.size()-1) + ">" + 
                                               ovl_ctg2.substr(ovl_ctg2.size()-1, 1);
                    string this_overlap_key2 = ovl_ctg2.substr(0, ovl_ctg2.size()-1) + ">" + 
                                               ovl_ctg2.substr(ovl_ctg2.size()-1, 1) + ">" +
                                               ovl_ctg.substr(0, ovl_ctg.size()-1)   + ">" + 
                                               ovl_ctg.substr(ovl_ctg.size()-1, 1);
                    if(this_ctg_relevant_overlaps.find( this_overlap_key1 ) == this_ctg_relevant_overlaps.end() &&
                       this_ctg_relevant_overlaps.find( this_overlap_key2 ) == this_ctg_relevant_overlaps.end()                    
                      )
                    {
                        this_ctg_relevant_overlaps.insert(std::pair<string, unsigned long>(this_overlap_key2, ovl_len2));
                    }
                    //
                    if(connected_ctg_extended.find(ovl_ctg2) == connected_ctg_extended.end())
                    {
                        connected_ctg_extended.insert(std::pair<string, int>(ovl_ctg2, 1));
                    }
                    iitr ++;
                }
            }else
            {
                // ovl_ctg as "-" => change direction of other contigs overlapping ovl_ctg
                map<string, NODE>::iterator citr_tmp2 = ctg_info.find(ovl_ctg.substr( 0, ovl_ctg.size()-1) );
                NODE tmp2_node   = (*citr_tmp2).second;
                //
                map<string, unsigned long>::iterator oitr;
                map<string, unsigned long>::iterator oitr_end;
                oitr     = tmp2_node.out.begin();
                oitr_end = tmp2_node.out.end();
                while(oitr != oitr_end)
                {
                    string ovl_ctg2        = (*oitr).first; // "ctg_id+/-"
                    string ovl_ctg2_dir    = ovl_ctg2.substr(ovl_ctg2.size()-1, 1); // "+" or "-"
                    if(ovl_ctg2_dir.compare("+")==0)
                    {
                        ovl_ctg2_dir = "-";
                    }else
                    {
                        ovl_ctg2_dir = "+";
                    }
                    // update direction
                    ovl_ctg2 = ovl_ctg2.substr(0, ovl_ctg2.size()-1) + ovl_ctg2_dir;
                    //
                    unsigned long ovl_len2 = (*oitr).second;
                    cout << "         : " << ovl_ctg << " <<<<<<<<<< " << ovl_ctg2 << " by " << ovl_len2 << " bp" << endl;
                    // collect overlap info                    
                    string this_overlap_key1 = ovl_ctg.substr(0, ovl_ctg.size()-1)   + ">" + 
                                               ovl_ctg.substr(ovl_ctg.size()-1, 1)   + ">" + 
                                               ovl_ctg2.substr(0, ovl_ctg2.size()-1) + ">" + 
                                               ovl_ctg2.substr(ovl_ctg2.size()-1, 1);
                    string this_overlap_key2 = ovl_ctg2.substr(0, ovl_ctg2.size()-1) + ">" + 
                                               ovl_ctg2.substr(ovl_ctg2.size()-1, 1) + ">" +
                                               ovl_ctg.substr(0, ovl_ctg.size()-1)   + ">" + 
                                               ovl_ctg.substr(ovl_ctg.size()-1, 1);
                    if(this_ctg_relevant_overlaps.find( this_overlap_key1 ) == this_ctg_relevant_overlaps.end() &&
                       this_ctg_relevant_overlaps.find( this_overlap_key2 ) == this_ctg_relevant_overlaps.end()                    
                      )
                    {
                        // need this_overlap_key2 here! always keep as "A+/-" -> "B+/-"
                        this_ctg_relevant_overlaps.insert(std::pair<string, unsigned long>(this_overlap_key2, ovl_len2));
                    }
                    //
                    if(connected_ctg_extended.find(ovl_ctg2) == connected_ctg_extended.end())
                    {
                        connected_ctg_extended.insert(std::pair<string, int>(ovl_ctg2, 1));
                    }
                    oitr ++;
                }        
            }
            //
            ccitr ++;
        }
        // overlaps related to this_ctg
        cout << "   unique: "  << endl;
        map<string, unsigned long>::iterator ovlitr;
        map<string, unsigned long>::iterator ovlitr_end;
        ovlitr     = this_ctg_relevant_overlaps.begin();
        ovlitr_end = this_ctg_relevant_overlaps.end();
        while(ovlitr != ovlitr_end)
        {
            cout << "         : " << (*ovlitr).first << " by " << (*ovlitr).second << " bp " << endl;
            //
            ovlitr ++;
        }
        // define clipping positions among these overlaps

        find_cliping_site(this_ctg_relevant_overlaps);








        //
        citr ++;
    }
    //
    return true;
}
bool find_cliping_site(map<string, unsigned long> this_ctg_relevant_overlaps)
{
    /*   STRATEGY:
                              |                      |
         -------------a-------p1------p2------>      |
                              |    <--|-------------p3-------b----------
              --------X-------|-------|------------->|
                              |-------|--------------|--------------------Y---->
                              |       |-----c--------|----->
                              |                      |
         s1: find the maximum overlap of length L between two contigs, supposing they are cgX and ctgY
         s2: for other contigs {a,b,c}, if overlapping with X or Y
                 if a overlap Y at p1 with given direction, clip a at [p1, a_size], or other direction: [1, p1] 
                 if b overlap X at p3 with given direction, clip b at [p3, b_size], or other direction: [1, p2]
                 if c overlap a at p2 with given direction, clip c at [1, L-(p2-p1+1)], or the other direction, [c_size - ( L - (p2-p1+1) ) + 1, c_size]
    */
    if(this_ctg_relevant_overlaps.size()==0)
    {
        return true;
    }
    cout << "   analysis of overlap: "  << endl;    
    // find max length
    string max_overlap_ids    = "";
    unsigned long max_overlap = 0;
    map<string, unsigned long>::iterator ovlitr;
    map<string, unsigned long>::iterator ovlitr_end;
    ovlitr     = this_ctg_relevant_overlaps.begin();
    ovlitr_end = this_ctg_relevant_overlaps.end();
    while(ovlitr != ovlitr_end)
    {
        if((*ovlitr).second > max_overlap)
        {
            max_overlap     = (*ovlitr).second;
            max_overlap_ids = (*ovlitr).first; // ctga_id + ">" + ctga_dir + ">" + ctgb_id + ">" + ctgb_dir
        }
        //
        ovlitr ++;
    }
    cout << "         : max overlaps as " << max_overlap << " bp; edge as " << max_overlap_ids << endl;
    vector<string> max_overlap_ctg_info = split_string(max_overlap_ids, '>');
    string ctgX_id  = max_overlap_ctg_info[0];
    string ctgX_dir = max_overlap_ctg_info[1];
    string ctgY_id  = max_overlap_ctg_info[2];
    string ctgY_dir = max_overlap_ctg_info[3];
    map<string, int> tmp_fixed_ctg;
    tmp_fixed_ctg.insert(std::pair<string, int>(ctgX_id, 1));
    tmp_fixed_ctg.insert(std::pair<string, int>(ctgY_id, 1));    
    // variable collecting clipping info
    map<string, unsigned long> final_clip_tmp; // <ctg_id + "+/-" + "l/r">, "ctg_id-r", meaning rc ctg and clipping at right end with given size
    // initialize with contig with max clip; just select the X
    final_clip_tmp.insert(std::pair<string, unsigned long>(ctgX_id+ctgX_dir+"r", max_overlap));
    // clip ctgX_id+ctgX_dir+"r" = clip ctgY_id+ctgY_dir+"l", that is, clip one of them is enough!
    if(equivalent_clip.find( ctgX_id+ctgX_dir+"r" ) == equivalent_clip.end())
    {
        equivalent_clip.insert(std::pair<string, unsigned long>(ctgX_id+ctgX_dir+"r", max_overlap));
    }else
    {
        if(equivalent_clip[ ctgX_id+ctgX_dir+"r" ] < max_overlap)
        {
            equivalent_clip[ ctgX_id+ctgX_dir+"r" ] = max_overlap;
        }
    } 
    cout << "   equivalent-map: " << ctgX_id+ctgX_dir+"r" << " updated" << endl;
    if(equivalent_clip.find( ctgY_id+ctgY_dir+"l" ) == equivalent_clip.end())
    {
        equivalent_clip.insert(std::pair<string, unsigned long>(ctgY_id+ctgY_dir+"l", max_overlap));
    }else
    {
        if(equivalent_clip[ ctgY_id+ctgY_dir+"l" ] < max_overlap)
        {
            equivalent_clip[ ctgY_id+ctgY_dir+"l" ] = max_overlap;
        }
    }
    cout << "   equivalent-map: " << ctgY_id+ctgY_dir+"l" << " updated" << endl;    
    if(ctgX_dir.compare("-") ==0)
    {
        if(equivalent_clip.find(ctgX_id+"+"+"l") == equivalent_clip.end())
        {
            equivalent_clip.insert(std::pair<string, unsigned long>(ctgX_id+"+"+"l", max_overlap));
        }else
        {
            if(equivalent_clip[ctgX_id+"+"+"l"] < max_overlap)
            {
                equivalent_clip[ctgX_id+"+"+"l"] = max_overlap;
            }
        }
        cout << "   equivalent-map: " << ctgX_id+"+"+"l" << " updated" << endl;    
    }else
    {
        if(equivalent_clip.find(ctgX_id+"-"+"l") == equivalent_clip.end())
        {
            equivalent_clip.insert(std::pair<string, unsigned long>(ctgX_id+"-"+"l", max_overlap));
        }else
        {
            if(equivalent_clip[ctgX_id+"-"+"l"] < max_overlap)
            {
                equivalent_clip[ctgX_id+"-"+"l"] = max_overlap;
            }
        }
        cout << "   equivalent-map: " << ctgX_id+"-"+"l" << " updated" << endl;                    
    }
    if(ctgY_dir.compare("-") ==0)
    {
        if(equivalent_clip.find(ctgY_id+"+"+"r") == equivalent_clip.end())
        {
            equivalent_clip.insert(std::pair<string, unsigned long>(ctgY_id+"+"+"r", max_overlap));
        }else
        {
            if(equivalent_clip[ctgY_id+"+"+"r"] < max_overlap)
            {
                equivalent_clip[ctgY_id+"+"+"r"] = max_overlap;
            }
        }
        cout << "   equivalent-map: " << ctgY_id+"+"+"r" << " updated" << endl;                            
    }else
    {
        if(equivalent_clip.find(ctgY_id+"-"+"r") == equivalent_clip.end())
        {
            equivalent_clip.insert(std::pair<string, unsigned long>(ctgY_id+"-"+"r", max_overlap));
        }else
        {
            if(equivalent_clip[ctgY_id+"-"+"r"] < max_overlap)
            {
                equivalent_clip[ctgY_id+"-"+"r"] = max_overlap;
            }
        } 
        cout << "   equivalent-map: " << ctgY_id+"-"+"r" << " updated" << endl;                                           
    }
    // clip those contigs overlapping any of max_overlap_ids but with smaller overlaps than max_overlap
    ovlitr     = this_ctg_relevant_overlaps.begin();
    ovlitr_end = this_ctg_relevant_overlaps.end();
    while(ovlitr != ovlitr_end)
    {
        if((*ovlitr).first.compare(max_overlap_ids) != 0)
        {
            //
            string this_overlap_ids              = (*ovlitr).first; // ctga_id + ">" + ctga_dir + ">" + ctgb_id + ">" + ctgb_dir
            unsigned long this_overlap_len       = (*ovlitr).second;
            vector<string> this_overlap_ctg_info = split_string(this_overlap_ids, '>');
            string ctga_id  = this_overlap_ctg_info[0];
            string ctga_dir = this_overlap_ctg_info[1];
            string ctgb_id  = this_overlap_ctg_info[2];
            string ctgb_dir = this_overlap_ctg_info[3];                
            //
            if(tmp_fixed_ctg.find(ctga_id)==tmp_fixed_ctg.end() &&
               tmp_fixed_ctg.find(ctgb_id)!=tmp_fixed_ctg.end()
              )
            {
                // clip ctga_id: if "-", rc and clip at right end; otherwise, clip at right end
                cout << "         : clip " << ctga_dir << ctga_id << " at the right end by " << this_overlap_len << endl;
                string tmp_key = ctga_id + ctga_dir + "r";
                if(final_clip_tmp.find(tmp_key) == final_clip_tmp.end() )
                {
                    final_clip_tmp.insert(std::pair<string, unsigned long>(tmp_key, this_overlap_len));
                }else
                {
                    if(final_clip_tmp[tmp_key] < this_overlap_len)
                    {
                        final_clip_tmp[tmp_key] = this_overlap_len;
                    }
                }
            }else
            if(tmp_fixed_ctg.find(ctga_id)!=tmp_fixed_ctg.end() &&
               tmp_fixed_ctg.find(ctgb_id)==tmp_fixed_ctg.end()
              )
            {
                // clip ctgb_id: if "-", rc and clip at left end; otherwise, clip at right end
                cout << "         : clip " << ctgb_dir << ctgb_id << " at the left end by " << this_overlap_len << endl;
                string tmp_key = ctgb_id + ctgb_dir + "l";
                if(final_clip_tmp.find(tmp_key) == final_clip_tmp.end() )
                {
                    final_clip_tmp.insert(std::pair<string, unsigned long>(tmp_key, this_overlap_len));
                }else
                {
                    if(final_clip_tmp[tmp_key] < this_overlap_len)
                    {
                        final_clip_tmp[tmp_key] = this_overlap_len;
                    }
                }                
            }else
            {
                cout << "         : need new rule for " << this_overlap_ids << endl;
                string tmp_keya  = ctga_id + ctga_dir + "r";
                string tmp_keyb  = ctgb_id + ctgb_dir + "l";
                string tmp_keya2 = ctga_id + ctga_dir + "l";
                string tmp_keyb2 = ctgb_id + ctgb_dir + "r";                
                if(final_clip_tmp.find(tmp_keya) != final_clip_tmp.end())
                {
                    if(final_clip_tmp[tmp_keya] < this_overlap_len)
                    {
                        final_clip_tmp[tmp_keya] = this_overlap_len;
                    }
                }else
                if(final_clip_tmp.find(tmp_keyb) != final_clip_tmp.end())
                {
                    if(final_clip_tmp[tmp_keyb] < this_overlap_len)
                    {
                        final_clip_tmp[tmp_keyb] = this_overlap_len;
                    }
                }else
                if(final_clip_tmp.find(tmp_keya2) != final_clip_tmp.end())
                {
                    if(final_clip_tmp[tmp_keya2] < this_overlap_len)
                    {
                        final_clip_tmp[tmp_keya2] = this_overlap_len;
                    }
                }else     
                if(final_clip_tmp.find(tmp_keyb2) != final_clip_tmp.end())
                {
                    if(final_clip_tmp[tmp_keyb2] < this_overlap_len)
                    {
                        final_clip_tmp[tmp_keyb2] = this_overlap_len;
                    }
                }else                           
                {
                    cout << "   Error: either " << tmp_keya  << " or " 
                                                << tmp_keyb  << " or "
                                                << tmp_keya2 << " or "
                                                << tmp_keyb2 << " expected. " << endl;
                    //return false;
                }
            }                
        }        
        //
        ovlitr ++;
    }   
    // check
    cout << "   final clipping: " << endl;    
    map<string, unsigned long>::iterator clitr     = final_clip_tmp.begin();
    map<string, unsigned long>::iterator clitr_end = final_clip_tmp.end();
    while(clitr != clitr_end)
    {
        string key1 = (*clitr).first; // ctg_id + "+/-" + "l/r"
        string dir1 = key1.substr( key1.size()-2, 1);        
        string end1 = key1.substr( key1.size()-1, 1);
        string key2 = "";
        if(dir1.compare("+")==0 && end1.compare("l")==0 )
        {
            key2 = key1.substr(0, key1.size()-2) + "-r";
        }else
        if(dir1.compare("+")==0 && end1.compare("r")==0 )
        {
            key2 = key1.substr(0, key1.size()-2) + "-l";
        }else
        if(dir1.compare("-")==0 && end1.compare("l")==0 )
        {
            key2 = key1.substr(0, key1.size()-2) + "+r";
        }else        
        if(dir1.compare("-")==0 && end1.compare("r")==0 )
        {
            key2 = key1.substr(0, key1.size()-2) + "+l";
        }else;
        //
        cout << "         : cc " << (*clitr).first << " by " << (*clitr).second << endl;
        if(equivalent_clip.find(key1) == equivalent_clip.end() &&
           equivalent_clip.find(key2) == equivalent_clip.end() )
        {            
            if(clip_all.find( key1 ) == clip_all.end() && 
               clip_all.find( key2 ) == clip_all.end() 
              )
            {
                clip_all.insert(std::pair<string, unsigned long>( key1, (*clitr).second));
            }else
            {
                if(clip_all.find( key1 ) != clip_all.end())
                {
                    clip_all[key1] = (*clitr).second;
                }else
                {
                    clip_all[key2] = (*clitr).second;
                }
            }
        }
        //
        clitr ++;
    }   
    //
    return true;
}
//
bool prune_gfa(map<string, bool> ctg_low_cov, map<string, NODE>* ctg_info)
{
    if((*ctg_info).size() == 0)
    {
        cout << "   Error: empty ctg info from gfa. " << endl;
        return false;
    }
    // prune raw graph and output a new graph
    map<string, NODE>::iterator citr;
    map<string, NODE>::iterator citr_end;
    citr     = (*ctg_info).begin();
    citr_end = (*ctg_info).end();
    while(citr != citr_end)
    {
        string this_ctg = (*citr).first; // the contig to be removed
        NODE tmp_node   = (*citr).second;
        if(ctg_low_cov.find(this_ctg) != ctg_low_cov.end() && ctg_low_cov[this_ctg]==false)
        {
            cout << "      check-pru: " << this_ctg 
                 << "\t"                << tmp_node.len << " bp to be removed from graph, "
                 << "\toverlapping "    << endl;            
            // remove all overlapping information related to this contig
            //
            map<string, unsigned long>::iterator iitr;
            map<string, unsigned long>::iterator iitr_end;
            iitr     = tmp_node.in.begin();
            iitr_end = tmp_node.in.end();
            while(iitr != iitr_end)
            {
                string ovl_ctg = (*iitr).first;
                ovl_ctg        = ovl_ctg.substr(0, ovl_ctg.size()-1);
                //
                if((*ctg_info).find(ovl_ctg) != (*ctg_info).end()
                  )
                {
                    if( (*ctg_info)[ovl_ctg].in.find(this_ctg+"+") !=  (*ctg_info)[ovl_ctg].in.end() )
                    {
                        (*ctg_info)[ovl_ctg].in.erase(this_ctg+"+");
                    }else
                    if( (*ctg_info)[ovl_ctg].in.find(this_ctg+"-") !=  (*ctg_info)[ovl_ctg].in.end() )
                    {
                        (*ctg_info)[ovl_ctg].in.erase(this_ctg+"-");
                    }else
                    if( (*ctg_info)[ovl_ctg].out.find(this_ctg+"+") !=  (*ctg_info)[ovl_ctg].out.end() )
                    {
                        (*ctg_info)[ovl_ctg].out.erase(this_ctg+"+");
                    }else
                    if( (*ctg_info)[ovl_ctg].out.find(this_ctg+"-") !=  (*ctg_info)[ovl_ctg].out.end() )
                    {
                        (*ctg_info)[ovl_ctg].out.erase(this_ctg+"-");
                    }else ;
                    cout << "      " << this_ctg << " removed from overlapping list of " << ovl_ctg << endl;
                }else
                {
                    cout << "      " << ovl_ctg  << " not found in ctg_info." << endl;
                }
                //
                iitr ++;
            }   
            //
            map<string, unsigned long>::iterator oitr;
            map<string, unsigned long>::iterator oitr_end;
            oitr     = tmp_node.out.begin();
            oitr_end = tmp_node.out.end();
            while(oitr != oitr_end)
            {
                string ovl_ctg = (*oitr).first;
                ovl_ctg        = ovl_ctg.substr(0, ovl_ctg.size()-1);                
                //
                if((*ctg_info).find(ovl_ctg) != (*ctg_info).end()
                  )
                {
                    if( (*ctg_info)[ovl_ctg].in.find(this_ctg+"+") !=  (*ctg_info)[ovl_ctg].in.end() )
                    {
                        (*ctg_info)[ovl_ctg].in.erase(this_ctg+"+");
                    }else
                    if( (*ctg_info)[ovl_ctg].in.find(this_ctg+"-") !=  (*ctg_info)[ovl_ctg].in.end() )
                    {
                        (*ctg_info)[ovl_ctg].in.erase(this_ctg+"-");
                    }else
                    if( (*ctg_info)[ovl_ctg].out.find(this_ctg+"+") !=  (*ctg_info)[ovl_ctg].out.end() )
                    {
                        (*ctg_info)[ovl_ctg].out.erase(this_ctg+"+");
                    }else
                    if( (*ctg_info)[ovl_ctg].out.find(this_ctg+"-") !=  (*ctg_info)[ovl_ctg].out.end() )
                    {
                        (*ctg_info)[ovl_ctg].out.erase(this_ctg+"-");
                    }else ;
                    cout << "      " << this_ctg << " removed from overlapping list of " << ovl_ctg << endl;
                }
                else
                {
                    cout << "      " << ovl_ctg  << " not found in ctg_info." << endl;
                }                
                //
                oitr ++;
            }                      
            //
            // (*ctg_info).erase(citr++);
            (*ctg_info)[this_ctg].rm_ctg = true;
        }
        //
        citr ++;
    }
    //
    return true;
}
//
bool read_cov(string cov_file, 
              double min_hap_cov, 
              map<string, map<unsigned long, double> >* ctg_cov, 
              map<string, bool>* ctg_low_cov)
{
    ifstream ifp;
    ifp.open(cov_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open coverage file " << cov_file << endl;
        return false;
    }
    string             last_ctg = "";
    bool          last_ctg_keep = false;
    unsigned long  last_ctg_len = 0;
    //
    string       this_ctg;
    unsigned long win_sta;
    double        win_cov;
    unsigned long ctg_len;    
    int           win_n;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 7)
        {
            cout << "   warning: skipping insufficient line " << line << endl;
            continue;
        }
        this_ctg = lineinfo[0];
        win_sta  = strtoul(lineinfo[1].c_str(), NULL, 0);
        win_cov  = atof(lineinfo[5].c_str());
        ctg_len  = strtoul(lineinfo[6].c_str(), NULL, 0);
        //
        if(last_ctg.size()==0 && (*ctg_cov).find(this_ctg) == (*ctg_cov).end())
        {
            // first contig
            map<unsigned long, double> tmp_cov;
            tmp_cov.insert(std::pair<unsigned long, double>(win_sta, win_cov));
            (*ctg_cov).insert(std::pair<string, map<unsigned long, double> >(this_ctg, tmp_cov));
            //
            last_ctg      = this_ctg;
            last_ctg_keep = false;
            last_ctg_len  = ctg_len;
            win_n         = 1;
        }else
        if((*ctg_cov).find(this_ctg) == (*ctg_cov).end())
        {
            // additional contigs
            map<unsigned long, double> tmp_cov;
            tmp_cov.insert(std::pair<unsigned long, double>(win_sta, win_cov));
            (*ctg_cov).insert(std::pair<string, map<unsigned long, double> >(this_ctg, tmp_cov));
            //
            (*ctg_low_cov).insert(std::pair<string, bool>(last_ctg, last_ctg_keep));
            if(last_ctg_keep == true)
            {
                cout << "      check-c: "   << last_ctg 
                     << " cov-kept\t"       << last_ctg_len 
                     << "\tbp\t"            << win_n 
                     << "\twinds"           << endl;
            }else
            {
                cout << "      check-c: "   << last_ctg 
                     << " cov-removed\t"    << last_ctg_len 
                     << "\tbp\t"            << win_n 
                     << "\twinds"           << endl;
            }
            //
            last_ctg      = this_ctg;
            last_ctg_keep = false;
            last_ctg_len  = ctg_len;       
            win_n         = 1;     
        }        
        else
        {
            (*ctg_cov)[this_ctg].insert(std::pair<unsigned long, double>(win_sta, win_cov));
            win_n ++;
        }
        //
        if(win_cov >= min_hap_cov)
        {
            last_ctg_keep = true;
        }
    }
    // last contig info 
    if(1)
    {
        // additional contigs
        map<unsigned long, double> tmp_cov;
        tmp_cov.insert(std::pair<unsigned long, double>(win_sta, win_cov));
        (*ctg_cov).insert(std::pair<string, map<unsigned long, double> >(this_ctg, tmp_cov));
        //
        (*ctg_low_cov).insert(std::pair<string, bool>(last_ctg, last_ctg_keep));
        if(last_ctg_keep == true)
        {
            cout << "      check-c: " << last_ctg << " cov-kept\t" << last_ctg_len << "\tbp\t"<< win_n << "\twinds" << endl;
        }else
        {
            cout << "      check-c: " << last_ctg << " cov-removed\t" << last_ctg_len << "\tbp\t"<< win_n << "\twinds" << endl;
        }   
    }
    //
    ifp.close();
    //
    return true;
}
//
bool read_gfa(string gfa_file, map<string, NODE>* ctg_info)
{
    ifstream ifp;
    ifp.open(gfa_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << gfa_file << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        //
        if(lineinfo[0].compare("S")==0)
        {
            string this_ctg = "";
            string this_seq = ""; 
            this_ctg = lineinfo[1];
            this_seq = lineinfo[2];
            assert( (*ctg_info).find(this_ctg) == (*ctg_info).end() );
            //
            NODE tmp_node;
            tmp_node.seq    = this_seq;
            tmp_node.len    = this_seq.size();
            tmp_node.pleft  = 0;                 // left  clipping position given '+unitig'
            tmp_node.pright = this_seq.size()+1; // right clipping position given '+unitig'              
            tmp_node.rm_ctg = false;             // initialized as keep the contig
            //
            (*ctg_info).insert(std::pair<string, NODE>(this_ctg, tmp_node));
        }else
        if(lineinfo[0].compare("L")==0 && lineinfo[5].compare("*")!=0)
        {
            /*  utg019078l	33437 = 20065+13372 = 1852+31585 = 21501+11936 = 12764+20673
		L   utg019078l  +   utg022897l  +   20065M  L1:i:13372
		L   utg019078l  +   utg035008l  -   1852M   L1:i:31585
		L   utg019078l  -   utg026750l  +   21501M  L1:i:11936
		L   utg019078l  -   utg011655l  +   12764M  L1:i:20673	

		                     +++++++++++++++utg022897l++++++++++++++++>:43111bp
		                                          <---utg035008l--:16885bp                         
		           ++++++++++++utg019078l++++++++++>:33437bp
		       <---utg026750l------------:25249bp
	   	      <--utg011655l---:17285bp
            */
            if(lineinfo[2].compare("+") == 0)
            {
                string this_ctg = lineinfo[1];  
                assert( (*ctg_info).find(this_ctg) != (*ctg_info).end() );                          
                vector<string> pos_info = split_string(lineinfo[6], ':');
                unsigned long overlap_pos = strtoul(pos_info[2].c_str(), NULL, 0); // from this pos to right end of focal contig to clip off
                if(overlap_pos < (*ctg_info)[this_ctg].pright)
                {
                    (*ctg_info)[this_ctg].pright = overlap_pos;                    
                }
                //
                string overlap_ctg        = lineinfo[3]; // ctg_id
                string overlap_direction  = lineinfo[4]; // +-
                string overlap_len_str    = lineinfo[5].substr(0, lineinfo[5].size()-1);
                unsigned long overlap_len = strtoul(overlap_len_str.c_str(), NULL, 0);
                (*ctg_info)[this_ctg].out.insert(std::pair<string, unsigned long>(overlap_ctg+overlap_direction, overlap_len));
            }else
            {
                // when the focus contig is "-", reverse and complement the other overlapping contigs, and update overlapping pos and size
                string this_ctg = lineinfo[1];       
                assert( (*ctg_info).find(this_ctg) != (*ctg_info).end() );                     
                vector<string> pos_info = split_string(lineinfo[6], ':');
                unsigned long overlap_pos = strtoul(pos_info[2].c_str(), NULL, 0); // from this pos to left end of focal contig to clip off
                overlap_pos = (*ctg_info)[this_ctg].len - overlap_pos;             // always keep focus ctg as "+", and update related into correspondingly, 
                								   //    to clip [1, overlap_pos] (if 1-based system)
                if(overlap_pos > (*ctg_info)[this_ctg].pleft)
                {
                    (*ctg_info)[this_ctg].pleft = overlap_pos;
                }
                //
                string overlap_ctg        = lineinfo[3]; // ctg_id
                string overlap_direction  = lineinfo[4]; // +-
                if(overlap_direction.compare("+") == 0)
                {
                    overlap_direction = "-";// the current focal contig is "-" to "+", overlapping one as "-" to "+"
                }else
                {
                    overlap_direction = "+";// the current focal contig is "-" to "+", overlapping one as "+" to "-"
                }
                string overlap_len_str    = lineinfo[5].substr(0, lineinfo[5].size()-1);
                unsigned long overlap_len = strtoul(overlap_len_str.c_str(), NULL, 0);
                (*ctg_info)[this_ctg].in.insert(std::pair<string, unsigned long>(overlap_ctg+overlap_direction, overlap_len));
            }
        }else 
        {
            continue;
        }       
    }
    ifp.close();
    return true;
}
































