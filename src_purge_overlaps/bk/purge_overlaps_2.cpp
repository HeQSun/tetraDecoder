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
map<string, unsigned long> clip_all;            // <"ctgX_id+-/lr", overlap length>
map<string, unsigned long> equivalent_clip;     // <"ctgX_id+-/lr", overlap length> and <"ctgY_id+-/lr", overlap length>: clip one and keep the other!
map<string, unsigned long> max_overlap_ids_all; // <ctga_id + ">" + ctga_dir + ">" + ctgb_id + ">" + ctgb_dir, overlap length>
map<string, NODE> ctg_info;                     // <ctg_id, NODE>
unsigned long total_gfa_ctg_size = 0;
unsigned long total_gfa_ctg_cnt  = 0;
//
bool read_gfa(string gfa_file); // => ctg_info: <ctg_id, NODE>, total_gfa_ctg_size and total_gfa_ctg_cnt
bool read_cov(string cov_file, 
              double min_hap_cov, 
              map<string, map<unsigned long, double> >* ctg_cov, 
              map<string, bool>* ctg_low_cov);
bool find_links(); // handling map<string, NODE> ctg_info
bool find_cliping_site(map<string, unsigned long> this_ctg_relevant_overlaps);
bool remove_equivalent_clip(map<string, unsigned long> max_overlap_ids_all);
bool clip_output_contigs(string clipped_fa_file, long min_ctg_len); // handling map<string, NODE> ctg_info
//
int main(int argc, char* argv[])
{
    if(argc < 6)
    {
        printf("\nFunction: clean overlaps between unitigs in .gfa; TODO consider coverage for cleaning.\n");
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl << endl;
        printf("Usage: purge_overlaps unitigs.gfa min_ctg_len cnv_winsize1000_step1000_hq.txt hap_cov min_hap_ratio\n\n");        
        cout << "      min_ctg_len is minimum size of clipped contigs to output. " 
             << endl
             << "      cnv_winsize1000_step1000_hq.txt is output of CNV_HQ_v3." 
             << endl
             << "      hap_cov is expected coverage on haplotigs. "
             << endl
             << "      min_hap_ratio is used to determined min_hap_cov - below this value could be clipped."
             << endl
             << endl;
        cout << "Note, coverage info has not been used yet!!!" << endl;
        exit(1);
    }
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    printf("\nUnitig end cleaning started on %s\n", asctime(tinfo));    
    // s0. get inputs
    string gfa_file     = argv[1];
    long min_ctg_len    = strtol(argv[2], NULL, 0);
    string cov_file     = argv[3];
    double hap_cov      = atof(argv[4]);         // not applied yet 20220811
    double min_hap_cov  = hap_cov*atof(argv[5]); // not applied yet 20220811
    cout << "   Info: minimum length of clipped contigs to output: " << min_ctg_len << endl;
    cout << "   Info: contigs with coverage < "                      << min_hap_cov
         << " for all windows will be purged directly. "             << endl;
    // s1. get initial graph
    cout << "   Info: reading graph from gfa file..." << endl;
    if(!read_gfa(gfa_file))
    {
        cout << "   Error: reading gfa file failed. Check " << gfa_file << endl;
        return false;
    }
    cout << "   Info: reading graph done." << endl;
    //
    find_links();
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
    // s3. clean equivalent clips
    cout << "   Info: " << clip_all.size() << " clips before cleaning" << endl;
    remove_equivalent_clip(max_overlap_ids_all);
    cout << "   Info: " << clip_all.size() << " clips after  cleaning" << endl;    
    // s4. get final clipping info: ctg_info: <ctg_id, NODE>
    final_clipped_off_size = 0;
    map<string, unsigned long>::iterator clitr     = clip_all.begin();
    map<string, unsigned long>::iterator clitr_end = clip_all.end();
    while(clitr != clitr_end)
    {
        string key_tmp = (*clitr).first; // ctg_id+/-l/r
        string ctg_id  = key_tmp.substr(0, key_tmp.size()-2);
        string ctg_dir = key_tmp.substr(key_tmp.size()-2, 1); // +/-
        string ctg_end = key_tmp.substr(key_tmp.size()-1, 1); // l/r
        assert(ctg_info.find(ctg_id) != ctg_info.end() );
        // convert to "+" for all clipping (=>no need to reverse and complement the contig sequence when clipping)
        if(ctg_dir.compare("-")==0 && ctg_end.compare("r")==0)
        {
            ctg_info[ctg_id].pleft = (*clitr).second;         
            cout << "Final: clip " << ctg_id                 << "+l " << (*clitr).second 
                 << " bp = [1,"    << ctg_info[ctg_id].pleft << "]"   << endl; 
        }else
        if(ctg_dir.compare("-")==0 && ctg_end.compare("l")==0)
        {
            ctg_info[ctg_id].pright = ctg_info[ctg_id].len - (*clitr).second + 1;        
            cout << "Final: clip " << ctg_id          << "+r " << (*clitr).second 
                 << " bp = ["      << ctg_info[ctg_id].pright  << "," << ctg_info[ctg_id].len << "]" << endl;            
        }else       
        if(ctg_dir.compare("+")==0 && ctg_end.compare("l")==0)         
        {
            ctg_info[ctg_id].pleft = (*clitr).second;        
            cout << "Final: clip " << key_tmp                << " " << (*clitr).second 
                 << " bp = [1,"    << ctg_info[ctg_id].pleft << "]" << endl;
        }else
        if(ctg_dir.compare("+")==0 && ctg_end.compare("r")==0)         
        {
            ctg_info[ctg_id].pright = ctg_info[ctg_id].len - (*clitr).second + 1;        
            cout << "Final: clip " << key_tmp                 << " " << (*clitr).second 
                 << " bp = ["      << ctg_info[ctg_id].pright << "," << ctg_info[ctg_id].len << "]" << endl;
        }
        final_clipped_off_size += (*clitr).second;
        //
        clitr ++;
    } 
    //
    cout << "   Info: contigs in map "                    << ctg_info.size() << endl;        
    cout << "   Info: total contig number........: "      << total_gfa_ctg_cnt                           << endl;
    cout << "       : total contig size observed.: "      << total_gfa_ctg_size                          << endl;
    cout << "       : total regional size clipped2: "     << final_clipped_off_size                      << endl;    
    cout << "       : total informative size kept2: "     << total_gfa_ctg_size - final_clipped_off_size << endl;    
    cout << "       : contigs requiring no clipping are not calculated in clipped2 and kept2. "          << endl;
    cout << "   Info: get potential clipping sites done." << endl; 
    // s5. output clipped contigs
    string clipped_fa_file = "";
    vector<string> gfainfo = split_string(gfa_file, '/');
    clipped_fa_file = "clipped_" + gfainfo[gfainfo.size()-1] + ".fa";
    if(!clip_output_contigs(clipped_fa_file, min_ctg_len))
    {
        cout << "   Error: failed in output clipped sequences. " << endl;
        return -1;
    }
    //
    time_t endftime;
    struct tm* endtinfo;
    time(&endftime);
    endtinfo = localtime(&endftime);
    printf("\nUnitig end cleaning successfully finished on %s\n", asctime(endtinfo));
    return 0;
}
//
bool clip_output_contigs(string clipped_fa_file, long min_ctg_len)
{
    // map<string, NODE> ctg_info....<ctg_id, NODE>
    ofstream ofp;
    ofp.open(clipped_fa_file.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open file " << clipped_fa_file << endl;
        return false;
    }
    unsigned long collected_contig_size = 0;
    unsigned long collected_contig_cnt  = 0;
    unsigned long discarded_ctg_cnt     = 0;
    unsigned long discarded_ctg_size    = 0;
    unsigned long clipped_size          = 0;
    map<string, NODE>::iterator citr;
    map<string, NODE>::iterator citr_end;
    citr     = ctg_info.begin();
    citr_end = ctg_info.end();
    while(citr != citr_end)
    {
        string this_ctg = (*citr).first;
        NODE  this_node = (*citr).second;
        //
        long        pleft  = (long)this_node.pleft;
        long        pright = (long)this_node.pright;
        long this_ctg_size = (long)this_node.len;
        // skip case: ..........pl[--------keep--------]pr......... => (pr-1) - (pl+1) + 1
        if(pright-1-pleft < min_ctg_len)
        {
            //cout << "   discard "                                  << this_ctg 
            //     << " as ctg_size - pleft - (ctg_size-pright+1)="  << this_ctg_size - pleft - (this_ctg_size-pright+1)
            //     << " < given cutoff "                             << min_ctg_len
            //     << " bp"                                          << endl; 
            discarded_ctg_cnt ++;
            discarded_ctg_size += this_ctg_size;
            citr ++;
            continue;
        }
        //
        string clipped_seq = this_node.seq.substr(0, pright-1); // clip right
        clipped_seq        = clipped_seq.substr(pleft);         // clip left
        ofp << ">" << this_ctg << "c" << endl;
        ofp << clipped_seq << endl;
        collected_contig_size += clipped_seq.size();
        collected_contig_cnt ++;
        //
        clipped_size += this_node.seq.size()-clipped_seq.size();
        //
        citr ++;
    }
    ofp.close();
    //
    cout << "   Info: "                                    << discarded_ctg_cnt 
         << " contigs discarded fully with total size of " << discarded_ctg_size << " bp " << endl;
    cout << "   Info: paritally clipped size "             << clipped_size       << " bp " << endl;
    cout << "   Info: "                                    << collected_contig_cnt 
         << " contigs collected with total size of "       << collected_contig_size 
         << " bp."                                         << endl;
    //
    return true;
}
// 
bool remove_equivalent_clip(map<string, unsigned long> max_overlap_ids_all)
{
    /*
        given: L   utg000161l  +   utg001792l  -   21838M  L1:i:130099,
               clip utg001792l+r 21838 bp is equivalent to clip utg000161l+r 21838 bp
               we keep the long contig not changed!
        ctg_info.............<ctg_id, NODE>
        clip_all.............<ctg_id + "+/-" + "l/r", overlap length> (global)
        max_overlap_ids_all..<ctga_id + ">" + ctga_dir + ">" + ctgb_id + ">" + ctgb_dir, overlap length>
    */
    map<string, unsigned long>::iterator moitr;
    map<string, unsigned long>::iterator moitr_end;
    moitr     = max_overlap_ids_all.begin();
    moitr_end = max_overlap_ids_all.end();
    while(moitr != moitr_end)
    {
        string        max_overlap_ids = (*moitr).first;
        unsigned long max_overlap     = (*moitr).second;
        vector<string> max_overlap_ctg_info = split_string(max_overlap_ids, '>');
        string ctga_id  = max_overlap_ctg_info[0];
        string ctga_dir = max_overlap_ctg_info[1];
        string ctgb_id  = max_overlap_ctg_info[2];
        string ctgb_dir = max_overlap_ctg_info[3];        
        //
        string keya1 = ctga_id + ctga_dir + "r";
        string keyb1 = ctgb_id + ctgb_dir + "l";
        assert(ctg_info.find(ctga_id) != ctg_info.end() );
        assert(ctg_info.find(ctgb_id) != ctg_info.end() );
        unsigned long ctga_size = ctg_info[ctga_id].len;
        unsigned long ctgb_size = ctg_info[ctgb_id].len;        
        //
        string keya2 = "";
        string keyb2 = "";
        if(ctga_dir.compare("+")==0)
        {
            keya2 = ctga_id + "-" + "l";
        }else
        {
            keya2 = ctga_id + "+" + "l";
        }
        if(ctgb_dir.compare("+")==0)
        {
            keyb2 = ctgb_id + "-" + "r";
        }else
        {
            keyb2 = ctgb_id + "+" + "r";
        }
        //
        if(clip_all.find(keya1) != clip_all.end() &&
           clip_all.find(keyb1) != clip_all.end())
        {
            if(clip_all[keya1] == clip_all[keyb1])
            {
                if(ctga_size > ctgb_size)
                {
                    clip_all.erase(keyb1);
                    cout << "   clean: " << keyb1 << " removed from clip list as ovl " << keya1 << endl;
                }else
                {
                    clip_all.erase(keya1);
                    cout << "   clean: " << keya1 << " removed from clip list as ovl " << keyb1 << endl;                    
                }
            }
        }else
        if(clip_all.find(keya1) != clip_all.end() &&
           clip_all.find(keyb2) != clip_all.end())
        {
            if(clip_all[keya1] == clip_all[keyb2])
            {
                if(ctga_size > ctgb_size)
                {            
                    clip_all.erase(keyb2);
                    cout << "   clean: " << keyb2 << " removed from clip list as ovl " << keya1 << endl;
                }else
                {
                    clip_all.erase(keya1);
                    cout << "   clean: " << keya1 << " removed from clip list as ovl " << keyb2 << endl;                    
                }
            }        
        }else
        if(clip_all.find(keya2) != clip_all.end() &&
           clip_all.find(keyb1) != clip_all.end())
        {
            if(clip_all[keya2] == clip_all[keyb1])
            {
                if(ctga_size > ctgb_size)
                {            
                    clip_all.erase(keyb1);
                    cout << "   clean: " << keyb1 << " removed from clip list as ovl " << keya2 << endl;
                }else
                {
                    clip_all.erase(keya2);
                    cout << "   clean: " << keya2 << " removed from clip list as ovl " << keyb1 << endl;                    
                }
            }           
        }else
        if(clip_all.find(keya2) != clip_all.end() &&
           clip_all.find(keyb2) != clip_all.end())
        {
            if(clip_all[keya2] == clip_all[keyb2])
            {
                if(ctga_size > ctgb_size)
                {            
                    clip_all.erase(keyb2);
                    cout << "   clean: " << keyb2 << " removed from clip list as ovl " << keya2 << endl;
                }else
                {
                    clip_all.erase(keya2);
                    cout << "   clean: " << keya2 << " removed from clip list as ovl " << keyb2 << endl;                    
                }
            }            
        }else;        
        //
        moitr ++;
    }    
    //
    return true;
}
//
bool find_links()
{
    // map<string, NODE> ctg_info
    // given a node, find linked-nodes and overlaps    
    map<string, NODE>::iterator citr;
    map<string, NODE>::iterator citr_end;
    citr     = ctg_info.begin();
    citr_end = ctg_info.end();
    while(citr != citr_end)
    {
        string this_ctg = (*citr).first; // the contig to be removed
        cout << endl;
        // collect all overlaps on the "+" end of this_ctg
        // < ctg1_id + ">" + ctg1_dir + ">" + ctg2_id + ">" + ctg2_dir, length >
        map<string, unsigned long> this_ctg_relevant_overlaps; 
        // find contigs as children of this_ctg, i.e. "this_ctg+" -> "child_ctg+/-"
        map<string, int> connected_ctg; // <"ctg_id+/-", 1>
        // check 1: "+r" end
        cout << "check-cls: " << this_ctg << " as focus +, as start node" << endl;        
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
        //
        cout << "   unique-1: "  << endl;
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
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        cout << "check-cls: " << this_ctg << " as focus +, as stop-node..." << endl;        
        connected_ctg.clear();
        this_ctg_relevant_overlaps.clear();
        //
        oitr     = tmp_node.in.begin();
        oitr_end = tmp_node.in.end();
        while(oitr != oitr_end)
        {
            string ovl_ctg        = (*oitr).first; // "ctg_id+/-"
            unsigned long ovl_len = (*oitr).second;
            cout << "         : " << ovl_ctg << " >>>>>>>>>> " << this_ctg << "+ by " << ovl_len << " bp" << endl;
            // collect overlap info
            string this_overlap_key = ovl_ctg.substr(0, ovl_ctg.size()-1) + ">" + 
                                      ovl_ctg.substr(ovl_ctg.size()-1, 1) + ">" +
                                      this_ctg                            + ">" + "+";
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
                map<string, unsigned long>::iterator oitr;
                map<string, unsigned long>::iterator oitr_end;
                oitr     = tmp2_node.out.begin();
                oitr_end = tmp2_node.out.end();
                while(oitr != oitr_end)
                {
                    string ovl_ctg2        = (*oitr).first; // "ctg_id+/-"
                    unsigned long ovl_len2 = (*oitr).second;
                    //cout << "         : " << ovl_ctg  << " overlapping " << ovl_ctg2 << " by " << ovl_len2 << " bp" << endl;
                    cout << "         : " << ovl_ctg << " >>>>>>>>>> "  << ovl_ctg2  << " by " << ovl_len2 << " bp" << endl;                    
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
                        this_ctg_relevant_overlaps.insert(std::pair<string, unsigned long>(this_overlap_key1, ovl_len2));
                    }
                    //
                    oitr ++;
                }
            }else
            {
                // ovl_ctg as "-" => change direction of other contigs overlapping ovl_ctg
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
                    unsigned long ovl_len2 = (*iitr).second;
                    cout << "         : " << ovl_ctg2 << " <<<<<<<<<< " << ovl_ctg << " by " << ovl_len2 << " bp" << endl;
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
                        this_ctg_relevant_overlaps.insert(std::pair<string, unsigned long>(this_overlap_key1, ovl_len2));
                    }
                    //
                    iitr ++;
                }        
            }
            //
            ccitr ++;
        }        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // overlaps related to this_ctg
        cout << "   unique-2: "  << endl;
        //map<string, unsigned long>::iterator ovlitr;
        //map<string, unsigned long>::iterator ovlitr_end;
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
    // collect max_overlap_ids
    max_overlap_ids_all.insert(std::pair<string, unsigned long>(max_overlap_ids, max_overlap));
    //
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
    cout << "   equivalent-map: " << ctgX_id+ctgX_dir+"r" << " added" << endl;
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
    cout << "   equivalent-map: " << ctgY_id+ctgY_dir+"l" << " added" << endl;    
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
        cout << "   equivalent-map: " << ctgX_id+"+"+"l" << " added" << endl;    
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
        cout << "   equivalent-map: " << ctgX_id+"-"+"l" << " added" << endl;                    
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
        cout << "   equivalent-map: " << ctgY_id+"+"+"r" << " added" << endl;                            
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
        cout << "   equivalent-map: " << ctgY_id+"-"+"r" << " added" << endl;                                           
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
                   cout << "         : need new rule for " << this_overlap_ids;                
                    cout << ": either " << tmp_keya  << " or " 
                                        << tmp_keyb  << " or "
                                        << tmp_keya2 << " or "
                                        << tmp_keyb2 << " expected. " << endl;
                    // such case may need complex rule for clipping! not fixed! 20220812
                }
            }                
        }        
        //
        ovlitr ++;
    }   
    // check
    cout << "   potential clipping: " << endl;    
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
        //if(equivalent_clip.find(key1) == equivalent_clip.end() &&
        //   equivalent_clip.find(key2) == equivalent_clip.end() )
        if(1)
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
    if(0)
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
bool read_gfa(string gfa_file)
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
            assert( ctg_info.find(this_ctg) == ctg_info.end() );
            //
            NODE tmp_node;
            tmp_node.seq    = this_seq;
            tmp_node.len    = this_seq.size();
            tmp_node.pleft  = 0;                 // left  clipping position given '+unitig'
            tmp_node.pright = this_seq.size()+1; // right clipping position given '+unitig'              
            tmp_node.rm_ctg = false;             // initialized as keep the contig
            //
            ctg_info.insert(std::pair<string, NODE>(this_ctg, tmp_node));
            //
            total_gfa_ctg_size += tmp_node.len;
            total_gfa_ctg_cnt  ++;
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
                assert( ctg_info.find(this_ctg) != ctg_info.end() );                          
                vector<string> pos_info = split_string(lineinfo[6], ':');
                unsigned long overlap_pos = strtoul(pos_info[2].c_str(), NULL, 0); // from this pos to right end of focal contig to clip off
                if(overlap_pos < ctg_info[this_ctg].pright)
                {
                    // ctg_info[this_ctg].pright = overlap_pos;                    
                }
                //
                string overlap_ctg        = lineinfo[3]; // ctg_id
                string overlap_direction  = lineinfo[4]; // +-
                string overlap_len_str    = lineinfo[5].substr(0, lineinfo[5].size()-1);
                unsigned long overlap_len = strtoul(overlap_len_str.c_str(), NULL, 0);
                ctg_info[this_ctg].out.insert(std::pair<string, unsigned long>(overlap_ctg+overlap_direction, overlap_len));
            }else
            {
                // when the focus contig is "-", reverse and complement the other overlapping contigs, and update overlapping pos and size
                string this_ctg = lineinfo[1];       
                assert( ctg_info.find(this_ctg) != ctg_info.end() );                     
                vector<string> pos_info = split_string(lineinfo[6], ':');
                unsigned long overlap_pos = strtoul(pos_info[2].c_str(), NULL, 0); // from this pos to left end of focal contig to clip off
                overlap_pos = ctg_info[this_ctg].len - overlap_pos;             // always keep focus ctg as "+", and update related into correspondingly, 
                								   //    to clip [1, overlap_pos] (if 1-based system)
                if(overlap_pos > ctg_info[this_ctg].pleft)
                {
                    // ctg_info[this_ctg].pleft = overlap_pos;
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
                ctg_info[this_ctg].in.insert(std::pair<string, unsigned long>(overlap_ctg+overlap_direction, overlap_len));
            }
        }else 
        {
            continue;
        }       
    }
    ifp.close();
    return true;
}
































