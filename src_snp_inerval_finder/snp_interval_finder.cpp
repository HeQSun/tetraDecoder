/*
    Given 
    
        a list of phased snp markers within a genomic region
        
    cluster snps into intervals.
    
    2022-09-23 Started: Hequan Sun, MPIPZ
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
//
struct PAMA
{
    string patallele;
    string matallele;
};
//
bool read_marker(string fmarker,
                 map<string, map<unsigned long, PAMA> >* mkr);
//
int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        cout << "\n   Function: given a list of snps, find their clusters of genomic regions. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "snp_interval_finder snp_list.txt gap_size out_prefix_str"
             << endl << endl;
        return 1;
    }
    double startT= clock();
    //
    string marker_file = (string)argv[1];
    int gap_size       = atoi(argv[2]);
    string out_prefix  = (string)argv[3];
    cout << "   Info: maximum gap allowed between snps of the same cluster is " << gap_size << " bp." << endl;
    // s1. read markers
    cout << "   Info: reading snp markers from " << marker_file << endl;
    map<string, map<unsigned long, PAMA> > snp_marker; // map<chr, map<position, {mat-allele, pat-allele} > >
    if(!read_marker(marker_file, &snp_marker) )
    {
        cout << "   Error: read snp marker failed. " << endl;
        return false;
    }      
    cout << "   Info: reading snp marker done. " << endl;
    // s2. find clustes of snps
    map<string, map<unsigned long, unsigned long> > snp_interval; // map<chr, map<interval_sta, interval_end> >    
    map<string, map<unsigned long, PAMA> >::iterator ctg_itr;
    map<string, map<unsigned long, PAMA> >::iterator ctg_itr_end;
    ctg_itr     = snp_marker.begin();
    ctg_itr_end = snp_marker.end();
    while(ctg_itr != ctg_itr_end)
    {
        string this_ctg = (*ctg_itr).first;
        map<unsigned long, PAMA> this_snp_list = (*ctg_itr).second;
        //
        map<unsigned long, PAMA>::iterator pos_itr;
        map<unsigned long, PAMA>::iterator pos_itr_end;
        pos_itr     = this_snp_list.begin();
        pos_itr_end = this_snp_list.end();
        // initialize
        unsigned long last_pos = (*pos_itr).first;
        map<unsigned long, unsigned long> tmp_interval;
        tmp_interval.insert(std::pair<unsigned long, unsigned long>(last_pos, last_pos));
        //
        pos_itr ++;        
        while(pos_itr != pos_itr_end)
        {
            unsigned long this_pos = (*pos_itr).first;
            if(this_pos - tmp_interval[last_pos] < gap_size)
            {
                tmp_interval[last_pos] = this_pos;
            }else
            {
                if(snp_interval.find(this_ctg) == snp_interval.end())
                {
                    snp_interval.insert(std::pair<string, map<unsigned long, unsigned long> >(this_ctg, tmp_interval));                    
                }else
                {
                    snp_interval[this_ctg].insert(std::pair<unsigned long, unsigned long>(last_pos, tmp_interval[last_pos]));                    
                }
                // new window
                last_pos = this_pos;
                tmp_interval.clear();                
                tmp_interval.insert(std::pair<unsigned long, unsigned long>(last_pos, last_pos));
            }
            //
            pos_itr ++;
        }
        // last window
        if(snp_interval.find(this_ctg) == snp_interval.end())
        {
            snp_interval.insert(std::pair<string, map<unsigned long, unsigned long> >(this_ctg, tmp_interval));                    
        }else
        {
            snp_interval[this_ctg].insert(std::pair<unsigned long, unsigned long>(last_pos, tmp_interval[last_pos]));                    
        }        
        //
        ctg_itr ++;
    }
    // s3. output clusters
    string ofilename = out_prefix + "_snp_clusters.txt";
    ofstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open output file " << ofilename << endl;
        return 1;
    }
    //
    long n_ctg_with_snps = snp_interval.size();
    long n_snp_interval  = 0;
    map<string, map<unsigned long, unsigned long> >::iterator ctg_itr2;
    map<string, map<unsigned long, unsigned long> >::iterator ctg_itr2_end;
    ctg_itr2     = snp_interval.begin();
    ctg_itr2_end = snp_interval.end();
    while(ctg_itr2 != ctg_itr2_end)
    {
        map<unsigned long, unsigned long>::iterator clu_itr;
        map<unsigned long, unsigned long>::iterator clu_itr_end;        
        clu_itr     = (*ctg_itr2).second.begin();
        clu_itr_end = (*ctg_itr2).second.end();
        while(clu_itr != clu_itr_end)
        {
            unsigned long region_size = (*clu_itr).second-(*clu_itr).first;
            ofp  << (*ctg_itr2).first << ":" 
                 << (*clu_itr).first  << "-" 
                 << (*clu_itr).second << "\t" 
                 << region_size
                 << endl;
            clu_itr ++;
        }
        n_snp_interval += (*ctg_itr2).second.size();
        //
        ctg_itr2 ++;
    }
    cout << "   Info: on " << n_ctg_with_snps << " ctgs, " << n_snp_interval << " clusters of snps found." << endl;
    //
    ofp.close();
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;
}
// read markers
bool read_marker(string fmarker,
                 map<string, map<unsigned long, PAMA> >* mkr)
{
    /*
      format: contig_id   position .   maternal paternal  optional
              utg000103lc 727428   .   A        T         0.2 RefCall .....
    */
    ifstream ifp;
    ifp.open(fmarker.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open marker file " << fmarker << endl;
        return false;
    }
    map<string, bool>::iterator cpitr;
    map<string, map<unsigned long, PAMA> >::iterator cmitr;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0]  =='#') continue;        
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 5) 
        {
            cout << "   Warning: unexpected line: " << line;
            continue;
        }
        string        thisctg   = lineinfo[0];
        unsigned long thispos   = strtoul(lineinfo[1].c_str(), NULL, 0);
        string        matallele = lineinfo[3]; // M - ref
        string        patallele = lineinfo[4]; // P - alt
        //
        PAMA tpama;
        tpama.matallele = matallele;
        tpama.patallele = patallele;
        //
        cmitr = (*mkr).find(thisctg);
        if(cmitr != (*mkr).end())
        {
            (*cmitr).second.insert(std::pair<unsigned long, PAMA>(thispos, tpama));
        }else
        {
            map<unsigned long, PAMA> posmkr;
            posmkr.insert(std::pair<unsigned long, PAMA>(thispos, tpama));
            (*mkr).insert(std::pair<string, map<unsigned long, PAMA> >(thisctg, posmkr));
        }
    }
    ifp.close();
    return true;
}
