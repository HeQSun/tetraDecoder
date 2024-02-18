/*  this function takes as input blast of a contig to ncbi nt database: ${ctg}_againt_ncbi_nt_nt_all.oblast 
	finds the subject feature with best highest identities/match

    given a ctg="the_contig_id", blast should be done with below:
    
    for i in {0..9}; do 
        blastn -query ${ctg}.fa -db /opt/share/blastdb/ncbi/nt.0${i} -out ${ctg}_againt_ncbi_nt_nt_0${i}.oblast -outfmt 0 -max_target_seqs 5
    done
    for i in {10..26}; do 
        blastn -query ${ctg}.fa -db /opt/share/blastdb/ncbi/nt.${i} -out ${ctg}_againt_ncbi_nt_nt_${i}.oblast -outfmt 0 -max_target_seqs 5
    done        
    cat ${ctg}_againt_ncbi_nt_nt_*.oblast > ${ctg}_againt_ncbi_nt_nt_all.oblast 
    
    #
    written by Hequan Sun, MPIPZ
    Email: sun@mpipz.mpg.de/sunhequan@gmail.com
    
*/
#include       <map>
#include    <string>
#include    <vector>
#include   <fstream>
#include   <sstream>
#include  <iostream>
#include <algorithm>
#include  <string.h>
#include  <stdlib.h>
#include  <assert.h>
#include <sys/stat.h>
#include   <dirent.h>
#include   <iomanip>      // std::setprecision
#include "split_string.h"

struct NODE{
    string        subject;
    unsigned long identity;
    unsigned long coverage;
};
//
using namespace std;
bool collect_contig_size(string                             ctgsize_file, 
                  map<string, unsigned long>*               contigsize);
//
int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        cout << "\nFunction: find best mathing subject sequence info with a blast output of a contig. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "Usage: get_organelle_ctg_stat ${ctg}_againt_ncbi_nt_nt_all.oblast ctg_sizes.txt"
             << endl << endl;
        return 1;
    }
    double startT= clock();
    // get input info 
    string blast_file    = (string)argv[1];
    string ctgsize_file  = (string)argv[2];
    // read size info
    cout << endl;
    cout << "   Info: read contig size info... " << endl;
    map<string, unsigned long> contigsize;
    if(!collect_contig_size(ctgsize_file, &contigsize))
    {
        cout << "   Error: failed in collecting contig size info." << endl;
        exit(1);
    } 
    cout << "   Info: read contig size info done. " << endl;       
    //
    cout << "   Info: read blast output of the contig..." << endl;
    ifstream ifp; 
    ifp.open(blast_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << blast_file << endl;
        return 1;
    }
    NODE tmpnode;
    tmpnode.subject    = ">nohit";
    tmpnode.identity   = 0;
    tmpnode.coverage   = 0;
    vector<NODE> best_identity_subject; // collecting 3 only 
    tmpnode.subject    = ">nohit1";  
    best_identity_subject.push_back(tmpnode);
    tmpnode.subject    = ">nohit2";      
    best_identity_subject.push_back(tmpnode);
    tmpnode.subject    = ">nohit3";      
    best_identity_subject.push_back(tmpnode);  
    string query=""; 
    string line("");
    getline(ifp, line);         
    while(ifp.good())
    {
        if(line.size()==0) 
        {
            getline(ifp, line);        
            continue;
        }
        if(line.find("Query= ") != std::string::npos)
        {
            vector<string> lineinfo = split_string(line, ' ');
            if(query.size() != 0)
            {
                assert(lineinfo[1].compare(query) == 0);
            }
            query = lineinfo[1];
        }
        if(line.find(">") != std::string::npos)
        {
            string subject = line;
            //
            string line2("");
            getline(ifp, line2);  // Identities = 1926/1961 (98%), Gaps = 20/1961 (1%)
            if(line2.find(">") != std::string::npos) 
            {
                line = line2; 
                continue;
            }  
            while(ifp.good() && line2.find(">")==std::string::npos)
            {
                getline(ifp, line2);  // Identities = 1926/1961 (98%), Gaps = 20/1961 (1%)
                if(line2.find("Identities =") != std::string::npos)
                {
                    //cout << line2 << endl;
                    vector<string> lineinfo = split_string(line2, ' ');
                    vector<string> ideninfo = split_string(lineinfo[2], '/');
                    unsigned long coverage  = strtoul(ideninfo[1].c_str(), NULL, 0);
                    unsigned long identity  = strtoul(ideninfo[0].c_str(), NULL, 0);
                    if(coverage > best_identity_subject[0].coverage)
                    {
                        best_identity_subject[2] = best_identity_subject[1];
                        best_identity_subject[1] = best_identity_subject[0];
                        best_identity_subject[0].subject  = subject;
                        best_identity_subject[0].coverage = coverage;
                        best_identity_subject[0].identity = identity;
                        //cout << "   check: max coverage updated as " << best_identity_subject[0].coverage << endl;                        
                    }else
                    if(coverage > best_identity_subject[1].coverage)
                    {
                        best_identity_subject[2] = best_identity_subject[1];
                        best_identity_subject[1].subject  = subject;
                        best_identity_subject[1].coverage = coverage;    
                        best_identity_subject[1].identity = identity;
                    }else 
                    if(coverage > best_identity_subject[2].coverage)
                    {
                        best_identity_subject[2].subject  = subject;
                        best_identity_subject[2].coverage = coverage;
                        best_identity_subject[2].identity = identity;
                    }else
                    {
                        ;
                    } 
                }  
            }
            if(line2.find(">") != std::string::npos) 
            {
                line = line2; 
            }                      
        }else
        {
            getline(ifp, line);                
        }
    }
    cout << endl;
    cout << "   Report: top 3 hits of contig " << query << " of size " << contigsize[query] << " bp includes " << endl;
    for(int i=0; i<3; i++)
    {
        cout << "           " 
             << query                                      << "\t"
             << contigsize[query]                          << "\t"
             << best_identity_subject[i].subject.substr(1) << "\t"
             << best_identity_subject[i].coverage          << "\t"
             << best_identity_subject[i].identity          << endl;
    }
    cout << "   Note: contig-id size-bp subject coverage-bp identity-bp; tab-separated. " << endl;
    cout << endl;
    //     
    ifp.close();
    cout << "   Info: read blast output of the contig done." << endl;
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    //
    return 0;
}
//
bool collect_contig_size(string ctgsize_file, map<string, unsigned long>* contigsize)
{
    cout << "   Info: reading contig size info from " << ctgsize_file << endl;
    ifstream ifp;
    ifp.open(ctgsize_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open contig size file " << ctgsize_file << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size() == 0) continue;
        if(line[0]=='#')     continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size()<2) continue;
        string        ctgid   = lineinfo[0];
        unsigned long ctgsize = strtoul(lineinfo[1].c_str(), NULL, 0);
        (*contigsize).insert(std::pair<string, unsigned long>(ctgid, ctgsize));
        //// cout << ctgid << "\t" << ctgsize << endl;
    }
    if((*contigsize).size() == 0)
    {
        cout << "   Erorr: no info collected from contig size file. " << endl;
        return false;
    }else
    {
        cout << "   Info: " << (*contigsize).size() << " contig size info collected. " << endl;
    }
    ifp.close();
    //
    return true;
}
