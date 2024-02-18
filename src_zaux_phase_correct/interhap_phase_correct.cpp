/* this tool  

       given A: "hapa.assembly,hapa.ctgs.fasta", 
             B: "hapb.assembly,hapb.ctgs.fasta", 
             and a list of contigs in A,
             
             move those list of A contigs to B,
             correct the assemblies as:
                A: "hapa.corrected.assembly,hapa.ctgs.corrected.fasta", 
                B: "hapb.corrected.assembly,hapb.ctgs.corrected.fasta"
         
   Written by Hequan Sun, MPIPZ
   Email: sunhequan@gmail.com/sun@mpipz.mpg.de
   Date.: 20230213
   
*/
#include        <stdio.h>
#include       <string.h>
#include       <stdlib.h>
#include        <fstream>
#include       <iostream>
#include        <iomanip>
#include        <sstream>
#include       <algorithm>
#include         <vector>
#include            <map>
#include         <time.h>
#include       <assert.h>
#include       <libgen.h>
#include "split_string.h"
//
bool read_assembly(string                 asm_file, 
                   vector<string>*        ctg_str, 
                   vector<string>*        ctg_str_original,                    
                   vector<int>*           ctg_int, 
                   vector<unsigned long>* ctg_size,
                   string*                ctg_tour);
bool read_moving_ctgs(string              ctg_list_file,
                   map<string, unsigned long>* ctg_moving_list,
                   vector<string>*        ctg_moving_order);
bool select_fasta(string                  genome_file, 
                   vector<string>         ctg_str_update,
                   string                 outprex, 
                   string*                ofilename);
//
int main(int argc, char* argv[])
{
    if(argc < 6)
    {
        printf("\nFunction: contig phase correction between two haplotypes.\n");
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl << endl;
        cout << "Usage: interhap_phase_correct 1.hapa.assembly 2.hapa.ctgs.fa 3.hapb.assembly 4.hapb.ctgs.fa 5.hapactg_to_hapb.list"
             << endl << endl;
        exit(1);
    }
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    printf("\nCombining tours started on %s\n", asctime(tinfo));
    //   
    // s0. read .assembly of haplotype a
    string hapa_asm = (string)argv[1];
    if(hapa_asm.find(".assembly") == std::string::npos)
    {
        cout << "   Error: assembly file must be with suffix .assembly. " << endl;
        return -1;
    }    
    string hapa_fas = (string)argv[2];    
    vector<string>        hapa_ctg_str;
    vector<string>        hapa_ctg_str_original;
    vector<int>           hapa_ctg_int;
    vector<unsigned long> hapa_ctg_size;
    string                hapa_ctg_tour;   
    if(!read_assembly(hapa_asm, 
                      &hapa_ctg_str, 
                      &hapa_ctg_str_original,
                      &hapa_ctg_int, 
                      &hapa_ctg_size,
                      &hapa_ctg_tour))
    {
        cout << "   Error: reading .assembly of file " << hapa_asm << " failed. " << endl;
        return -1;        
    }
    // s1. read .assembly of haplotype b
    string hapb_asm = (string)argv[3];
    if(hapb_asm.find(".assembly") == std::string::npos)
    {
        cout << "   Error: assembly file must be with suffix .assembly. " << endl;
        return -1;
    }
    string hapb_fas = (string)argv[4];
    vector<string>        hapb_ctg_str;
    vector<string>        hapb_ctg_str_original;    
    vector<int>           hapb_ctg_int;
    vector<unsigned long> hapb_ctg_size;
    string                hapb_ctg_tour;   
    if(!read_assembly(hapb_asm, 
                      &hapb_ctg_str, 
                      &hapb_ctg_str_original,
                      &hapb_ctg_int, 
                      &hapb_ctg_size,
                      &hapb_ctg_tour))
    {
        cout << "   Error: reading .assembly of file " << hapb_asm << " failed. " << endl;
        return -1;
    }
    //
    // s2. read ctg ids of haplotype a to be moved to haplotype b
    string ctg_list_file = (string)argv[5];
    map<string, unsigned long> ctg_moving_list;
    vector<string>             ctg_moving_order;
    if(!read_moving_ctgs(ctg_list_file,
                        &ctg_moving_list,
                        &ctg_moving_order))
    {
        return -1;
    }
    // s3. update .assembly of haplotype a
    cout << "   info: update .assembly of haplotype a.." << endl;
    vector<string>        hapa_ctg_str_update;
    vector<int>           hapa_ctg_int_update;
    vector<unsigned long> hapa_ctg_size_update;
    map<int, int>         hapa_old_to_new_id;
    int ctg_int_decrease = 0;
    for(int ci = 0; ci < hapa_ctg_str.size(); ci ++)
    {
        string this_ctg = hapa_ctg_str[ci];
        //
        if(ctg_moving_list.find(this_ctg) != ctg_moving_list.end() )
        {
            ctg_int_decrease ++;
            ctg_moving_list[ this_ctg ] = hapa_ctg_size[ci]; // collect size
        }else
        {
            hapa_ctg_str_update.push_back(  this_ctg );      // ctgs updated fasta to collect
            hapa_ctg_int_update.push_back(  hapa_ctg_int[ci] - ctg_int_decrease );
            hapa_ctg_size_update.push_back( hapa_ctg_size[ci]);
            hapa_old_to_new_id.insert(std::pair<int, int>(hapa_ctg_int[ci], hapa_ctg_int[ci] - ctg_int_decrease));
        }
    }
    // output update .assembly of haplotype a
    vector<string> ofileinfo = split_string(basename((char*)hapa_asm.c_str()), '.'); // ".assembly"
    string ofilename = "";
    for(int oi=0; oi<ofileinfo.size()-1; oi++)   
    {
        if(ofilename.size() > 0)
        {
            ofilename += ".";
        }
        ofilename += ofileinfo[ oi ];
    }
    ofilename += ".updated.assembly";
    cout << "   info: update .assembly for haplotype a: " << ofilename << endl;
    ofstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    for(int ao = 0; ao < hapa_ctg_str_update.size(); ao ++)
    {
        ofp << ">"
            << hapa_ctg_str_update[ ao ]  << " "
            << hapa_ctg_int_update[ ao ]  << " "
            << hapa_ctg_size_update[ ao ] << "\n";
    }
    // output update .tour of haplotype a
    bool first_id = true;
    vector<string> tourinfo = split_string(hapa_ctg_tour, ' ');
    for(int ti = 0; ti < tourinfo.size(); ti ++)
    {
        int this_ctg_int = atoi( tourinfo[ti].c_str() );
        if(this_ctg_int > 0)
        {
            if(hapa_old_to_new_id.find(this_ctg_int) != hapa_old_to_new_id.end())
            {
                if(first_id)
                {
                    ofp << hapa_old_to_new_id[this_ctg_int];
                    first_id = false;                
                }else
                {
                    ofp << " " << hapa_old_to_new_id[this_ctg_int];
                }
            }
        }else
        {
            if(hapa_old_to_new_id.find(-1*this_ctg_int) != hapa_old_to_new_id.end())
            {
                if(first_id)
                {
                    ofp << -1*hapa_old_to_new_id[-1*this_ctg_int];
                    first_id = false;                
                }else
                {
                    ofp << " " << -1*hapa_old_to_new_id[-1*this_ctg_int];
                }
            }        
        }
    }
    ofp << "\n";      
    ofp.close();
    // update .fasta of haplotype a
    string hapa_ofilename;
    if(!select_fasta(hapa_fas, hapa_ctg_str_update, "updated", &hapa_ofilename))
    {
        cout << "   Error: failed in updating fasta of haplotype a." << endl;
        return -1;
    }    
    // s4. update .assembly of haplotype b
    cout << "   info: update .assembly of haplotype b.." << endl;    
    vector<string>        hapb_ctg_str_update  = hapb_ctg_str;
    vector<int>           hapb_ctg_int_update  = hapb_ctg_int;
    vector<unsigned long> hapb_ctg_size_update = hapb_ctg_size;
    string                hapb_ctg_tour_update = hapb_ctg_tour;
    vector<string>::iterator ncitr;
    vector<string>::iterator ncitr_end;
    ncitr     = ctg_moving_order.begin();
    ncitr_end = ctg_moving_order.end();
    while(ncitr != ncitr_end)
    {
        string this_mov_ctg = *ncitr;
        //
        hapb_ctg_str_update.push_back(  this_mov_ctg ); // ctgs updated fasta to collect
        hapb_ctg_int_update.push_back(  hapb_ctg_int_update.size() + 1 );
        hapb_ctg_size_update.push_back( ctg_moving_list[ *ncitr ] );
        //
        std::stringstream ss;
        ss.str() = "";
        ss << hapb_ctg_int_update.size();
        hapb_ctg_tour += " ";
        hapb_ctg_tour += ss.str();
        //
        ncitr ++;
    }
    // output update .assembly of haplotype b
    ofileinfo = split_string(basename((char*)hapb_asm.c_str()), '.'); // ".assembly"
    ofilename = "";
    for(int oi=0; oi<ofileinfo.size()-1; oi++)   
    {
        if(ofilename.size() > 0)
        {
            ofilename += ".";
        }
        ofilename += ofileinfo[ oi ];
    }
    ofilename += ".updated.assembly";
    cout << "   info: update .assembly for haplotype b: " << ofilename << endl;    
    // ofstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    for(int ao = 0; ao < hapb_ctg_str_update.size(); ao ++)
    {
        ofp << ">"
            << hapb_ctg_str_update[ ao ]  << " "
            << hapb_ctg_int_update[ ao ]  << " "
            << hapb_ctg_size_update[ ao ] << "\n";
    }
    // output update .tour of haplotype b
    first_id = true;
    tourinfo = split_string(hapb_ctg_tour, ' ');
    for(int ti = 0; ti < tourinfo.size(); ti ++)
    {
        int this_ctg_int = atoi( tourinfo[ti].c_str() );
        if(1)
        {
            if(first_id)
            {
                ofp << this_ctg_int;
                first_id = false;
            }else
            {
                ofp << " " << this_ctg_int;
            }
        }
    }  
    ofp << "\n";  
    ofp.close();    
    // update .fasta of haplotype b    
    string move_ofilename;
    if(!select_fasta(hapa_fas, ctg_moving_order, "tomove", &move_ofilename))
    {
        cout << "   Error: failed in selected fasta of haplotype a to move." << endl;
        return -1;
    }    
    // prepare output file of haplotype b
    vector<string> split_genomefilename = split_string(hapb_fas, '/');    
    std::stringstream ss;
    ss.str("");
    int replace_len = 0;
    size_t pos;
    string outprex = "updated";
    string hapb_ofilename = split_genomefilename[split_genomefilename.size()-1];
    if( hapb_ofilename.find(".fasta") != std::string::npos)
    {
        ss << "." << outprex << ".fasta";
        replace_len = 6;
        pos = hapb_ofilename.find(".fasta");
    }
    else
    if( hapb_ofilename.find(".fas") != std::string::npos)
    {
        ss << "." << outprex << ".fas";
        replace_len = 4;
        pos = hapb_ofilename.find(".fas");
    }
    else
    if( hapb_ofilename.find(".fa") != std::string::npos)
    {
        ss << "." << outprex << ".fa";
        replace_len = 3;
        pos = hapb_ofilename.find(".fa");
    }
    else
    {
        cout << "   Error: input file is not with .fa, .fas or .fasta suffix. " << endl;
        return false;
    }    
    hapb_ofilename.replace(pos, pos+replace_len-1, ss.str());    
    cout << "   hapb_ofilename = " << hapb_ofilename << endl;
    std::stringstream sscmd;
    sscmd.str("");
    sscmd << "awk \'FNR==1{print \"\"}1\' " << hapb_fas << " " << move_ofilename << " | grep . > " << hapb_ofilename;            
    system( sscmd.str().c_str() );
    // cleaning
    system( "rm *tomove*.fasta" );
    system( "mkdir updated" );
    system( "mv *updated*.* ./updated" ); 
    //
    time_t endftime;
    struct tm* endtinfo;
    time(&endftime);
    endtinfo = localtime(&endftime);
    cout << "\n\nOne-seq .assembly successfully finished on " << asctime(endtinfo) << endl;
    return 0;
}
//
bool select_fasta(string genome_file, vector<string> ctg_str_update, string outprex, string* ofilename)
{
    // remove fragment issue 
    map<string, int> ctgs_to_copy;
    for(int i = 0; i < ctg_str_update.size(); i ++)
    {
        vector<string> this_ctg_info = split_string(ctg_str_update[i], ':');
        ctgs_to_copy.insert(std::pair<string, int>(this_ctg_info[0], 1));
    }
    cout << "   Info: " << ctg_str_update.size() << " fragments from " << ctgs_to_copy.size() << " contigs. " << endl;
    // prepare output file
    vector<string> split_genomefilename = split_string(genome_file, '/');
    
    std::stringstream ss;
    ss.str("");
    int replace_len = 0;
    size_t pos;
    string this_ofilename = split_genomefilename[split_genomefilename.size()-1];
    if( this_ofilename.find(".fasta") != std::string::npos)
    {
        ss << "." << outprex << ".fasta";
        replace_len = 6;    
        pos = this_ofilename.find(".fasta");
    }
    else
    if( this_ofilename.find(".fas") != std::string::npos)
    {
        ss << "." << outprex << ".fas";
        replace_len = 4;
        pos = this_ofilename.find(".fas");
    }
    else
    if( this_ofilename.find(".fa") != std::string::npos)
    {
        ss << "." << outprex << ".fa";
        replace_len = 3;
        pos = this_ofilename.find(".fa");
    }
    else
    {
        cout << "   Error: input file is not with .fa, .fas or .fasta suffix. " << endl;
        return false;
    }
    this_ofilename.replace(pos, pos+replace_len-1, ss.str());
    std::ofstream ofp;
    ofp.open(this_ofilename.c_str(), std::ios::out);
    if(!ofp.is_open())
    {
        cout << "   Error: cannot open " << this_ofilename << " to write selected sequences. " << endl;
        return false;
    }    
    *ofilename = this_ofilename;
    // open input fasta to read sequences
    std::ifstream infp (genome_file.c_str());
    if(!infp.is_open())
    {
        printf("Cannot open file \'%s\' to read sequences in fasta format. Exited.\n", genome_file.c_str());
        exit(1);
    }
    printf("   Reading sequence info from file:\t%s...\n", genome_file.c_str());    
    // traverse input fasta file
    int total_seq_num      = 0;
    int total_seq_selected = 0;
    std::string line("");
    getline(infp, line);
    while(infp.good())
    {
        if(line.find(">") != std::string::npos)
        {
            string seq_name = line;
            total_seq_num ++;
            string seq("");
            int seq_line_id = 0;
            while(infp.good())
            {
                std::string chrseq("");                
                getline(infp, chrseq);
                seq_line_id ++;
                line    = chrseq;
                if(line.find(">")!=std::string::npos) break; // next seq
                if(seq_line_id > 1) seq += '\n';
                seq    += chrseq;
            }
            if( ctgs_to_copy.find(seq_name.substr(1) ) != ctgs_to_copy.end() )
            {
                total_seq_selected ++;
                if(total_seq_selected >= 2)
                ofp << endl;
                ofp << seq_name << endl;
                ofp << seq;
            }
        }
    }
    cout << "   Info: total number of sequences in give fasta: " << total_seq_num      << endl;
    cout << "         total number of sequences selected:      " << total_seq_selected << endl;
    infp.close();
    ofp.close();
    return true;
}
//
bool read_assembly(string                 asm_file, 
                   vector<string>*        ctg_str, 
                   vector<string>*        ctg_str_original, 
                   vector<int>*           ctg_int, 
                   vector<unsigned long>* ctg_size,
                   string*                ctg_tour)
{
    ifstream ifp;
    ifp.open(asm_file.c_str(), ios::in);
    if(!ifp.good() )
    {
        cout << "   Error: cannot open file " << asm_file << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        if(line.find(">")!=std::string::npos)
        {
            // ctg info: e.g., ">hap46_utg000064l 61 25224" or ">hap46_utg000003l:::fragment_1 57 2557730"
            vector<string> lineinfo = split_string(line.substr(1), ' ');
            string         this_ctg_str  = lineinfo[0];
            vector<string> this_ctg_str_info = split_string(this_ctg_str, ':');
            int            this_ctg_int  = atoi(lineinfo[1].c_str());
            unsigned long  this_ctg_size = strtoul(lineinfo[2].c_str(), NULL, 0);
            (*ctg_str).push_back(this_ctg_str);
            (*ctg_int).push_back(this_ctg_int);
            (*ctg_size).push_back(this_ctg_size);
        }else
        {
            // tour info
            if((*ctg_tour).size() > 0)
            {
                (*ctg_tour) += " ";
            }
            (*ctg_tour) += line;            
        }
    }
    ifp.close();
    cout << "   info: " << (*ctg_str).size() << " ctgs read from assembly " << asm_file << endl;
    cout << "   tour: " << (*ctg_tour)       << endl;
    return true;
}                   
//
bool read_moving_ctgs(string                   ctg_list_file,
                   map<string, unsigned long>* ctg_moving_list,
                   vector<string>*             ctg_moving_order)
{
    ifstream snfp;
    snfp.open(ctg_list_file.c_str() );
    if(!snfp.is_open())
    {
        cout << "   Error: cannot open ctg name file: " << ctg_list_file << endl;
        return false;
    }
    while(snfp.good())
    {
        string line("");
        getline(snfp, line);
        if(line.size()==0) continue;        
        // >hap1_utg000023l 85 125374
        vector<string> lineinfo = split_string(line, ' ');
        if(lineinfo[0][0]=='>')
        { 
            lineinfo[0] = lineinfo[0].substr(1);
        }
        (*ctg_moving_list).insert(std::pair<string, unsigned long>(lineinfo[0], 0));
        (*ctg_moving_order).push_back(lineinfo[0]);
    }
    cout << "   info: " << (*ctg_moving_list).size() << " contigs to be moved. " << endl;
    snfp.close();    
    //
    return true;
}





























