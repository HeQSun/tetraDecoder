/*
    Given 
    
        >hap10_utg000009l 1 1910294
	>hap10_utg000011l 2 1560245
	>hap10_utg000039l 3 1590523
	>hap10_utg000018l 4 1066164
	
	1 2 3 4
	1 -2 3 4
    
    find out where direction of ctg is different for contigs over 1 Mb.       
    
    2023-02-27: Hequan Sun, MPIPZ
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
int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout << "\n   Function: .assembly file with last row as hic scaffolding and ragoo scaffolding, " << endl
             << "               check if contigs over 1 Mb showing consistent strand info."  << endl;
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\n   Usage: "
             << "compare_scaf_ragoo_hic updated.ragoo.assembly"
             << endl << endl;
        return 1;
    }
    double startT= clock();
    // read ref group:contigs
    string asmfile = (string)argv[1];
    ifstream ifp;
    ifp.open(asmfile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open fil " << asmfile << endl;
    }
    vector<string> ctg_str_list;         // ctg string id
    vector<string> ctg_num_list;         // ctg number id
    vector<unsigned long> ctg_size_list; // ctg size
    int i = 0;
    string hic_scaffold   = "";
    string ragoo_scaffold = "";
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size() == 0 || line[0]=='#') continue;
        // cout << line << endl;
        if(line.find(">") != std::string::npos)
        {
            vector<string> lineinfo = split_string(line, ' ');
            string ctg_str = lineinfo[0].substr(1);
            string ctg_num = lineinfo[1];
            unsigned long ctg_size = strtoul(lineinfo[2].c_str(), NULL, 0);
            ctg_str_list.push_back(ctg_str);
            ctg_num_list.push_back(ctg_num);
            ctg_size_list.push_back(ctg_size);
            // cout << ctg_str << " " << ctg_num << " " <<  ctg_size << " collected. " << endl;
        }else
        {
            if(hic_scaffold.size() == 0)
            {
                hic_scaffold   = line;
            }else
            {
                ragoo_scaffold = line;
            }
        }
    }
    // cout << "   Hi-cc scaf: " << hic_scaffold   << endl;
    // cout << "   Ragoo scaf: " << ragoo_scaffold << endl;
    //
    map<string, string> strand; // ctg hic_strand:ragoo_strand
    vector<string> hic_scaf = split_string(hic_scaffold, ' ');
    for(int i = 0; i < hic_scaf.size(); i ++)
    {
        string this_str_num = hic_scaf[i];
        string this_str_lab = "";
        if(this_str_num.find("-") != std::string::npos)
        {
            this_str_num = this_str_num.substr(1);
            this_str_lab = "-";            
        }else
        {
            this_str_num = this_str_num.substr(0);
            this_str_lab = "+";             
        }
        strand.insert(std::pair<string, string>(this_str_num, this_str_lab));
    }
    vector<string> ragoo_scaf = split_string(ragoo_scaffold, ' ');
    for(int i = 0; i < ragoo_scaf.size(); i ++)
    {
        string this_str_num = ragoo_scaf[i];
        string this_str_lab = "";
        if(this_str_num.find("-") != std::string::npos)
        {
            this_str_num = this_str_num.substr(1);
            this_str_lab = "-";            
        }else
        {
            this_str_num = this_str_num.substr(0);
            this_str_lab = "+";             
        }
        if(strand.find(this_str_num) == strand.end() )
        {
            cout << "   warning: " << this_str_num << " not found in hic scaffolding." << endl;
        }else
        {
            strand[this_str_num] += this_str_lab;
        }
    }
    //
    map<string, string>::iterator nitr;
    map<string, string>::iterator nitr_end;
    nitr     = strand.begin();
    nitr_end = strand.end();
    while(nitr != nitr_end)
    {
        string this_str_num = (*nitr).first;
        string this_str_lab = (*nitr).second;
        //
        int str_location = atoi(this_str_num.c_str());
        if(ctg_size_list[str_location-1] > 500000)
        {
            //if(this_str_lab.compare("++") != 0 && this_str_lab.compare("--") != 0)
            {
                cout << "   check: " << ctg_str_list[str_location-1] 
                     << " "          << this_str_num 
                     << " "          << ctg_size_list[str_location-1] 
                     << " "          << this_str_lab;
               if(this_str_lab.compare("++") != 0 && this_str_lab.compare("--") != 0)                     
               {
                   cout << "!!!!!!!!!!!!!!!" << endl;
               }else
               {
                   cout << endl;
               }
            }
        }
        //
        nitr ++;
    }
    //
    ifp.close();
    //
    double finishT= clock();
    // cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;    
    return 0;
}

