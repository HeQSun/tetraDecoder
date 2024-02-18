/* this function  
       given a lists of markers 
       check how many variants can be derived from these cultivars.
   
   Written by Hequan Sun, MPIPZ
   Email: sunhequan@gmail.com/sun@mpipz.mpg.de
   Date.: 20211226
     
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
int main(int argc, char* argv[])
{
    // g++ var_summer.cpp split_string.cpp -O3 -o var_summer
    if(argc < 5)
    {
        printf("\nFunction: given snp markers of multiple cultivars, check how many variants there are in total.\n");
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl << endl;        
        printf("Usage: var_summer common-path output_list label var_set.txt \n\n");
        cout << "Info: common-path: path to resequencing analysis of cultivars. " << endl;
        cout << "Info: label......: A,B,C,...U"                                   << endl;
        cout << "Info: var_set.txt: the list of variants to be analyzed. "        << endl;
        cout << "Info: common-path+label => /path/to/sample/quality_variant_iQ32_simple.txt" << endl;
        cout << "Info: variaton file should follow format: 1.CHR tab 2.POS tab 3.REF tab 4.ALT" << endl;
        cout << "Info: output_list=0, only output count variations; otherwise, output counts and variations. " 
             << endl
             << endl;
        cout << "Example: var_summer \"common-path-p5312_\" 0 A,B,C "
             << "/shore_consensus_het_local_mkdup/ConsensusAnalysis/quality_variant_iQ32_simple.txt" 
             << endl
             << "         checking quality_variant_iQ32_simple.txt of samples A,B,C under common path. "
             << endl
             << endl;
        exit(1);
    }
    // get files
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    printf("\nCounting variations started on %s\n", asctime(tinfo));    
    // step 1: read allele of parent 1: parent2 alleles
    string common_path        = (string)argv[1];
    int output_list           = atoi(argv[2]);
    string label              = (string)argv[3]; // A/B/C/...U
    string var_set            = (string)argv[4];
    vector<string> sampleinfo = split_string(label, ',');
    if(sampleinfo.size() < 2)
    {
        cout << "   Waring: it seems only one sample provided, no need to continue. Exited. " << endl;
        return 0;
    }
    //
    cout << "   Info: variation reading starting..." << endl;
    map<string, int> all_var; // <key=chr\tpos\talt, val=count>
    unsigned long observed_var_cnt = 0;
    vector<string>::iterator vitr;
    vector<string>::iterator vitr_end;
    vitr     = sampleinfo.begin();
    vitr_end = sampleinfo.end();
    while(vitr != vitr_end)
    {
        string this_var_file = common_path + "_" + *vitr + "/" + var_set;        
        //
        ifstream ifp;
        ifp.open(this_var_file.c_str(), ios::in);
        if(!ifp.good())
        {
           cout << "   Error: cannot open var file " << this_var_file << endl;
        }
        cout << "   Info: checking var file for "    << *vitr << endl;
        unsigned long observed_var_cnt_this_file = 0;
        while(ifp.good())
        {
            string line("");
            getline(ifp, line);
            if(line.size()==0 || line[0]=='#') continue;
            vector<string> lineinfo = split_string(line, '\t');
            if(lineinfo.size() < 4)
            {
                cout << "   Warning: skipping line with insufficient info: " << line << endl;
                continue;
            }
            string this_key = lineinfo[0] + "\t" + lineinfo[1] + "\t" + lineinfo[3];
            //cout << "   check: line = " << line << ", key=" << this_key << endl;
            if(all_var.find(this_key) == all_var.end() )
            {
                all_var.insert(std::pair<string, int>(this_key, 1));
            }else
            {
                all_var[this_key] += 1;
            }
            observed_var_cnt ++;
            observed_var_cnt_this_file ++;
            if(observed_var_cnt_this_file%5000000==0)
            {
                cout << "   Info: " << observed_var_cnt_this_file << "th observed variations." << endl;
            }
        }
        cout << "   Info: total number of variations observed in current file: " << observed_var_cnt_this_file 
             << ", accumulated unique variations: "                              << all_var.size() 
             << endl;
        //
        vitr ++;
    }
    cout << "   Info: variation reading done." << endl;
    if(output_list != 0)
    {
        cout << "   Info: output list of variations..." << endl;
        /* TODO */
        cout << "   Info: output list of variations done." << endl;        
    }else
    {
        cout << "   Info: output list of variations was not asked. " << endl;
    } 
    cout << "   Info: total number of unique variations: " << all_var.size()  
         << ", from observed: "                            << observed_var_cnt 
         << std::fixed          << std::setprecision(4)  
         << ", uniqueness: "    << all_var.size()*1.0 / observed_var_cnt 
         << ", from cultivars " << label
         << endl;
    //
    time_t endftime;
    struct tm* endtinfo;
    time(&endftime);
    endtinfo = localtime(&endftime);
    printf("\nCounting variations successfully finished on %s\n", asctime(endtinfo));
    return 0;
}
//
