/* this tool  

       given a .assembly file, reverse the tour of contigs.
         
   Written by Hequan Sun, MPIPZ
   Email: sunhequan@gmail.com/sun@mpipz.mpg.de
   Date.: 20230130
   
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
    if(argc < 2)
    {
        printf("\nFunction: given a .assembly file, reverse the tour of contigs (last row).\n");
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl << endl;
        cout << "Usage: reverse xxxx.assembly " << endl;
        cout << "Note, only works with one chr tour of contigs (i.e., merging will be applied for separate chrs) . "
             << endl << endl;
        exit(1);
    }
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    printf("\nUnitig end cleaning started on %s\n", asctime(tinfo));
    // s0. get inputs
    string asm_file     = (string)argv[1];
    //
    vector<string> contigs;
    string tour_string="";
    ifstream ifp;
    ifp.open(asm_file.c_str(), ios::in);
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        if(line.find(">") != std::string::npos) 
        {
            // collect contigs 
            contigs.push_back(line);
            continue;
        }
        if(line.size() > 0)
        {
            if(tour_string.size() > 0)
            {
                tour_string += " ";
            }
            tour_string += line;
        }
    }
    ifp.close();
    cout << "   Info: "          << contigs.size() << " contings collected. " << endl;
    cout << "   Tour (merged, observed): " << endl << tour_string             << endl;
    //
    cout << "   Tour (merged, reversed): " << endl;
    string tour_string_rev = "";
    bool first_ctg = true;
    vector<string> this_tour = split_string(tour_string, ' ');
    for(int ti = this_tour.size()-1; ti>=0; ti --)
    {
        string this_ctg_id = this_tour[ti];
        if(this_ctg_id[0]=='-')
        {
            this_ctg_id = this_ctg_id.substr(1);
        }else
        {
            this_ctg_id = "-" + this_ctg_id;
        }
        if(first_ctg)
        {
            cout << this_ctg_id;
            tour_string_rev += this_ctg_id;
            first_ctg = false;
        }else
        {
            cout << " " << this_ctg_id;
            tour_string_rev += " " + this_ctg_id;            
        }
    }
    // output reversed .assembly
    vector<string> ipathinfo     = split_string(asm_file, '/');
    vector<string> ofilenameinfo = split_string(ipathinfo[ipathinfo.size() - 1], '.');
    string ofilename     = "";
    for(int oi = 0; oi < ofilenameinfo.size()-1; oi++)
    {
        if(oi > 0)
        {
            ofilename += ".";
        }
        ofilename += ofilenameinfo[oi];        
    }
    ofilename += ".rev.assembly\0";
    ofstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open output file " << ofilename << endl;
    }
    for(int ci = 0 ; ci < contigs.size(); ci ++)
    {
        ofp << contigs[ci] << endl;
    }
    ofp << tour_string_rev << endl;
    ofp.close();
    //
    time_t endftime;
    struct tm* endtinfo;
    time(&endftime);
    endtinfo = localtime(&endftime);
    cout << "\n\nReversing .assembly successfully finished on " << asctime(endtinfo) << endl;
    return 0;
}































