/* this function, given bam_depth.txt,
                  get sharing time of haplotype intervals.
   Written by Hequan Sun, MPIPZ/LMU
   Email: sunhequan@gmail.com/sun@mpipz.mpg.de
   Date : 20230430

*/
#include        <stdio.h>
#include       <string.h>
#include         <string>
#include            <map>
#include       <stdlib.h>
#include       <iostream>
#include        <sstream>
#include        <fstream>
#include         <time.h>
#include       <assert.h>
#include        <iomanip>
#include     <sys/stat.h>
#include       <dirent.h>
#include      <algorithm>
#include "split_string.h"
//
using namespace std;
struct NODE{
  unsigned long end; // end of haplotype interval
  int           cnt; // shared times of haplotype interval
};
//
int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        cout << endl;
        cout << "   Function: given bam_depth.txt: "                                      << endl;
        cout << "             check how many times of an interval is shared by hap chrs." << endl;
        const char *buildString = __DATE__ ", " __TIME__ "";
        cout << "             (version 1.0 - compiled on " << buildString << ")"     << endl<< endl;
        cout << "   Usage: get_interval_sharing_time bam_depth.txt out_prefix"<< endl;
        cout << endl;
        exit(1);
    }
    string dfile      = (string)argv[1];
    string out_prefix = (string)argv[2];
    cout << "   info: given depth file " << dfile   << endl;
    cout << "   info: output will be labeled with " << out_prefix << endl;
    //
    map<unsigned long, NODE> interval_sharing_cnt;
    ifstream ifp;
    ifp.open(dfile.c_str(), ios::in);
    if(!ifp.good())
    {
      cout << "   Error: cannot open depth file " << dfile << endl;
      return -1;
    }
    string line("");
    while(ifp.good())
    {
       getline(ifp, line);
       if(line.size()==0 || line[0]=='#') continue;
       break;
    }
    // cout << "   check: " << line << endl;
    // initialize the map
    vector<string> lineinfo = split_string(line, '\t');
    string        this_hap  = lineinfo[0];
    unsigned long this_pos  = strtoul(lineinfo[1].c_str(), NULL, 0);
    unsigned long this_cnt  = strtoul(lineinfo[2].c_str(), NULL, 0);
    NODE tmpnode;
    tmpnode.end = this_pos;
    tmpnode.cnt = this_cnt;
    unsigned long last_key  = this_pos;
    interval_sharing_cnt.insert(std::pair<unsigned long, NODE>(last_key, tmpnode));
    //
    while( ifp.good() )
    {
      string line("");
      getline(ifp, line);
      if(line.size()==0 || line[0]=='#') continue;
      vector<string> lineinfo = split_string(line, '\t');
      if(lineinfo.size() < 3)
      {
        cout << "   warning: skipping insufficient line " << line << endl;
        continue;
      }
      if(lineinfo[0].compare(this_hap) != 0 )
      {
        cout << "   Error: this code only works with one haplotype, but "
             << this_hap << " and " << lineinfo[0] << " found in " << dfile << endl;
        cout << "          pls reformat. " << endl;
        return -1;
      }
      //
      this_pos  = strtoul(lineinfo[1].c_str(), NULL, 0);
      this_cnt  = strtoul(lineinfo[2].c_str(), NULL, 0);
      //cout << "this_pos=" << this_pos << endl;
      //cout << "this_cnt=" << this_cnt << endl;
      //cout << "interval_sharing_cnt[last_key].cnt=" << interval_sharing_cnt[last_key].cnt << endl;
      //cout << "interval_sharing_cnt[last_key].end=" << interval_sharing_cnt[last_key].end << endl;
      if( interval_sharing_cnt[last_key].cnt     == this_cnt &&
          interval_sharing_cnt[last_key].end + 1 == this_pos
        )
      {
        interval_sharing_cnt[last_key].end = this_pos;
        //cout << "      continuing..." << endl;
      }else
      {
        tmpnode.end = this_pos;
        tmpnode.cnt = this_cnt;
        last_key  = this_pos;
        interval_sharing_cnt.insert(std::pair<unsigned long, NODE>(last_key, tmpnode));
        // cout << "      new cases; last_key updated as " << last_key << " with cnt " << tmpnode.cnt << endl;
      }
    }
    ifp.close();
    // output haplotype sharing intervals and the number of haplotypes sharing them
    string ofilename = out_prefix + "_IBD_sharing_stat.txt";
    ofstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    if(!ofp.good())
    {
      cout << "   Error: cannot open file " << ofilename << endl;
      return -1;
    }
    ofp << "#hap_chr\tinterval_sta\tinterval_end\tinterval_length\tshared_time-1" << endl;
    map<unsigned long, NODE>::iterator pitr;
    map<unsigned long, NODE>::iterator pitr_end;
    pitr     = interval_sharing_cnt.begin();
    pitr_end = interval_sharing_cnt.end();
    while(pitr != pitr_end)
    {
      ofp  << this_hap           << "\t"
           << (*pitr).first      << "\t"
           << (*pitr).second.end << "\t"
           << (*pitr).second.end - (*pitr).first + 1 << "\t"
           << (*pitr).second.cnt << endl;//note, this is coverage, shared time needs to increase by 1.
      pitr ++;
    }
    ofp.close();
    //
    return 0;
}
