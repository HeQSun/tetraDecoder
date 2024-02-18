/* this tool

       given A: "hapa.chr.fasta and chr sta end",
             B: "hapb.ctg.fasta, hapb.ctg.assembly, hapb.ctg.agp, and a insertion site",

             extract sequence from hapa.chr.fasta within [sta, end] => an additional contig x
             insert contig x into hapb.ctg.fasta
             insert contig x into hapb.ctg.assembly according to closest gap to given insertion site.

   Written by Hequan Sun, MPIPZ
   Email: sunhequan@gmail.com/sun@mpipz.mpg.de
   Date.: 20230321
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
#include "extract_substr_chr.h"
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
bool insert_seq_to_fasta(string Bctgfas, string extractedfas, string outprefix);
bool find_insertion_site(string Bctgagp, long Bip, string* target_ctg);
bool update_assembly(string Bctgasm, string target_ctg, string subseq_name, long subseq_size, string outprefix);
//
int main(int argc, char* argv[])
{
    if(argc < 10)
    {
        printf("\nFunction: extract sequence from hapA and insert it to hapB and update .assembly of hapB.\n");
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl << endl;
        cout << "Usage: interhap_gap_filling 1.Achr.fasta 2.A-chr 3.A-sta 4.A-end 5.Bctgs.fa 6.Bctgs.assembly 7.Bctgs.gap 8.B.insertion.pos 9.outprefix"
             << endl << endl;
        cout << "Caution! In .agp file, chromosome ids for insertion/selection must start with \"chr\"" << endl << endl;
        exit(1);
    }
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    printf("\nExtract subseq from hap A to fill gaps in hap B started on %s\n", asctime(tinfo));
    // s1. get options
    string Achrfas = (string)argv[1]; // chr level
    string Achr    = (string)argv[2];
    long   Asta    = atol(argv[3]);
    long   Aend    = atol(argv[4]);
    string Bctgfas = (string)argv[5]; // ctg level
    string Bctgasm = (string)argv[6];
    string Bctgagp = (string)argv[7];
    long   Bip     = atol(argv[8]);   // potential insertion position of ctg from A into B.assembly
    string outprefix = (string)argv[9];
    // s2. extract subsequence from the given Achr of A
    string subseq_name = "";
    long   subseq_size = 0;
    if(!extract_substr_chr(Achrfas, Achr, Asta, Aend, &subseq_name, &subseq_size))
    {
      cout << "   Error: failed in extracting subsequence from chr fasta: " << Achrfas << endl;
      return -1;
    }
    // s3. insert extracted sequence into fasta of hap A
    if(!insert_seq_to_fasta(Bctgfas, (string)"subseq_tmp.fasta", outprefix))
    {
      cout << "   Error: failed in inserting new sequence into " << Bctgfas << endl;
      return -1;
    }
    // s4. read B.gap and find proper insertion position in .tour line
    cout << "   Info: given potential insertion position: " << Bip << endl;
    string target_ctg;
    if(!find_insertion_site(Bctgagp, Bip, &target_ctg) )
    {
      cout << "   Error: failed in finding proper insertion site. " << endl;
      return -1;
    }
    cout << "   Info: subseq will be inserted after " << target_ctg << endl;
    // s5. read B.assembly and insert new subseq as contig
    if(!update_assembly(Bctgasm, target_ctg, subseq_name, subseq_size, outprefix) )
    {
      cout << "   Info: failed in updating assembly: " << Bctgasm << endl;
      return -1;
    }
    //
    time_t endftime;
    struct tm* endtinfo;
    time(&endftime);
    endtinfo = localtime(&endftime);
    cout << "\nGap filling in B.apg successfully finished on " << asctime(endtinfo) << endl;
    return 0;
}
//
bool update_assembly(string Bctgasm, string target_ctg, string subseq_name, long subseq_size, string outprefix)
{
  string ofilename = "res_" + outprefix + ".assembly";
  ofstream ofp(ofilename.c_str(), ios::out);
  if(!ofp.good())
  {
    cout << "   Error: cannot open output file " << ofilename << endl;
    return false;
  }
  //
  ifstream ifp;
  ifp.open(Bctgasm.c_str(), ios::in);
  if(!ifp.good())
  {
    cout << "   Error: cannot open file " << Bctgasm << endl;
    return false;
  }
  int last_ctg_cnt   = 0;
  int target_ctg_cnt = 0;
  int subseq_cnt    = -1;
  while(ifp.good())
  {
    string line("");
    getline(ifp, line);
    if(line.size()==0 || line[0]=='#') continue;
    vector<string> lineinfo = split_string(line, ' ');
    // >hap4_utg000009l 13 270977
    // -37 -1 76 4 -6 7 -8 -9 99 -51
    // 84 85 87 86 -90 88 -91 79 -105 -89
    if(line.find(">") != std::string::npos)
    {
      string ctg_id    = lineinfo[0].substr(1);
      int    ctg_cnt   = atoi(lineinfo[1].c_str());
      long   ctg_size  = atol(lineinfo[2].c_str());
      last_ctg_cnt     = ctg_cnt;
      if(ctg_id.compare(target_ctg) == 0)
      {
        target_ctg_cnt = ctg_cnt;
      }
      ofp << line << endl;
    }else
    {
      if(subseq_cnt < 0)
      {
        subseq_cnt = last_ctg_cnt + 1;
        ofp << ">" << subseq_name << " " << subseq_cnt << " " << subseq_size << endl;
      }
      // update the first .tour line
      for(int ti = 0; ti < lineinfo.size(); ti ++)
      {
        int this_ctg_cnt = atoi(lineinfo[ti].c_str());
        if(this_ctg_cnt < 0)
        {
          this_ctg_cnt = this_ctg_cnt*(-1);
        }
        // output
        if(this_ctg_cnt == target_ctg_cnt)
        {
          // if this is the place to insert extracted subseq
          if(ti == 0)
          {
            ofp << lineinfo[ti] << " " << subseq_cnt;
          }else
          {
            ofp << " ";
            ofp << lineinfo[ti] << " " << subseq_cnt;
          }
          cout << "   Info: "        << subseq_name  << " " << subseq_cnt
               << " inserted after " << target_ctg   << " " << this_ctg_cnt << endl;
        }else
        {
          // other cases
          if(ti == 0)
          {
            ofp << lineinfo[ti];
          }else
          {
            ofp << " ";
            ofp << lineinfo[ti];
          }
        }
      }
      ofp << endl;
    }
  }
  ifp.close();
  ofp.close();
  cout << "   Info: new assembly created as " << ofilename << endl;
  //
  return true;
}
//
bool find_insertion_site(string Bctgagp, long Bip, string* target_ctg)
{
  /* find the best position (after target_ctg) to insert the new sequence in .tour */
  ifstream ifp;
  ifp.open(Bctgagp.c_str(), ios::in);
  if(!ifp.good())
  {
    cout << "   Error: cannot open file " << Bctgagp << endl;
    return false;
  }
  long dist = 999999999;
  (*target_ctg) = "";
  long target_ctg_end = -1;
  while(ifp.good())
  {
    string line("");
    getline(ifp, line);
    if(line.size()==0 || line[0]=='#') continue;
    if(line.find("\tU\t")!=std::string::npos || line.find("chr")==std::string::npos) continue;
    // chr01_hap4      79949013        81779195        207     W       hap1_utg000020l 1       1830183 -
    // cout << line << endl;
    vector<string> lineinfo = split_string(line, '\t');
    if(lineinfo.size() < 9)
    {
      cout << "   warning: skipping line with insufficient info: " << line << endl;
      continue;
    }
    long ctg_end_at_chr = atol(lineinfo[2].c_str());
    string ctg_id = lineinfo[5];
    if(ctg_end_at_chr<=Bip && Bip-ctg_end_at_chr<dist)
    {
      (*target_ctg)  = ctg_id;
      dist           = Bip-ctg_end_at_chr;
      target_ctg_end = ctg_end_at_chr;
    }else
    if(ctg_end_at_chr>Bip && ctg_end_at_chr-Bip<dist)
    {
      (*target_ctg)  = ctg_id;
      dist           = ctg_end_at_chr - Bip;
      target_ctg_end = ctg_end_at_chr;
    }else
    {
      ;
    }
  }
  cout << "   Info: insertion will happen after chr-position " << target_ctg_end
       << " = end of contig "                                  << (*target_ctg)
       << endl;
  ifp.close();
  return true;
}
//
bool insert_seq_to_fasta(string Bctgfas, string extractedfas, string outprefix)
{
  string ofilename = "res_" + outprefix + ".fasta";
  ofstream ofp(ofilename.c_str(), ios::out);
  if(!ofp.good())
  {
    cout << "   Error: cannot open output file " << ofilename << endl;
    return false;
  }
  vector<string> fasmap;
  fasmap.push_back(Bctgfas);
  fasmap.push_back(extractedfas);
  for(int fi=0; fi<fasmap.size(); fi++)
  {
    string this_fas = fasmap[fi];
    ifstream ifp;
    ifp.open(this_fas.c_str(), ios::in);
    if(!ifp.good())
    {
      cout << "   Error: cannot open file " << this_fas << endl;
      return false;
    }
    while(ifp.good())
    {
      string line("");
      getline(ifp, line);
      if(line.size()==0 || line[0]=='#') continue;
      ofp << line << endl;
     }
    ifp.close();
  }
  //
  ofp.close();
  cout << "   Info: new fasta file created as " << ofilename << endl << endl;
  return true;
}
