/* this function, given syri.out between n pairs of hap-chrs,
                  get the variations in the pair-aligned regions, and analyze how IBDs are shared.

                  v4: including syn and nonsyn regions - labelled in output files.

   Written by Hequan Sun, MPIPZ/LMU
   Email: sunhequan@gmail.com/sun@mpipz.mpg.de
   Date : 20230319
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
//
struct VAR{
    string         rlg;             // ref lg id
    unsigned long rsta;             // ref lg loci sta
    unsigned long rend;             // ref lg loci end
    string         qlg;             // qry lg
    unsigned long qsta;             // qry lg loci sta
    unsigned long qend;             // qry lg loci end
    map<string, string> var_list;   // list of variations in [sta, end]
};
struct IBD{
    unsigned long sta;              // start of this IBD interval
    unsigned long end;              // end   of this IBD interval
    int           cnt;              // region shared by how many haplotype chrs.
    unsigned long sta_origin;       // start of this IBD interval
    unsigned long end_origin;       // end   of this IBD interval
    map<string, int> hap;           // haps sharing this IBD
};
string report_indels = "0";         // "0": do not report indels vs "1": report indels
//
bool read_so_list(string                   so_file_list,
                  vector<string>*          so_list);
bool read_syri_out(string                  so_file,
                  string                   report_indels,
                  int                      IBD_max,
                  string                   out_dir,
                  string                   out_prefix,
                  vector<string>*          ibd_files);
//
int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        cout << endl;
        cout << "   Function: given syri.out between two hap-chrs: "                        << endl;
        cout << "             1. get variations in syntenic & rearranged regions of 2 haps" << endl;
        cout << "             2. select potential IBD regions. "                            << endl;
        const char *buildString = __DATE__ ", " __TIME__ "";
        cout << "             (version 1.0 - compiled on " << buildString << ")"     << endl<< endl;
        cout << "   Usage: get_vars_within_region 1.syri.out.list 2.indels 3.IBD-max"<< endl;
        cout << "          note: if indels==1, will consider indels within the intervals."  << endl;
        cout << "          note: the larger the value of IBD-max, the more likely it is IBD"<< endl;
        cout << endl;
        exit(1);
    }
    // get options:
    string so_file_list = (string)argv[1]; // list of syri.out files
    report_indels       = (string)argv[2]; // 1: report indels; 0: iganore indels
    int    IBD_max      =   atoi(argv[3]); // the larger the value, the more likely it is IBS or IBD
    string out_prefix;                     //
    // s1. read syri.out filenames
    vector<string> so_list;
    if(!read_so_list(so_file_list.c_str(), &so_list))
    {
      cout << "   Error: failed in reading the list of syri.out files. " << endl;
      return -1;
    }
    // s2. prepare output folder
    string out_dir = "result_stat_etc";
    DIR* dir = opendir( out_dir.c_str() );
    if (dir)
    {
        /* Directory exists. */
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        const int dir_err = mkdir( out_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << out_dir << endl;
            return -1;
        }
    }
    else ;
    // s3. find statistics of variations in syn regions and potential IBD regions.
    vector<string> ibd_files;
    for(int fi=0; fi<so_list.size(); fi++)
    {
        string this_so_file = so_list[fi]; // 1.samplei 2.lgii 3.hapii 4.samplej 5.lgjj 6.hapjj 7.syri.out
        vector<string> fileinfo = split_string(this_so_file, ' ');
        if(fileinfo.size()<7)
        {
          cout << "   warning: insufficient file info: " << this_so_file << endl;
        }
        out_prefix = fileinfo[0] + "_" + fileinfo[1] + "_" + fileinfo[2] + "_" +
                     fileinfo[3] + "_" + fileinfo[4] + "_" + fileinfo[5] ;
        cout << "   info: reading " << fileinfo[6] << "..." << endl;
        //
        if(!read_syri_out(fileinfo[6], report_indels, IBD_max, out_dir, out_prefix, &ibd_files))
        {
          cout << "   Error: failed in reading syri.out file " << fileinfo[6] << endl;
          return -1;
        }
        cout << "   info: reading " << fileinfo[6] << " done." << endl;
    }
    //
    return 0;
}
//
bool read_so_list(string so_file_list, vector<string>* so_list)
{
  ifstream ifp;
  ifp.open(so_file_list.c_str(), ios::in);
  if(!ifp.good())
  {
    cout << "   Error: cannot open file of list of syri.out files. " << endl;
    return false;
  }
  while(ifp.good())
  {
    string line("");
    getline(ifp, line);
    if(line.size()==0) continue;
    if(line[0]=='#') continue;
    (*so_list).push_back(line);
  }
  ifp.close();
  cout << "   info: " << (*so_list).size() << " files collected." << endl;
  return true;
}
//
bool read_syri_out(string so_file,
                   string report_indels,
                   int    IBD_max,
                   string out_dir,
                   string out_prefix,
                   vector<string>* ibd_files)
{
  ifstream ifp;
  // out_prefix: e.g., A_chr06_21_O_chr06_22 = refsample_reflg_refhap_qrysample_qrylg_qryhap
  vector<string> out_flags = split_string(out_prefix, '_');
  string refsample = out_flags[0];
  string qrysample = out_flags[3];
  // 1st round: read syn regions
  ifp.open(so_file.c_str(), ios::in);
  if(!ifp.good())
  {
    cout << "   Error: cannot open file " << so_file << endl;
    return false;
  }
  map<string, VAR> syn_var; // < "SYNi/TRANj/INVk", {sample_rlg, rsta, rend, {var}, qlg, qsta, qend} >
  vector<string> syn_order;
  while(ifp.good())
  {
    string line("");
    getline(ifp, line);
    if(line.size()==0 || line[0]=='#') continue;
    if(line.find("SYN")   == std::string::npos &&
       line.find("TRANS") == std::string::npos &&
       line.find("INV")   == std::string::npos
      )
    {
      continue;
    }
    /*
      chr06_hap21 624290  627316  -   -   chr06_hap22 219806  222828  SYN1    -   SYN -
      chr06_hap21 624290  627316  -   -   chr06_hap22 219806  222828  SYNAL1  SYN1    SYNAL   -
      chr06_hap21 624297  624297  C   A   chr06_hap22 219813  219813  SNP2038 SYN1    SNP -
      chr06_hap21 624299  624299  A   G   chr06_hap22 219815  219815  SNP2039 SYN1    SNP -
      chr06_hap21 624307  624307  T   A   chr06_hap22 219823  219823  SNP2040 SYN1    SNP -
    */
    vector<string> lineinfo = split_string(line, '\t');
    if(lineinfo[8].find("AL") != std::string::npos) continue;
    if(lineinfo.size() < 11)
    {
      cout << "   warning: skipping insufficient line " << line << endl;
      continue;
    }
    if( (lineinfo[8].find("SYN")!=std::string::npos   && lineinfo[10].compare("SYN")  ==0) ||
        (lineinfo[8].find("TRANS")!=std::string::npos && lineinfo[10].compare("TRANS")==0) ||
        (lineinfo[8].find("INVTR")!=std::string::npos && lineinfo[10].compare("INVTR")==0) ||
        (lineinfo[8].find("INV")!=std::string::npos   && lineinfo[10].compare("INV")  ==0)
      )
    {
      // cout << "   check: " << lineinfo[8] << ": " << line << endl;
      string this_syn = lineinfo[8];
      VAR tmpvar;
      tmpvar.rlg  = refsample + "_" + lineinfo[0];
      tmpvar.rsta = strtoul(lineinfo[1].c_str(),NULL,0);
      tmpvar.rend = strtoul(lineinfo[2].c_str(),NULL,0);
      tmpvar.qlg  = qrysample + "_" + lineinfo[5];
      tmpvar.qsta = strtoul(lineinfo[6].c_str(),NULL,0);
      tmpvar.qend = strtoul(lineinfo[7].c_str(),NULL,0);
      syn_var.insert(std::pair<string, VAR>(this_syn, tmpvar));
      syn_order.push_back(this_syn);
    }
  }
  cout << "   info: " << syn_var.size() << " syntenic regions collected. " << endl;
  ifp.close();
  // 2nd round: read variations
  ifp.open(so_file.c_str(), ios::in);
  while(ifp.good())
  {
    string line("");
    getline(ifp, line);
    if(line.size()==0 || line[0]=='#') continue;
    if(line.find("SYN")   == std::string::npos &&
       line.find("TRANS") == std::string::npos &&
       line.find("INV")   == std::string::npos
      )
    {
      continue;
    }
    vector<string> lineinfo = split_string(line, '\t');
    if(line.find("AL") != std::string::npos) continue;
    if(lineinfo.size() < 11)
    {
      cout << "   warning: skipping insufficient line " << line << endl;
      continue;
    }
    if(lineinfo[10].compare("SNP")!=0 &&
       lineinfo[10].compare("INS")!=0 &&
       lineinfo[10].compare("DEL")!=0)
    {
      continue;
    }
    // cout << line << endl;
    if(report_indels.compare("0")==0 && lineinfo[10].compare("SNP")!=0) continue; // skipping indels
    string this_syn = lineinfo[9];
    // cout << this_syn << endl;
    assert( syn_var.find(this_syn) != syn_var.end() );
    string var_sta = lineinfo[1];
    string var_end = lineinfo[2];
    syn_var[this_syn].var_list.insert(std::pair<string, string>(var_sta, var_end));
  }
  ifp.close();
  // output number of variations for each syn/rearranged region into the directory
  string ofilename = out_dir + "/" + out_prefix + "_var_density_in_paired_regions.txt"; // density of variations
  ofstream ofp;
  ofp.open(ofilename.c_str(), ios::out);
  if(!ofp.good())
  {
    cout << "   Error: cannot open file " << ofilename << " to write var density. " << endl;
    return false;
  }
  string ofilename2 = out_dir + "/" + out_prefix + "_var_density_defined_IBD.txt"; // potential IBD regions
  ofstream ofp2;
  ofp2.open(ofilename2.c_str(), ios::out);
  if(!ofp2.good())
  {
    cout << "   Error: cannot open file " << ofilename2 << " to write var defined IBD. " << endl;
    return false;
  }
  string var_type="snp";
  if(report_indels.compare("1")==0)
  {
    var_type="sdi";
  }
  vector<string>::iterator ositr;
  vector<string>::iterator ositr_end;
  ositr     = syn_order.begin();
  ositr_end = syn_order.end();
  while(ositr != ositr_end)
  {
    string this_syn = *ositr;
    map<string, VAR>::iterator sitr = syn_var.find(this_syn);
    VAR    this_var = (*sitr).second;
    //
    double region_size        = this_var.rend - this_var.rsta + 1;     //
    double region_var_num     = this_var.var_list.size();
    double region_var_density = region_size/(region_var_num+0.000001); // 1 snp per x bp
    double var_unit = 1;
    if(region_var_num == 0 )
    {
      var_unit           = 0;
      region_var_density = region_size;
    } // out_flags
    string region_label="SYN";
    if(this_syn.find("SYN")==std::string::npos)
    {
        region_label="SV";
    }
    ofp  << this_syn
         << " " << this_var.rlg
         << " " << this_var.rsta  << " "  << this_var.rend
         << " " << this_var.qlg
         << " " << this_var.qsta  << " "  << this_var.qend
         << " " << region_var_num << " "     << var_type
         << " " << region_size    << " bp: " << var_unit
         << " " << var_type       << " per " << region_var_density    << " bp " << region_label << endl;
    if(region_var_density > IBD_max)
    {
      ofp2 << this_syn
           << " " << this_var.rlg
           << " " << this_var.rsta  << " "  << this_var.rend
           << " " << this_var.qlg
           << " " << this_var.qsta  << " "  << this_var.qend
           << " " << region_var_num << " "     << var_type
           << " " << region_size    << " bp: " << var_unit
           << " " << var_type       << " per " << region_var_density    << " bp " << region_label <<  endl;
    }
    //
    ositr ++;
  }
  ofp.close();
  (*ibd_files).push_back(ofilename2);
  cout << "   info: coordinates of IBDs       collected in " << ofilename2 << endl;
  cout << "   info: statistics  of variations collected in " << ofilename  << endl;
  //
  return true;
}
