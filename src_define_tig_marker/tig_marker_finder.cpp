/* this function, given 

	inputs:
		cnv_winsize10000_step10000_hq_merged_vs_hifi.txt -- need last row as hifi depth - 9 columns
		seq-depth: hap-high, dip-high, trip-high, tetra-high
	
	Output contig-regional-markers from 
		haplotig/diplotig/triplotig/tetraplotig
	
	Note: some contigs may generate multiple regional markers 
	      e.g., ctg_1_50000 -> 
	            ctg_1_50000_mkr_hap_1_20000
	            ctg_1_50000_mkr_tetrap_20001_50000

   Hequan Sun, MPIPZ
   Email: sunhequan@gmail.com/sun@mpipz.mpg.de
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
#include    <math.h>
#include   <iomanip>  // std::setprecision
#include "split_string.h"

string find_win_type(double  winDepth, 
                     double* depCut);
vector<string> find_win_type_v2(double  winDepth, 
                     double* depCut,
                     double  winDepth_hifi,
                     double* hifiDepCut);                     
void update_type_size(double avgHap,
                     string last_type, 
                     unsigned long last_sta,
                     unsigned long last_end,
                     unsigned long last_dep,
                     unsigned long* hapSize, 
                     unsigned long* dipSize, 
                     unsigned long* tripSize,
                     unsigned long* tetrapSize,
                     unsigned long* repSize);
using namespace std;
int main(int argc, char* argv[])
{
    if(argc < 6)
    {
        // g++ tig_marker_finder.cpp split_string.cpp -O3 -o tig_marker_finder_v2
        cout << "\nFunction: define tig markers according to sequencing depth. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")"                         << endl;
        cout << "\nUsage: tig_marker_finder_v2 cnv_winsize10000_step10000_hq_merged_vs_hifi.txt hapDep hifiDep new_winsize output_prefix"
             << endl << endl;
        cout << "         cnv_winsize10000_step10000_hq_merged_vs_hifi.txt is with 9 columns. " << endl;
        cout << "         hapDep  gives average illu depth of being haplotigs; max hap/dip/trip/tetrap will be inferred. "<< endl;
        cout << "         hifiDep gives average hifi depth of being haplotigs; max hap/dip/trip/tetrap will be inferred. "<< endl;        
        cout << "         new_winsize gives the new window size to create markers. "  << endl;
        cout << "         output_prefix gives a flag of output file."         << endl << endl;
        return 1;
    }
    double startT= clock();
    cout << endl;
    cout << "   Info: starting merging hap/dip/trip/tetrap windows based on depth..." << endl;
    unsigned long hapSize[2]    = {0, 0}; // {observed-ctg-size, represented-genome-size}
    unsigned long dipSize[2]    = {0, 0};
    unsigned long tripSize[2]   = {0, 0};
    unsigned long tetrapSize[2] = {0, 0};
    unsigned long repSize[2]    = {0, 0};
    // step 1. get and check input from options
    string wfile = (string)argv[1];
    ifstream ifp;
    ifp.open(wfile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open window-depth file " << wfile << endl;
        return 1;
    }
    // hap=116 , haphigh=182 , diphigh=289 , triphigh=398 , tetrahigh=507 
    string dstr   = (string)argv[2];
    double avgHap = atof(argv[2]);
    double depCut[4];
    for(int ii = 0; ii < 4; ii ++)
    {
        // below is theoretical way
        // depCut[ii] = ceil(avgHap * (ii + 1) + avgHap / 2);
    }
    /* however, for example, as the number of hap-windows >> that of dip-windows, we allowed a bit more on hap.
       similar to other pairwise comparisons. */
    /*
    depCut[0] = 182; depCut[0] = ceil(avgHap * 1 + avgHap/2)+8;
    depCut[1] = 289; depCut[1] = ceil(avgHap * 2 + avgHap/2)-1;
    depCut[2] = 398; depCut[2] = ceil(avgHap * 3 + avgHap/2)-8;
    depCut[3] = 507; depCut[3] = ceil(avgHap * 4 + avgHap/2)-15;
    */    
    // 20200629
    depCut[0] = 182; depCut[0] = ceil(avgHap * 1 + avgHap/2);
    depCut[1] = 289; depCut[1] = ceil(avgHap * 2 + avgHap/2);
    depCut[2] = 398; depCut[2] = ceil(avgHap * 3 + avgHap/2);
    depCut[3] = 507; depCut[3] = ceil(avgHap * 4 + avgHap/2);        
    cout << "   Info: given haplotig average depth of hifi+illu"<< avgHap       << "x:"                     << endl;
    cout << "   Info:    haplotig-markers defined with depth [" << 0            << ", " << depCut[0] << "]" << endl;
    cout << "   Info:    diplotig-markers defined with depth [" << depCut[0]+1  << ","  << depCut[1] << "]" << endl;
    cout << "   Info:   triplotig-markers defined with depth [" << depCut[1]+1  << ","  << depCut[2] << "]" << endl;
    cout << "   Info: tetraplotig-markers defined with depth [" << depCut[2]+1  << ","  << depCut[3] << "]" << endl;
    cout << "   Info:    replotig-markers defined with depth [" << depCut[3]+1  << ","  << "INF"     << "]" << endl;
    // new 20200909
    string hifidstr   = (string)argv[3];
    double avgHifiHap = atof(argv[3]);
    double hifiDepCut[4];   
    // 
    hifiDepCut[0] = ceil(avgHifiHap * 1 + avgHifiHap/2);
    hifiDepCut[1] = ceil(avgHifiHap * 2 + avgHifiHap/2);
    hifiDepCut[2] = ceil(avgHifiHap * 3 + avgHifiHap/2);
    hifiDepCut[3] = ceil(avgHifiHap * 4 + avgHifiHap/2);    
    //    
    cout << "   Info: given haplotig average hifi depth of "    << avgHifiHap       << "x:"                         << endl;
    cout << "   Info:    haplotig-markers defined with depth [" << 0                << ", " << hifiDepCut[0] << "]" << endl;
    cout << "   Info:    diplotig-markers defined with depth [" << hifiDepCut[0]+1  << ","  << hifiDepCut[1] << "]" << endl;
    cout << "   Info:   triplotig-markers defined with depth [" << hifiDepCut[1]+1  << ","  << hifiDepCut[2] << "]" << endl;
    cout << "   Info: tetraplotig-markers defined with depth [" << hifiDepCut[2]+1  << ","  << hifiDepCut[3] << "]" << endl;
    cout << "   Info:    replotig-markers defined with depth [" << hifiDepCut[3]+1  << ","  << "INF"         << "]" << endl;   
    cout << "   Info: these hifi-depth values would work with hifi+illu depth values to determine marker types. "   << endl; 
    //
    unsigned long new_winsize = strtoul(argv[4], NULL, 0);
    cout << "   Info: new window size to create markers is " << new_winsize << " bp. " << endl;
    //
    string oflag = (string)argv[5];
    // step 2. prepare output
    vector<string> wfinfo = split_string(wfile, '/');
    size_t pos            = wfinfo[wfinfo.size()-1].find(".txt");
    string ofilename      = wfinfo[wfinfo.size()-1].substr(0, pos) + "_markers_"+ oflag + ".txt";
    ofstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot prepare output file " << ofilename << endl;
        return 1;
    }
    ofp  << "#header: ctg\tstart\tend\tavg_depth\tnum_win_merged\tctg_size\tmarker_type" << endl;
    cout << "   Info: intermediate output will be collected in file: " << ofilename << endl;
    // step 3. define ctg-markers during traversing the input file
    bool          first_win = true;
    int           count_win = 0;
    string        last_ctg  = "";
    unsigned long last_sta  = 0;
    unsigned long last_end  = 0;
    double        last_dep  = 0;
    unsigned long last_size = 0;
    string        last_type = ""; // haplotig, diplotig, triplotig, tetraplotig.
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 9) continue;
        // get info of this window
        string         this_ctg       = lineinfo[0];
        unsigned long  this_sta       = strtoul(lineinfo[1].c_str(), NULL, 0);
        unsigned long  this_end       = strtoul(lineinfo[2].c_str(), NULL, 0);
        double         this_dep       = round(atof(lineinfo[5].c_str())); // real depth -- non-scaled hifi+illu
        unsigned long  this_size      = strtoul(lineinfo[6].c_str(), NULL, 0);
        double         this_dep_hifi  = round(atof(lineinfo[8].c_str())); // real depth -- non-scaled only hifi        
      //string         this_type      = find_win_type(this_dep, depCut);
        vector<string> this_type      = find_win_type_v2(this_dep, depCut, this_dep_hifi, hifiDepCut); // new 20200909 
        //
        if(first_win)
        {
            count_win = 1;
            last_ctg  = this_ctg;
            last_sta  = this_sta;
            last_end  = this_end;
            last_dep  = this_dep;
            last_size = this_size;
            last_type = this_type[0];
            first_win = false;
            continue;
        }
        //
        if(0)
        cout << "   check: "
             << this_ctg           << "\t"
             << this_sta           << "\t"
             << this_end           << "\t"
             << "illu:"            << this_type[0] << "\t"
             << "hifi:"            << this_type[1]
             << endl;             
        // output or merge
        if( (this_ctg.compare(last_ctg) != 0) ||
            (this_ctg.compare(last_ctg) == 0 && (this_type[0].compare(last_type) != 0 && this_type[1].compare(last_type) != 0) )
          )
        {
            // output
            ofp << last_ctg           << "\t"
                << last_sta           << "\t"
                << last_end           << "\t"
                << last_dep/count_win << "\t"
                << count_win          << "\t"
                << last_size          << "\t"
                << last_type          << endl;
            update_type_size(avgHap,
                             last_type, 
                             last_sta,
                             last_end,
                             last_dep/count_win,
                             hapSize,
                             dipSize,
                             tripSize,
                             tetrapSize,
                             repSize);
            // update completely
            count_win = 1;
            last_ctg  = this_ctg;
            last_sta  = this_sta;
            last_end  = this_end;
            last_dep  = this_dep;
            last_size = this_size;
            last_type = this_type[0];
        }else
        {
            // same contig, same type
            // update partially
            last_end   = this_end;
            last_dep  += this_dep;
            count_win += 1;
        }
    }
    /// output last regional marker
    // output
    ofp << last_ctg           << "\t"
        << last_sta           << "\t"
        << last_end           << "\t"
        << last_dep/count_win << "\t"
        << count_win          << "\t"
        << last_size          << "\t"
        << last_type          << endl;
    update_type_size(avgHap,
                     last_type, 
                     last_sta,
                     last_end,
                     last_dep/count_win,
                     hapSize,
                     dipSize,
                     tripSize,
                     tetrapSize,
                     repSize);
    // close files
    ifp.close();
    ofp.close();
    //
    cout << "   Info: merging hap/dip/trip/tetrap windows based on depth done. " << endl;
    cout << "   Summary: " << endl;
    cout << "         haplotig-marker size    s1=" << *hapSize   << " bp; GS += "   << *hapSize    * 1<< " bp; " << endl
         << "         diplotig-marker size    s2=" << *dipSize   << "  bp; GS += "  << *dipSize    * 2<< "  bp; "<< endl
         << "         triplotig-marker size   s3=" << *tripSize  << "  bp; GS += "  << *tripSize   * 3<< "  bp; "<< endl
         << "         tetraplotig-marker size s4=" << *tetrapSize<< "   bp; GS += " << *tetrapSize * 4<< "  bp; "<< endl
         << "         repeat-marker size      s5=" << *repSize   << "    bp; " << endl
         << "         these lead to inferred genome size: "
         << "s1*1 + s2*2 + s3*3 + s4*4 + s5*(dep/avgHap) = " 
         <<           (*(hapSize+1)) + (*(dipSize+1)) + (*(tripSize+1)) + (*(tetrapSize+1)) + (*(repSize+1)) 
         << " bp." << endl;
    cout << "   Note: due to merged repeat windows, there might be some diff in est compared with R version"
         << " bp, \n         while unique regions are highly consistent. " 
         << endl;
    cout << "         Both estimations can be used. " << endl << endl;
    // recreate markers with new_winsize
    cout << "   Info: creating markers with window size " << new_winsize << " bp..."     << endl;
    ifstream ifp2;
    ifp2.open(ofilename.c_str(), ios::in);
    if(!ifp2.good())
    {
        cout << "   Error: cannot open file " << ofilename << endl;
        return 1;
    }
    // prepare final marker list
    pos                   = ofilename.find(".txt");
    unsigned long tmpsize = new_winsize/1000;
    std::stringstream sso;
    sso.str(""); 
    sso << ofilename.substr(0, pos) << "_wsize" << tmpsize << "kb_final.txt";
    ofstream ofp2;
    ofp2.open(sso.str().c_str(), ios::out);
    if(!ofp2.good())
    {
        cout << "   Error: cannot prepare final marker file " << sso.str() << endl;
        return 1;
    }
    unsigned long total_marker_num  = 0;
    unsigned long total_marker_size = 0;
    cout << "   Info: final output will be collected in file: " << sso.str() << endl;
    while(ifp2.good())
    {
        string line("");
        getline(ifp2, line);
        if(line.size()==0 || line[0]=='#') continue;        
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 7)
        {
            cout << "   Warning: skipped line without enough: " << endl;
            continue;
        }
        string        ctgid     = lineinfo[0];
        unsigned long start     = strtoul(lineinfo[1].c_str(), NULL, 0);
        unsigned long end       = strtoul(lineinfo[2].c_str(), NULL, 0);
        unsigned long this_size = end - start + 1;
        if(this_size <= new_winsize)
        {
            ofp2 << line << endl;
            total_marker_num  += 1;
            total_marker_size += end - start + 1;            
        }else
        {
            // this means, if last window is smaller than new_winsize, merge it to the last second
            int new_win_num = round(this_size*1.0 / new_winsize);
            for(int mi = 0; mi < new_win_num; mi ++)
            {
                unsigned long mkr_sta = start + new_winsize * mi;
                unsigned long mkr_end = start + new_winsize * (mi+1) - 1;
                if(mi == new_win_num - 1)
                {
                    mkr_end = end;
                }
                ofp2 << ctgid   << "\t"
                     << mkr_sta << "\t"
                     << mkr_end << "\t"
                     << lineinfo[3] << "\t"
                     << lineinfo[4] << "\t"
                     << lineinfo[5] << "\t"
                     << lineinfo[6] << endl;
                total_marker_num  += 1;
                total_marker_size += mkr_end - mkr_sta + 1;
            }
        }
    }
    ifp2.close(); 
    ofp2.close();   
    cout << "   Info: " << total_marker_num  << " markers generated, reflecting (assembly) size " 
                        << total_marker_size << " bp. " << endl;
    cout << "   Info: creating markers with window size " << new_winsize << " bp done. " << endl;
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    //
    return 0;
}
/*
#header: ctg	start	end	avg_depth	num_win_merged	ctg_size	marker_type
utg000001l_pilon	1	2264094	115.379	227	2264094	hap
utg000002l_pilon	1	11150000	117.444	1115	11812137	hap
utg000002l_pilon	11150001	11160000	230	1	11812137	dip
utg000002l_pilon	11160001	11812137	117.939	66	11812137	hap
utg000003l_pilon	1	2847836	118.565	285	2847836	hap
utg000004l_pilon	1	414462	112.762	42	414462	hap

*/
//
void update_type_size(double avgHap,
                      string last_type, 
                      unsigned long last_sta,
                      unsigned long last_end,
                      unsigned long last_dep,
                      unsigned long* hapSize, 
                      unsigned long* dipSize, 
                      unsigned long* tripSize,
                      unsigned long* tetrapSize,
                      unsigned long* repSize)
{
    // collect size 
    if(last_type.compare("hap") == 0)
    {
        *(hapSize+0)    += (last_end - last_sta + 1);
        *(hapSize+1)    += (last_end - last_sta + 1)*1;
    }else
    if(last_type.compare("dip") == 0)
    {
        *(dipSize+0)    += (last_end - last_sta + 1);
        *(dipSize+1)    += (last_end - last_sta + 1)*2;
    }else
    if(last_type.compare("trip") == 0)
    {
        *(tripSize+0)   += (last_end - last_sta + 1);
        *(tripSize+1)   += (last_end - last_sta + 1)*3;
    }else
    if(last_type.compare("tetrap") == 0)
    {
        *(tetrapSize+0) += (last_end - last_sta + 1);
        *(tetrapSize+1) += (last_end - last_sta + 1)*4;
    }else
    if(last_type.compare("rep") == 0)
    {
        double r         = round(last_dep/avgHap);
        *(repSize+0)    += (last_end - last_sta + 1);
        *(repSize+1)    += (last_end - last_sta + 1) * r;
    }
}
//
vector<string> find_win_type_v2(double  winDepth, 
                                double* depCut,
                                double  winDepth_hifi,                                
                                double* hifiDepCut)
{
    // determine type of a window according to its illu depth and only hifi depth
    vector<string> type; // haplotig, diplotig, triplotig, tetraplotig.
    type.push_back("type1");
    type.push_back("type2");
    if(winDepth <= depCut[0])
    {
        type[0] = "hap";
    }else
    if(winDepth >= depCut[0] + 1 && winDepth <= depCut[1])
    {
        type[0] = "dip";
    }else
    if(winDepth >= depCut[1] + 1 && winDepth <= depCut[2])
    {
        type[0] = "trip";
    }else
    if(winDepth >= depCut[2] + 1 && winDepth <= depCut[3])
    {
        type[0] = "tetrap";
    }else
    {
        type[0] = "rep"; // repeats
    }
    //
    if(winDepth_hifi <= hifiDepCut[0])
    {
        type[1] = "hap";
    }else
    if(winDepth_hifi >= hifiDepCut[0] + 1 && winDepth_hifi <= hifiDepCut[1])
    {
        type[1] = "dip";
    }else
    if(winDepth_hifi >= hifiDepCut[1] + 1 && winDepth_hifi <= hifiDepCut[2])
    {
        type[1] = "trip";
    }else
    if(winDepth_hifi >= hifiDepCut[2] + 1 && winDepth_hifi <= hifiDepCut[3])
    {
        type[1] = "tetrap";
    }else
    {
        type[1] = "rep"; // repeats
    }    
    //
    return type;
}                        
//
string find_win_type(double winDepth, double* depCut)
{
    // determine type of a window according to its depth 
    string type = ""; // haplotig, diplotig, triplotig, tetraplotig.
    if(winDepth <= depCut[0])
    {
        type = "hap";
    }else
    if(winDepth >= depCut[0] + 1 && winDepth <= depCut[1])
    {
        type = "dip";
    }else
    if(winDepth >= depCut[1] + 1 && winDepth <= depCut[2])
    {
        type = "trip";
    }else
    if(winDepth >= depCut[2] + 1 && winDepth <= depCut[3])
    {
        type = "tetrap";
    }else
    {
        type = "rep"; // repeats
    }
    return (string)type;
}
/*
ctg		wstart	wend	ignore	nor_avg_dep	raw_avg_dep	ctg_size
utg000001l	1	10000	1	62	80.1722	2264093
utg000001l	10001	20000	1	88	114.566	2264093
utg000001l	20001	30000	1	113	147.37	2264093
utg000001l	30001	40000	1	116	151.274	2264093
utg000001l	40001	50000	1	105	136.007	2264093
utg000001l	50001	60000	1	110	142.969	2264093
utg000001l	60001	70000	1	114	147.903	2264093
utg000001l	70001	80000	1	107	138.673	2264093
utg000001l	80001	90000	1	106	138.026	2264093
utg000001l	90001	100000	1	112	145.654	2264093
*/

