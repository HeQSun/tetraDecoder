/* HiC-based binning (for polyploid species) - aiming to achieve haplotype-resolved assembly: 
#
FUNCTION:
   analyzes Hi-C contacts at contig markers
   cluster contigs into p*n groups, where p is ploidy and n is the haploid number (Tetraploid potato: p=4, n=1)
   classify contig markers to linkage groups, i.e., long reads would be classified into linkage groups.	
INPUTT:
   1.read alignments at contigs markers .....................................: local.bam
   2.minimum correlation to build contact graph of haplotig markers..........: cor=0.55
   3.contig sizes............................................................: ctgsize=ctgs.ctgsizes
   4.expected number of linkage groups.......................................: nLG=4
   6.minimum size of haplotig markers to build the backbone of linkage groups: min_hapctg_size=100000
OUTPUT: 
   s1:
   contact matrix of haplotigs...........: final_contig_contact_matrix.txt
   s2:
   inter-haplotig cor edges..............: s3_genotype_haplotig_GT_similarity_matrix.dot
   raw clustering with cor edges.........: s3_genotype_haplotig_GT_similarity_matrix_subclusters_raw.dot
   intra-raw-clustering breaking up......: s3_genotype_haplotig_GT_similarity_matrix_subclusters_raw_advanced.dot
   inter-cluster paired-haplotig cor.....: s3_genotype_haplotig_GT_similarity_matrix_subclusters_corr_check_torm.txt		
   1st-level merging.....................: s3_genotype_haplotig_GT_similarity_matrix_subclusters_raw_advanced_merge.dot
   final clusters after 2nd level merging: s3_genotype_haplotig_GT_similarity_matrix_subclusters_final.dot		
   final homolougous linkage groups......: s3_res_homologous_linkgage_groups.dot
   s3:
   marker with limited Hi-C contact......: s4p1_all_raw_unlinked_markers_with_limited_hic_support.bed
   raw tetraplotig markers not groupped..: s4p2_all_raw_unlinked_tetraplotig_markers.bed
   raw grouped markers...................: s4p3_all_raw_linked_hap_dip_triplotig_markers.bed
   raw replotig markers not groupped.....: s4p4_all_raw_rep_markers.bed
   tmp on re-linking non-grouped-markers.: s4p5_refine_grouping_tmp.txt
   final grouping of all possible markers: s4p6_refine_grouping_final_contig_markers.txt
#
Written by Hequan Sun, MPIPZ (Germany), 2022
Email: sun@mpipz.mpg.de/sunhequan@gmail.com
#
*/
#include      <algorithm>
#include            <map>
#include         <string>
#include         <vector>
#include        <fstream>
#include        <sstream>
#include       <iostream>
#include        <iomanip>  // std::setprecision
#include     <sys/stat.h>
#include       <string.h>
#include       <stdlib.h>
#include       <assert.h>
#include         <math.h>
#include       <dirent.h>
#include "split_string.h"
using namespace std;
//
struct NODE 
{
    unsigned long             end; // window end
    string                   type; // hap/dip/trip/tetrap/rep
    unsigned long        win_size; // size of this marker
    unsigned long        ctg_size; // size of the contig having this marker
    vector<string>           sctg; // id of linked subject contig
    vector<double>          qscor; // contact value linked query and subject contigs
};
//
double min_contact        =  0.50;  // correlation of genotype sequences at two markers defined by sequencing depth 
int    N_cluster          =    4;   // number of expected linkage groups: 1 * iploidy
int    iploidy            =    4;   // ploidy level; default 4 for tetraploid
double haplotig_win_ratio =   0.95; // haplotig win ratio along contigs; used to select "pure"-haplotig markers
unsigned long min_hapctg_size = 50000; // min sizes of contigs to build backbone of chr-groups
//
bool create_folder(string folder);
bool read_target_ctg(string ctg_file, 
                     map<string, unsigned long>* target_contig_size);
bool read_win_marker_ctg(string win_marker_file, 
                         map<string, unsigned long> target_contig_size, 
                         map<string, map<unsigned long, NODE> >* win_marker);
bool categorize_ctg(map<string, map<unsigned long, NODE> >  win_marker, 
                    double haplotig_win_ratio,
                    unsigned long min_hapctg_size,
                    map<string, map<unsigned long, NODE> >* hap_win_marker,
                    map<string, map<unsigned long, NODE> >* other_win_marker);
int decipher_cigar(string cigar, vector<char>* operation, vector<int>* count);                    
string get_R12(int thisflag_int);
bool output_raw_links(map<string, unsigned long> hic_cross_link, 
                      map<string, unsigned long> target_contig_size,
                      string                     tmpfolder,
                      string                     this_label);                    
//
int main(int argc, char* argv[])
{
    if(argc < 7)
    {
        cout << "\nFunction: analyze Hi-C contacts of tig markers, cluster markers and bin long reads. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl
             << "\nUsage: samtools view -f 0x2 -F 3840 x_merged_chrxx_extract.bam"
             << " |"
             << " hic_binning - contact_cutoff target_contig_size.txt window_markers.txt min_hapctg_size out_prefix_str"
             << endl 
             << endl
             << "      -f 0x2              -- only properly paired read alignments. "
             << endl
             << "      -F 3840             -- only primary alignments. "
             << endl
             << "      *contact_cutoff     -- minimum Hi-C contact to build the contact graph of contigs [0.5]. "  
             << endl
             << "      *window_markers.txt -- window markers defined with sequencing coverage. "              
             << endl
             << "      *contig_size.txt    -- list of contig sizes (with two tab-separated columns). " 
             << endl
             << "      *min_hapctg_size    -- minimum contig size to select haplotig markers for linkage clusterings [50000]."                              
             << endl
             << "      *output_prefix      -- flag of output file." << endl 
             << endl;
        cout << "   Note-S3: to filter    dot: dot_filter_v2 xxx.dot 0.75"               << endl;
        cout << "   Note-S3: to visualize dot: circo -Tpdf xxx.dot > xxx_circo.pdf, or " << endl;
        cout << "                                fdp -Tpdf xxx.dot > xxx_fdp.pdf . "     << endl;
        return 1;
    }
    double startT= clock();
    //
    int cii = 0;
    cout << "   CMDline: ";
    while(cii < argc)
    {
        cout << argv[cii] << " ";
        cii ++;
    }
    cout << endl;
    //
    min_contact = atof(argv[2]);
    cout << "   Info: contig contact graph will be built with minimum contact cutoff of: " << min_contact     << endl;
    string ctg_file = (string)argv[3];  
    cout << "   Info: target contig list file: "                                           << ctg_file        << endl;    
    string win_marker_file = string(argv[4]);
    cout << "   Info: coverage-defined window marker file: "                               << win_marker_file << endl;
    min_hapctg_size = strtoul(argv[5], NULL, 0);
    cout << "   Info: min size in bp of haplotigs to build backbone of linkage groups "    << min_hapctg_size << endl; 
    string out_prefix = (string)argv[6];    
    cout << "   Info: output will be labeled with "                                        << out_prefix      << endl;
    // step s1: read target contig and size info
    cout << "   Info: step 1 - read contig size info... " << endl;
    map<string, unsigned long> target_contig_size; // <target_ctg, size>
    if(!read_target_ctg(ctg_file, &target_contig_size))
    {
        cout << "   Error: failed in collecting target contig and size info." << endl;
        return 1;
    }
    cout << "   Info: step 1 - read contig size info done. " << endl;    
    cout << endl;        
    // step s2: read window marker info along target contigs
    cout << "   Info: step 2 - read info on window markers... " << endl;    
    map<string, map<unsigned long, NODE> > win_marker;
    if(!read_win_marker_ctg(win_marker_file, target_contig_size, &win_marker) )
    {
        cout << "   Error: failed in reading coverage-defined window markers. " << endl;
        return 1;
    }
    cout << "   Info: step 2 - read info on window markers done. " << endl;        
    cout << endl;        
    // step s3: select long haplotigs
    cout << "   Info: step 3 - select contigs over " << min_hapctg_size 
         << " bp with "                              << haplotig_win_ratio  
         << " hap-windows. "                         << endl;        
    map<string, map<unsigned long, NODE> > hap_win_marker;
    map<string, map<unsigned long, NODE> > other_win_marker;
    if(!categorize_ctg(win_marker, haplotig_win_ratio, min_hapctg_size, &hap_win_marker, &other_win_marker))
    {
        cout << "   Error: failed in categorizing contigs. " << endl;
        return 1;
    }
    cout << "   Info: step 3- select long haplotigs done. "  << endl;
    cout << endl;        
    // step s4: traverse the bam file and collect Hi-C links
    cout << "   Info: step 4 - traverse bam file from samtools view pipe to get Hi-C link info... " << endl;        
    map<string, unsigned long> ctg_cross_link_long_hap; // main: <ctg1#ctg2, hic_link_read_pair_cnt_on_long_haplotigs>
    map<string, unsigned long> ctg_cross_link_other;    // main: <ctg1#ctg2, hic_link_read_pair_cnt_other_case>    
    unsigned long numraw                  = 0;
    unsigned long num_pair_skip           = 0;             
    unsigned long num_pair_inter_long_hap = 0; 
    unsigned long num_pair_other          = 0;    
    unsigned long lenNoCG  = 0; // no cigar string unmapped - length
    unsigned long numNoCG  = 0; // no cigar string unmapped - number
    string exampleNoCG     = "";             
    std::string read_a; // first  row
    std::string read_b; // second row
    while (std::getline(std::cin, read_a)) 
    {
        if(read_a.size()==0) continue;
        if(read_a[0] == '@') continue;
        numraw += 1;         
        // get second read
        std::getline(std::cin, read_b);
        if(read_b.size()==0)
        {
            cout << "   Warning: you have read_a = " << read_a << " but read_b does not exist. " << endl;
            num_pair_skip ++;
            continue;
        }
        numraw += 1;         
        //
        if(numraw%1000000 == 0)
        {
            cout << "   Info: " << numraw << "th aligm..." << endl;
        }
        //
        vector<string> read_a_info = split_string(read_a, '\t');
        vector<string> read_b_info = split_string(read_b, '\t');
        if(read_a_info.size()<11 || read_b_info.size()<11)
        {
            cout << "   Warning: insufficient read_a info, skipped: " << read_a << endl;
            cout << "   Warning: insufficient read_b info, skipped: " << read_b << endl;
            num_pair_skip ++;
            continue;
        }      
        //
        string readname_a   = read_a_info[0]; // a read name
        string thisflag_a   = read_a_info[1]; // a flag value
        string this_a_cigar = read_a_info[5]; // a cigar string        
        string this_a_seq   = read_a_info[9]; // a read sequence         
        string readname_b   = read_b_info[0]; // b read name
        string thisflag_b   = read_b_info[1]; // b flag value   
        string this_b_cigar = read_b_info[5]; // b cigar string        
        string this_b_seq   = read_b_info[9]; // a read sequence
        int hexflag_a = strtol(thisflag_a.c_str(), NULL, 0); // a flag value 
        int hexflag_b = strtol(thisflag_b.c_str(), NULL, 0); // b flag value    
        // special case 1: cigar string as star: not aligned
        if(this_a_cigar.compare("*")==0 || this_b_cigar.compare("*")==0) 
        {
            //
            if(exampleNoCG.size()==0)
            {
                exampleNoCG = read_a;
                cout << endl
                     << "   Warning: there are paired-alignments without explicit CIGAR string, skipped., e.g.: " 
                     << read_a 
                     << endl
                     << read_b
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;
            }            
            numNoCG += 2; // individual reads; read pair number = numNoCG/2
            lenNoCG  += this_a_seq.size();
            lenNoCG  += this_b_seq.size();
            num_pair_skip ++;
            //
            continue;
        }        
        if(readname_a.compare(readname_b) != 0)
        {
            cout << "   Error: probably unsorted bam provided? Please sort it with read name. " << endl;
            cout << read_a << endl;
            cout << read_b << endl;
            return 1;
        }
        // read a: covered len by alignment
        vector<char>  operation_a;
        vector<int>   count_a;        
        int           covlen_a  = decipher_cigar(this_a_cigar, &operation_a, &count_a); 
        // read b: covered len by alignment
        vector<char>  operation_b;
        vector<int>   count_b;        
        int           covlen_b  = decipher_cigar(this_b_cigar, &operation_b, &count_b);        
        // check if read is R1 or R2
        string read_a_R12 = get_R12( hexflag_a );
        string read_b_R12 = get_R12( hexflag_b );
        if(0)
        cout << "   check: read a: " << read_a_R12 << ", and read b: " << read_b_R12 << endl;
        if(!(read_a_R12.compare("R1")==0 && read_b_R12.compare("R2")==0) &&
           !(read_a_R12.compare("R2")==0 && read_b_R12.compare("R1")==0) )
        {
            cout << "   Warning: unclear read pair info: " << endl;
            cout << "      read a = " << read_a << endl;
            cout << "      read b = " << read_b << endl;
            num_pair_skip ++;
            continue;
        }
        string this_a_ctg   = read_a_info[2];  // a reference contig name
        string this_b_ctg   = read_b_info[2];  // b reference contig name
        //
        if( this_a_ctg.compare(this_b_ctg) != 0 &&
            hap_win_marker.find(this_a_ctg) != hap_win_marker.end() &&
            hap_win_marker.find(this_b_ctg) != hap_win_marker.end()             
          )
        {
            // read pair mapped to different contigs and both in the long haplotig list 
            num_pair_inter_long_hap ++;
            string key1 = this_a_ctg + "#" + this_b_ctg;
            string key2 = this_b_ctg + "#" + this_a_ctg;
            if(ctg_cross_link_long_hap.find(key1) == ctg_cross_link_long_hap.end() &&
               ctg_cross_link_long_hap.find(key2) == ctg_cross_link_long_hap.end() )
            {
                ctg_cross_link_long_hap.insert(std::pair<string, unsigned long>(key1, 1));
            }else
            if(ctg_cross_link_long_hap.find(key1) != ctg_cross_link_long_hap.end() &&
               ctg_cross_link_long_hap.find(key2) == ctg_cross_link_long_hap.end() )            
            {
                ctg_cross_link_long_hap[key1] += 1;
            }else
            if(ctg_cross_link_long_hap.find(key1) == ctg_cross_link_long_hap.end() &&
               ctg_cross_link_long_hap.find(key2) != ctg_cross_link_long_hap.end() )            
            {
                ctg_cross_link_long_hap[key2] += 1;
            }else ;
        }else
        { 
            // read pairs mapped to other contigs
            num_pair_other ++;
            string key1 = this_a_ctg + "#" + this_b_ctg;
            string key2 = this_b_ctg + "#" + this_a_ctg;
            if(ctg_cross_link_other.find(key1) == ctg_cross_link_other.end() &&
               ctg_cross_link_other.find(key2) == ctg_cross_link_other.end() )
            {
                ctg_cross_link_other.insert(std::pair<string, unsigned long>(key1, 1));
            }else
            if(ctg_cross_link_other.find(key1) != ctg_cross_link_other.end() &&
               ctg_cross_link_other.find(key2) == ctg_cross_link_other.end() )            
            {
                ctg_cross_link_other[key1] += 1;
            }else
            if(ctg_cross_link_other.find(key1) == ctg_cross_link_other.end() &&
               ctg_cross_link_other.find(key2) != ctg_cross_link_other.end() )            
            {
                ctg_cross_link_other[key2] += 1;
            }else ;            
        }
    }
    //
    unsigned long num_pair_kept = num_pair_inter_long_hap + num_pair_other;
    cout << "   Info: summary on "                               << numraw                    
         << " read alignment lines checked, among these, "       << endl
         << "         "                                          << num_pair_skip            
         << " read pairs skipped (/aligned to other contigs), "  << endl
         << "         "                                          << num_pair_kept            
         << " collected into given group: "                      << endl
         << "         "                                          << num_pair_inter_long_hap  
         << " R1-R2 aligned to two long haplotigs, over "        << endl
         << "         "                                          << min_hapctg_size 
         << " bp with "                                          << haplotig_win_ratio
         << " as haplotig windows; "                             << endl
         << "         "                                          << num_pair_other 
         << " R1-R2 aligned to other cases. "                    << endl;
    cout << "   Info: step 4 - traverse bam file to get Hi-C link  done. " << endl;
    //
    cout << "   Info: step 5 - output raw Hi-C link info (for checking purpose)... " << endl;    
    string tmpfolder_s5 = "s5_"+out_prefix+"_raw_tig_marker_cross_link_count";
    if(!create_folder(tmpfolder_s5))
    {
        return 1;
    } 
    if(!output_raw_links(ctg_cross_link_other,
                         target_contig_size,
                         tmpfolder_s5, 
                         "other") )
    {
        cout << "   Error: cannot output raw Hi-C link info for other cases. " << endl;
        return 1;
    }  
    if(!output_raw_links(ctg_cross_link_long_hap,
                         target_contig_size,
                         tmpfolder_s5, 
                         "hap") )
    {
        cout << "   Error: cannot output raw Hi-C link info for long haplotigs. " << endl;
        return 1;
    }          
    cout << "   Info: step 5 - output raw Hi-C link info (for checking purpose) done. " << endl;        
    cout << endl;    
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    //
    return 0;
}
//
bool output_raw_links(map<string, unsigned long> hic_cross_link, 
                      map<string, unsigned long> target_contig_size,
                      string                     tmpfolder,
                      string                     this_label)
{
    /* functon: output raw Hi-C links for checking purpose */
    if(hic_cross_link.size() == 0)
    {
        cout << "   Error: size of hic_cross_link is 0 - no cross-linking found between contigs? " << endl;
        return false;
    }
    string siminfo = tmpfolder + "/s5_cross_link_" + this_label + ".dot\0"; // 
    ofstream sofp;
    sofp.open(siminfo.c_str(), ios::out);
    if(!sofp.good())
    {
        cout   << "   Error: cannot open file " << siminfo << endl;
        return false;
    }
    // un-directed graph
    sofp << "/* Hi-C contact graph of contigs */" << endl;
    sofp << "graph\tGraph_1 {" << endl;    
    map<string, unsigned long>::iterator litr;
    map<string, unsigned long>::iterator litr_end;
    litr     = hic_cross_link.begin();
    litr_end = hic_cross_link.end();
    while(litr != litr_end)
    {
        string this_key = (*litr).first;
        unsigned long this_cnt = (*litr).second;
        vector<string> keyinfo = split_string(this_key, '#');
        assert(target_contig_size.find(keyinfo[0]) != target_contig_size.end() );
        assert(target_contig_size.find(keyinfo[1]) != target_contig_size.end() );        
        sofp << "\t"
             << keyinfo[0] << " -- " << keyinfo[1] << " "            
             << "[color=gold, penwidth=1, arrowsize=1, label=" << this_cnt << "];" 
             << " /* "
             << target_contig_size[ keyinfo[0] ]
             << ", "
             << target_contig_size[ keyinfo[1] ]
             << " */"  
             << endl;    
        //        
        litr ++;
    }
    sofp << "}" << endl; 
    //
    sofp.close();    
    return true;
}
//
bool categorize_ctg(map<string, map<unsigned long, NODE> >  win_marker, 
                    double haplotig_win_ratio,
                    unsigned long min_hapctg_size,
                    map<string, map<unsigned long, NODE> >* hap_win_marker,
                    map<string, map<unsigned long, NODE> >* other_win_marker)
{
    /* function: select long "pure" haplotigs: haplotig_win_ratio are hap-windows */
    if(win_marker.size() == 0)
    {
        cout << "   Error: win_marker is unexpectedly empty." << endl;
        return false;
    }
    unsigned long total_ctg_size     = 0;    
    unsigned long total_ctg_size_hap = 0;
    map<string, map<unsigned long, NODE> >::iterator raw_ctg_itr;
    map<string, map<unsigned long, NODE> >::iterator raw_ctg_itr_end;
    raw_ctg_itr     = win_marker.begin();
    raw_ctg_itr_end = win_marker.end();
    while(raw_ctg_itr != raw_ctg_itr_end)
    {
        string this_ctg = (*raw_ctg_itr).first;
        map<unsigned long, NODE> tmp_node = (*raw_ctg_itr).second;
        map<unsigned long, NODE>::iterator win_itr;
        map<unsigned long, NODE>::iterator win_itr_end;
        win_itr     = tmp_node.begin();
        win_itr_end = tmp_node.end();
        unsigned long hap_size = 0;
        unsigned long this_ctg_size = 0;
        while(win_itr != win_itr_end)
        {
            string this_type = (*win_itr).second.type;
            if(this_type.compare("hap") == 0)
            {
                hap_size += (*win_itr).second.win_size;
            }
            if(this_ctg_size == 0)
            {
                this_ctg_size = (*win_itr).second.ctg_size;
            }
            win_itr ++;
        }
        assert(this_ctg_size > 0);
        total_ctg_size += this_ctg_size;
        double this_hap_ratio = hap_size*1.0 / this_ctg_size;
        cout << "   check: " << this_ctg << " with hap_size = " << hap_size << ", ratio = " << this_hap_ratio << endl;        
        if(this_hap_ratio>=haplotig_win_ratio && this_ctg_size>=min_hapctg_size)
        {
            (*hap_win_marker).insert(std::pair<string, map<unsigned long, NODE> >(this_ctg, tmp_node));
            total_ctg_size_hap += this_ctg_size;
        }else
        {
            (*other_win_marker).insert(std::pair<string, map<unsigned long, NODE> >(this_ctg, tmp_node));            
        }
        //
        raw_ctg_itr ++;
    }
    //
    cout << "   Info: "                          << win_marker.size() 
         << " contigs with a total size of "     << total_ctg_size 
         << " bp checked, among these: "         << endl;
    cout << "         "                          << (*hap_win_marker).size()  
         << " contigs with a total size of "     << total_ctg_size_hap
         << " bp selected as over "              << min_hapctg_size
         << " bp and with a hap win ratio over " << haplotig_win_ratio
         << ". "                                 << endl;
    //
    return true;
}                    
//
bool read_win_marker_ctg(string win_marker_file, 
                         map<string, unsigned long> target_contig_size, 
                         map<string, map<unsigned long, NODE> >* win_marker)
{
    /*
        target_contig_size..... <target_ctg_id, size>
        win_marker............. <ctg, <win_sta, {win_end, hap/dip/trip/tetrap/rep, size_ctg, size_window, ...}> >
        NODE = {
		    unsigned long             end; // window end
		    string                   type; // hap/dip/trip/tetrap/rep
		    unsigned long        win_size; // size of this marker
		    unsigned long        ctg_size; // size of the contig having this marker
		    vector<string>           sctg; // id of linked subject contig
		    vector<double>          qscor; // correlation value linked query and subject contigs        
   	        }
   	Note, here only contigs in given list of target_contig_size will be collected
    */
    ifstream ifp;
    ifp.open(win_marker_file.c_str(), ios::in);
    if(!ifp.good())
    {
        return false;
    }
    int num_win = 0;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        /*
            utg000027lc	1	50000	174.013	156	1551343	hap
	    utg000027lc	50001	100000	174.013	156	1551343	hap
	    utg000027lc	100001	150000	174.013	156	1551343	hap
        */
        vector<string> lineinfo = split_string(line, '\t');
        string this_ctg         = lineinfo[0];
        if(target_contig_size.find(this_ctg) != target_contig_size.end())
        {
            num_win ++;
            // only collect info in given contig list.
            unsigned long this_win_sta  = strtoul(lineinfo[1].c_str(), NULL, 0);
            unsigned long this_win_end  = strtoul(lineinfo[2].c_str(), NULL, 0);
            unsigned long this_ctg_size = strtoul(lineinfo[5].c_str(), NULL, 0);
            string this_type            = lineinfo[6];
            NODE tmp_node;
            tmp_node.end      = this_win_end;
            tmp_node.type     = this_type;
            tmp_node.win_size = this_win_end - this_win_sta + 1;
            tmp_node.ctg_size = this_ctg_size;
            tmp_node.sctg.clear();
            tmp_node.qscor.clear();
            if( (*win_marker).find(this_ctg) == (*win_marker).end() )
            {
                map<unsigned long, NODE> tmp_win;
                tmp_win.insert(std::pair<unsigned long, NODE>(this_win_sta, tmp_node) );
                (*win_marker).insert(std::pair<string, map<unsigned long, NODE> >(this_ctg, tmp_win));
            }else
            {
                (*win_marker)[this_ctg].insert(std::pair<unsigned long, NODE>(this_win_sta, tmp_node));
            }
        }
    }
    ifp.close();
    //
    cout << "   Info: " << target_contig_size.size() << " target contigs given, " << endl;
    cout << "   Info: " << (*win_marker).size()      << " ctgs collected with "   <<  num_win << " windows. " << endl;
    //
    return true;
}
//
bool read_target_ctg(string ctg_file, map<string, unsigned long>* target_contig_size)
{
    /* read target contigs and their size */
    fstream ifp;
    ifp.open(ctg_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open ctg file " << ctg_file << endl;
        return false;
    }
    //
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 2)
        {
            cout << "   Warning: skipping line with insufficient info: " << line << endl;
            continue;
        }
        string this_ctg_id = lineinfo[0];
        unsigned long this_ctg_size = strtoul(lineinfo[1].c_str(), NULL, 0);
        if( (*target_contig_size).find(this_ctg_id) != (*target_contig_size).end() )
        {
            if( (*target_contig_size)[this_ctg_id] != this_ctg_size)
            {
                cout << "   Warning: different sizes with the same contig id, please check: " << this_ctg_id << endl;
            }else
            {
                cout << "   Warning: repeated contig id: " << this_ctg_id << endl;
            }
        }else
        {
            (*target_contig_size).insert(std::pair<string, unsigned long>(this_ctg_id, this_ctg_size));
        }
    }
    //
    return true;
}
//
string get_R12(int thisflag_int)
{
    /* function: check if read is R1 or R2*/
    string read_R12 = "RU";
    if((thisflag_int & 0x40) == 0x40 && (thisflag_int & 0x80) == 0)
    {
        read_R12 = "R1";
    }
    else
    if((thisflag_int & 0x40) == 0    && (thisflag_int & 0x80) == 0x80)
    {
        read_R12 = "R2";
    }
    else
    {
        read_R12 = "RU";
        cout << "   Warning: unclear read info: "    << endl;
        cout << "      flag = "      << thisflag_int << endl;
        cout << "      flag & 0x40 == 0x40 = " << (thisflag_int & 0x40) << endl;
        cout << "      flag & 0x80 == 0    = " << (thisflag_int & 0x80) << endl;            
    } 
    //
    return read_R12;
}
//
int decipher_cigar(string cigar, vector<char>* operation, vector<int>* count)
{
    /*
	https://samtools.github.io/hts-specs/SAMv1.pdf: 
        //
	Op	BAM	Description						Consumes_query	Consumes_reference	
	M	0	alignment match (can be a sequence match or mismatch)	yes		yes	
	I	1	insertion to the reference				yes		no	
	D	2	deletion from the reference				no		yes	
	N	3	skipped region from the reference			no		yes	
	S	4	soft clipping (clipped sequences present in	SEQ)	yes		no	
	H	5	hard clipping (clipped sequences NOT present in	SEQ)	no		no	
	P	6	padding (silent deletion from padded reference)		no		no	
	=	7	sequence match						yes		yes	
	X	8	sequence mismatch					yes		yes  	
	NOTE:
        // CIGAR alphabet: [0-9]+[MIDNSHPX=] and *: 23S23M1D8M1I16M1D29M38H 		
        // Sum of lengths of the MIS=X operations shall equal the length of SEQ in bam 
        // H can only be present as the first and/or last operation.
        // S may only have H operations between them and the ends of the CIGAR string
        // or mRNA-to-genome alignment, an N operation represents an intron.  	  
    */
    char *cstr = (char*)cigar.c_str(); // 
    string numstr("");
    int covlen = 0;
    for (int i=0; i<cigar.size(); i++)
    {
        if(cstr[i]>='0' && cstr[i]<='9')
        {
            numstr += cstr[i];
        }
        else
        if(cstr[i] == 'M' || cstr[i] == 'I' || cstr[i] == 'D' || 
           cstr[i] == 'N' || 
           cstr[i] == 'S' || cstr[i] == 'H' || 
           cstr[i] == 'P' || 
           cstr[i] == 'X' || cstr[i] == '=')
        {
            (*operation).push_back(cstr[i]);    
            (*count).push_back(atoi(numstr.c_str()));
            //                    
            if(cstr[i] == 'M' || cstr[i] == 'D' || cstr[i] == 'N' || cstr[i] == 'X' || cstr[i] == '=')
            {
                covlen += atoi(numstr.c_str()); // positions of ref covered by read
            }
            numstr.clear();
        }
        if(cstr[i] == 'X' || cstr[i] == '=' || cstr[i] == 'P' || cstr[i] == 'N') 
        {            
            cout << "   Warning: you have special operation \'" << cstr[i] << "\' in Cigar: " << cigar << endl;
        }else
        if(cstr[i] == 'H') 
        {            
            cout << "   Warning: hard clipping occurred: \'"    << cstr[i] << "\' in Cigar: " << cigar << endl;
        }        
    }
    //
    return covlen;
}
//
bool create_folder(string folder)
{
    DIR* dir = opendir(folder.c_str());
    if(dir)
    {
        /* Directory exists. */
        closedir(dir);
        return true;
    }
    else if (ENOENT == errno)
    {
        /* Directory does not exist. */
        const int dir_err = mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << folder << endl;
           return false;
        }
    }
    else;
    return true;
}
