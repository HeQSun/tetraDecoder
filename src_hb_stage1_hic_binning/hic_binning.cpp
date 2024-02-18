/* HiC-based haplotype separation for polyploid species - aiming to help achieve haplotype-resolved assembly............
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	FUNCTION:
	   analyzes Hi-C contacts at contig markers
	   cluster contigs into p*n groups, where p is ploidy and n is the haploid number (Tetraploid: p=4, n=1)
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	INPUT:
	   1.read alignments at contigs markers .....................................: local.bam
	   2.minimum hic-contact to build contact graph of haplotig markers..........: contact_cutoff=70
	   3.contig sizes............................................................: target_contig_size.txt
	   4.window markers defined with coverage....................................: window_markers.txt
	   5.minimum size of haplotig markers to build the backbone of linkage groups: min_hapctg_size=500000
	   6.allelic table defined with gmap.........................................: allelic.table.txt
	   7.alignment of focal contig assembly to a haploid reference genome........: align_to_ref_paf_file
	   8.initial haplotig markers to be used for backbone clustering.............: haplotig_win_ratio 
	   9.output labelling on the linkage group...................................: out_prefix_str	   
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	   
	OUTPUT: 
	   s6,s7: these are some intermediate data for the purpose of checking
	       raw clustering with cor edges........: s6_haplotig_hic_contact_matrix_subclusters_raw.dot
	       1st-level merging....................: s6_haplotig_hic_contact_matrix_subclusters_raw_merged.dot
	       extended with long haplotigs.........: s7_haplotig_hic_contact_matrix_subclusters_hap_extended.dot   
	   s8: this is the final result!
	       window markers grouped by Hi-C.......: s8_grouping_window_markers_refined_1st.txt	       
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        IDEA and THOUGHTS: 
           1. Tetraplotig goes to all four; triplotig goes to the opposite three haps of allelic-and-assigned haplotigs
           2. Clean hic matrix: 
              if ctg1 is allelic to ctg2, ctg1 and ctg3 in backbone clusters, 
              the link between ctg2 and ctg3 should be removed!    
           3. How do I make use of hic links among non-haplotigs? 
              E.g, a diplotig highly linked to a triplotig, indicating they have at least one hap in common? 
           4. Inter-group correction? Currently not used! 20221212              
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
   Written by Hequan Sun, MPIPZ (Germany), 2022.........................................................................
   Email: sun@mpipz.mpg.de/sunhequan@gmail.com..........................................................................
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
#include "group_refinement.h"
using namespace std;
//
struct NODE 
{
    unsigned long         end;         // window end
    string               type;         // hap/dip/trip/tetrap/rep
    unsigned long    win_size;         // size of this marker
    unsigned long    ctg_size;         // size of the contig having this marker
    map<int, double> grp_link_val;     // best link between this marker and each hap-extended-group
    map<int, string> grp_link_ctg;     // best ctg linking this marker and each hap-extended-group    
    map<int, double> grp_link_val_sum; // total   ctg linking this marker and each hap-extended-group    
    map<int, double> grp_link_val_cnt; // cnt  of ctg linking this marker and each hap-extended-group    
    map<int, double> grp_link_val_size;// size of ctg linking this marker and each hap-extended-group    
};
/*
struct WHIC 
{
  //unsigned long sta1; // marker-window sta of ctg1
    unsigned long end1; // marker-window end of ctg1
    string        typ1; // type of window of ctg1: hap/dip/...
  //unsigned long sta2; // marker-window sta of ctg2
    unsigned long end2; // marker-window end of ctg2 
    string        typ2; // type of window of ctg2: hap/dip/...
    unsigned long hcnt; // hic-link-cnt between {ctg1:sta1-end1} and {ctg2:sta2-end2}    
};
struct NONALLELIE
{
 // unsigned long sta;    // sta of window marker
    unsigned long end;    // end of window marker
    map<int, double> grp; // <non_allelic_grp_id, Hi-C_score>
};
*/
struct CONTACT
{
    double LL; // left  of ctg1 has contact with left  of ctg2
    double LR; // left  of ctg1 has contact with right of ctg2
    double RL; // right of ctg1 has contact with left  of ctg2
    double RR; // right of ctg1 has contact with right of ctg2
};
struct ENDREGION
{
    unsigned long ctg_left_sta;
    unsigned long ctg_left_end;
    unsigned long ctg_right_sta;
    unsigned long ctg_right_end;    
};
// globals
double min_contact             = 70;     // hic contact of two markers!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double min_contact_i           = 3.3;    // normalized_reads_scale / 30000
int    N_cluster               = 4;      // number of clusters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
int    N_cluster_final         = 4;      // number of expected linkage groups: 1 * iploidy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
int    iploidy                 = 4;      // ploidy level; default 4 for tetraploid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double max_allelic_ratio       = 0.01;   // a ctg is allelic to a group, if it overlaps by this ratio with the group
double haplotig_win_ratio      = 1.0000; // haplotig win ratio along contigs; used to select "pure"-haplotig markers!!!!
unsigned long min_hapctg_size  = 100000; // min size of haplotigs to build backbone of chr-groups!!!!!!!!!!!!!!!!!!!!!!!
unsigned long min_end_region   = 500000; // min size of end-regions of haplotigs for initializing hic contact !!!!!!!!!!
map<string, CONTACT> ctg_cross_link_long_hap2; // main: < ctg1-ctg2, {LL, LR, RL, RR} >!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
map<string, CONTACT> ctg_cross_link_other2;    // main: < ctg1-ctg2, {LL, LR, RL, RR} >!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
map<string, map<unsigned long, NODE> > win_marker;      // all window markers along given ctgs!!!!!!!!!!!!!!!!!!!!!!!!!!
map<string, map<unsigned long, NODE> > hap_win_marker;  // long and pure haplotigs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
map<string, map<unsigned long, NODE> > other_win_marker;// short haplotigs and non-pure/non-haplotigs!!!!!!!!!!!!!!!!!!!
map<string, double> ctg_hap_ratio;       // ratio of haplotig-windows along the contigs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
map<string, map<unsigned long, map<unsigned long, WHIC> > > hic_matrix;//<ctg1-ctg2,<sta1,<sta2,{end1,end2,hic-cnt}> > >
map<string, string> allelic_map;         // allelic between contigs <ctg1-ctg2, "0"> - gmap/gene based infer!!!!!!!!!!!!
map<string, string> allelic_map_seq;     // allelic between contigs <ctg1-ctg2, "0"> - ref-align based infer!!!!!!!!!!!!
map<string, unsigned long> allelic_map_seq_size;// allelic overlap size between contigs <ctg1-ctg2, overlap_size>!!!!!!!
// these four variables are used to locate a region of ref chr that a contig was aligned to.!!!!!!!!!!!!!!!!!!!!!!!!!!!!
unsigned long resolution_region = 20000000;   // 20 Mb regions along chr!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
int                  region_cnt = 50;         // default: 50 regions of 20 Mb along chr!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double   normalized_reads_scale = min_end_region/5.0; // cnt of hi-c links of two contigs will be scaled with the value!
map<unsigned long, unsigned long> chr_region; // divide a chr region into [1,20 Mb],...,[180.000001, 200 Mb]!!!!!!!!!!!!
map<string, int> ctg_at_chr_region; // ctg has been aligned to one of the int=1,.......,int=9 region!!!!!!!!!!!!!!!!!!!!
bool               verbose = false; // if true, output many intermediate info
//
bool create_folder(string                                  folder);
bool read_target_ctg(string                                ctg_file,
                  map<string, unsigned long>*              target_contig_size);
bool read_win_marker_ctg(string                            win_marker_file,
                  map<string, unsigned long>               target_contig_size,
                  map<string, map<unsigned long, NODE> >*  win_marker);
bool categorize_ctg(map<string, map<unsigned long, NODE> > win_marker,
                  double                                   haplotig_win_ratio,
                  unsigned long                            min_hapctg_size,
                  map<string, map<unsigned long, NODE> >*  hap_win_marker,
                  map<string, map<unsigned long, NODE> >*  other_win_marker,
                  map<string, double>*                     ctg_hap_ratio);
bool get_allelic_ctg(string                                allelic_file);
bool get_allelic_ctg_align_based(string                    align_file, 
                  map<string, unsigned long>               target_contig_size);
bool clean_alignment(map<unsigned long, unsigned long>     tmp_align,
                  string                                   ctg,
                  map<unsigned long, unsigned long>*       cleaned_align); // no gap allowed when merging intervals
bool fuzzy_clean_alignment(map<unsigned long, unsigned long>  tmp_align,
                  unsigned long                            this_contig_size,
                  map<unsigned long, unsigned long>*       cleaned_align); // allowed gap... when merging intervals
void print_alignments(map<unsigned long, unsigned long>    tmp_align);
unsigned long calculate_overlap(map<unsigned long, unsigned long> tmp_align_1, 
                  map<unsigned long, unsigned long>        tmp_align_2);   
int locate_ctg_at_chr(map<unsigned long, unsigned long>    tmp_align, 
                  map<unsigned long, unsigned long>        chr_region,
                  string                                   ctg);           // => ctg_at_chr_region<ctg, region0-9>
int decipher_cigar(string cigar, vector<char>* operation,  vector<int>* count);
string get_R12(int                                         thisflag_int);
void update_links(string                                   this_a_ctg,
                  string                                   this_b_ctg,
                  string                                   contact_type, // LL, LR, RL, RR
                  unsigned long*                           n_pair_inter_long_hap,
                  unsigned long*                           n_pair_other );
bool output_raw_links(map<string, CONTACT>                 hic_cross_link,
                  map<string, unsigned long>               target_contig_size,
                  string                                   tmpfolder,
                  string                                   this_label,
                  map<string, ENDREGION>                   ctg_end_region);
bool cluster_hap_links(map<string, CONTACT>                hic_cross_link, 
                  map<string, unsigned long>               target_contig_size,
                  string                                   tmpfolder,
                  string                                   this_label,
                  map<string, ENDREGION>                   ctg_end_region,
                  map<int, vector<int> >*                  merging_record,
                  map<int, int>*                           cluster_id_upd,
                  map<int, map<string, double> >*          cluster_edge_upd,
                  map<int, map<string, int> >*             cluster_vertex_upd,
                  map<int, unsigned long>*                 cluster_vertex_size_upd);      
bool get_clusters(map<string, int>                         vertex,
                  map<string, vector<string> >             edge,
                  map<string, double>                      cor,
                  map<int, map<string, double> >*          cluster_edge,
                  string                                   tmpfolder);
bool merge_clusters_v2(map<string, double>                    cor_any,
                  map<int, map<string, double> >           cluster_edge,
                  map<int, map<string, int> >              cluster_vertex,   
                  map<int, unsigned long>                  cluster_vertex_size, 
                  map<string, unsigned long>               target_contig_size,
                  map<string, string>                      allelic_map, 
                  map<string, string>                      allelic_map_seq,
                  map<string, unsigned long>               allelic_map_seq_size,
                  string                                   tmpfolder,
                  map<int, vector<int> >*                  merging_record,
                  map<int, int>*                           cluster_id_upd,                  
                  map<int, map<string, double> >*          cluster_edge_upd,
                  map<int, map<string, int> >*             cluster_vertex_upd,   
                  map<int, unsigned long>*                 cluster_vertex_size_upd);                  
int check_allelic(map<string, int>                         smallest_cluster_vertex,
                  map<string, int>                         larger_cluster_vertex);
bool read_hic_binning_lg(string                            hic_group_file, 
                  map<string, int>*                        ctg2_hic_group,
                  map<int, map<string, int> >*             hic_group_2ctg);
bool add_vet_to_group(string                               vet, 
                  int                                      hic_grpid_old,
                  int                                      hic_grpid_new,
                  map<int, map<string, int> >*             hic_group_2ctg,
                  map<string, int>*                        ctg2_hic_group);
bool integrate_other_contigs(map<string, unsigned long>    target_contig_size,
                  map<string, ENDREGION>                   ctg_end_region,
                  string                                   tmpfolder,
                  map<int, vector<int> >                   merging_record,                             
                  map<int, int>*                           cluster_id_upd,
                  map<int, map<string, double> >*          cluster_edge_upd,
                  map<int, map<string, int> >*             cluster_vertex_upd,
                  map<int, unsigned long>*                 cluster_vertex_size_upd,
                  map<string, int>*                        ctg2_hic_group,
                  map<int, map<string, int> >*             hic_group_2ctg);
bool integrate_other_contigs_v2(map<string, unsigned long> target_contig_size,
                  map<string, ENDREGION>                   ctg_end_region,
                  string                                   tmpfolder,
                  map<int, vector<int> >                   merging_record,                             
                  map<int, int>*                           cluster_id_upd,
                  map<int, map<string, double> >*          cluster_edge_upd,
                  map<int, map<string, int> >*             cluster_vertex_upd,
                  map<int, unsigned long>*                 cluster_vertex_size_upd,
                  map<string, int>*                        ctg2_hic_group,
                  map<int, map<string, int> >*             hic_group_2ctg);
bool update_hic_matrix(unsigned long                       read1_sta, 
                  unsigned long                            read2_sta, 
                  string                                   ctg1, 
                  string                                   ctg2); // involve global win_marker and hic_matrix
bool write_hic_matrix(string                               tmpfolder);
bool analyze_hic_matrix(map<string, int>                   ctg2_hic_group,
                  map<int, map<string, int> >              hic_group_2ctg,
                  map<string, map<int, map<string, unsigned long> > > ctg_allelic_group,
                  map<string, unsigned long>               target_contig_size,
                  map<string, map<unsigned long, NONALLELIE> >* ctg_win_clear_allelic_group,
                  string                                   out_prefix,
                  string                                   out_folder); // involve global win_marker and hic_matrix
bool define_group_with_allelic(map<string, int>            ctg2_hic_group,
                  map<int, map<string, int> >              hic_group_2ctg,
                  map<string, string>                      allelic_map, 
                  map<string, string>                      allelic_map_seq,
                  map<string, unsigned long>               target_contig_size,
                  map<string, map<int, map<string, unsigned long> > >* ctg_allelic_group);
bool get_options(int                                       argc, 
                  char*                                    argv[], 
                  double*                                  min_contact, 
                  double*                                  min_contact_i,
                  double*                                  haplotig_win_ratio,                  
                  unsigned long*                           min_hapctg_size,                  
                  double*                                  max_allelic_ratio,
                  string*                                  ctg_file,
                  string*                                  win_marker_file,
                  string*                                  allelic_file,
                  string*                                  align_file,
                  string*                                  out_prefix);            
//
int main(int argc, char* argv[])
{
    if(argc < 12)
    {
        cout << "\nFunction: analyze Hi-C contacts among contig markers and phase them for binning long reads. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl
             << "\nUsage: samtools view -F 3840 x_merged_chrxx_extract.bam"
             << " |"
             << " hic_binning [options]"
             << endl 
             << endl
             << "      --min-hic           FLOAT    minimum Hi-C contact to build the contact graph of contigs [70]. "  
             << endl
             << "                                   note, this should be based on depth in Hi-C sequencing. "
             << endl
             << "      --min-hic-i         FLOAT    minimum Hi-C contact to integrate contigs into groups  [3.33]."
             << endl               
             << "      --hap-win-r         FLOAT    a contig with r\% hap-windows will be selected to build initial clusters. "
             << endl  
             << "      --min-hap-size      INT      minimum contig size to select haplotig markers for linkage clustering [100000]."     
             << endl
             << "                                   note, contact at ends of paired-contigs also calculated at this resolution. "                                      
             << endl                        
             << "      --max-allelic-r     FLOAT    ratio of a ctg overlapping a group larger this would be considered as allelic [0.01]."
             << endl          
             << "      --ctg               FILE     list of contig sizes (with two tab-separated columns). " 
             << endl
             << "      --win-marker        FILE     window markers defined with sequencing coverage. "              
             << endl
             << "      --gmap-allelic      FILE     allelic conditions among contigs defined with gmap tool. "
             << endl
             << "      --align-paf         FILE     alignment of contigs to a reference genome. "
             << endl
             << "      -o                  FILE     flag of output file, e.g., chr1, chr2, chr3, ...." 
             << endl 
             << endl;
        cout << "   Note-1: -F 3840 will keep only primary alignments. "                << endl;             
        cout << "   Note-2: to visualize dot: circo -Tpdf xxx.dot > xxx_circo.pdf, or " << endl;
        cout << "                             fdp -Tpdf xxx.dot > xxx_fdp.pdf. "        << endl;
        cout << endl;        
        //   contact_cutoff target_contig_size.txt window_markers.txt min_hapctg_size allelic.table.txt align_file haplotig_win_ratio out_prefix_str
        return 1;
    }
    double startT= clock();
    cout << endl;
    // step s-1: get options
    min_contact;
    haplotig_win_ratio;    
    min_hapctg_size;   
    max_allelic_ratio; 
    string ctg_file ;  
    string win_marker_file;
    string allelic_file;
    string align_file;
    string out_prefix;    
    if(!get_options(argc, argv, 
                         &min_contact, 
                         &min_contact_i,
                         &haplotig_win_ratio,
                         &min_hapctg_size, 
                         &max_allelic_ratio,
                         &ctg_file,
                         &win_marker_file,
                         &allelic_file,
                         &align_file,
                         &out_prefix) )
    {
        cout << "   Info: failed in reading options. " << endl;
        return 1;
    }
    // steo s0. set up chr_region - caution on assuming chr size < 1000 Mb
    for(int ri=0; ri < region_cnt; ri++)
    {
        chr_region.insert(std::pair<unsigned long, unsigned long>(resolution_region*ri+1, resolution_region*(ri+1)));
    }
    // step s1: read target contig and size info
    cout << endl;
    cout << "   Info: step 1 - read contig size info... " << endl;
    map<string, unsigned long> target_contig_size; // <target_ctg, size>
    if(!read_target_ctg(ctg_file, &target_contig_size))
    {
        cout << "   Error: failed in collecting target contig and size info." << endl;
        return 1;
    }
    cout << "         " << target_contig_size.size() << " contigs mapped to current reference chr. " << endl;
    cout << "   Info: step 1 - read contig size info done. " << endl;    
    cout << endl;        
    // step s2: read window marker info along target contigs
    cout << "   Info: step 2 - read info on window markers... " << endl;
    if(!read_win_marker_ctg(win_marker_file, target_contig_size, &win_marker) )
    {
        cout << "   Error: failed in reading coverage-defined window markers. " << endl;
        return 1;
    }
    cout << "   Info: step 2 - read info on window markers done. " << endl;        
    cout << endl;        
    // step s3: select long haplotigs: set up gloabl variables: hap_win_marker and other_win_marker
    cout << "   Info: step 3 - select contigs over " << min_hapctg_size 
         << " bp with "                              << haplotig_win_ratio * 100 
         << "\% hap-windows. "                       << endl;        
    if(!categorize_ctg(win_marker, 
                       haplotig_win_ratio, 
                       min_hapctg_size, 
                       &hap_win_marker, 
                       &other_win_marker,
                       &ctg_hap_ratio))
    {        
        cout << "   Error: failed in categorizing contigs. " << endl;
        return 1;
    }
    cout << "   Info: step 3 - select long haplotigs done. " << endl;
    cout << endl;            
    // step s4. read allelic contigs 
    cout << "   Info: step 4 - reading allelic relationship of contigs..." << endl; // => allelic_map=<"ctg1-ctg2", "0">
    if(!get_allelic_ctg(allelic_file))
    {        
        cout << "   Error: failed in reading allelic map. "  << endl;
        return 1;
    }
    if(!get_allelic_ctg_align_based(align_file, target_contig_size)) // => allelic_map_seq=<"ctg1-ctg2", "0">
    {
        cout << "   Error: failed in creating alignment based allelic map. " << endl;
        return 1;
    }
    cout << "   Info: step 4 - reading allelic relationship of contigs done." << endl;
    cout << endl;
    // step s5: traverse the bam file and collect Hi-C links
    cout << "   Info: step 5 - traverse bam file from samtools view pipe to get Hi-C link info... " << endl;              
    map<string, ENDREGION> ctg_end_region;              // main: < ctg, {ctg_left_sta/end, ctg_right_sta/end} >
    unsigned long numraw                   = 0;
    unsigned long num_pair_skip            = 0;
    unsigned long num_pair_inter_long_hap  = 0;
    unsigned long num_pair_other           = 0;
    unsigned long num_pair_intra           = 0;
    unsigned long num_pair_allelic         = 0;    
    unsigned long num_pair_not_at_ctg_ends = 0;
    unsigned long lenNoCG  = 0;  // no cigar string unmapped - length
    unsigned long numNoCG  = 0;  // no cigar string unmapped - number
    string exampleNoCG     = "";
    string last_read_name  = ""; // avoid repeating counting of the same read pair on the same pair of contigs.        
    string last_ctg_pair   = "";
    std::string read_a; // first row
    while (std::getline(std::cin, read_a)) 
    {
        if(read_a.size()==0) continue;
        if(read_a[0] == '@') continue;
        numraw += 1;         
        //
        if(numraw%1000000 == 0)
        {
            cout << "   Info: " << numraw << "th aligm..." << endl;
        }
        //
        vector<string> read_a_info = split_string(read_a, '\t');
        if(read_a_info.size()<11)
        {
            cout << "   Warning: insufficient read_a info, skipped: " << read_a << endl;
            num_pair_skip ++;
            continue;
        }      
        //
        string readname_a   = read_a_info[0]; // a read name
        string thisflag_a   = read_a_info[1]; // a flag value
        string this_a_cigar = read_a_info[5]; // a cigar string        
        string this_a_seq   = read_a_info[9]; // a read sequence         
        int hexflag_a = strtol(thisflag_a.c_str(), NULL, 0); // a flag value 
        // special case 1: cigar string as star: not aligned
        if(this_a_cigar.compare("*")==0 ) 
        {
            //
            if(exampleNoCG.size()==0)
            {
                exampleNoCG = read_a;
                cout << endl
                     << "   Warning: there are paired-alignments without explicit CIGAR string, skipped., e.g.: " 
                     << read_a 
                     << endl;
                cout << "   Info: you will get a total number of such alignments in the end of program. " 
                     << endl;
            }            
            numNoCG += 1; // individual reads; read pair number = numNoCG/2
            lenNoCG  += this_a_seq.size();
            num_pair_skip ++;
            continue;
        }        
        // read a: covered len by alignment
        vector<char>  operation_a;
        vector<int>   count_a;        
        int           covlen_a  = decipher_cigar(this_a_cigar, &operation_a, &count_a);  
        // check if read is R1 or R2
        string read_a_R12 = get_R12( hexflag_a );
        if(0)
        cout << "   check: read a: " << read_a_R12 << endl;
        string this_a_ctg   = read_a_info[2];  // a reference contig name
        string this_b_ctg   = read_a_info[6];  // b reference contig name
        string key1 = this_a_ctg + "-" + this_b_ctg;
        string key2 = this_b_ctg + "-" + this_a_ctg;          
        // skip read pairs mapped within a contig.
        if(this_b_ctg.compare("=") == 0)
        {
            num_pair_intra ++;
            num_pair_skip ++;            
            this_b_ctg = this_a_ctg;
            continue;
        }
        // if the contigs show allelic relationship, skip the links - gene-based
        if(allelic_map.find(key1)!=allelic_map.end() || allelic_map.find(key2)!=allelic_map.end() )
        {
            num_pair_allelic ++;
            num_pair_skip ++;
            continue;
        }
        // if the contigs show allelic relationship, skip the links - sequence-alignment-based
        if(allelic_map_seq.find(key1)!=allelic_map_seq.end() || allelic_map_seq.find(key2)!=allelic_map_seq.end() )
        {
            num_pair_allelic ++;
            num_pair_skip ++;
            continue;
        }        
        // only when "this_a_ctg!=this_b_ctg", the followings needed!
        assert( this_a_ctg.compare(this_b_ctg) != 0 );
        unsigned long this_a_pos = strtoul(read_a_info[3].c_str(), NULL, 0); // align pos for regional link checking
        unsigned long this_b_pos = strtoul(read_a_info[7].c_str(), NULL, 0); // align pos for regional link checking
        assert(target_contig_size.find(this_a_ctg) != target_contig_size.end());
        assert(target_contig_size.find(this_b_ctg) != target_contig_size.end());        
        unsigned long length_a_ctg = target_contig_size[this_a_ctg];
        unsigned long length_b_ctg = target_contig_size[this_b_ctg];
        // build up hic contact matrix for all window markers !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(!update_hic_matrix(this_a_pos, this_b_pos, this_a_ctg, this_b_ctg))
        {
            cout << "   Error: failed in adding hi-c links." << endl;
            return 1;
        }
        // left/right region of contig a
        unsigned long left_sta_a_ctg = 1;
        unsigned long left_end_a_ctg = min_end_region;
        if(left_end_a_ctg > length_a_ctg)
        {
            left_end_a_ctg = length_a_ctg;
        }
        unsigned long right_sta_a_ctg = 1;
        unsigned long right_end_a_ctg = length_a_ctg;        
        if(length_a_ctg > min_end_region)
        {
            right_sta_a_ctg = length_a_ctg - min_end_region + 1;
        }
        if(ctg_end_region.find(this_a_ctg) == ctg_end_region.end())
        {
            ENDREGION tmp_endregion;
            tmp_endregion.ctg_left_sta  = left_sta_a_ctg;
            tmp_endregion.ctg_left_end  = left_end_a_ctg;
            tmp_endregion.ctg_right_sta = right_sta_a_ctg;
            tmp_endregion.ctg_right_end = right_end_a_ctg;               
            ctg_end_region.insert(std::pair<string, ENDREGION>(this_a_ctg, tmp_endregion));
        }       
        // left/right region of contig b
        unsigned long left_sta_b_ctg = 1;
        unsigned long left_end_b_ctg = min_end_region;
        if(left_end_b_ctg > length_b_ctg)
        {
            left_end_b_ctg = length_b_ctg;
        }
        unsigned long right_sta_b_ctg = 1;
        unsigned long right_end_b_ctg = length_b_ctg;        
        if(length_b_ctg > min_end_region)
        {
            right_sta_b_ctg = length_b_ctg - min_end_region + 1;
        } 
        if(ctg_end_region.find(this_b_ctg) == ctg_end_region.end())
        {
            ENDREGION tmp_endregion;
            tmp_endregion.ctg_left_sta  = left_sta_b_ctg;
            tmp_endregion.ctg_left_end  = left_end_b_ctg;
            tmp_endregion.ctg_right_sta = right_sta_b_ctg;
            tmp_endregion.ctg_right_end = right_end_b_ctg;               
            ctg_end_region.insert(std::pair<string, ENDREGION>(this_b_ctg, tmp_endregion));
        }      
        // skip redundant link: if R1 checked and R1-R2 refer to same pair of contigs, no need to check R2 again!
        if(last_read_name.compare(readname_a) == 0 && 
           (last_ctg_pair.compare(key1)==0 || last_ctg_pair.compare(key2)==0 )
          )
        {
            num_pair_skip ++;        
            continue;            
        }else
        {
            last_read_name = readname_a;
            last_ctg_pair  = key1;
        }
        // check regional links: this_a_pos this_b_pos 
        // [left_sta_a_ctg, left_end_a_ctg] [right_sta_a_ctg, right_end_a_ctg]
        // [left_sta_b_ctg, left_end_b_ctg] [right_sta_b_ctg, right_end_b_ctg]
        bool collected = false;
        unsigned long inter_ll = 0;
        unsigned long inter_lr = 0;
        unsigned long inter_rl = 0;
        unsigned long inter_rr = 0;
        unsigned long other_ll = 0;
        unsigned long other_lr = 0;
        unsigned long other_rl = 0;
        unsigned long other_rr = 0;
        if( left_sta_a_ctg<=this_a_pos && this_a_pos<=left_end_a_ctg &&
            left_sta_b_ctg<=this_b_pos && this_b_pos<=left_end_b_ctg
          )
        {
            // LL - left region of contig-a contact left region of contig-b 
            string contact_type = "LL";
            update_links(this_a_ctg, 
                         this_b_ctg,
                         contact_type,                   // LL, LR, RL, RR
                         &inter_ll,
                         &other_ll                         
                        );     
            collected = true;       
        }
        //
        if( left_sta_a_ctg<=this_a_pos  && this_a_pos<=left_end_a_ctg &&
            right_sta_b_ctg<=this_b_pos && this_b_pos<=right_end_b_ctg
          )
        {
            // LR - left region of contig-a contact right region of contig-b 
            string contact_type = "LR";
            update_links(this_a_ctg, 
                         this_b_ctg,
                         contact_type,                   // LL, LR, RL, RR
                         &inter_lr,
                         &other_lr                         
                        );     
            collected = true;                                    
        }  
        //
        if( right_sta_a_ctg<=this_a_pos && this_a_pos<=right_end_a_ctg &&
            left_sta_b_ctg<=this_b_pos && this_b_pos<=left_end_b_ctg
          )
        {
            // RL - right region of contig-a contact left region of contig-b 
            string contact_type = "RL";
            update_links(this_a_ctg, 
                         this_b_ctg,
                         contact_type,                   // LL, LR, RL, RR
                         &inter_rl,
                         &other_rl                         
                        ); 
            collected = true;                                        
        }   
        //
        if( right_sta_a_ctg<=this_a_pos && this_a_pos<=right_end_a_ctg &&
            right_sta_b_ctg<=this_b_pos && this_b_pos<=right_end_b_ctg
          )
        {
            // RR - right region of contig-a contact right region of contig-b 
            string contact_type = "RR";
            update_links(this_a_ctg, 
                         this_b_ctg,
                         contact_type,                   // LL, LR, RL, RR
                         &inter_rr,
                         &other_rr
                        ); 
            collected = true;                                        
        }
        if(inter_ll + inter_lr + inter_rl + inter_rr > 0)
        {
            num_pair_inter_long_hap ++;
        }
        if(other_ll + other_lr + other_rl + other_rr > 0)
        {
            num_pair_other ++;
        }
        //
        if(!collected)   
        {
            // not in regions of interests!
            num_pair_not_at_ctg_ends ++;
        }
    }
    //
    unsigned long num_pair_kept = num_pair_inter_long_hap + num_pair_other;
    cout << "   Info: summary on a0="                            << numraw                    
         << " read alignment lines checked, among these, "       << endl
         << "         a1="                                       << num_pair_skip            
         << " skipped (including a1.1="                          << num_pair_intra
         << " intra-contig links and a1.2="                      << num_pair_allelic
         << " allelic links), "                                  << endl
         << "         a2="                                       << num_pair_not_at_ctg_ends
         << " not in end-regions of contigs, "                   << endl
         << "         a3="                                       << num_pair_kept            
         << " collected into given group, among these: "         << endl
         << "            a4="                                    << num_pair_inter_long_hap  
         << " R1-R2 cross-aligned to two long haplotigs, over "  << min_hapctg_size 
         << " bp with "                                          << haplotig_win_ratio * 100 << "%"
         << " as haplotig windows; "                             << endl
         << "            a5="                                    << num_pair_other 
         << " R1-R2 aligned as other cases. "                    << endl;
    cout << "         Note, a0=a1+a2+a3, where a3=a4+a5."        << endl; 
    cout << "   Info: step 5 - traverse bam file to get Hi-C link  done. " << endl;
    cout << endl;  
    // step s6. output raw hic links among haplotigs and other cases   
    cout << "   Info: step 6 - output raw Hi-C link info (for checking purpose)... "    << endl;
    string tmpfolder_s6 = "s6_"+out_prefix+"_raw_tig_marker_cross_link_count";
    if(!create_folder(tmpfolder_s6))
    {
        return 1;
    } 
    if(!output_raw_links(ctg_cross_link_other2,
                         target_contig_size,
                         tmpfolder_s6, 
                         "other",
                         ctg_end_region) )
    {
        cout << "   Error: cannot output raw Hi-C link info for other cases. "          << endl;
        return 1;
    }
    if(!output_raw_links(ctg_cross_link_long_hap2,
                         target_contig_size,
                         tmpfolder_s6, 
                         "hap",
                         ctg_end_region) )
    {
        cout << "   Error: cannot output raw Hi-C link info for long haplotigs. "       << endl;
        return 1;
    }
    cout << "   Info: step 6 - output raw Hi-C link info (for checking purpose) done. " << endl;
    cout << endl;
    // step s7. cluster haplotigs
    cout << "   Info: step 7 - initial clustering of haplotigs using Hi-C link... "     << endl;
    map<int, vector<int> >         merging_record;          // <clusterid_old, <merged-subcluster-ids> >    
    map<int, int>                  cluster_id_upd;          // updated <clusterid_old, clusterid_new>    
    map<int, map<string, double> > cluster_edge_upd;        // <clusterid_old, <ctg1-ctg2, hic.contact.value > >
    map<int, map<string, int> >    cluster_vertex_upd;      // <clusterid_old, <ctg, 1> >  
    map<int, unsigned long>        cluster_vertex_size_upd; // <clusterid_old, total_ctg_length >   
    if(!cluster_hap_links(ctg_cross_link_long_hap2, 
                         target_contig_size,
                         tmpfolder_s6,
                         "hap",
                         ctg_end_region,
                         &merging_record,
                         &cluster_id_upd,
                         &cluster_edge_upd,
                         &cluster_vertex_upd,
                         &cluster_vertex_size_upd ) )
    {
        cout << "   Error: failed in haplotig clustering. " << endl;
        return 1;
    }
    cout << "   Info: step 7 - initial clustering of haplotigs using Hi-C link done. " << endl;
    cout << endl;
    //
    cout << "   Info: step 8 - re-collecting backbone-clusters of haplotigs ..."       << endl;
    cout << "       : re-collect output of function: merge_clusters ..."               << endl;
    string hic_group_file = tmpfolder_s6 + "/s6_haplotig_hic_contact_matrix_subclusters_raw_merged.dot\0";
    map<string, int>            ctg2_hic_group; // <ctg_id, clusterid_old>
    map<int, map<string, int> > hic_group_2ctg; // <clusterid_old, <ctg_id, clusterid_new> >    
    if(!read_hic_binning_lg(hic_group_file, &ctg2_hic_group, &hic_group_2ctg) )
    {
        cout << "   Error: failed in collecting backbone-clusters of haplotigs. "  << endl;
        return 1;
    }
    if(ctg2_hic_group.size() == 0)
    {
        cout << "   Error: no contigs forming clusters with given settings. "      << endl;
        return 1;
    }
    cout << "   Info: step 8 - re-collecting backbone-clusters of haplotigs done." << endl;
    cout << endl;
    //
    cout << "   Info: step 9 - integrating short haplotigs..."                     << endl;
    if(!integrate_other_contigs_v2(target_contig_size,
                                ctg_end_region,
                                tmpfolder_s6,
                                merging_record,                             
                                &cluster_id_upd,
                                &cluster_edge_upd,
                                &cluster_vertex_upd,
                                &cluster_vertex_size_upd,
                                &ctg2_hic_group,
                                &hic_group_2ctg))
    {
        cout << "   Error: failed in integrating short/non-pure haplotigs and non-haplotigs. " << endl;
        return 1;
    }
    cout << "   Info: step 9 - integrating short haplotigs done."                  << endl;
    cout << endl;
    
    cout << "   Info: step 10 - assigning window markers to haplotyp groups..."    << endl;    
    // try using confidently clustered contigs for assigning window markers!
    string new_hic_group_file = tmpfolder_s6 + "/s7_haplotig_hic_contact_matrix_subclusters_hap_extended.dot\0";
    // string new_hic_group_file = tmpfolder_s6 + "/s6_haplotig_hic_contact_matrix_subclusters_raw_merged.dot\0";
    ctg2_hic_group.clear(); // <ctg_id, clusterid_old>
    hic_group_2ctg.clear(); // <clusterid_old, <ctg_id, clusterid_new> >    
    if(!read_hic_binning_lg(new_hic_group_file, &ctg2_hic_group, &hic_group_2ctg) )
    {
        cout << "   Error: failed in collecting extended backbone-clusters of haplotigs. " << endl;
        return 1;
    }                 
    /* <ctg, <allelic_grp, <allelic_ctg, allelic_ctg_size > > >
       ctg has a lower chance to be assigned to allelic_grp, where the prob can be defined by number of allelic ctgs.
    */
    map<string, map<int, map<string, unsigned long> > > ctg_allelic_group; 
    if(!define_group_with_allelic(ctg2_hic_group, 
                                  hic_group_2ctg, 
                                  allelic_map, 
                                  allelic_map_seq, 
                                  target_contig_size,
                                  &ctg_allelic_group))
    {
        // using the info below: 
        // map<string, string> allelic_map;       // allelic relationship between contigs <ctg1-ctg2, "0"> - gmap/gene based infer!
        // map<string, string> allelic_map_seq;   // allelic relationship between contigs <ctg1-ctg2, "0"> - ref-align based infer!
        // ctg2_hic_group.clear(); // <ctg_id, clusterid_old>
        // hic_group_2ctg.clear(); // <clusterid_old, <ctg_id, clusterid_new> >  
        cout << "   Error: failed in defining grouping based on allelic info. " << endl;
        return false; 
    }                              
    // below gives:  "tmpfolder_s6/s8_grouping_window_markers.txt"
    map<string, map<unsigned long, NONALLELIE> > ctg_win_clear_allelic_group;
    if( !analyze_hic_matrix(ctg2_hic_group, 
                            hic_group_2ctg, 
                            ctg_allelic_group, 
                            target_contig_size, 
                            &ctg_win_clear_allelic_group,
                            out_prefix, 
                            tmpfolder_s6) )
    {
        cout << "   Error: failed in assigning markers not in backbone clusters. " << endl;
        return false;
    }
    // check hic matrix
    if(!write_hic_matrix(tmpfolder_s6))
    {
        cout << "   Error: failed in writing out hic matrix. " << endl;
        return 1;
    }
    cout << "   Info: step 10 - assigning window markers to haplotyp groups done." << endl;
    cout << endl;
    //
    cout << "   Info: step 11 - refining hic-based window marker grouping..."      << endl;
    string raw_hb_grouping_file = tmpfolder_s6 + "/s8_grouping_window_markers.txt\0";    
    if(!group_refinement( raw_hb_grouping_file, 
                          tmpfolder_s6, 
                          hic_matrix, 
                          target_contig_size, 
                          ctg_hap_ratio,
                          ctg_win_clear_allelic_group,
                          ctg2_hic_group,
                          hic_group_2ctg,
                          allelic_map, 
                          allelic_map_seq,
                          allelic_map_seq_size,
                          normalized_reads_scale) )
    {
        cout << "   Error: failed in refining hic-based linkage grouping. "        << endl;
        return 1;
    }
    cout << "   Info: step 11 - refining hic-based window marker grouping done."   << endl;    
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    //
    return 0;
}
bool integrate_other_contigs_v2(map<string, unsigned long>   target_contig_size,
                             map<string, ENDREGION>          ctg_end_region,
                             string                          tmpfolder,
                             map<int, vector<int> >          merging_record,                             
                             map<int, int>*                  cluster_id_upd,
                             map<int, map<string, double> >* cluster_edge_upd,
                             map<int, map<string, int> >*    cluster_vertex_upd,
                             map<int, unsigned long>*        cluster_vertex_size_upd,
                             map<string, int>*               ctg2_hic_group,
                             map<int, map<string, int> >*    hic_group_2ctg)
{
    /*
        Function: integrate short haplotigs into current backbone-clusters formed by long haplotigs.
        # local
        target_contig_size       // <ctg_id, size>        
        ctg_end_region           // <ctg_id, {ctg_left_sta/end, ctg_right_sta/end} >        
        cluster_id_upd;          // updated <clusterid_old, clusterid_new>    
        merging_record;          // <clusterid_old, <merged-subcluster-old-ids> >                 
        cluster_edge_upd;        // <clusterid_old, <ctg1-ctg2, hic.contact.value > >
        cluster_vertex_upd;      // <clusterid_old, <ctg_id, 1> >  
        cluster_vertex_size_upd; // <clusterid_old, total_ctg_length >   
        ctg2_hic_group           //                 <ctg_id, clusterid_old>
        hic_group_2ctg           // <clusterid_old, <ctg_id, clusterid_new> > 
        # global
        ctg_cross_link_long_hap2 // < ctg1-ctg2, {LL, LR, RL, RR} >
        ctg_cross_link_other2    // < ctg1-ctg2, {LL, LR, RL, RR} >      
        hap_win_marker           // <ctg_id, <win_sta, {win_end, hap/dip/trip/tetrap/rep, size_ctg, size_window, ...}> >
                                 // hap_win_marker is long and pure haplotigs                                         
        other_win_marker         // includes short or non-pure haplotigs and non-haplotigs
        allelic_map              // allelic relationship between contigs <ctg1-ctg2, "0"> - gmap/gene based inferrence
        allelic_map_seq          // allelic relationship between contigs <ctg1-ctg2, "0"> - ref-align based inferrence        
        allelic_map_seq_size     // allelic overlap size between contigs <ctg1-ctg2, overlap_size>                                 
    */
    ////////////////////////////////////// Part I: integrate long pure haplotigs ///////////////////////////////////////
    for(int roundi = 0; roundi < 1; roundi ++)
    {
        map<string, int> just_inserted; // just inserted would not be used as a bait
        map<string, CONTACT>::iterator hitr;
        map<string, CONTACT>::iterator hitr_end;
        hitr     = ctg_cross_link_long_hap2.begin();
        hitr_end = ctg_cross_link_long_hap2.end();
        while(hitr != hitr_end)
        {
            string this_key = (*hitr).first;
            vector<string> keyinfo = split_string(this_key, '-');
            string ctg_1 = keyinfo[0];
            string ctg_2 = keyinfo[1];
            for(int ci = 0; ci < 2; ci ++)
            {
                string this_ctg = keyinfo[ci];
                if( (*ctg2_hic_group).find(this_ctg) == (*ctg2_hic_group).end() )
                {
                    // this ctg is not in any existing cluster ~ need potential cluster
                    // pp1: traverse existing cluster and find one with the highest hi-c link
                    map<int, double> sum_link_val; // <clusterid_old, sum-link-value>
                    map<int, double> sum_link_cnt; // <clusterid_old, cnt-link-value>                
                    map<int, double> max_link_val; // <clusterid_old, max-link-value>
                    map<int, string> max_link_ctg; // <clusterid_old, max-link-ctg>
                    map<int, int>           allelic_cnt;  // <clusterid_old, cnt-of-known_ctg-being-allelic-to-this_ctg>
                    map<int, unsigned long> allelic_ctg_size;        // <clusterid_old, size-of-ctg-involving-allelic>
                    map<int, unsigned long> allelic_ctg_size_overlap;// <clusterid_old, size-of-overlap-involving-allelic>                    
                    map<int, double>        allelic_ctg_size_overlap_ratio;// <clusterid_old, ratio-of-overlap-involving-allelic>                                        
                    //
                    map<int, map<string, int> >::iterator citr;
                    map<int, map<string, int> >::iterator citr_end;
                    citr     = (*cluster_vertex_upd).begin();
                    citr_end = (*cluster_vertex_upd).end();
                    while(citr != citr_end)
                    {
                        int this_old_group_id           = (*citr).first;
                        map<string, int> this_old_group = (*citr).second;
                        // initialize 
                        if(allelic_cnt.find( this_old_group_id ) == allelic_cnt.end() )
                        {
                            allelic_cnt.insert(  std::pair<int, int>   (this_old_group_id, 0));
                            max_link_val.insert( std::pair<int, double>(this_old_group_id, 0));
                            max_link_ctg.insert( std::pair<int, string>(this_old_group_id, ""));
                            sum_link_val.insert( std::pair<int, double>(this_old_group_id, 0));
                            sum_link_cnt.insert( std::pair<int, double>(this_old_group_id, 0));
                            allelic_ctg_size.insert( std::pair<int, unsigned long>(this_old_group_id, 0));
                            allelic_ctg_size_overlap.insert( std::pair<int, unsigned long>(this_old_group_id, 0));
                            allelic_ctg_size_overlap_ratio.insert(std::pair<int, double>(this_old_group_id, 0));
                        }
                        // pp2: traverse each vertex in this cluster 
                        map<string, int>::iterator vitr;
                        map<string, int>::iterator vitr_end;
                        vitr     = this_old_group.begin();
                        vitr_end = this_old_group.end();
                        while(vitr != vitr_end)
                        {                        
                            string known_ctg = (*vitr).first;
                            /* this makes clustering worse, turned off 2021102
                            if( just_inserted.find(known_ctg) != just_inserted.end() )
                            {
                                // a contig inserted into the group just now should not be used as a bait
                                vitr ++;
                                continue;
                            }
                            */
                            string key1      = this_ctg  + "-" + known_ctg;
                            string key2      = known_ctg + "-" + this_ctg;
                            // check allelic condition - same effect if turned on, turned off for now 20221102.
                            if(allelic_map.find(key1)     != allelic_map.end()     || 
                               allelic_map.find(key2)     != allelic_map.end()     ||
                               allelic_map_seq.find(key1) != allelic_map_seq.end() || 
                               allelic_map_seq.find(key2) != allelic_map_seq.end()  
                               )
                            {
                                allelic_cnt[ this_old_group_id ] += 1;
                                // vitr ++;
                                // continue; 
                                // calculate allelic ctg overlap size 
                                if(allelic_map_seq_size.find(key1) != allelic_map_seq_size.end())
                                {
                                    allelic_ctg_size_overlap[this_old_group_id] += allelic_map_seq_size[key1];
                                    allelic_ctg_size[this_old_group_id] += (target_contig_size[this_ctg] + target_contig_size[known_ctg]);
                                }else
                                if(allelic_map_seq_size.find(key2) != allelic_map_seq_size.end())
                                {
                                    allelic_ctg_size_overlap[this_old_group_id] += allelic_map_seq_size[key2];
                                    allelic_ctg_size[this_old_group_id] += (target_contig_size[this_ctg] + target_contig_size[known_ctg]);                                    
                                }else ;   
                                // calculate allelic ctg size                                                             
                            }
                            // check hic link
                            ////////////////////////////////////////////////////////////////////////////////////////////      
                            map<unsigned long, map<unsigned long, WHIC> > ctg_pair_hic;
                            string key_found = key1;                        
                            if( hic_matrix.find(key1) != hic_matrix.end() )
                            {
                                ctg_pair_hic = hic_matrix[key1];
                            }else
                            if( hic_matrix.find(key2) != hic_matrix.end() )
                            {
                                ctg_pair_hic = hic_matrix[key2];    
                                key_found    = key2;                    
                            }else ;
                            //
                            map<unsigned long, map<unsigned long, WHIC> >::iterator hitr;
                            map<unsigned long, map<unsigned long, WHIC> >::iterator hitr_end;
                            hitr     = ctg_pair_hic.begin();
                            hitr_end = ctg_pair_hic.end();
                            while(hitr != hitr_end)
                            {
                                unsigned long sta1 = (*hitr).first;
                                map<unsigned long, WHIC> tmp_hic_val = (*hitr).second;
                                //
                                map<unsigned long, WHIC>::iterator whitr;
                                map<unsigned long, WHIC>::iterator whitr_end;
                                whitr     = tmp_hic_val.begin();
                                whitr_end = tmp_hic_val.end();
                                while(whitr != whitr_end)
                                {   
                                    unsigned long sta2 = (*whitr).first;                                                         
                                    WHIC tmp_whic = (*whitr).second;
                                    unsigned long ctg1_win_size = tmp_whic.end1 - sta1 + 1;
                                    unsigned long ctg2_win_size = tmp_whic.end2 - sta2 + 1;
                                    double normalized_reads = tmp_whic.hcnt*1.0 / (ctg1_win_size+ctg2_win_size) * normalized_reads_scale;                                
                                    sum_link_val[this_old_group_id] += normalized_reads;
                                    sum_link_cnt[this_old_group_id] ++;
                                    //
                                    if( normalized_reads >=  max_link_val[ this_old_group_id ])
                                    {
                                        max_link_val[ this_old_group_id ] = normalized_reads;
                                        max_link_ctg[ this_old_group_id ] = known_ctg;
                                    }
                                    //
                                    whitr ++;
                                }
                                //
                                hitr ++;
                            }                                                                                                                                                
                            ////////////////////////////////////////////////////////////////////////////////////////////
                            vitr ++;
                        }
                        //
                        if(allelic_ctg_size[this_old_group_id] == 0) allelic_ctg_size[this_old_group_id] = 1;
                        allelic_ctg_size_overlap_ratio[ this_old_group_id ] = 
                        allelic_ctg_size_overlap[this_old_group_id] * 2.0 / allelic_ctg_size[this_old_group_id];
                        //
                        citr ++;
                    }
                    // determine best group to insert 
                    int    best_old_group_id  = -1;
                    double best_max_link_val  =  0;
                    string best_max_link_ctg  = "";
                    map<int, double>::iterator lvitr;
                    map<int, double>::iterator lvitr_end;
                    lvitr     = max_link_val.begin();
                    lvitr_end = max_link_val.end();
                    while(lvitr != lvitr_end)
                    {
                        int this_old_group_id = (*lvitr).first;
                        // TODO need an option on threshold of "5"
                        if( allelic_cnt[ this_old_group_id ] > 5 && allelic_ctg_size_overlap_ratio[ this_old_group_id ]> max_allelic_ratio) 
                        {
                            /*
                            cout << "       : ---- found "    << allelic_cnt[ this_old_group_id ] 
                                 << " allelic cases between " << this_ctg 
                                 << " and cluster "           << this_old_group_id
                                 << endl;
                            */
                            lvitr ++;
                            continue;
                        }
                        //
                        if(sum_link_cnt[this_old_group_id] == 0) sum_link_cnt[this_old_group_id] ++;
                        double ctg_grp_avg_hic = sum_link_val[this_old_group_id] / sum_link_cnt[this_old_group_id];
                        /*
                        cout << "   check: round " << roundi 
                             << ": ctg "           << this_ctg 
                             << " linking grp "    << this_old_group_id 
                             << " with avg hic: "  << ctg_grp_avg_hic
                             << ", with max hic: " << max_link_val[ this_old_group_id ]
                             << " to ctg "         << max_link_ctg[ this_old_group_id ]
                             << endl;
                        */
                        if( ctg_grp_avg_hic > best_max_link_val)
                        {
                            best_max_link_val = ctg_grp_avg_hic;
                            best_max_link_ctg = max_link_ctg[ this_old_group_id ];
                            best_old_group_id = this_old_group_id;                        
                        }
                        lvitr ++;
                    }
                    // CAUTION on setting here! TODO setting up option                
                    // if(best_max_link_val > normalized_reads_scale / 30000 )
                    if(best_max_link_val > min_contact_i)
                    {
                        // update edge, vertex at best group
                        cout << "       : " << this_ctg << " will be inserted into cluster " << best_old_group_id << endl;
                        //
                        just_inserted.insert(std::pair<string, int>(this_ctg, 1));
                        (*ctg2_hic_group).insert(std::pair<string, int>(this_ctg, best_old_group_id));
                        // update edge info 
                        cout << "       : --------- cluster " 
                             << best_old_group_id 
                             << " with " 
                             << (*cluster_edge_upd)[ best_old_group_id ].size() 
                             << " edges updated as ";
                        string best_connect_edge = this_ctg + "-" + best_max_link_ctg;
                        // add the "bridging" edge
                        (*cluster_edge_upd)[ best_old_group_id ].insert(std::pair<string, double>(best_connect_edge, best_max_link_val));
                        cout << (*cluster_edge_upd)[ best_old_group_id ].size() 
                             << " edges. " 
                             << endl;
                        // update vertex info 
                        cout << "       : --------- cluster " 
                             <<  best_old_group_id  
                             << " with " 
                             << (*cluster_vertex_upd)[ best_old_group_id ].size() 
                             << " vertices updated as ";            
                        (*cluster_vertex_upd)[ best_old_group_id ].insert(std::pair<string, int>(this_ctg, 1));
                        cout << (*cluster_vertex_upd)[ best_old_group_id ].size() 
                             << " vertices. " 
                             << endl;    
                        // update cluster size info 
                        cout << "       : --------- cluster " 
                             <<  best_old_group_id  
                             << " with size of " 
                             << (*cluster_vertex_size_upd)[ best_old_group_id ] 
                             << " bp updated as ";              
                        assert( target_contig_size.find(this_ctg) != target_contig_size.end() );
                        (*cluster_vertex_size_upd)[ best_old_group_id ] += target_contig_size[ this_ctg ];
                        cout << (*cluster_vertex_size_upd)[ best_old_group_id ] 
                             << " bp. " 
                             << endl;                                                                                
                    }else
                    {
                        // cout << "       : warning on " << this_ctg << " cannot be inserted into backbone-clusters."<<endl;
                    }                     
                }
            }
            //
            hitr ++;
        } 
    }          
    // output subclusters: this is final, merged from raw subclusters thus equals expected number of numbers.
    string subclusterinfo = tmpfolder + "/s7_haplotig_hic_contact_matrix_subclusters_hap_extended.dot\0"; // 
    ofstream subcluster_ofp;
    subcluster_ofp.open(subclusterinfo.c_str(), ios::out);
    if(!subcluster_ofp.good())
    {
        cout   << "   Error: cannot open file " << subclusterinfo << endl;
        return false;
    } 
    subcluster_ofp << "/* Here are the merged subclusters of contigs */" << endl;
    subcluster_ofp << "graph\tGraph_1 {"                                 << endl;
    // to collect all clusters-specific vertex
    map<int, map<string, double> >::iterator clitr;
    map<int, map<string, double> >::iterator clitr_end;
    clitr     = (*cluster_edge_upd).begin();
    clitr_end = (*cluster_edge_upd).end();
    unsigned long total_contig_size = 0;
    unsigned long total_contig_numb = 0;
    int count_ci = 0;
    while(clitr != clitr_end)
    {
        count_ci ++;
        subcluster_ofp << "\tsubgraph cluster_" << (*clitr).first << " {" << endl;
        // check merging record
        map<int, vector<int> >::iterator mergeitr;
        mergeitr = merging_record.find( (*clitr).first );
        if(mergeitr == merging_record.end() )
        {
            subcluster_ofp << "\t/* no merging related to this cluster */ " << endl;
        }else
        {
            vector<int> merged_from = merging_record[ (*clitr).first ];
            subcluster_ofp << "\t/* merged with subclusters: ";
            for(int jj=0; jj<merged_from.size(); jj++)
            {
                if(jj>0)
                {
                    subcluster_ofp << ", ";
                }
                subcluster_ofp << merged_from[jj];
            }
            subcluster_ofp << " */" << endl;
        }
        //
        unsigned long cluster_ctg_size = 0;
        map<string, int> subcluster_vertex;
        //
        map<string, double>::iterator eitr;
        map<string, double>::iterator eitr_end;
        eitr     = (*clitr).second.begin();
        eitr_end = (*clitr).second.end();
        while(eitr != eitr_end)
        {
            vector<string> edgeinfo = split_string((*eitr).first, '-');
            subcluster_ofp << "\t"
            << edgeinfo[0] << " -- " << edgeinfo[1]
            << " "            << "[color=gold, fontcolor=gold, penwidth=1, label=" << (*eitr).second << "];" 
            << " /* cluster " << (*clitr).first << " */"
            << endl;  
            //
            if(subcluster_vertex.find(edgeinfo[0]) == subcluster_vertex.end() )
            {
                subcluster_vertex.insert(std::pair<string, int>(edgeinfo[0], 1));
                assert(target_contig_size.find(edgeinfo[0]) != target_contig_size.end());
                cluster_ctg_size += target_contig_size[edgeinfo[0]];
            }
            if(subcluster_vertex.find(edgeinfo[1]) == subcluster_vertex.end() )
            {
                subcluster_vertex.insert(std::pair<string, int>(edgeinfo[1], 1));
                assert(target_contig_size.find(edgeinfo[1]) != target_contig_size.end());
                cluster_ctg_size += target_contig_size[edgeinfo[1]];
            }
            //            
            eitr ++;
        }
        // output nodes/contigs 
        map<string, int>::iterator vvitr;
        map<string, int>::iterator vvitr_end;
        vvitr     = subcluster_vertex.begin();
        vvitr_end = subcluster_vertex.end();
        while(vvitr != vvitr_end)
        {
            if(ctg_at_chr_region.find((*vvitr).first ) != ctg_at_chr_region.end())
            {
                subcluster_ofp << "\t" 
                           << (*vvitr).first 
                           << " [color=orangered, style=filled, fillcolor=orangered, fontcolor=white];"
                           << " /* " 
                           << ctg_at_chr_region[ (*vvitr).first ]
                           << " */"
                           << endl;
            }else            
            {
                subcluster_ofp << "\t" 
                           << (*vvitr).first 
                           << " [color=orangered, style=filled, fillcolor=orangered, fontcolor=white];"
                           << " /* " 
                           << " -1 "
                           << " */"
                           << endl;                
            }            
            vvitr ++;
        }
        // sub cluster id 
        (*cluster_id_upd).insert(std::pair<int, int>((*clitr).first, count_ci)); // <clusterid_old, clusterid_new> 
        subcluster_ofp << "\tlabel=\"G_" << count_ci << "\";" << endl;
        subcluster_ofp << "\tfontsize=90;"            << endl; // Make title stand out by giving a large font size
        subcluster_ofp << "\tfontcolor=orangered;"    << endl;
        subcluster_ofp << "\tcolor=gray;"             << endl;
        //
        subcluster_ofp << "\t/* "          << subcluster_vertex.size() << " contigs with total size of " 
                       << cluster_ctg_size << " bp */"              << endl;
        subcluster_ofp << "\t}" << endl;    
        // collect info for later finding "homologous haplotype-specific LGs"
        //
        total_contig_numb += subcluster_vertex.size();
        total_contig_size += cluster_ctg_size;  
        //
        clitr ++;
    }
    //
    subcluster_ofp << "}" << endl;
    subcluster_ofp.close();    
    //
    cout << "       : "           << total_contig_numb 
         << " contigs of "        << total_contig_size 
         << " bp clustered into " << (*cluster_edge_upd).size() 
         << " clusters. "         << endl;               
    //
    return true;
}
//
bool cluster_hap_links(map<string, CONTACT>            hic_cross_link, 
                       map<string, unsigned long>      target_contig_size,
                       string                          tmpfolder,
                       string                          this_label,
                       map<string, ENDREGION>          ctg_end_region,
                       map<int, vector<int> >*         merging_record,
                       map<int, int>*                  cluster_id_upd,
                       map<int, map<string, double> >* cluster_edge_upd,
                       map<int, map<string, int> >*    cluster_vertex_upd,
                       map<int, unsigned long>*        cluster_vertex_size_upd)
{
    /* functon: cluster haplotigs using Hi-C links */
    if(hic_cross_link.size() == 0)
    {
        cout << "   Error: size of hic_cross_link is 0 - no cross-linking found between contigs? " << endl;
        return false;
    }
    map<string, int>               vertex_all;          // <ctg, 1>
    map<string, vector<string> >   edge_all;            // <ctg1,     <ctgs-connecting-ctg1 > >
    map<string, double>            cor_high;            // <ctg1-ctg2, hic.contact.value >  - this includes > min_contact
    map<string, double>            cor_any;             // <ctg1-ctg2, hic.contact.value >  - this includes everything      
    map<string, int>               unique_contigs;
    map<int, map<string, double> > cluster_edge;        // <clusterid_old, <ctg1-ctg2, hic.contact.value > >
    map<int, map<string, int> >    cluster_vertex;      // <clusterid_old, <ctg, 1 > >    
    map<int, unsigned long>        cluster_vertex_size; // <clusterid_old, total_ctg_length >   
    // while(cluster_edge.size() < N_cluster)
    if(1)
    {
        cout << "   Info: clustering with minimum hi-c contact of " << min_contact  << endl;
        cout << "         only edges with contact > this value will be kept here. " << endl;
        vertex_all.clear();
        edge_all.clear();
        cor_high.clear();
        cor_any.clear();
        unique_contigs.clear();
        cluster_edge.clear();
        cluster_vertex.clear();
        cluster_vertex_size.clear();        
        // select hic links among haplotigs above a certain threshold: min_contact
        map<string, CONTACT>::iterator litr;
        map<string, CONTACT>::iterator litr_end;
        litr     = hic_cross_link.begin();
        litr_end = hic_cross_link.end();
        while(litr != litr_end)
        {
            string this_key        = (*litr).first;
            CONTACT this_contact   = (*litr).second;
            vector<string> keyinfo = split_string(this_key, '-');
            assert(ctg_end_region.find(keyinfo[0]) != ctg_end_region.end() );
            assert(ctg_end_region.find(keyinfo[1]) != ctg_end_region.end() );
            unsigned long this_cnt = this_contact.LL;        
            string contact_type    = "LL";
            unsigned long ctg1_end_size;     
            unsigned long ctg2_end_size; 
            ctg1_end_size    = ctg_end_region[ keyinfo[0] ].ctg_left_end - ctg_end_region[ keyinfo[0] ].ctg_left_sta + 1;
            ctg2_end_size    = ctg_end_region[ keyinfo[1] ].ctg_left_end - ctg_end_region[ keyinfo[1] ].ctg_left_sta + 1;                  
            if(this_contact.LR > this_cnt)
            {
                this_cnt     = this_contact.LR;
                contact_type = "LR"; 
                ctg1_end_size    = ctg_end_region[ keyinfo[0] ].ctg_left_end  - ctg_end_region[ keyinfo[0] ].ctg_left_sta + 1;
                ctg2_end_size    = ctg_end_region[ keyinfo[1] ].ctg_right_end - ctg_end_region[ keyinfo[1] ].ctg_right_sta + 1; 
            }
            if(this_contact.RL > this_cnt)
            {
                this_cnt     = this_contact.RL;
                contact_type = "RL";    
                ctg1_end_size    = ctg_end_region[ keyinfo[0] ].ctg_right_end - ctg_end_region[ keyinfo[0] ].ctg_right_sta + 1;
                ctg2_end_size    = ctg_end_region[ keyinfo[1] ].ctg_left_end  - ctg_end_region[ keyinfo[1] ].ctg_left_sta + 1;                           
            }
            if(this_contact.RR > this_cnt)
            {
                this_cnt     = this_contact.RR;
                contact_type = "RR";             
                ctg1_end_size    = ctg_end_region[ keyinfo[0] ].ctg_right_end - ctg_end_region[ keyinfo[0] ].ctg_right_sta + 1;
                ctg2_end_size    = ctg_end_region[ keyinfo[1] ].ctg_right_end - ctg_end_region[ keyinfo[1] ].ctg_right_sta + 1;               
            }
            string vet1  = keyinfo[0];
            string vet2  = keyinfo[1];
            string edge1 = vet1 + "-" + vet2;
            string edge2 = vet2 + "-" + vet1;         
            double normalized_reads = this_cnt*1.0 / (ctg1_end_size+ctg2_end_size) * normalized_reads_scale; // reads per 100 kb      
            if(this_label.compare("hap")==0 && normalized_reads > min_contact) // another level of haplotig marker selection
            {
                // collect vertex <ctg, 1>
                if( vertex_all.find(vet1) == vertex_all.end() )
                {
                    vertex_all.insert(std::pair<string, int>(vet1, 1));
                }
                if( vertex_all.find(vet2) == vertex_all.end() )
                {
                    vertex_all.insert(std::pair<string, int>(vet2, 1));
                }  
                // collect edge <ctg, <ctgs> >
                if( edge_all.find(vet1) == edge_all.end() )
                {
                    vector<string> tmpedge;
                    tmpedge.push_back(vet2);
                    edge_all.insert(std::pair<string, vector<string> >(vet1, tmpedge) );
                }else
                {
                    edge_all[vet1].push_back(vet2);
                }
                if(edge_all.find(vet2) == edge_all.end() )
                {
                    vector<string> tmpedge;
                    tmpedge.push_back(vet1);
                    edge_all.insert(std::pair<string, vector<string> >(vet2, tmpedge) );
                }else
                {
                    edge_all[vet2].push_back(vet1);
                }  
                // collect correlation: <ctg1-ctg2, cor>             
                if(cor_high.find(edge1) == cor_high.end() )
                {
                    cor_high.insert(std::pair<string, double>(edge1, normalized_reads) );
                }else
                {
                    // should never happen
                    cout << "   Warning: edge " << edge1 << " seen more than 1 time. " << endl;
                }
            }
            // cor_any
            if(cor_any.find(edge1) == cor_any.end() )
            {
                cor_any.insert(std::pair<string, double>(edge1, normalized_reads) );
            }           
            //
            if( unique_contigs.find(keyinfo[0]) == unique_contigs.end() )    
            {
                unique_contigs.insert(std::pair<string, int>(keyinfo[0], 1));
            }
            if( unique_contigs.find(keyinfo[1]) == unique_contigs.end() )    
            {
                unique_contigs.insert(std::pair<string, int>(keyinfo[1], 1));
            }        
            //        
            litr ++;
        }
        cout << "       : " << unique_contigs.size() << " raw unique contigs checked. " << endl;    
        // ss1: initial clustering of haplotigs        
        if(this_label.compare("hap")==0)
        {
            cout << "       : " << vertex_all.size()     << " vertices selected with "
                 <<                cor_high.size()        << " edges. " 
                 << endl;  
            // get sub-clusters = backbone of potential linkage groups
            if( !get_clusters(vertex_all, edge_all, cor_high, &cluster_edge, tmpfolder) )
            {
                return false;
            }
            // output subclusters: this is raw, can be with more linkage groups than expectd, thus needs further merging.
            string subclusterinfo = tmpfolder + "/s6_haplotig_hic_contact_matrix_subclusters_raw.dot\0"; // 
            ofstream subcluster_ofp;
            subcluster_ofp.open(subclusterinfo.c_str(), ios::out);
            if(!subcluster_ofp.good())
            {
                cout   << "   Error: cannot open file " << subclusterinfo << endl;
                return false;
            } 
            subcluster_ofp << "/* Here are the raw subclusters of selected haplotigs */" << endl;
            subcluster_ofp << "graph\tGraph_1 {"                              << endl;
            // to collect all clusters-specific vertex 
            map<int, map<string, double> >::iterator clitr;
            map<int, map<string, double> >::iterator clitr_end;
            clitr     = cluster_edge.begin();
            clitr_end = cluster_edge.end();
            unsigned long total_contig_size = 0;
            unsigned long total_contig_numb = 0;
            int count_ci = 0;
            while(clitr != clitr_end)
            {
                count_ci ++;
                subcluster_ofp << "\tsubgraph cluster_" << (*clitr).first << " {" << endl;
                unsigned long cluster_ctg_size = 0;    
                map<string, int> subcluster_vertex;
                //
                map<string, double>::iterator eitr;
                map<string, double>::iterator eitr_end;
                eitr     = (*clitr).second.begin();
                eitr_end = (*clitr).second.end();
                while(eitr != eitr_end)
                {
                    vector<string> edgeinfo = split_string((*eitr).first, '-');
                    subcluster_ofp << "\t"
                    << edgeinfo[0] << " -- " << edgeinfo[1]
                    << " "            << "[color=gold, penwidth=1, arrowsize=1, label=" << (*eitr).second << "];" 
                    << " /* cluster " << (*clitr).first << " */"                    
                    << endl;
                    if(subcluster_vertex.find(edgeinfo[0]) == subcluster_vertex.end() )
                    {
                        subcluster_vertex.insert(std::pair<string, int>(edgeinfo[0], 1));
                        assert(target_contig_size.find(edgeinfo[0]) != target_contig_size.end());
                        cluster_ctg_size += target_contig_size[edgeinfo[0]];
                    }
                    if(subcluster_vertex.find(edgeinfo[1]) == subcluster_vertex.end() )
                    {
                        subcluster_vertex.insert(std::pair<string, int>(edgeinfo[1], 1));
                        assert(target_contig_size.find(edgeinfo[1]) != target_contig_size.end());
                        cluster_ctg_size += target_contig_size[edgeinfo[1]];
                    }         
                    eitr ++;
                }
                // output nodes/contigs 
                map<string, int>::iterator vvitr;
                map<string, int>::iterator vvitr_end;
                vvitr     = subcluster_vertex.begin();
                vvitr_end = subcluster_vertex.end();
                while(vvitr != vvitr_end)
                {
                    if(ctg_at_chr_region.find((*vvitr).first ) != ctg_at_chr_region.end())
                    {
                        int max_region_i = ctg_at_chr_region[ (*vvitr).first ];
                        subcluster_ofp << "\t"
                                   << (*vvitr).first 
                                   << " [color=orangered, style=filled, fillcolor=orangered, fontcolor=white];"
                                   << " /* " 
                                   << max_region_i
                                   << " : "
                                   << (max_region_i+0)*resolution_region+1
                                   << "~"                   
                                   << (max_region_i+1)*resolution_region
                                   << " */"
                                   << endl;   
                    }
                    else            
                    {
                        subcluster_ofp << "\t" 
                                   << (*vvitr).first 
                                   << " [color=orangered, style=filled, fillcolor=orangered, fontcolor=white];"
                                   << " /* " 
                                   << " -1 "
                                   << " */"
                                   << endl;                
                    }                                
                    vvitr ++;
                }
                // sub cluster id 
                subcluster_ofp << "\tlabel=\"G_" << count_ci << "\";" << endl;
                subcluster_ofp << "\tfontsize=90;"            << endl; // Make title stand out by giving a large font size
                subcluster_ofp << "\tfontcolor=orangered;"    << endl;
                subcluster_ofp << "\tcolor=gray;"             << endl;                
                //
                subcluster_ofp << "\t/* "          << subcluster_vertex.size() << " contigs with total size of " 
                               << cluster_ctg_size << " bp */"              << endl;
                subcluster_ofp << "\t}" << endl;     
                //
                total_contig_numb += subcluster_vertex.size();
                total_contig_size += cluster_ctg_size;  
                //
                cluster_vertex.insert(std::pair<int, map<string, int> >( (*clitr).first, subcluster_vertex) );
                cluster_vertex_size.insert(std::pair<int, unsigned long>( (*clitr).first, cluster_ctg_size) );
                //
                clitr ++;
            }
            //
            subcluster_ofp << "}" << endl;
            subcluster_ofp.close();    
            //
            cout << "   Info: "           << total_contig_numb 
                 << " contigs of "        << total_contig_size 
                 << " bp clustered into " << cluster_edge.size() 
                 << " clusters. "         << endl;          
        }
      //if(cluster_edge.size() < N_cluster)
        if(0)
        {
            min_contact += 5;
            cout << "   Info: less clusters than expected number of "     << N_cluster   << endl
                 << "         re-clustering with increased hi-c contact " << min_contact << endl
                 << " by 5 would be performed. " << endl;
        }
    }
    cout << "       : " << cor_any.size()  << " edges collected, among these, " << endl
         << "         " << cor_high.size() << " edges showing contact > "       << min_contact << endl;
    //      
    if(!merge_clusters_v2(cor_any,
                       cluster_edge,
                       cluster_vertex,
                       cluster_vertex_size,
                       target_contig_size,
                       allelic_map, 
                       allelic_map_seq,
                       allelic_map_seq_size,
                       tmpfolder,
                       merging_record,                       
                       cluster_id_upd,
                       cluster_edge_upd,
                       cluster_vertex_upd,
                       cluster_vertex_size_upd))
    {
        cout << "   Error: failed in merging clusters. " << endl;
        return false;
    }
    //  
    return true;
}
bool merge_clusters_v2(map<string, double>        cor_any,
                  map<int, map<string, double> >  cluster_edge,
                  map<int, map<string, int> >     cluster_vertex,   
                  map<int, unsigned long>         cluster_vertex_size, 
                  map<string, unsigned long>      target_contig_size,
                  map<string, string>             allelic_map, 
                  map<string, string>             allelic_map_seq,
                  map<string, unsigned long>      allelic_map_seq_size,
                  string                          tmpfolder,
                  map<int, vector<int> >*         merging_record,
                  map<int, int>*                  cluster_id_upd,                  
                  map<int, map<string, double> >* cluster_edge_upd,
                  map<int, map<string, int> >*    cluster_vertex_upd,   
                  map<int, unsigned long>*        cluster_vertex_size_upd)
{
    /*
        Function: find best linked clusters (independent of cluster sizes), merge until getting expected cluste number.
        Inputs:
         map<string, double>             cor_any;                // <ctg1-ctg2, hic.contact.value >    
         map<int, map<string, double> >  cluster_edge;           // <clusterid_old, <ctg1-ctg2, hic.contact.value > >
         map<int, map<string, int> >     cluster_vertex;         // <clusterid_old, <ctg, 1> >    
         map<int, unsigned long>         cluster_vertex_size;    // <clusterid_old, total_ctg_length >               
         map<string, unsigned long>      target_contig_size,     // <ctg_id, size>    
         map<string, string>             allelic_map,            // <ctg1-ctg2, "" >
         map<string, string>             allelic_map_seq,        // <ctg1-ctg2, "" >
         map<string, unsigned long>      allelic_map_seq_size;   // <ctg1-ctg2, overlap_size>
         map<int, vector<int> >*         merging_record          // <clusterid_old, <merged-subcluster-old_ids> >         
         map<int, int>*                  cluster_id_upd          // <clusterid_old, clusterid_new> ....................updated
         map<int, map<string, double> >* cluster_edge_upd        // <clusterid_old, <ctg1-ctg2, hic.contact.value > >..updated
         map<int, map<string, int> >*    cluster_vertex_upd      // <clusterid_old, <ctg, 1> > ........................updated
         map<int, unsigned long>*        cluster_vertex_size_upd // <clusterid_old, total_ctg_length >.................updated
         #
         map<string, map<unsigned long, map<unsigned long, WHIC> > > hic_matrix; < ctg1-ctg2, <sta1, <sta2, {end1,typ1,end2,typ2,hic-cnt} > > >
             WHIC
   	       unsigned long     end1; // marker-window end of ctg1
    	       string            typ1; // type of window of ctg1: hap/dip/...
	       unsigned long     end2; // marker-window end of ctg2 
	       string            typ2; // type of window of ctg2: hap/dip/...
	       unsigned long     hcnt; // hic-link-cnt between {ctg1:sta1-end1} and {ctg2:sta2-end2}
         //
    */
    /* Function: check the hic links between raw clusters. 
               hb_group_2ctg........: <old_group_id, <ctg_id, new_group_id> >
    */
    map<string, string> allelic_map_merged;
    allelic_map_merged.insert(allelic_map.begin(), allelic_map.end() );
    allelic_map_merged.insert(allelic_map_seq.begin(), allelic_map_seq.end() );            
    //
    int last_cluster_num = cluster_vertex.size();
    int this_cluster_num = cluster_vertex.size();    
    while(cluster_vertex.size() > N_cluster)
    {
        cout << "   Info: ---- number of current clusters " << cluster_vertex.size()
             << " > expected cluster number of "            << N_cluster
             << endl;
        int    to_merge_grp1_id     = -1;
        int    to_merge_grp2_id     = -1;
        double to_merge_hic_max     = 0; // larger  => higher priority to merge two groups
        double to_merge_allelic_cnt = 0; // smaller => higher priority to merge two groups
        map<int, map<int, double> > best_linked_hic_val; // <grp1, <grp2, best_hic_val__ctg1-ctg2> >
        map<int, map<int, string> > best_linked_hic_key; // <grp1, <grp2, best_hic_val> >
        //
        map<int, map<int, double> > grp_contact; // <grp1, <grp2, sum_grp12_hic > >
        //
        map<int, map<string, int> >::iterator gitr1;
        map<int, map<string, int> >::iterator gitr_end;
        gitr1    = cluster_vertex.begin();
        gitr_end = cluster_vertex.end();
        while(gitr1 != gitr_end)
        {
            // get cluster 1
            int group1_id = (*gitr1).first;
            // cout << "   check: old group1 id = " << group1_id << endl;
            map<string, int> group1 = (*gitr1).second;
            // initialize result - a
            map<int, double> sum_grp12_hic;
            grp_contact.insert(std::pair<int, map<int, double> >(group1_id, sum_grp12_hic));
            map<int, double> best_grp12_hic;
            best_linked_hic_val.insert(std::pair<int, map<int, double> >(group1_id, best_grp12_hic) );
            map<int, string> best_grp12_hic_key;
            best_linked_hic_key.insert(std::pair<int, map<int, string> >(group1_id, best_grp12_hic_key) );
            //
            map<int, map<string, int> >::iterator gitr2 = gitr1;
            gitr2 ++;        
            while(gitr2 != gitr_end)
            {
                // get cluster 2
                int group2_id = (*gitr2).first;
                // cout << "        : old group2 id = " << group2_id << endl;
                map<string, int> group2 = (*gitr2).second;
                // initialize result - b
                grp_contact[group1_id].insert(std::pair<int, double>(group2_id, 0));
                int individual_cnts = 0;
                best_linked_hic_val[group1_id].insert(std::pair<int, double>(group2_id, 0));
                best_linked_hic_key[group1_id].insert(std::pair<int, string>(group2_id, ""));
                // initialize allelic cnt between the current two groups
                int           allelic_ctg_cnt_group12        = 0;            
                unsigned long allelic_ctg_size               = 0;
                unsigned long allelic_ctg_size_overlap       = 0;
                double        allelic_ctg_size_overlap_ratio = 0; 
                // check conctact between ctg1's and ctg2's
                map<string, int>::iterator citr1;
                map<string, int>::iterator citr1_end;
                citr1     = group1.begin();
                citr1_end = group1.end();
                while(citr1 != citr1_end)
                {
                    string ctg1 = (*citr1).first;
                    //
                    map<string, int>::iterator citr2;
                    map<string, int>::iterator citr2_end;
                    citr2     = group2.begin();
                    citr2_end = group2.end();
                    while(citr2 != citr2_end)
                    {
                        string ctg2 = (*citr2).first;
                        //
                        string key1 = ctg1 + "-" + ctg2;
                        string key2 = ctg2 + "-" + ctg1;          
                        map<unsigned long, map<unsigned long, WHIC> > ctg_pair_hic;
                        string key_found = key1;                        
                        if( hic_matrix.find(key1) != hic_matrix.end() )
                        {
                            ctg_pair_hic = hic_matrix[key1];
                        }else
                        if( hic_matrix.find(key2) != hic_matrix.end() )
                        {
                            ctg_pair_hic = hic_matrix[key2];    
                            key_found    = key2;                    
                        }else ;
                        // check if allelic 
                        if(allelic_map_merged.find(key1) != allelic_map_merged.end() ||
                           allelic_map_merged.find(key2) != allelic_map_merged.end()
                          )
                        {
                            allelic_ctg_cnt_group12 ++;
                        }
                        // check ctg size and overlap size 
                        if(allelic_map_seq_size.find(key1) != allelic_map_seq_size.end() )
                        {
                            allelic_ctg_size_overlap += allelic_map_seq_size[key1]; 
                            allelic_ctg_size += (target_contig_size[ctg1] + target_contig_size[ctg2]);
                        }else
                        if(allelic_map_seq_size.find(key2) != allelic_map_seq_size.end() )
                        {
                            allelic_ctg_size_overlap += allelic_map_seq_size[key1];
                            allelic_ctg_size += (target_contig_size[ctg1] + target_contig_size[ctg2]);                            
                        }
                        //
                        map<unsigned long, map<unsigned long, WHIC> >::iterator hitr;
                        map<unsigned long, map<unsigned long, WHIC> >::iterator hitr_end;
                        hitr     = ctg_pair_hic.begin();
                        hitr_end = ctg_pair_hic.end();
                        while(hitr != hitr_end)
                        {
                            unsigned long sta1 = (*hitr).first;
                            map<unsigned long, WHIC> tmp_hic_val = (*hitr).second;
                            //
                            map<unsigned long, WHIC>::iterator whitr;
                            map<unsigned long, WHIC>::iterator whitr_end;
                            whitr     = tmp_hic_val.begin();
                            whitr_end = tmp_hic_val.end();
                            while(whitr != whitr_end)
                            {   
                                unsigned long sta2 = (*whitr).first;                                                         
                                WHIC tmp_whic = (*whitr).second;
                                unsigned long ctg1_win_size = tmp_whic.end1 - sta1 + 1;
                                unsigned long ctg2_win_size = tmp_whic.end2 - sta2 + 1;
                                double normalized_reads = tmp_whic.hcnt*1.0 / (ctg1_win_size+ctg2_win_size) * normalized_reads_scale;                                
                                grp_contact[group1_id][group2_id] += normalized_reads;
                                individual_cnts ++;
                                //
                                if( normalized_reads >= best_linked_hic_val[group1_id][group2_id])
                                {
                                    best_linked_hic_val[group1_id][group2_id] = normalized_reads;
                                    best_linked_hic_key[group1_id][group2_id] = key_found;
                                }
                                //
                                whitr ++;
                            }
                            //
                            hitr ++;
                        }
                        //
                        citr2 ++;
                    }
                    //
                    citr1 ++;
                } 
                // here 'allelic_ctg_size' collected sizes of two contigs, so the size of 'overlap*2.0' as below
                if(allelic_ctg_size == 0) allelic_ctg_size = 1;
                allelic_ctg_size_overlap_ratio = allelic_ctg_size_overlap * 2.0 / allelic_ctg_size;
                //
                if(individual_cnts == 0) individual_cnts ++;
                cout << "       : hic of group " << group1_id 
                     <<                      "-" << group2_id 
                     <<                     ": " << grp_contact[group1_id][group2_id] 
                     <<            ", avg hic: " << grp_contact[group1_id][group2_id] / individual_cnts
                     <<           ", best hic: " << best_linked_hic_key[group1_id][group2_id]
                     <<                     ": " << best_linked_hic_val[group1_id][group2_id]
                     <<     ", #allelic pairs: " << allelic_ctg_cnt_group12
                     <<         ", overlap bp: " << allelic_ctg_size_overlap
                     <<              ", ratio: " << allelic_ctg_size_overlap_ratio;
                if(allelic_ctg_size_overlap_ratio < max_allelic_ratio)                 
                {
                    cout <<                "." << endl;
                }else
                {
                    // when the number of allelic contigs is large, the two clusters might be allelic groups.
                    cout << " (low priority)." << endl; 
                }
                /* TODO: set on options on threshold
                         maybe with both ratio of allelic pairs among all pairs and 
                                           cnt of allelic pairs.
                */
                if(grp_contact[group1_id][group2_id] / individual_cnts > to_merge_hic_max && 
                   (allelic_ctg_size_overlap_ratio < max_allelic_ratio) 
                  )
                {
                    to_merge_hic_max     = grp_contact[group1_id][group2_id] / individual_cnts;
                    to_merge_grp1_id     = group1_id;
                    to_merge_grp2_id     = group2_id;
                    to_merge_allelic_cnt = allelic_ctg_cnt_group12;
                }
                //
                gitr2 ++;
            }
            //
            gitr1 ++;
        }
        // merge 
        if(to_merge_grp1_id != -1 && to_merge_grp2_id != -1 && to_merge_grp1_id != to_merge_grp2_id)
        {
            cout << "       : ---- cluster "         << to_merge_grp1_id 
                 << " will be merged with cluster "  << to_merge_grp2_id 
                 << ", with avg hic "                << to_merge_hic_max
                 << " and allelic pair cnt "         << to_merge_allelic_cnt
                 << endl;        
            // update edge info 
            cout << "       : --------- cluster " 
                 << to_merge_grp1_id 
                 << " with " 
                 << cluster_edge[to_merge_grp1_id].size() 
                 << " edges updated as ";
            cluster_edge[to_merge_grp1_id].insert(cluster_edge[to_merge_grp2_id].begin(), 
                                                  cluster_edge[to_merge_grp2_id].end());
            // add the "bridging" edge
            cluster_edge[to_merge_grp1_id].insert(std::pair<string, double>(best_linked_hic_key[to_merge_grp1_id][to_merge_grp2_id], 
                                                                            best_linked_hic_val[to_merge_grp1_id][to_merge_grp2_id] ) );
            cout << cluster_edge[to_merge_grp1_id].size() 
                 << " edges. " 
                 << endl;
            // remove edges of smallest group
            cluster_edge.erase(to_merge_grp2_id);
            // update vertex info 
            cout << "       : --------- cluster " 
                 << to_merge_grp1_id 
                 << " with " 
                 << cluster_vertex[to_merge_grp1_id].size() 
                 << " vertices updated as ";            
            cluster_vertex[to_merge_grp1_id].insert(cluster_vertex[ to_merge_grp2_id ].begin(),
                                                    cluster_vertex[ to_merge_grp2_id ].end());
            cout << cluster_vertex[to_merge_grp1_id].size() 
                 << " vertices. " 
                 << endl;    
            // remove vertices of smallest group 
            cluster_vertex.erase(to_merge_grp2_id);    
            // update cluster size info 
            cout << "       : --------- cluster " 
                 << to_merge_grp1_id 
                 << " with size of " 
                 << cluster_vertex_size[to_merge_grp1_id] 
                 << " bp updated as ";              
            cluster_vertex_size[to_merge_grp1_id] += cluster_vertex_size[to_merge_grp2_id];
            cout << cluster_vertex_size[to_merge_grp1_id] 
                 << " bp. " 
                 << endl;
            // remove sizes of smallest group
            cluster_vertex_size.erase(to_merge_grp2_id);                                    
            //
            map<string, int> group2 = cluster_vertex[ to_merge_grp2_id ];
            cluster_vertex[to_merge_grp1_id].insert( group2.begin(), group2.end() ); 
            cluster_vertex.erase( to_merge_grp2_id );
            cout << "   Info: group-"              << to_merge_grp2_id 
                 << " has been merged into group-" << to_merge_grp1_id 
                 << ", and thus removed list."     << endl;
            //
            cluster_vertex_size.erase( to_merge_grp2_id );
            //
            if((*merging_record).find( to_merge_grp1_id ) != (*merging_record).end() && 
               (*merging_record).find( to_merge_grp2_id ) != (*merging_record).end() 
              )
            {
                (*merging_record)[ to_merge_grp1_id ].push_back( to_merge_grp2_id );
                //
                (*merging_record)[ to_merge_grp1_id ].insert((*merging_record)[ to_merge_grp1_id ].begin(),
                                                             (*merging_record)[ to_merge_grp2_id ].begin(),
                                                             (*merging_record)[ to_merge_grp2_id ].end() );
            }else            
            if((*merging_record).find( to_merge_grp1_id ) != (*merging_record).end() && 
               (*merging_record).find( to_merge_grp2_id ) == (*merging_record).end() 
              )
            {
                (*merging_record)[ to_merge_grp1_id ].push_back( to_merge_grp2_id );
            }else
            if((*merging_record).find( to_merge_grp1_id ) == (*merging_record).end() && 
               (*merging_record).find( to_merge_grp2_id ) != (*merging_record).end() 
              )
            {
                vector<int> tmp_merged_grp;
                tmp_merged_grp.push_back( to_merge_grp2_id );
                (*merging_record).insert( std::pair<int, vector<int> >(to_merge_grp1_id, tmp_merged_grp) );
                //
                (*merging_record)[ to_merge_grp1_id ].insert((*merging_record)[ to_merge_grp1_id ].begin(),
                                                             (*merging_record)[ to_merge_grp2_id ].begin(),
                                                             (*merging_record)[ to_merge_grp2_id ].end() );
            }else
            {
                vector<int> tmp_merged_grp;
                tmp_merged_grp.push_back( to_merge_grp2_id );
                (*merging_record).insert( std::pair<int, vector<int> >(to_merge_grp1_id, tmp_merged_grp) );
            }
        }else
        {
            cout << "   warning: cannot find a proper target cluster to merge current clusters, no action." 
                 << endl;
        } 
        //
        this_cluster_num = cluster_vertex.size();  
        if( last_cluster_num > this_cluster_num )
        {
            last_cluster_num = this_cluster_num;
        }else
        {
            cout << "   Warning: cluster number cannot be merged into the expected value of " << N_cluster 
                 << "." << endl;
            break;
        }
    }
    //
    if(1)
    {
        map<int, map<string, int> >::iterator gitr;
        map<int, map<string, int> >::iterator gitr_end;
        gitr     = cluster_vertex.begin();
        gitr_end = cluster_vertex.end();
        while(gitr != gitr_end)
        {
            int this_old_group_id = (*gitr).first;
            map<string, int> this_group = (*gitr).second;
            cout << "       : merged old group id = " << this_old_group_id 
                 << ": "                         << this_group.size() 
                 << " ctgs: "                    << endl;
            map<string, int>::iterator citr;
            map<string, int>::iterator citr_end;
            citr     = this_group.begin();
            citr_end = this_group.end();
            while(citr != citr_end)
            {
                string this_ctg = (*citr).first;
                int    this_new_group_id = (*citr).second;
                cout << "        : " << this_ctg << endl;
                //
                citr ++;
            }
            //
            gitr ++;
        }
    }
    //
    if(cluster_vertex.size() < N_cluster)
    {
        cout << "   warning: number of clusters is only " << cluster_vertex.size() << ", smaller than " << N_cluster
             << "." << endl;
        cout << "            you should increase initial Hi-C score to get more clusters. " << endl;
    }
    cout << "   Info: merging clusters done." << endl;    
    // output subclusters: this is final, merged from raw subclusters thus equals expected number of numbers.
    string subclusterinfo = tmpfolder + "/s6_haplotig_hic_contact_matrix_subclusters_raw_merged.dot\0"; // 
    ofstream subcluster_ofp;
    subcluster_ofp.open(subclusterinfo.c_str(), ios::out);
    if(!subcluster_ofp.good())
    {
        cout   << "   Error: cannot open file " << subclusterinfo << endl;
        return false;
    } 
    subcluster_ofp << "/* Here are the merged subclusters of contigs */" << endl;
    subcluster_ofp << "graph\tGraph_1 {"                           << endl;
    // to collect all clusters-specific vertex
    map<int, map<string, double> >::iterator clitr;
    map<int, map<string, double> >::iterator clitr_end;
    clitr     = cluster_edge.begin();
    clitr_end = cluster_edge.end();
    unsigned long total_contig_size = 0;
    unsigned long total_contig_numb = 0;
    int count_ci = 0;
    while(clitr != clitr_end)
    {
        count_ci ++;
        subcluster_ofp << "\tsubgraph cluster_" << (*clitr).first << " {" << endl;
        // check merging record
        map<int, vector<int> >::iterator mergeitr;
        mergeitr = (*merging_record).find( (*clitr).first );
        if(mergeitr == (*merging_record).end() )
        {
            subcluster_ofp << "\t/* no merging related to this cluster */ " << endl;
        }else
        {
            vector<int> merged_from = (*merging_record)[ (*clitr).first ];
            subcluster_ofp << "\t/* merged with subclusters: ";
            for(int jj=0; jj<merged_from.size(); jj++)
            {
                if(jj>0)
                {
                    subcluster_ofp << ", ";
                }
                subcluster_ofp << merged_from[jj];
            }
            subcluster_ofp << " */" << endl;
        }
        //
        unsigned long cluster_ctg_size = 0;
        map<string, int> subcluster_vertex;
        //
        map<string, double>::iterator eitr;
        map<string, double>::iterator eitr_end;
        eitr     = (*clitr).second.begin();
        eitr_end = (*clitr).second.end();
        while(eitr != eitr_end)
        {
            vector<string> edgeinfo = split_string((*eitr).first, '-');
            subcluster_ofp << "\t"
            << edgeinfo[0] << " -- " << edgeinfo[1]
            << " "            << "[color=gold, fontcolor=gold, penwidth=1, label=" << (*eitr).second << "];" 
            << " /* cluster " << (*clitr).first << " */"
            << endl;  
            //
            if(subcluster_vertex.find(edgeinfo[0]) == subcluster_vertex.end() )
            {
                subcluster_vertex.insert(std::pair<string, int>(edgeinfo[0], 1));
                assert(target_contig_size.find(edgeinfo[0]) != target_contig_size.end());
                cluster_ctg_size += target_contig_size[edgeinfo[0]];
            }
            if(subcluster_vertex.find(edgeinfo[1]) == subcluster_vertex.end() )
            {
                subcluster_vertex.insert(std::pair<string, int>(edgeinfo[1], 1));
                assert(target_contig_size.find(edgeinfo[1]) != target_contig_size.end());
                cluster_ctg_size += target_contig_size[edgeinfo[1]];
            }
            //            
            eitr ++;
        }
        // output nodes/contigs 
        map<string, int>::iterator vvitr;
        map<string, int>::iterator vvitr_end;
        vvitr     = subcluster_vertex.begin();
        vvitr_end = subcluster_vertex.end();
        while(vvitr != vvitr_end)
        {
            if(ctg_at_chr_region.find((*vvitr).first ) != ctg_at_chr_region.end())
            {
                subcluster_ofp << "\t" 
                               << (*vvitr).first 
                               << " [color=orangered, style=filled, fillcolor=orangered, fontcolor=white];"
                               << " /* " 
                               << ctg_at_chr_region[ (*vvitr).first ]
                               << " */"
                               << endl;
            }else
            {
                subcluster_ofp << "\t" 
                               << (*vvitr).first 
                               << " [color=orangered, style=filled, fillcolor=orangered, fontcolor=white];"
                               << " /* " 
                               << -1
                               << " */"
                               << endl;            
            }
            vvitr ++;
        }
        // sub cluster id 
        (*cluster_id_upd).insert(std::pair<int, int>((*clitr).first, count_ci)); // <clusterid_old, clusterid_new> 
        subcluster_ofp << "\tlabel=\"G_" << count_ci << "\";" << endl;
        subcluster_ofp << "\tfontsize=90;"            << endl; // Make title stand out by giving a large font size
        subcluster_ofp << "\tfontcolor=orangered;"    << endl;
        subcluster_ofp << "\tcolor=gray;"             << endl;
        //
        subcluster_ofp << "\t/* "          << subcluster_vertex.size() << " contigs with total size of " 
                       << cluster_ctg_size << " bp */"              << endl;
        subcluster_ofp << "\t}" << endl;    
        // collect info for later finding "homologous haplotype-specific LGs"
        //
        total_contig_numb += subcluster_vertex.size();
        total_contig_size += cluster_ctg_size;  
        //
        clitr ++;
    }
    //
    subcluster_ofp << "}" << endl;
    subcluster_ofp.close();    
    //
    cout << "   Info: "         << total_contig_numb << " contigs of " << total_contig_size << " bp clustered into "
         << cluster_edge.size() << " clusters. "     << endl; 
    // collect result, including (*cluster_id_upd) updated in the process of output clusters
    (*cluster_edge_upd).insert( cluster_edge.begin(), cluster_edge.end() );
    (*cluster_vertex_upd).insert( cluster_vertex.begin(), cluster_vertex.end() );
    (*cluster_vertex_size_upd).insert( cluster_vertex_size.begin(), cluster_vertex_size.end() );    
    //
    return true;
}             
//
int check_allelic(map<string, int> smallest_cluster_vertex, 
                  map<string, int> larger_cluster_vertex)
{
    /* function: check if contigs in small list having allelic relationship with those in large list, with 
       map<string, string> allelic_map --- <ctg1-ctg2, "0">
    */
    int allelic_cnt = 0;
    map<string, int>::iterator sitr;
    map<string, int>::iterator sitr_end;
    sitr     = smallest_cluster_vertex.begin();
    sitr_end = smallest_cluster_vertex.end();    
    while(sitr != sitr_end)
    {
        string ctg_s = (*sitr).first;
        //
        map<string, int>::iterator litr;
        map<string, int>::iterator litr_end;
        litr     = larger_cluster_vertex.begin();
        litr_end = larger_cluster_vertex.end();
        while(litr != litr_end)
        {
            string ctg_l = (*litr).first;
            string key1  = ctg_s + "-" + ctg_l;
            string key2  = ctg_l + "-" + ctg_s;
            if(allelic_map.find(key1)!=allelic_map.end() || allelic_map.find(key2)!=allelic_map.end() ||
               allelic_map_seq.find(key1)!=allelic_map_seq.end() || allelic_map_seq.find(key2)!=allelic_map_seq.end()  )
            {
                allelic_cnt ++;
            }
            //
            litr ++;
        }
        //
        sitr++;
    }    
    return allelic_cnt;
}
//
bool output_raw_links(map<string, CONTACT>       hic_cross_link, 
                      map<string, unsigned long> target_contig_size,
                      string                     tmpfolder,
                      string                     this_label,
                      map<string, ENDREGION>     ctg_end_region)
{
    /* functon: output raw Hi-C links for checking purpose */
    if(hic_cross_link.size() == 0)
    {
        cout << "   Error: size of hic_cross_link is 0 - no cross-linking found between contigs? " << endl;
        return false;
    }
    map<string, int> unique_contigs;
    //
    string siminfo = tmpfolder + "/s6_cross_link_" + this_label + "_raw_full.dot\0"; // 
    ofstream sofp;
    sofp.open(siminfo.c_str(), ios::out);
    if(!sofp.good())
    {
        cout   << "   Error: cannot open file " << siminfo << endl;
        return false;
    }
    // un-directed graph
    sofp << "/* Hi-C contact graph of contigs: label=this_cnt/(ctg1_win_size+ctg2_win_size) * 100000 */" << endl;
    sofp << "graph\tGraph_1 {" << endl;
    map<string, CONTACT>::iterator litr;
    map<string, CONTACT>::iterator litr_end;
    litr     = hic_cross_link.begin();
    litr_end = hic_cross_link.end();
    while(litr != litr_end)
    {
        string this_key        = (*litr).first;
        CONTACT this_contact   = (*litr).second;
        vector<string> keyinfo = split_string(this_key, '-');
        assert(ctg_end_region.find(keyinfo[0]) != ctg_end_region.end() );
        assert(ctg_end_region.find(keyinfo[1]) != ctg_end_region.end() );
        unsigned long this_cnt = this_contact.LL;        
        string contact_type    = "LL";
        unsigned long ctg1_size;     
        unsigned long ctg2_size; 
        ctg1_size    = ctg_end_region[ keyinfo[0] ].ctg_left_end - ctg_end_region[ keyinfo[0] ].ctg_left_sta + 1;
        ctg2_size    = ctg_end_region[ keyinfo[1] ].ctg_left_end - ctg_end_region[ keyinfo[1] ].ctg_left_sta + 1;                  
        if(this_contact.LR > this_cnt)
        {
            this_cnt     = this_contact.LR;
            contact_type = "LR"; 
            ctg1_size    = ctg_end_region[ keyinfo[0] ].ctg_left_end  - ctg_end_region[ keyinfo[0] ].ctg_left_sta + 1;
            ctg2_size    = ctg_end_region[ keyinfo[1] ].ctg_right_end - ctg_end_region[ keyinfo[1] ].ctg_right_sta + 1; 
        }
        if(this_contact.RL > this_cnt)
        {
            this_cnt     = this_contact.RL;
            contact_type = "RL";    
            ctg1_size    = ctg_end_region[ keyinfo[0] ].ctg_right_end - ctg_end_region[ keyinfo[0] ].ctg_right_sta + 1;
            ctg2_size    = ctg_end_region[ keyinfo[1] ].ctg_left_end  - ctg_end_region[ keyinfo[1] ].ctg_left_sta + 1;                           
        }
        if(this_contact.RR > this_cnt)
        {
            this_cnt     = this_contact.RR;
            contact_type = "RR";             
            ctg1_size    = ctg_end_region[ keyinfo[0] ].ctg_right_end - ctg_end_region[ keyinfo[0] ].ctg_right_sta + 1;
            ctg2_size    = ctg_end_region[ keyinfo[1] ].ctg_right_end - ctg_end_region[ keyinfo[1] ].ctg_right_sta + 1;               
        }    
        string vet1  = keyinfo[0];
        string vet2  = keyinfo[1];
        string edge1 = vet1 + "-" + vet2;
        string edge2 = vet2 + "-" + vet1;         
        double normalized_reads = this_cnt*1.0 / (ctg1_size+ctg2_size) * normalized_reads_scale;  // reads per 100 kb  
        sofp << "\t"
             << keyinfo[0] << " -- " << keyinfo[1] << " "            
             << "[color=gold, fontcolor=gold, penwidth=1, arrowsize=1, label=" << normalized_reads << "];" 
             << " /* "
             << contact_type 
             << " "
             << target_contig_size[ keyinfo[0] ]
             << " "
             << target_contig_size[ keyinfo[1] ]             
             << " "
             << ctg1_size
             << " "
             << ctg2_size
             << " "
             << this_cnt
             << " */"  
             << endl;
        //
        if( unique_contigs.find(keyinfo[0]) == unique_contigs.end() )    
        {
            unique_contigs.insert(std::pair<string, int>(keyinfo[0], 1));
        }
        if( unique_contigs.find(keyinfo[1]) == unique_contigs.end() )    
        {
            unique_contigs.insert(std::pair<string, int>(keyinfo[1], 1));
        }        
        //        
        litr ++;
    }
    // TODO: add edge, example:" 	utg000007l_pilon [color=orangered, style=filled, fillcolor=orangered, fontcolor=white];"
    sofp << "}" << endl; 
    sofp.close();  
    cout << "       : " << unique_contigs.size() << " unique contigs in " << siminfo << endl;
    //  
    return true;
}
//
bool get_clusters(map<string, int> vertex,
                  map<string, vector<string> > edge,
                  map<string, double> cor,
                  map<int, map<string, double> >* cluster_edge,
                  string                          tmpfolder
                 )
{
    // function: traverse ctg1-ctg2 graph, connected vertices will be collected as one group = potential LGs.
    //
    // vertex, edge and correlation
    // map<string, int> vertex;           // <ctg, 1>
    // map<string, vector<string> > edge; // <ctg, <ctgs> >
    // map<string, double> cor;           // <ctg1-ctg2, hic.contact.value >
    // map<int, map<string, double> > cluster_edge;   // <clusterid, <ctg1-ctg2, hic.contact.value> >
    map<int, map<string, int> > cluster_vertex; // <clusterid, <ctg, 1> >
    int cluster_id = -1;
    bool check_on  = false;
    map<string, int> visited_vet;
    while(cor.size() > 0)
    {        
        // get an <contig1-contig2, cor>
        map<string, double>::iterator citr = cor.begin();
        string         this_edge    = (*citr).first; // contig1-contig2
        vector<string> this_vetinfo = split_string(this_edge, '-');
        // initialize a cluster id
        cluster_id ++;
        if(check_on) cout << "   check: building cluster " << cluster_id << endl;
        // initialize edge of this cluster
        map<string, double> tmpcor; // <contig1-contig2, cor>
        tmpcor.insert(std::pair<string, double>( (*citr).first, (*citr).second) );
        (*cluster_edge).insert(std::pair<int, map<string, double> >(cluster_id, tmpcor) );
        // initialize vertex of this cluster
        map<string, int> tmpvet;
        tmpvet.insert(std::pair<string, int>(this_vetinfo[0], 1));
        tmpvet.insert(std::pair<string, int>(this_vetinfo[1], 1));
        cluster_vertex.insert(std::pair<int, map<string, int> >(cluster_id, tmpvet));
        // visited vertex
        visited_vet.insert(std::pair<string, int>(this_vetinfo[0], 1));
        visited_vet.insert(std::pair<string, int>(this_vetinfo[1], 1));
        // remove this <contig1-contig2, cor>
        cor.erase(citr);
        //
        vector<string> vet_queue;
        vet_queue.push_back(this_vetinfo[0]);
        vet_queue.push_back(this_vetinfo[1]);
        while(vet_queue.size() > 0)
        {
            string this_vet = *(vet_queue.begin());
            if(check_on) cout << "   check: this_vet: " << this_vet << endl;            
            // remove this vertex from queue
            vet_queue.erase( vet_queue.begin() );
            // find all edges involving this vertex 
            map<string, vector<string> >::iterator eitr;
            eitr = edge.find(this_vet);
            vector<string> connected_vertex = (*eitr).second;
            if(connected_vertex.size() > 0)
            {
                vector<string>::iterator vitr;
                vector<string>::iterator vitr_end;
                vitr     = connected_vertex.begin();
                vitr_end = connected_vertex.end();
                while(vitr != vitr_end)
                {
                    // update vertex 
                    cluster_vertex[cluster_id].insert(std::pair<string, int>(*vitr, 1));
                    map<string, int>::iterator existence_itr = visited_vet.find(*vitr);
                    if(existence_itr == visited_vet.end())
                    {
                        vet_queue.push_back(*vitr);
                    }
                    // update edge with cor
                    string tmp_edge1 = this_vet + "-" + *vitr; 
                    string tmp_edge2 = *vitr    + "-" + this_vet;                                         
                    citr = cor.find(tmp_edge1);
                    if(citr == cor.end())
                    {
                        citr = cor.find(tmp_edge2);
                    }
                    if(citr != cor.end() )
                    {
                        string         tmp_edge    = (*citr).first; // contig1-contig2
                        vector<string> tmp_vetinfo = split_string(tmp_edge, '-');
                        // update edge
                        (*cluster_edge)[cluster_id].insert(std::pair<string, double>( (*citr).first, (*citr).second) );
                        // update visited vertex 
                        visited_vet.insert(std::pair<string, int>(tmp_vetinfo[0], 1));
                        visited_vet.insert(std::pair<string, int>(tmp_vetinfo[1], 1));                        
                        // remove this <contig1-contig2, cor>
                        cor.erase(citr);
                    }
                    //
                    vitr ++;
                }
            }
        }
    }    
    //
    cout << "       : " << (*cluster_edge).size() << " clusters found. " << endl;
    int N_cluster = 4;    
    if((*cluster_edge).size() != N_cluster)
    {
        string hic_sum_file = out_folder + "/s8_grouping_window_markers_refined_1st_group_hic_stats.txt\0"; 
        ofstream ofp;
        ofp.open(hic_sum_file.c_str(), ios::out);
        if(!ofp.good())
        {
            cout << "   Error: failed in opening file for writing Hi-C stats: " << hic_sum_file << endl;
            return false;
        }
        cout << "   warning: the number of clusters is " << (*cluster_edge).size()
             << ", not the expected value of "           << N_cluster
             << ". you should tune initial Hi-C score to get an expected number of clusters." 
             << endl;
        for(int gi = 1; gi <= N_cluster; gi ++)
        {
            for(int gj = gi; gj <= N_cluster; gj ++)
            {
                ofp << "group\t" << gi << "\t" << gj << "\t0\t" << "avg\t0" << endl; 
            }
        }
        ofp.close();
        return false; 
    }
    //
    return true;
}
bool analyze_hic_matrix(map<string, int>            ctg2_hic_group,
                        map<int, map<string, int> > hic_group_2ctg,
                        map<string, map<int, map<string, unsigned long> > > ctg_allelic_group,
                        map<string, unsigned long>  target_contig_size,
                        map<string, map<unsigned long, NONALLELIE> >* ctg_win_clear_allelic_group,
                        string                      out_prefix,
                        string                      out_folder)
{
    /* Function: for each marker not phased, assign it to 1~4 initial backbone-clusters based on marker type 
          ctg2_hic_group.............: < ctg_id, old_group_id> 
          hic_group_2ctg.............: < clusterid_old, <ctg_id, clusterid_new> >           
          ctg_allelic_group..........: <ctg, <allelic_grp_id, <allelic_ctg, allelic_ctg_size > > >
          out_prefix.................: lable linkage group, chr1, chr2, ....
          win_marker.................: < ctg_id, <win_sta, {win_end, type, win_size, ctg_size, grp_link_val, grp_link_ctg, grp_link_val_sum, grp_link_val_cnt, grp_link_val_size} > >
          hic_matrix.................: < ctg1-ctg2, <sta1, <sta2, {end1,typ1,end2,typ2,hic-cnt} > > >          
          ctg_win_clear_allelic_group: <ctg, <sta, <end, {grp1/grp2..., hic1/hic2/...} > > >
             WHIC
   	       unsigned long     end1; // marker-window end of ctg1
    	       string            typ1; // type of window of ctg1: hap/dip/...
	       unsigned long     end2; // marker-window end of ctg2 
	       string            typ2; // type of window of ctg2: hap/dip/...
	       unsigned long     hcnt; // hic-link-cnt between {ctg1:sta1-end1} and {ctg2:sta2-end2}
	     NODE
               unsigned long     end;  // window end
               string            type; // hap/dip/trip/tetrap/rep
               unsigned long     win_size; // size of this marker
               unsigned long     ctg_size; // size of the contig having this marker
               map<int, double> grp_link_val;     // the best   link between this marker and each hap-extended-group
               map<int, string> grp_link_ctg;     // the best    ctg linking this marker and each hap-extended-group
               map<int, double> grp_link_val_sum; // the total   ctg linking this marker and each hap-extended-group
               map<int, double> grp_link_val_cnt; // the cnt of  ctg linking this marker and each hap-extended-group
               map<int, double> grp_link_val_size;// the size of ctg linking this marker and each hap-extended-group
	     NONALLELIE
               unsigned long end;    // end of window marker
               map<int, double> grp; // <non_allelic_grp_id, Hi-C_score>
    */
    //
    map<string, map<unsigned long, map<unsigned long, WHIC> > >::iterator kitr;
    map<string, map<unsigned long, map<unsigned long, WHIC> > >::iterator kitr_end;
    kitr     = hic_matrix.begin();
    kitr_end = hic_matrix.end();
    while(kitr != kitr_end)
    {
        string this_key = (*kitr).first;
        vector<string> keyinfo = split_string(this_key, '-');
        string ctg1 = keyinfo[0];
        string ctg2 = keyinfo[1];
        // if one ctg is in backbone cluster and the other is not, check if any marker of the other ctg can be assigned.
        if( ctg2_hic_group.find(ctg1) != ctg2_hic_group.end() && ctg2_hic_group.find(ctg2)==ctg2_hic_group.end() )
        {                 
            // find old and new group id of backbone cluster
            // cout << "   check: " << ctg1 << "\t" << ctg2 << ", the former in backbone. " << endl;
            int old_group_id = ctg2_hic_group[ctg1];
            assert(hic_group_2ctg.find(old_group_id) != hic_group_2ctg.end() );
            assert(hic_group_2ctg[old_group_id].find(ctg1) != hic_group_2ctg[old_group_id].end() );
            int new_group_id = hic_group_2ctg[old_group_id][ctg1];
            //
            map<unsigned long, map<unsigned long, WHIC> > tmp1_pos = (*kitr).second;
            //
            map<unsigned long, map<unsigned long, WHIC> >::iterator p1_itr;
            map<unsigned long, map<unsigned long, WHIC> >::iterator p1_itr_end;
            p1_itr     = tmp1_pos.begin();
            p1_itr_end = tmp1_pos.end();
            while(p1_itr != p1_itr_end)
            {
                unsigned long sta1                = (*p1_itr).first;
                map<unsigned long, WHIC> tmp2_pos = (*p1_itr).second;
                //
                map<unsigned long, WHIC>::iterator p2_itr;
                map<unsigned long, WHIC>::iterator p2_itr_end;
                p2_itr     = tmp2_pos.begin();
                p2_itr_end = tmp2_pos.end();
                while(p2_itr != p2_itr_end)
                {
                    unsigned long sta2 = (*p2_itr).first;
                    WHIC tmp_whic      = (*p2_itr).second;
                    unsigned long ctg1_win_size = tmp_whic.end1 - sta1 + 1;
                    unsigned long ctg2_win_size = tmp_whic.end2 - sta2 + 1;
                    double normalized_reads = tmp_whic.hcnt*1.0 / (ctg1_win_size+ctg2_win_size) * normalized_reads_scale;                
                    //
                    assert( win_marker.find(ctg2) != win_marker.end() ); 
                    assert( win_marker[ctg2].find(sta2) != win_marker[ctg2].end() );
                    if( win_marker[ctg2][sta2].grp_link_val.find(new_group_id) == win_marker[ctg2][sta2].grp_link_val.end() )
                    {
                        // initialize best link
                        win_marker[ctg2][sta2].grp_link_val.insert(std::pair<int, double>(new_group_id, normalized_reads));
                        win_marker[ctg2][sta2].grp_link_ctg.insert(std::pair<int, string>(new_group_id, ctg1));
                        // initialize total link 
                        win_marker[ctg2][sta2].grp_link_val_sum.insert(std::pair<int, double>(new_group_id, normalized_reads));
                        win_marker[ctg2][sta2].grp_link_val_cnt.insert(std::pair<int, double>(new_group_id, 1)); 
                        win_marker[ctg2][sta2].grp_link_val_size.insert(std::pair<int, double>(new_group_id, ctg1_win_size+ctg2_win_size));
                    }else
                    {
                        if( win_marker[ctg2][sta2].grp_link_val[new_group_id] < normalized_reads)
                        {
                            // update best link
                            win_marker[ctg2][sta2].grp_link_val[new_group_id] = normalized_reads;
                            win_marker[ctg2][sta2].grp_link_ctg[new_group_id] = ctg1;
                        }                        
                        // update total link 
                        win_marker[ctg2][sta2].grp_link_val_sum[new_group_id]  += normalized_reads;
                        win_marker[ctg2][sta2].grp_link_val_cnt[new_group_id]  += 1;
                        win_marker[ctg2][sta2].grp_link_val_size[new_group_id] += ctg1_win_size+ctg2_win_size;                        
                    }
                    //
                    p2_itr ++;
                }
                //
                p1_itr ++;
            }
        }else
        if( ctg2_hic_group.find(ctg1) == ctg2_hic_group.end() && ctg2_hic_group.find(ctg2)!=ctg2_hic_group.end()  )
        {
            // find old and new group id of backbone cluster
            // cout << "   check: " << ctg1 << "\t" << ctg2 << ", the latter in backbone. " << endl;            
            int old_group_id = ctg2_hic_group[ctg2];
            assert(hic_group_2ctg.find(old_group_id) != hic_group_2ctg.end() );
            assert(hic_group_2ctg[old_group_id].find(ctg2) != hic_group_2ctg[old_group_id].end() );
            int new_group_id = hic_group_2ctg[old_group_id][ctg2];
            //
            map<unsigned long, map<unsigned long, WHIC> > tmp1_pos = (*kitr).second;
            //
            map<unsigned long, map<unsigned long, WHIC> >::iterator p1_itr;
            map<unsigned long, map<unsigned long, WHIC> >::iterator p1_itr_end;
            p1_itr     = tmp1_pos.begin();
            p1_itr_end = tmp1_pos.end();
            while(p1_itr != p1_itr_end)
            {
                unsigned long sta1                = (*p1_itr).first;
                map<unsigned long, WHIC> tmp2_pos = (*p1_itr).second;
                //
                map<unsigned long, WHIC>::iterator p2_itr;
                map<unsigned long, WHIC>::iterator p2_itr_end;
                p2_itr     = tmp2_pos.begin();
                p2_itr_end = tmp2_pos.end();
                while(p2_itr != p2_itr_end)
                {
                    unsigned long sta2 = (*p2_itr).first;
                    WHIC tmp_whic      = (*p2_itr).second;
                    unsigned long ctg1_win_size = tmp_whic.end1 - sta1 + 1;
                    unsigned long ctg2_win_size = tmp_whic.end2 - sta2 + 1;
                    double normalized_reads = tmp_whic.hcnt*1.0 / (ctg1_win_size+ctg2_win_size) * normalized_reads_scale;                
                    //
                    assert( win_marker.find(ctg1) != win_marker.end() ); 
                    assert( win_marker[ctg1].find(sta1) != win_marker[ctg1].end() );
                    if( win_marker[ctg1][sta1].grp_link_val.find(new_group_id) == win_marker[ctg1][sta1].grp_link_val.end() )
                    {
                        // initialize best link
                        win_marker[ctg1][sta1].grp_link_val.insert(std::pair<int, double>(new_group_id, normalized_reads));
                        win_marker[ctg1][sta1].grp_link_ctg.insert(std::pair<int, string>(new_group_id, ctg2));
                        // initialize total link 
                        win_marker[ctg1][sta1].grp_link_val_sum.insert(std::pair<int, double>(new_group_id, normalized_reads));
                        win_marker[ctg1][sta1].grp_link_val_cnt.insert(std::pair<int, double>(new_group_id, 1));
                        win_marker[ctg1][sta1].grp_link_val_size.insert(std::pair<int, double>(new_group_id, ctg1_win_size+ctg2_win_size));                        
                    }else
                    {
                        if( win_marker[ctg1][sta1].grp_link_val[new_group_id] < normalized_reads)
                        {
                            // update best link
                            win_marker[ctg1][sta1].grp_link_val[new_group_id] = normalized_reads;
                            win_marker[ctg1][sta1].grp_link_ctg[new_group_id] = ctg2;                   
                        }
                        // update total link
                        win_marker[ctg1][sta1].grp_link_val_sum[new_group_id]  += normalized_reads;
                        win_marker[ctg1][sta1].grp_link_val_cnt[new_group_id]  += 1;
                        win_marker[ctg1][sta1].grp_link_val_size[new_group_id] += ctg1_win_size+ctg2_win_size;                                 
                    }
                    //
                    p2_itr ++;
                }
                //
                p1_itr ++;
            }            
        }else ;
        //
        kitr ++;
    }
    // output assignment of window markers
    /* Function: output hic contact among window markers */
    string out_file = out_folder + "/s8_grouping_window_markers.txt\0";
    ofstream ofp;
    ofp.open(out_file.c_str(), ios::out);
    if(!ofp.good())
    {
        cout   << "   Error: cannot open file " << out_file << endl;
        return false;
    }    
    //
    map<string, map<unsigned long, NODE> >::iterator mitr;
    map<string, map<unsigned long, NODE> >::iterator mitr_end;
    mitr     = win_marker.begin();
    mitr_end = win_marker.end();
    while(mitr != mitr_end)
    {
        string this_ctg = (*mitr).first;
        cout << "   checky: phasing windows of " << this_ctg << endl;
        /* check allelic contigs in groups .... */
        map<int, map<string, unsigned long> > allelic_grps;
        if(ctg_allelic_group.find(this_ctg) != ctg_allelic_group.end() )
        {
            allelic_grps = ctg_allelic_group[this_ctg];
        }
        // sort the linkage groups according to the number of allelic groups: grps with more ctgs with less priority
        vector<int>    final_lgs;      // <1st_allelic_lg,     2nd_allelic_lg,     ...>
        vector<double> final_ctg_cnt;  // <1st_allelic_lg_hic, 2nd_allelic_lg_hic, ...>
        vector<unsigned long> final_ctg_size; // <1st_allelic_lg_size,2nd_allelic_lg_size, ...>
        vector<unsigned long> final_ctg_size_overlap; // <1st_allelic_lg_overlap_size,2nd_allelic_lg_overlap_size, ...>
        vector<double> final_ctg_size_overlap_ratio; // <1st_allelic_lg_overlap_ratio,2nd_allelic_lg_overlap_rato, ...>        
        map<int, map<string, unsigned long> >::iterator hitr; // <allelic_grp_id, <allelic_ctg, allelic_ctg_size > >
        map<int, map<string, unsigned long> >::iterator hitr_end;
        hitr     = allelic_grps.begin();
        hitr_end = allelic_grps.end();
        while(hitr != hitr_end)
        {
            int linked_lg = (*hitr).first;
            int allelic_ctg_num = (*hitr).second.size();
            unsigned long allelic_ctg_size = 0;
            unsigned long allelic_ctg_size_overlap = 0;            
            double  allelic_ctg_size_overlap_ratio = 0;
            map<string, unsigned long>::iterator sitr;
            map<string, unsigned long>::iterator sitr_end;
            sitr     = (*hitr).second.begin();
            sitr_end = (*hitr).second.end();
            while(sitr != sitr_end)
            {
                allelic_ctg_size += (*sitr).second;
                string allelic_ctg_id = (*sitr).first;
                string allelic_key1 = this_ctg       + "-" + allelic_ctg_id;
                string allelic_key2 = allelic_ctg_id + "-" + this_ctg;
                if( allelic_map_seq_size.find( allelic_key1 ) != allelic_map_seq_size.end() )
                {
                    allelic_ctg_size_overlap += allelic_map_seq_size[ allelic_key1 ];
                }else
                if( allelic_map_seq_size.find( allelic_key2 ) != allelic_map_seq_size.end() )
                {
                    allelic_ctg_size_overlap += allelic_map_seq_size[ allelic_key2 ];
                }else
                {
                    allelic_ctg_size_overlap += 0;
                }
                //
                sitr ++;
            }
            // how much of this_ctg itself is allelic to any given group!
            allelic_ctg_size_overlap_ratio = allelic_ctg_size_overlap * 1.0 / (0 + target_contig_size[this_ctg]);            
            cout << "   | num of allelic ctgs in LG-"<< linked_lg 
                 << ", cnt: "                        << allelic_ctg_num 
                 << ", size: "                       << allelic_ctg_size 
                 << ", overlap size: "               << allelic_ctg_size_overlap 
                 << ", overlap ratio: "              << allelic_ctg_size_overlap_ratio << endl;            
            if(allelic_ctg_size_overlap_ratio < max_allelic_ratio) 
            {
                // if overlap too short, remove allelic relationship
                allelic_ctg_num = 0;
                allelic_ctg_size_overlap = 0;
            }
            //
            if(final_lgs.size() == 0 && allelic_ctg_num > 0)
            {
                final_lgs.push_back(           linked_lg );
                final_ctg_cnt.push_back( allelic_ctg_num );
                final_ctg_size.push_back( allelic_ctg_size );
                final_ctg_size_overlap.push_back( allelic_ctg_size_overlap );
                final_ctg_size_overlap_ratio.push_back(allelic_ctg_size_overlap_ratio);
            }else
            if(allelic_ctg_num > 0)            
            {
                bool tmp_inserted = false;
                vector<int>::iterator                tmp_lg_itr = final_lgs.begin();
                vector<double>::iterator            tmp_hic_itr = final_ctg_cnt.begin();
                vector<unsigned long>::iterator    tmp_size_itr = final_ctg_size.begin();
                vector<unsigned long>::iterator tmp_size_ol_itr = final_ctg_size_overlap.begin();                
                vector<double>::iterator      tmp_size_ol_r_itr = final_ctg_size_overlap_ratio.begin();                                
                for(int ii = 0; ii < final_lgs.size(); ii ++)
                {
                    if( final_ctg_cnt[ii] == allelic_ctg_num )
                    {
                        if(final_ctg_size[ii] <= allelic_ctg_size)
                        {
                            final_lgs.insert(tmp_lg_itr, linked_lg);
                            final_ctg_cnt.insert(tmp_hic_itr, allelic_ctg_num);
                            final_ctg_size.insert(tmp_size_itr, allelic_ctg_size);
                            final_ctg_size_overlap.insert(tmp_size_ol_itr, allelic_ctg_size_overlap);
                            final_ctg_size_overlap_ratio.insert(tmp_size_ol_r_itr, allelic_ctg_size_overlap_ratio);                            
                            tmp_inserted = true;
                            break;
                        }
                    }else                
                    if( final_ctg_cnt[ii] <  allelic_ctg_num )
                    {
                        final_lgs.insert(tmp_lg_itr, linked_lg);
                        final_ctg_cnt.insert(tmp_hic_itr, allelic_ctg_num);
                        final_ctg_size.insert(tmp_size_itr, allelic_ctg_size);
                        final_ctg_size_overlap.insert(tmp_size_ol_itr, allelic_ctg_size_overlap);                        
                        final_ctg_size_overlap_ratio.insert(tmp_size_ol_r_itr, allelic_ctg_size_overlap_ratio);                                                    
                        tmp_inserted = true;
                        break;
                    }
                    tmp_lg_itr ++;
                    tmp_hic_itr ++;
                    tmp_size_itr ++;
                    tmp_size_ol_itr ++;
                    tmp_size_ol_r_itr ++;
                }
                if(!tmp_inserted)
                {
                    final_lgs.insert(tmp_lg_itr, linked_lg);
                    final_ctg_cnt.insert(tmp_hic_itr, allelic_ctg_num);
                    final_ctg_size.insert(tmp_size_itr, allelic_ctg_size);   
                    final_ctg_size_overlap.insert(tmp_size_ol_itr, allelic_ctg_size_overlap);                                     
                    final_ctg_size_overlap_ratio.insert(tmp_size_ol_r_itr, allelic_ctg_size_overlap_ratio);
                }
            }else
            {
                ; // no action if allelic_ctg_num is 0
            }
            //
            hitr ++;
        }
        // print to screen: larger final_ctg_cnt[ii] => lower prob to join this_ctg into the group
        for(int ii = 0; ii < final_lgs.size(); ii ++)
        {
            cout << "   * allelic: LG-"  << final_lgs[ii] 
                 << " with #ctgs = "     << final_ctg_cnt[ii]
                 << " with bp = "        << final_ctg_size[ii]
                 << " with ovelap bp = " << final_ctg_size_overlap[ii]
                 << " with overlap r = " << final_ctg_size_overlap_ratio[ii]
                 << endl;
        } 
        if(final_lgs.size() > 0)
        {       
            cout << endl;
        }else
        {
            cout << "   * no allelic LG found. " << endl;
        }
        //
        /* check allelic contigs in groups done */
        if( ctg2_hic_group.find(this_ctg) == ctg2_hic_group.end() )
        {
            map<unsigned long, NODE> tmp_win = (*mitr).second;
            map<unsigned long, NODE>::iterator witr;
            map<unsigned long, NODE>::iterator witr_end;
            witr     = tmp_win.begin();
            witr_end = tmp_win.end();
            while(witr != witr_end)
            {
                unsigned long win_sta = (*witr).first;
                NODE tmp_node = (*witr).second;
                unsigned long win_end = tmp_node.end;
                string win_type = tmp_node.type;
                cout << "   "        << this_ctg 
                     << "\t"         << win_sta 
                     << "\t"         << win_end 
                     << "\t"         << win_type 
                     << ":"          << endl;                    
                /* sort candidate linkage groups ...  */
                vector<int>    wind_lgs; // <1st_best_lg,     2nd_best_lg,     ...>
                vector<double> wind_hic; // <1st_best_lg_hic, 2nd_best_lg_hic, ...>
                map<int, double>::iterator hitr;
                map<int, double>::iterator hitr_end;
                hitr     = tmp_node.grp_link_val.begin();
                hitr_end = tmp_node.grp_link_val.end();
                while(hitr != hitr_end)
                {
                    int    linked_lg  = (*hitr).first;
                    if(tmp_node.grp_link_val_cnt[linked_lg] == 0) tmp_node.grp_link_val_cnt[linked_lg] = 1;
                    if(tmp_node.grp_link_val_size[linked_lg] == 0) tmp_node.grp_link_val_size[linked_lg] = 1;
                    double linked_val = tmp_node.grp_link_val_sum[linked_lg] / tmp_node.grp_link_val_cnt[linked_lg];                    
                    cout << "      | hic-link score with LG-"<< linked_lg << ": "<< linked_val << endl;
                    //
                    if(wind_lgs.size() == 0)
                    {
                        wind_lgs.push_back( linked_lg  );
                        wind_hic.push_back( linked_val );
                    }else
                    {
                        bool tmp_inserted = false;
                        vector<int>::iterator    tmp_lg_itr  = wind_lgs.begin();
                        vector<double>::iterator tmp_hic_itr = wind_hic.begin();
                        for(int ii = 0; ii < wind_lgs.size(); ii ++)
                        {
                            if( wind_hic[ii] <= linked_val )
                            {
                                wind_lgs.insert(tmp_lg_itr, linked_lg);
                                wind_hic.insert(tmp_hic_itr, linked_val);
                                tmp_inserted = true;
                                break;
                            }
                            tmp_lg_itr ++;
                            tmp_hic_itr ++;                    
                        }
                        if(!tmp_inserted)
                        {
                            wind_lgs.insert(tmp_lg_itr,  linked_lg);
                            wind_hic.insert(tmp_hic_itr, linked_val);                    
                        }
                    }
                    //
                    hitr ++;
                }
                // print info to screen
                for(int ii = 0; ii < wind_lgs.size(); ii ++)
                {
                    cout << "      * candi-grp ordering before cleaning: LG-" << wind_lgs[ii] << " with score " << wind_hic[ii] << endl;
                }        
                if(wind_lgs.size() > 0)
                {
                    cout << "      * ----" << endl;   
                }else
                {
                    cout << "      * no candi-grp ordering before cleaning." << endl;
                    cout << "      * ----" << endl;
                }
                // clean candidate group with allelic groups 
                int num_to_keep = -1;
                if( win_type.compare("hap") == 0 )
                {
                    num_to_keep = 1;            
                }else
                if( win_type.compare("dip") == 0 )
                {
                    num_to_keep = 2;            
                }else
                if( win_type.compare("trip") == 0 )
                {
                    num_to_keep = 3;            
                }else
                if( win_type.compare("tetrap") == 0 )
                {
                    num_to_keep = 4;
                }else num_to_keep = -1;     
                //        
                for(int ii = 0; ii < final_lgs.size(); ii ++)
                {
                    if(wind_lgs.size() > num_to_keep)
                    {
                        for(int ri = 0; ri < wind_lgs.size(); ri ++)
                        {
                            if(wind_lgs[ri] == final_lgs[ii])
                            {
                                // record clear allelic groups
                                if( (*ctg_win_clear_allelic_group).find(this_ctg) == (*ctg_win_clear_allelic_group).end() )
                                {
                                    NONALLELIE tmp_nonallelic;
                                    tmp_nonallelic.end = win_end;
                                    tmp_nonallelic.grp.insert(std::pair<int, double>( wind_lgs[ri], wind_hic[ri] ));
                                    map<unsigned long, NONALLELIE> tmp_ctg_win;
                                    tmp_ctg_win.insert(std::pair<unsigned long, NONALLELIE>(win_sta, tmp_nonallelic));
                                    (*ctg_win_clear_allelic_group).insert(
                                        std::pair<string, map<unsigned long, NONALLELIE> >(this_ctg, tmp_ctg_win));
                                }else
                                {
                                    if( (*ctg_win_clear_allelic_group)[this_ctg].find(win_sta) == 
                                        (*ctg_win_clear_allelic_group)[this_ctg].end() )
                                    {
                                        NONALLELIE tmp_nonallelic;
                                        tmp_nonallelic.end = win_end;
                                        tmp_nonallelic.grp.insert(std::pair<int, double>( wind_lgs[ri], wind_hic[ri] ));
                                        (*ctg_win_clear_allelic_group)[this_ctg].insert(
                                            std::pair<unsigned long, NONALLELIE>(win_sta, tmp_nonallelic));
                                    }else
                                    {
                                        (*ctg_win_clear_allelic_group)[this_ctg][win_sta].grp.insert(
                                            std::pair<int, double>( wind_lgs[ri], wind_hic[ri] ));
                                    }
                                }                 
                                wind_lgs.erase( wind_lgs.begin() + ri );
                                wind_hic.erase( wind_hic.begin() + ri );                                              
                                break;
                            }                            
                        }
                    }
                }      
                // print to screen
                for(int ii = 0; ii < wind_lgs.size(); ii ++)
                {
                    cout << "      * candi-grp ordering after cleaning: LG-" << wind_lgs[ii] 
                         << " with score "                                   << wind_hic[ii] 
                         << endl;
                    //
                    if(ii < num_to_keep)
                    {
                        int    tmp_grp_link_id  = wind_lgs[ii];
                        double tmp_grp_link_val = wind_hic[ii];
                        if(tmp_node.grp_link_val_cnt[tmp_grp_link_id] == 0) 
                            tmp_node.grp_link_val_cnt[tmp_grp_link_id] = 1;
                        if(tmp_node.grp_link_val_size[tmp_grp_link_id] == 0) 
                            tmp_node.grp_link_val_size[tmp_grp_link_id] = 1;
                        /*
                        double tmp_grp_link_val_avg = tmp_node.grp_link_val_sum[tmp_grp_link_id] / 
                                                      tmp_node.grp_link_val_cnt[tmp_grp_link_id] / 
                                                      tmp_node.grp_link_val_size[tmp_grp_link_id]; 
                        */
                        double tmp_grp_link_val_avg = tmp_node.grp_link_val_sum[tmp_grp_link_id] / 
                                                      tmp_node.grp_link_val_cnt[tmp_grp_link_id];                        
                        //
                        double best_link_val = tmp_grp_link_val_avg;
                        int    best_link_grp = tmp_grp_link_id; 
                        string best_link_ctg = tmp_node.grp_link_ctg[ tmp_grp_link_id ];                                        
                        ofp  << this_ctg 
                             << "\t"         << win_sta 
                             << "\t"         << win_end 
                             << "\t"         << win_type                                                  
                             << "\t"         << wind_lgs[ii]  
                             << "\t"         << out_prefix // chr                                                
                             << "\t"         << best_link_ctg
                             << "\t"         << best_link_val
                             << "\tassign"   << ii+1
                             << endl; 
                    }                      
                }
                if(wind_lgs.size() == 0)
                {
                    cout << "      * no cleaning on candi-grp ordering." << endl;
                }
                if(wind_lgs.size() <  num_to_keep)
                {
                    int not_found_cnt = num_to_keep - wind_lgs.size();
                    for(int nfi = 0; nfi < not_found_cnt; nfi ++)
                    {
                        int    best_link_grp = -1;
                        double best_link_val = 0;
                        string best_link_ctg = "NA";
                        ofp  << this_ctg 
                             << "\t"          << win_sta 
                             << "\t"          << win_end 
                             << "\t"          << win_type                                                  
                             << "\t"          << best_link_grp  
                             << "\t"          << out_prefix // Chr                                                 
                             << "\t"          << best_link_ctg
                             << "\t"          << best_link_val
                             << "\tto_refine" << endl;                       
                    }
                }        
                cout << endl;                                                            
                /* sort candidate linkage groups done */                
                // CANNOT LINK MARKERS RELATED TO REPEATS!
                if( win_type.compare("rep") == 0 )
                {
                    int    best_link_grp = -1;
                    double best_link_val = 0;
                    string best_link_ctg = "NA";
                    ofp  << this_ctg 
                         << "\t"         << win_sta 
                         << "\t"         << win_end 
                         << "\t"         << win_type                                                  
                         << "\t"         << best_link_grp  
                         << "\t"         << out_prefix // Chr                                                 
                         << "\t"         << best_link_ctg
                         << "\t"         << best_link_val
                         << "\tassignno" << endl;
                    // remove this best
                    tmp_node.grp_link_val.erase(best_link_grp);                
                    tmp_node.grp_link_ctg.erase(best_link_grp);
                }
                //
                witr ++;
            }
        }else
        {
            // output contigs in backbone clusters 
            map<unsigned long, NODE> tmp_win = (*mitr).second;
            map<unsigned long, NODE>::iterator witr;
            map<unsigned long, NODE>::iterator witr_end;
            witr     = tmp_win.begin();
            witr_end = tmp_win.end();
            while(witr != witr_end)
            {
                unsigned long win_sta = (*witr).first;
                NODE tmp_node = (*witr).second;
                unsigned long win_end = tmp_node.end;
                string win_type = tmp_node.type; 
                //
                int old_group_id = ctg2_hic_group[this_ctg];
                assert(hic_group_2ctg.find(old_group_id) != hic_group_2ctg.end() );
                assert(hic_group_2ctg[old_group_id].find(this_ctg) != hic_group_2ctg[old_group_id].end() );
                int new_group_id = hic_group_2ctg[old_group_id][this_ctg];                
                //
                ofp  << this_ctg 
                     << "\t"         << win_sta 
                     << "\t"         << win_end 
                     << "\t"         << win_type                         
                     << "\t"         << new_group_id 
                     << "\t"         << out_prefix // Chr
                     << "\t"         << this_ctg
                     << "\t"         << 111
                     << "\tassignbb" << endl;  // from backbone cluster
                //
                witr ++;
            }            
        }
        mitr ++;
    }
    //
    ofp.close();
    cout << "       : window markers grouping result has been collected in " << out_file << endl;
    //
    return true;
}
bool define_group_with_allelic(map<string, int>          ctg2_hic_group,
                  map<int, map<string, int> >            hic_group_2ctg,
                  map<string, string>                    allelic_map, 
                  map<string, string>                    allelic_map_seq,
                  map<string, unsigned long>             target_contig_size,
                  map<string, map<int, map<string, unsigned long> > >* ctg_allelic_group)
{
    /* Function: for a ctgx not in any backbone clusters, if it has allelic contigs in some backbone cluster, 
                 this ctgx should not be assigned to those backbone clusters!
                 output:
                 map<string, vector<int> >                allelic_defined_group //<ctgx,<potential groups {1,2,3,4}> >
    */
    map<string, string> allelic_map_merged;
    allelic_map_merged.insert(allelic_map.begin(), allelic_map.end());
    allelic_map_merged.insert(allelic_map_seq.begin(), allelic_map_seq.end());
    cout << "       : "                              << allelic_map_seq.size() 
         << " allelic cases merged with "            << allelic_map.size() 
         << " allelic cases into "                   << allelic_map_merged.size()
         << " allelic cases (common ones removed). " << endl;
    // <this_ctg, <allelic_grp, <allelic_ctg, allelic_ctg_size > > >
    map<string, map<int, map<string, unsigned long> > > ctg_allelic_group_tmp; 
    map<string, string>::iterator citr;
    map<string, string>::iterator citr_end;
    citr     = allelic_map_merged.begin();
    citr_end = allelic_map_merged.end(); 
    while(citr != citr_end)
    {
        string key  = (*citr).first;
        vector<string> keyinfo = split_string(key, '-');
        string ctg1 = keyinfo[0];
        string ctg2 = keyinfo[1];
        //
        if( ctg2_hic_group.find(ctg1)!=ctg2_hic_group.end() && ctg2_hic_group.find(ctg2)==ctg2_hic_group.end() )
        {
            // ctg1     in backbone, but ctg2 not in backbone; ctg2 to be assignd
            int old_group_id = ctg2_hic_group[ctg1];
            int new_group_id = hic_group_2ctg[ old_group_id ][ctg1];            
            if(0)
            cout << "   check: allelic "        << ctg1 
                 << " vs "                      << ctg2
                 << ", where "                  << ctg1 
                 << " in backbone cluster "     << new_group_id
                 << ", so "                     << ctg2
                 << " could not go to cluster " << new_group_id
                 << endl;
            if(ctg_allelic_group_tmp.find(ctg2) == ctg_allelic_group_tmp.end() )
            {
                map<string, unsigned long> tmp_ctg_map;
                tmp_ctg_map.insert(std::pair<string, unsigned long>(ctg1, target_contig_size[ctg1] ));
                map<int, map<string, unsigned long> > tmp_grp_map;
                tmp_grp_map.insert(std::pair<int, map<string, unsigned long> >(new_group_id, tmp_ctg_map));
                ctg_allelic_group_tmp.insert(std::pair<string, map<int, map<string, unsigned long> > >(ctg2, tmp_grp_map));
            }else
            {
                if(ctg_allelic_group_tmp[ctg2].find(new_group_id) == ctg_allelic_group_tmp[ctg2].end())
                {
                    map<string, unsigned long> tmp_ctg_map;
                    tmp_ctg_map.insert(std::pair<string, unsigned long>(ctg1, target_contig_size[ctg1] ));
                    ctg_allelic_group_tmp[ctg2].insert(std::pair<int, map<string, unsigned long> >(new_group_id, tmp_ctg_map));                    
                }else
                {
                    if(ctg_allelic_group_tmp[ctg2][new_group_id].find(ctg1) == ctg_allelic_group_tmp[ctg2][new_group_id].end() )
                    {
                        ctg_allelic_group_tmp[ctg2][new_group_id].insert(std::pair<string, unsigned long>(ctg1, target_contig_size[ctg1] ));
                    }
                }
            }
        }else
        if( ctg2_hic_group.find(ctg1)==ctg2_hic_group.end() && ctg2_hic_group.find(ctg2)!=ctg2_hic_group.end() )
        {
            // ctg1 not in backbone, but ctg2     in backbone; ctg1 to be assigned
            int old_group_id = ctg2_hic_group[ctg2];
            int new_group_id = hic_group_2ctg[ old_group_id ][ctg2];
            if(0)            
            cout << "   check: allelic "        << ctg1 
                 << " vs "                      << ctg2
                 << ", where "                  << ctg2 
                 << " in backbone cluster "     << new_group_id
                 << ", so "                     << ctg1
                 << " could not go to cluster " << new_group_id
                 << endl;  
            if(ctg_allelic_group_tmp.find(ctg1) == ctg_allelic_group_tmp.end() )
            {
                map<string, unsigned long> tmp_ctg_map;
                tmp_ctg_map.insert(std::pair<string, unsigned long>(ctg2, target_contig_size[ctg2] ));
                map<int, map<string, unsigned long> > tmp_grp_map;
                tmp_grp_map.insert(std::pair<int, map<string, unsigned long> >(new_group_id, tmp_ctg_map));
                ctg_allelic_group_tmp.insert(std::pair<string, map<int, map<string, unsigned long> > >(ctg1, tmp_grp_map));
            }else
            {
                if(ctg_allelic_group_tmp[ctg1].find(new_group_id) == ctg_allelic_group_tmp[ctg1].end())
                {
                    map<string, unsigned long> tmp_ctg_map;
                    tmp_ctg_map.insert(std::pair<string, unsigned long>(ctg2, target_contig_size[ctg2] ));
                    ctg_allelic_group_tmp[ctg1].insert(std::pair<int, map<string, unsigned long> >(new_group_id, tmp_ctg_map));                    
                }else
                {
                    if(ctg_allelic_group_tmp[ctg1][new_group_id].find(ctg2) == ctg_allelic_group_tmp[ctg1][new_group_id].end() )
                    {
                        ctg_allelic_group_tmp[ctg1][new_group_id].insert(std::pair<string, unsigned long>(ctg2, target_contig_size[ctg2] ));
                    }
                }
            }                           
        }else ;
        //
        citr ++;
    } 
    //
    bool check_out = true;
    if(check_out)
    {
        map<string, map<int, map<string, unsigned long> > >::iterator pdgitr;
        map<string, map<int, map<string, unsigned long> > >::iterator pdgitr_end;
        pdgitr     = ctg_allelic_group_tmp.begin();
        pdgitr_end = ctg_allelic_group_tmp.end();
        while(pdgitr != pdgitr_end)
        {
            string this_ctg = (*pdgitr).first;
            cout << "       : checkx " << this_ctg << " allelic to backbone cluster(s): " << endl;
            map<int, map<string, unsigned long> > tmp_grp_map = (*pdgitr).second;
            map<int, map<string, unsigned long> >::iterator gitr;
            map<int, map<string, unsigned long> >::iterator gitr_end;
            gitr     = tmp_grp_map.begin();
            gitr_end = tmp_grp_map.end();
            while(gitr != gitr_end)
            {
                int this_allelic_grp = (*gitr).first;
                map<string, unsigned long> tmp_ctg_map = (*gitr).second;
                cout << "           : cluster "  << this_allelic_grp << " with " << tmp_ctg_map.size() << " ctgs, e.g., ";
                map<string, unsigned long>::iterator citr;
                map<string, unsigned long>::iterator citr_end;
                citr     = tmp_ctg_map.begin();
                citr_end = tmp_ctg_map.end();
                int i = 0;
                unsigned long total_allelic_size = 0;
                while(citr != citr_end)
                {
                    if(i<3) // number of example allelic ctgs to print
                    {
                        cout << (*citr).first << ",";
                    }
                    total_allelic_size += (*citr).second;
                    //
                    citr ++;
                    i ++;
                }
                cout << ", total allelic ctg bp = " << total_allelic_size << endl;
                //
                gitr ++;
            }
            //
            pdgitr ++;
        }
    }   
    //
    (*ctg_allelic_group) = ctg_allelic_group_tmp; // <ctg, <allelic_grp, <allelic_ctg, allelic_ctg_size > > >
    //
    return true;
}
//
bool write_hic_matrix(string tmpfolder)
{
    /* Function: output hic contact among window markers */
    string out_file = tmpfolder + "/zadd_hic_matrix.txt\0"; // 
    ofstream ofp;
    ofp.open(out_file.c_str(), ios::out);
    if(!ofp.good())
    {
        cout   << "   Error: cannot open file " << out_file << endl;
        return false;
    }
    //
    map<string, map<unsigned long, map<unsigned long, WHIC> > >::iterator kitr;
    map<string, map<unsigned long, map<unsigned long, WHIC> > >::iterator kitr_end;
    kitr     = hic_matrix.begin();
    kitr_end = hic_matrix.end();
    while(kitr != kitr_end)
    {
        string this_key = (*kitr).first;
        vector<string> keyinfo = split_string(this_key, '-');
        string ctg1 = keyinfo[0];
        string ctg2 = keyinfo[1];
        map<unsigned long, map<unsigned long, WHIC> > tmp1_pos = (*kitr).second;
        //
        map<unsigned long, map<unsigned long, WHIC> >::iterator p1_itr;
        map<unsigned long, map<unsigned long, WHIC> >::iterator p1_itr_end;
        p1_itr     = tmp1_pos.begin();
        p1_itr_end = tmp1_pos.end();
        while(p1_itr != p1_itr_end)
        {
            unsigned long sta1                = (*p1_itr).first;
            map<unsigned long, WHIC> tmp2_pos = (*p1_itr).second;
            //
            map<unsigned long, WHIC>::iterator p2_itr;
            map<unsigned long, WHIC>::iterator p2_itr_end;
            p2_itr     = tmp2_pos.begin();
            p2_itr_end = tmp2_pos.end();
            while(p2_itr != p2_itr_end)
            {
                unsigned long sta2 = (*p2_itr).first;
                WHIC tmp_whic      = (*p2_itr).second;
                unsigned long ctg1_win_size = tmp_whic.end1 - sta1 + 1;
                unsigned long ctg2_win_size = tmp_whic.end2 - sta2 + 1;
                double normalized_reads = tmp_whic.hcnt*1.0 / (ctg1_win_size+ctg2_win_size) * normalized_reads_scale;                
                ofp << ctg1             << "\t" 
                    << sta1             << "\t" 
                    << tmp_whic.end1    << "\t" 
                    << tmp_whic.typ1    << "\t"
                    << ctg2             << "\t"
                    << sta2             << "\t" 
                    << tmp_whic.end2    << "\t"
                    << tmp_whic.typ2    << "\t"
                    << normalized_reads << endl;                  
                //
                p2_itr ++;
            }
            //
            p1_itr ++;
        }
        //
        kitr ++;
    }
    cout << "       : hic matrix wrote out to file " << out_file << endl;
    //
    ofp.close();
    return true;
}
//
bool update_hic_matrix(unsigned long read1_sta, unsigned long read2_sta, string ctg1, string ctg2)
{
    /* Function: update the hic links between ctg1 window and ctg2 window
           read1_sta........... where read 1 is aligned along ctg1
           read2_sta........... where read 2 is aligned along ctg2
           map<string, map<unsigned long, NODE> > win_marker; < ctg, <win_sta, {win_end, type, win_size, ctg_size, grp_link_val, grp_link_ctg} > >
           map<string, map<unsigned long, map<unsigned long, WHIC> > > hic_matrix; < ctg1-ctg2, <sta1, <sta2, {end1,typ1,end2,typ2,hic-cnt} > > >
             WHIC
   	       unsigned long     end1; // marker-window end of ctg1
    	       string            typ1; // type of window of ctg1: hap/dip/...
	       unsigned long     end2; // marker-window end of ctg2 
	       string            typ2; // type of window of ctg2: hap/dip/...
	       unsigned long     hcnt; // hic-link-cnt between {ctg1:sta1-end1} and {ctg2:sta2-end2}	       
	     NODE
               unsigned long     end; // window end
               string            type; // hap/dip/trip/tetrap/rep
               unsigned long     win_size; // size of this marker
               unsigned long     ctg_size; // size of the contig having this marker
    */
    if( win_marker.find(ctg1) == win_marker.end() )
    {
        cout << "   Error: " << ctg1 << " not found in window marker list. " << endl;
        return false;
    }
    if( win_marker.find(ctg2) == win_marker.end() )
    {
        cout << "   Error: " << ctg2 << " not found in window marker list. " << endl;
        return false;
    }
    //
    unsigned long win1_sta = 0;
    unsigned long win1_end = 0;
    string        win1_typ = "";    
    unsigned long win2_sta = 0;
    unsigned long win2_end = 0;
    string        win2_typ = "";
    unsigned long hcnt     = 0;
    // locate marker windows
    map<unsigned long, NODE>::iterator witr;
    map<unsigned long, NODE>::iterator witr_end;    
    map<unsigned long, NODE> ctg1_win = win_marker[ctg1];    
    witr     = ctg1_win.begin();
    witr_end = ctg1_win.end();
    while(witr != witr_end)
    {
        win1_sta = (*witr).first;
        win1_end = (*witr).second.end;
        if(win1_sta<=read1_sta && read1_sta<=win1_end)
        {
            win1_typ = (*witr).second.type;
            break;
        }        
        //
        witr ++;
    }
    if(win1_typ.size()==0)
    {
        cout << "   Error: no window marker found harbouring position " << read1_sta << " along " << ctg1 << endl;
        return false;
    }
    map<unsigned long, NODE> ctg2_win = win_marker[ctg2];
    witr     = ctg2_win.begin();
    witr_end = ctg2_win.end();
    while(witr != witr_end)
    {
        win2_sta = (*witr).first;
        win2_end = (*witr).second.end;
        if(win2_sta<=read2_sta && read2_sta<=win2_end)
        {
            win2_typ = (*witr).second.type;
            break;
        }
        //
        witr ++;
    }
    if(win2_typ.size()==0)
    {
        cout << "   Error: no window marker found harbouring position " << read2_sta << " along " << ctg2 << endl;
        return false;
    } 
    if(0)
    cout << "   checkkk: " 
         << read1_sta << "\t"
         << read2_sta << "\t"
         << ctg1      << "\t" 
         << win1_sta  << "\t" 
         << win1_end  << "\t" 
         << win1_typ  << "\t"
         << ctg2      << "\t"
         << win2_sta  << "\t" 
         << win2_end  << "\t"
         << win2_typ  << endl;
    // map<ctg1-ctg2, map<sta1, map<sta2, WHIC> > >
    string key1 = ctg1 + "-" + ctg2;
    string key2 = ctg2 + "-" + ctg1;
    if( hic_matrix.find(key1)==hic_matrix.end() && hic_matrix.find(key2)==hic_matrix.end() )
    {
        WHIC tmp_whic;
        tmp_whic.end1 = win1_end;
        tmp_whic.typ1 = win1_typ;
        tmp_whic.end2 = win2_end;
        tmp_whic.typ2 = win2_typ;
        tmp_whic.hcnt = 1;
        map<unsigned long, WHIC> tmp2_pos;
        tmp2_pos.insert(std::pair<unsigned long, WHIC>(win2_sta, tmp_whic) );
        map<unsigned long, map<unsigned long, WHIC> > tmp1_pos;
        tmp1_pos.insert(std::pair<unsigned long, map<unsigned long, WHIC> >(win1_sta, tmp2_pos));
        hic_matrix.insert(std::pair<string, map<unsigned long, map<unsigned long, WHIC> > >(key1, tmp1_pos));
    }else
    if( hic_matrix.find(key1)!=hic_matrix.end() )
    {
        if( hic_matrix[key1].find(win1_sta) == hic_matrix[key1].end() )
        {
            WHIC tmp_whic;
            tmp_whic.end1 = win1_end;
            tmp_whic.typ1 = win1_typ;
            tmp_whic.end2 = win2_end;
            tmp_whic.typ2 = win2_typ;
            tmp_whic.hcnt = 1;  
            map<unsigned long, WHIC> tmp2_pos;
            tmp2_pos.insert(std::pair<unsigned long, WHIC>(win2_sta, tmp_whic) );
            hic_matrix[key1].insert(std::pair<unsigned long, map<unsigned long, WHIC> >(win1_sta, tmp2_pos));
        }else
        {
            if( hic_matrix[key1][win1_sta].find(win2_sta) == hic_matrix[key1][win1_sta].end() )
            {
                WHIC tmp_whic;
                tmp_whic.end1 = win1_end;
                tmp_whic.typ1 = win1_typ;
                tmp_whic.end2 = win2_end;
                tmp_whic.typ2 = win2_typ;
                tmp_whic.hcnt = 1;                  
                hic_matrix[key1][win1_sta].insert(std::pair<unsigned long, WHIC>(win2_sta, tmp_whic) );                
            }else
            {
                hic_matrix[key1][win1_sta][win2_sta].hcnt += 1;
            }
        }
    }else
    if( hic_matrix.find(key2)!=hic_matrix.end() )
    {
        /* swap ctg2-ctg1 as ctg1-ctg2 */
        string tmp_ctg = ctg1;
        ctg2           = ctg1;
        ctg1           = tmp_ctg;
        unsigned long tmp_sta = win1_sta;
        win1_sta              = win2_sta;
        win2_sta              = tmp_sta;
        //
        if( hic_matrix[key2].find(win1_sta) == hic_matrix[key2].end() )
        {
            WHIC tmp_whic;
            tmp_whic.end1 = win2_end; // swapped
            tmp_whic.typ1 = win2_typ; // swapped
            tmp_whic.end2 = win1_end; // swapped
            tmp_whic.typ2 = win1_typ; // swapped
            tmp_whic.hcnt = 1;  
            map<unsigned long, WHIC> tmp2_pos;
            tmp2_pos.insert(std::pair<unsigned long, WHIC>(win2_sta, tmp_whic) );
            hic_matrix[key2].insert(std::pair<unsigned long, map<unsigned long, WHIC> >(win1_sta, tmp2_pos));
        }else
        {
            if( hic_matrix[key2][win1_sta].find(win2_sta) == hic_matrix[key2][win1_sta].end() )
            {
                WHIC tmp_whic;
                tmp_whic.end1 = win2_end; // swapped
                tmp_whic.typ1 = win2_typ; // swapped
                tmp_whic.end2 = win1_end; // swapped 
                tmp_whic.typ2 = win1_typ; // swapped
                tmp_whic.hcnt = 1;                  
                hic_matrix[key2][win1_sta].insert(std::pair<unsigned long, WHIC>(win2_sta, tmp_whic) );                
            }else
            {
                hic_matrix[key2][win1_sta][win2_sta].hcnt += 1;
            }
        }
    } else;
    //
    return true;
}
//
void update_links(string this_a_ctg, 
                  string this_b_ctg,
                  string contact_type,                                 // LL, LR, RL, RR
                  unsigned long* n_pair_inter_long_hap,
                  unsigned long* n_pair_other
                  )
{
    /*
	map<string, CONTACT> ctg_cross_link_long_hap2;           // global: < ctg1-ctg2, {LL, LR, RL, RR} >
	map<string, CONTACT> ctg_cross_link_other2;              // global: < ctg1-ctg2, {LL, LR, RL, RR} >      
	map<string, map<unsigned long, NODE> > hap_win_marker;   // global
	map<string, map<unsigned long, NODE> > other_win_marker; // global
    */
    string key1 = this_a_ctg + "-" + this_b_ctg;
    string key2 = this_b_ctg + "-" + this_a_ctg;
    //
    if( hap_win_marker.find(this_a_ctg) != hap_win_marker.end() &&
        hap_win_marker.find(this_b_ctg) != hap_win_marker.end()         
      )
    {
        // read pair mapped to different contigs and both in the long haplotig list 
        (*n_pair_inter_long_hap) ++;
        if((ctg_cross_link_long_hap2).find(key1) == (ctg_cross_link_long_hap2).end() &&
           (ctg_cross_link_long_hap2).find(key2) == (ctg_cross_link_long_hap2).end() )
        {
            CONTACT tmp_contact;
            tmp_contact.LL = 0;
            tmp_contact.LR = 0;
            tmp_contact.RL = 0;
            tmp_contact.RR = 0;
            if(contact_type.compare("LL"))
            {
                tmp_contact.LL = 1;
            }else
            if(contact_type.compare("LR"))
            {
                tmp_contact.LR = 1;
            }else
            if(contact_type.compare("RL"))
            {
                tmp_contact.RL = 1;
            }else
            if(contact_type.compare("RR"))
            {
                tmp_contact.RR = 1;
            }else ;                                                
            (ctg_cross_link_long_hap2).insert(std::pair<string, CONTACT>(key1, tmp_contact));
        }else
        if((ctg_cross_link_long_hap2).find(key1) != (ctg_cross_link_long_hap2).end() &&
           (ctg_cross_link_long_hap2).find(key2) == (ctg_cross_link_long_hap2).end() )        
        {
            if(contact_type.compare("LL"))
            {
                (ctg_cross_link_long_hap2)[key1].LL += 1;
            }else
            if(contact_type.compare("LR"))
            {
                (ctg_cross_link_long_hap2)[key1].LR += 1;
            }else
            if(contact_type.compare("RL"))
            {
                (ctg_cross_link_long_hap2)[key1].RL += 1;
            }else
            if(contact_type.compare("RR"))
            {
                (ctg_cross_link_long_hap2)[key1].RR += 1;
            }else ;
        }else
        if((ctg_cross_link_long_hap2).find(key1) == (ctg_cross_link_long_hap2).end() &&
           (ctg_cross_link_long_hap2).find(key2) != (ctg_cross_link_long_hap2).end() )        
        {
            if(contact_type.compare("LL"))
            {
                (ctg_cross_link_long_hap2)[key2].LL += 1;
            }else
            if(contact_type.compare("LR"))
            {
                (ctg_cross_link_long_hap2)[key2].LR += 1;
            }else
            if(contact_type.compare("RL"))
            {
                (ctg_cross_link_long_hap2)[key2].RL += 1;
            }else
            if(contact_type.compare("RR"))
            {
                (ctg_cross_link_long_hap2)[key2].RR += 1;
            }else ;
        }else ;
    }else
    { 
        // read pairs mapped to other contigs
        (*n_pair_other) ++;
        if((ctg_cross_link_other2).find(key1) == (ctg_cross_link_other2).end() &&
           (ctg_cross_link_other2).find(key2) == (ctg_cross_link_other2).end() )
        {
            CONTACT tmp_contact;
            tmp_contact.LL = 0;
            tmp_contact.LR = 0;
            tmp_contact.RL = 0;
            tmp_contact.RR = 0;
            if(contact_type.compare("LL"))
            {
                tmp_contact.LL = 1;
            }else
            if(contact_type.compare("LR"))
            {
                tmp_contact.LR = 1;
            }else
            if(contact_type.compare("RL"))
            {
                tmp_contact.RL = 1;
            }else
            if(contact_type.compare("RR"))
            {
                tmp_contact.RR = 1;
            }else ;                                                
            (ctg_cross_link_other2).insert(std::pair<string, CONTACT>(key1, tmp_contact));
        }else
        if((ctg_cross_link_other2).find(key1) != (ctg_cross_link_other2).end() &&
           (ctg_cross_link_other2).find(key2) == (ctg_cross_link_other2).end() )        
        {
            if(contact_type.compare("LL"))
            {
                (ctg_cross_link_other2)[key1].LL += 1;
            }else
            if(contact_type.compare("LR"))
            {
                (ctg_cross_link_other2)[key1].LR += 1;
            }else
            if(contact_type.compare("RL"))
            {
                (ctg_cross_link_other2)[key1].RL += 1;
            }else
            if(contact_type.compare("RR"))
            {
                (ctg_cross_link_other2)[key1].RR += 1;
            }else ;        
        }else
        if((ctg_cross_link_other2).find(key1) == (ctg_cross_link_other2).end() &&
           (ctg_cross_link_other2).find(key2) != (ctg_cross_link_other2).end() )        
        {
            if(contact_type.compare("LL"))
            {
                (ctg_cross_link_other2)[key2].LL += 1;
            }else
            if(contact_type.compare("LR"))
            {
                (ctg_cross_link_other2)[key2].LR += 1;
            }else
            if(contact_type.compare("RL"))
            {
                (ctg_cross_link_other2)[key2].RL += 1;
            }else
            if(contact_type.compare("RR"))
            {
                (ctg_cross_link_other2)[key2].RR += 1;
            }else ;
        }else ;    
    }
}
//
bool get_allelic_ctg_align_based(string      align_file, 
                  map<string, unsigned long> target_contig_size)
{
    ifstream ifp;
    ifp.open(align_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << align_file << endl;
        return false;
    }
    // s1. read contig aligned to coordinates of ref
    map<string, map<unsigned long, unsigned long> > align_to_ref; 
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        /*  0.ctg_id          1.ctg_size  2.ctg_sta  3.ctg_end  4.strand  5.ref_seq  6.ref_size  7.ref_sta  8.ref_end...
            utg000353l_pilon  155552      3754       20136      +         chr01      88591686    189268     205431   ...
        */
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 9)
        {
            cout << "   Warning: skip line with insufficient info: " << line << endl;
            continue;
        }
        string ctg_id = lineinfo[0];
        if(target_contig_size.find(ctg_id) == target_contig_size.end())
        {
            // contig was not assigned to the current linkage group, skip it!
            continue;
        }
        unsigned long ctg_sta = strtoul(lineinfo[2].c_str(), NULL, 0);
        unsigned long ctg_end = strtoul(lineinfo[3].c_str(), NULL, 0);
        unsigned long ref_sta = strtoul(lineinfo[7].c_str(), NULL, 0);
        unsigned long ref_end = strtoul(lineinfo[8].c_str(), NULL, 0);
        if(ref_end - ref_sta > 3000)
        {
            if( align_to_ref.find(ctg_id) == align_to_ref.end() )
            {
                map<unsigned long, unsigned long> tmp_align;
                tmp_align.insert(std::pair<unsigned long, unsigned long>(ref_sta, ref_end));
                align_to_ref.insert(std::pair<string, map<unsigned long, unsigned long> >(ctg_id, tmp_align));
            }else
            {
                align_to_ref[ctg_id].insert( std::pair<unsigned long, unsigned long>(ref_sta, ref_end) );
            }
        }
    }
    ifp.close();
    // s2. clean redundant intervals
    cout << "       : cleaning alignments... " << endl;    
    map<string, map<unsigned long, unsigned long> >::iterator citr;
    map<string, map<unsigned long, unsigned long> >::iterator citr_end;    
    citr     = align_to_ref.begin();
    citr_end = align_to_ref.end();
    map<string, map<unsigned long, unsigned long> > align_to_ref_clean;
    while(citr != citr_end)
    {
        string ctg = (*citr).first;
        map<unsigned long, unsigned long> tmp_align = (*citr).second;        
        map<unsigned long, unsigned long> tmp_align_cleaned;
        assert( target_contig_size.find(ctg) != target_contig_size.end() );
        //if(!fuzzy_clean_alignment(tmp_align, target_contig_size[ctg], &tmp_align_cleaned))
        if(!clean_alignment(tmp_align, ctg, &tmp_align_cleaned))        
        {
            cout << "   Error: failed in cleaning alignments. " << endl;
            return false;
        }
        align_to_ref_clean.insert(std::pair<string, map<unsigned long, unsigned long> >(ctg, tmp_align_cleaned));
        // get chr-region where a ctg is aligned to
        if(ctg_at_chr_region.find(ctg) == ctg_at_chr_region.end() )
        {
            locate_ctg_at_chr(tmp_align_cleaned, chr_region, ctg);            
        }       
        //
        citr ++;
    }
    cout << "       : cleaning alignments done. " << endl;            
    // s3. analyze allelic relationship among contigs
    map<string, int> ctg_tmp;
    map<string, map<unsigned long, unsigned long> >::iterator citr_1st;
    citr_1st = align_to_ref_clean.begin();
    citr_end = align_to_ref_clean.end();
    while(citr_1st != citr_end)
    {
        string ctg_1 = (*citr_1st).first;
        map<unsigned long, unsigned long> tmp_align_1 = (*citr_1st).second;
        assert( target_contig_size.find(ctg_1) != target_contig_size.end() );        
        //
        map<string, map<unsigned long, unsigned long> >::iterator citr_2nd = citr_1st;
        citr_2nd ++;
        while(citr_2nd != citr_end)
        {
            string ctg_2 = (*citr_2nd).first;
            map<unsigned long, unsigned long> tmp_align_2 = (*citr_2nd).second;            
            assert( target_contig_size.find(ctg_2) != target_contig_size.end() );
            // calculate the degree of overlap between ctg_1 and ctg_2
            unsigned long overlap_size = calculate_overlap(tmp_align_1, tmp_align_2);
            if(overlap_size > 10000) // at least 10 kb overlap -  // TODO: set up option on threshold
            if(overlap_size*1.0 > target_contig_size[ctg_1]*0.1 || 
               overlap_size*1.0 > target_contig_size[ctg_2]*0.1 ) // TODO: set up option on threshold
            {
                if(verbose)
                cout << "       : "       << ctg_1 
                     << " and "           << ctg_2 
                     << " overlap by "    << overlap_size 
                     << " bp, with size " << target_contig_size[ctg_1]
                     << " vs "            << target_contig_size[ctg_2]
                     << endl;
                string key1  = ctg_1 + "-" + ctg_2;
                string key2  = ctg_2 + "-" + ctg_1;
                if( allelic_map_seq.find(key1) == allelic_map_seq.end() )
                {
                    allelic_map_seq.insert(std::pair<string, string>(key1, "0"));
                    allelic_map_seq_size.insert(std::pair<string, unsigned long>(key1, overlap_size));
                }
                if( allelic_map_seq.find(key2) == allelic_map_seq.end() )
                {
                    allelic_map_seq.insert(std::pair<string, string>(key2, "0"));
                    allelic_map_seq_size.insert(std::pair<string, unsigned long>(key1, overlap_size));                    
                }
                //
                if(ctg_tmp.find(ctg_1) == ctg_tmp.end())
                {
                    ctg_tmp.insert(std::pair<string, int>(ctg_1, 0));
                }
                if(ctg_tmp.find(ctg_2) == ctg_tmp.end())
                {
                    ctg_tmp.insert(std::pair<string, int>(ctg_2, 0));
                }
            }         
            //
            citr_2nd ++;
        }
        //
        citr_1st ++;
    }
    //
    cout << "   Info: "                            << allelic_map_seq.size() 
         << " allelic relationship collected for " << ctg_tmp.size()
         << " contigs - alignment based. "         << endl;    
    //
    return true;
}
//
int locate_ctg_at_chr(map<unsigned long, unsigned long> tmp_align, 
                      map<unsigned long, unsigned long> chr_region,
                      string ctg)
{
    /* function: given alignments of a contig to chr, roughly locate which region (int:0-9) of the chr it belongs to 
                    tmp_align..................<ctg_align_sta, ctg_align_end>
                    chr_region.................<reg_sta, reg_end>
                 => ctg_at_chr_region..........global <ctg, region0-9>                 
    */
    int res = 0;
    bool this_verbose = false;    
    // initialize overlapping size at each region as 0
    map<int, unsigned long> reg_overlap_size; // <region0-9, bp_overlapping_ctg>
    for(int ri = 0; ri < region_cnt; ri ++)
    {
        reg_overlap_size.insert(std::pair<int, unsigned long>(ri, 0) );
    }
    //
    map<unsigned long, unsigned long>::iterator pitr_ctg;
    map<unsigned long, unsigned long>::iterator pitr_ctg_end;
    pitr_ctg     = tmp_align.begin();
    pitr_ctg_end = tmp_align.end();
    while(pitr_ctg != pitr_ctg_end)
    {
        unsigned long sta_ctg = (*pitr_ctg).first;
        unsigned long end_ctg = (*pitr_ctg).second;
        //
        map<unsigned long, unsigned long>::iterator pitr_reg;
        map<unsigned long, unsigned long>::iterator pitr_reg_end;
        pitr_reg      = chr_region.begin();
        pitr_reg_end  = chr_region.end();
        int region_i  = -1;
        while(pitr_reg != pitr_reg_end)
        {
            region_i  ++;
            unsigned long sta_reg = (*pitr_reg).first;
            unsigned long end_reg = (*pitr_reg).second;            
            //
            if(end_reg < sta_ctg)
            {
                pitr_reg ++;
                continue;
            }
            if(sta_reg > end_ctg)
            {
                break;
            }
            //
            if(sta_reg<=sta_ctg && 
               sta_ctg<=end_reg && end_reg<=end_ctg)
            {
                /*
                             s_ctg---------------e_ctg
                     s_reg------------e_reg
                */
                reg_overlap_size[ region_i ] += (end_reg - sta_ctg + 1);                 
                if(this_verbose)    
                cout << "   locate:  " << sta_ctg << "~" << end_ctg << " ovl " 
                                       << sta_reg << "~" << end_reg << " by " 
                     << end_reg - sta_ctg + 1     << endl;               
            }else
            if(sta_reg<=sta_ctg && end_reg>end_ctg)
            {
                /*
                             s_ctg---------------e_ctg
                     s_reg----------------------------e_reg
                */
                reg_overlap_size[ region_i ] += (end_ctg - sta_ctg + 1);
                if(this_verbose)    
                cout << "   locate:  " << sta_ctg << "~" << end_ctg << " ovl " 
                                       << sta_reg << "~" << end_reg << " by " 
                     << end_ctg - sta_ctg + 1     << endl;              
            }else      
            if(sta_reg>sta_ctg && end_reg<=end_ctg)
            {
                /*
                             s_ctg---------------e_ctg
                                     s_reg----e_reg
                */
                reg_overlap_size[ region_i ] += (end_reg - sta_reg + 1);      
                if(this_verbose)    
                cout << "   locate:  " << sta_ctg << "~" << end_ctg << " ovl " 
                                       << sta_reg << "~" << end_reg << " by " 
                     << end_reg - sta_reg + 1     << endl;             
            }else     
            if(sta_ctg<sta_reg && sta_reg<end_ctg && 
               end_reg>end_ctg)
            {
                /*
                             s_ctg---------------e_ctg
                                     s_reg------------e_reg
                */
                reg_overlap_size[ region_i ] += (end_ctg - sta_reg + 1);      
                if(this_verbose)
                cout << "   locate:  " << sta_ctg << "~" << end_ctg << " ovl " 
                                       << sta_reg << "~" << end_reg << " by " 
                     << end_ctg - sta_reg + 1     << endl;
            }else ;                     
            //
            pitr_reg ++;
        }
        //
        pitr_ctg ++;
    }
    // find a region with max overlap with current contig 
    unsigned long max_ovl = 0;
    int max_region_i = 0;
    for(int ri = 0; ri < region_cnt; ri ++)
    {
        if(reg_overlap_size[ri] > max_ovl)
        {
            max_ovl      = reg_overlap_size[ri];
            max_region_i = ri;
        }
    }
    ctg_at_chr_region.insert(std::pair<string, int>(ctg, max_region_i));
    if(this_verbose)
    cout << "   check: "          << ctg 
         << " located at region " << max_region_i 
         << " = "                 << (max_region_i+0)*resolution_region+1
         << "~"                   << (max_region_i+1)*resolution_region
         << endl;
    //    
    return res;
}
//
unsigned long calculate_overlap(map<unsigned long, unsigned long> tmp_align_1, 
                                map<unsigned long, unsigned long> tmp_align_2)
{
    /* function: calculate overlapping size of two sets of alignments of ctg 1 and ctg 2 to the same ref */
    unsigned long overlap_size = 0;
    bool this_verbose = false;
    //
    map<unsigned long, unsigned long>::iterator pitr_1;
    map<unsigned long, unsigned long>::iterator pitr_1_end;
    pitr_1     = tmp_align_1.begin();
    pitr_1_end = tmp_align_1.end();
    while(pitr_1 != pitr_1_end)
    {
        unsigned long sta_1 = (*pitr_1).first;
        unsigned long end_1 = (*pitr_1).second;
        //
        map<unsigned long, unsigned long>::iterator pitr_2;
        map<unsigned long, unsigned long>::iterator pitr_2_end;
        pitr_2     = tmp_align_2.begin();
        pitr_2_end = tmp_align_2.end();
        while(pitr_2 != pitr_2_end)
        {
            unsigned long sta_2 = (*pitr_2).first;
            unsigned long end_2 = (*pitr_2).second;            
            //
            if(end_2 < sta_1)
            {
                pitr_2 ++;
                continue;
            }
            if(sta_2 > end_1)
            {
                break;
            }
            //
            if(sta_2<=sta_1 && 
               sta_1<=end_2 && end_2<=end_1)
            {
                /*
                             s1---------------e1
                     s2------------e2
                */
                overlap_size += (end_2 - sta_1 + 1);                 
                if(this_verbose)    
                cout << "   check:  " << sta_1 << "~" << end_1 << " ovl " 
                                      << sta_2 << "~" << end_2 << " by " 
                     << end_2 - sta_1 + 1      << endl;               
            }else
            if(sta_2<=sta_1 && end_2>end_1)
            {
                /*
                             s1---------------e1
                     s2----------------------------e2
                */
                overlap_size += (end_1 - sta_1 + 1);
                if(this_verbose)    
                cout << "   check:  " << sta_1 << "~" << end_1 << " ovl " 
                                      << sta_2 << "~" << end_2 << " by " 
                     << end_1 - sta_1 + 1      << endl;              
            }else      
            if(sta_2>sta_1 && end_2<=end_1)
            {
                /*
                             s1---------------e1
                                     s2----e2
                */
                overlap_size += (end_2 - sta_2 + 1);      
                if(this_verbose)    
                cout << "   check:  " << sta_1 << "~" << end_1 << " ovl " 
                                      << sta_2 << "~" << end_2 << " by " 
                     << end_2 - sta_2 + 1      << endl;             
            }else     
            if(sta_1<sta_2 && sta_2<end_1 && 
               end_2>end_1)
            {
                /*
                             s1---------------e1
                                     s2------------e2
                */
                overlap_size += (end_1 - sta_2 + 1);      
                if(this_verbose)
                cout << "   check:  " << sta_1 << "~" << end_1 << " ovl " 
                                      << sta_2 << "~" << end_2 << " by " 
                     << end_1 - sta_2 + 1      << endl;
            }else ;                     
            //
            pitr_2 ++;
        }
        //
        pitr_1 ++;
    }
    //
    return overlap_size;
}
//
bool fuzzy_clean_alignment(map<unsigned long, unsigned long>  tmp_align,
                           unsigned long this_contig_size,
                           map<unsigned long, unsigned long>* cleaned_align)
{
    /* function: remove highly overlapping alignments
                 this_contig_size*0.1 gap size would be allowed to merge "overlapping" intervals 
    */    
    if(tmp_align.size() == 0)
    {
        cout << "   Error: empty alignment " << endl;
        return false;
    }
    map<unsigned long, unsigned long>::iterator pitr_last;    
    map<unsigned long, unsigned long>::iterator pitr_1;
    map<unsigned long, unsigned long>::iterator pitr_1_end;
    pitr_1     = tmp_align.begin();
    pitr_1_end = tmp_align.end();
    (*cleaned_align).insert( std::pair<unsigned long, unsigned long>( (*pitr_1).first, (*pitr_1).second ) ); // 1st 
    pitr_1 ++; // 2nd and afterwards
    while(pitr_1 != pitr_1_end)
    {        
        pitr_last = (*cleaned_align).end();
        pitr_last --;
        unsigned long sta_last = (*pitr_last).first;
        unsigned long end_last = (*pitr_last).second;
        //
        unsigned long sta_next = (*pitr_1).first;
        unsigned long end_next = (*pitr_1).second;
        //
        if( end_last + this_contig_size*0.1 < sta_next)
        {
            // non-overlapping interval
            (*cleaned_align).insert(std::pair<unsigned long, unsigned long>(sta_next, end_next));
        }else
        if( end_last < end_next)
        {
            // overlapping interval 
            (*cleaned_align)[ sta_last ] = end_next; 
        }else ;
        //
        pitr_1 ++;
    }
    //
    if(0)
    {
        cout << "   check: " << tmp_align.size() << " alignments before cleaning: "   << endl;
        // print_alignments(tmp_align);
        cout << "   check: fuzzy level " << this_contig_size*0.1 << " bp  = 0.1*ctg size " << endl;        
        cout << "   check: " << (*cleaned_align).size() << " alignments after fuzzy cleaning: " << endl;
        print_alignments(*cleaned_align);
    }
    //
    return true;
}
//
bool clean_alignment(map<unsigned long, unsigned long>  tmp_align,
                     string ctg,
                     map<unsigned long, unsigned long>* cleaned_align)
{
    /* function: remove highly overlapping alignments */
    if(tmp_align.size() == 0)
    {
        cout << "   Error: empty alignment " << endl;
        return false;
    }
    map<unsigned long, unsigned long>::iterator pitr_last;    
    map<unsigned long, unsigned long>::iterator pitr_1;
    map<unsigned long, unsigned long>::iterator pitr_1_end;
    pitr_1     = tmp_align.begin();
    pitr_1_end = tmp_align.end();
    (*cleaned_align).insert( std::pair<unsigned long, unsigned long>( (*pitr_1).first, (*pitr_1).second ) ); // 1st 
    pitr_1 ++; // 2nd and afterwards
    while(pitr_1 != pitr_1_end)
    {        
        pitr_last = (*cleaned_align).end();
        pitr_last --;
        unsigned long sta_last = (*pitr_last).first;
        unsigned long end_last = (*pitr_last).second;
        //
        unsigned long sta_next = (*pitr_1).first;
        unsigned long end_next = (*pitr_1).second;
        //
        if( end_last < sta_next)
        {
            // non-overlapping interval
            (*cleaned_align).insert(std::pair<unsigned long, unsigned long>(sta_next, end_next));
        }else
        if( end_last < end_next)
        {
            // overlapping interval 
            (*cleaned_align)[ sta_last ] = end_next; 
        }else ;
        //
        pitr_1 ++;
    }
    //
    if(0)
    {
        cout << "   check: " << ctg              << endl; 
        cout << "        : " << tmp_align.size() << " alignments before cleaning: " << endl;
        print_alignments(tmp_align);
        cout << "        : " << (*cleaned_align).size() << " alignments after  cleaning: " << endl;
        print_alignments(*cleaned_align);
    }
    //
    return true;
}
void print_alignments(map<unsigned long, unsigned long> tmp_align)
{
    map<unsigned long, unsigned long>::iterator pitr;
    map<unsigned long, unsigned long>::iterator pitr_end;
    pitr     = tmp_align.begin();
    pitr_end = tmp_align.end();
    while(pitr != pitr_end)
    {
        cout << "        :" << (*pitr).first << "\t" << (*pitr).second << endl;
        pitr ++;
    }
}                
//
bool get_allelic_ctg(string allelic_file)
{
    /*
         map<string, string> allelic_map ---- gloabl <"ctg1-ctg2", "0">
    */
    ifstream ifp;
    ifp.open(allelic_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open allelic table file " << allelic_file << endl;
        return false;
    }
    map<string, int> ctg_tmp;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        /*
            chr01   432509      utg000314l_pilon
            chr01   19202957    utg000007l_pilon    utg000148l_pilon    utg000336l_pilon        
            chr01   46480655    utg000193l_pilon    utg000079l_pilon
        */
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() <= 3) continue;
        for(int i=2; i<lineinfo.size()-1; i++)
        {
            for(int j=i+1; j<lineinfo.size(); j++)
            {
                string ctg1 = lineinfo[i];
                string ctg2 = lineinfo[j];
                string key1  = ctg1 + "-" + ctg2;
                string key2  = ctg2 + "-" + ctg1;
                if( allelic_map.find(key1) == allelic_map.end() )
                {
                    allelic_map.insert(std::pair<string, string>(key1, "0"));
                }
                if( allelic_map.find(key2) == allelic_map.end() )
                {
                    allelic_map.insert(std::pair<string, string>(key2, "0"));
                }
                //
                if(ctg_tmp.find(ctg1) == ctg_tmp.end())
                {
                    ctg_tmp.insert(std::pair<string, int>(ctg1, 0));
                }
                if(ctg_tmp.find(ctg2) == ctg_tmp.end())
                {
                    ctg_tmp.insert(std::pair<string, int>(ctg2, 0));
                }                
            }
        }
    }
    ifp.close();
    //
    cout << "   Info: "                            << allelic_map.size() 
         << " allelic relationship collected for " << ctg_tmp.size()
         << " contigs - gmap/gene based. "         << endl;
    //
    return true;
}
//
bool categorize_ctg(map<string, map<unsigned long, NODE> >  win_marker, 
                    double haplotig_win_ratio,
                    unsigned long min_hapctg_size,
                    map<string, map<unsigned long, NODE> >* hap_win_marker,
                    map<string, map<unsigned long, NODE> >* other_win_marker,
                    map<string, double>*                    ctg_hap_ratio)
{
    /* function: select long "pure" haplotigs: haplotig_win_ratio are hap-windows 
                 map<string, double> ctg_hap_ratio; giving ratio of hap-windows for all contigs
    */
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
        cout << "       : " << this_ctg << " hap_size = " << hap_size << " bp, ratio = " << this_hap_ratio << endl;        
        if(this_hap_ratio>=haplotig_win_ratio && this_ctg_size>=min_hapctg_size)
        {
            (*hap_win_marker).insert(std::pair<string, map<unsigned long, NODE> >(this_ctg, tmp_node));
            total_ctg_size_hap += this_ctg_size;
        }else
        {
            (*other_win_marker).insert(std::pair<string, map<unsigned long, NODE> >(this_ctg, tmp_node));            
        }
        (*ctg_hap_ratio).insert(std::pair<string, double>(this_ctg, this_hap_ratio));
        //
        raw_ctg_itr ++;
    }
    //
    cout << "   Info: "                        << win_marker.size() 
         << " contigs with a total size of "   << total_ctg_size 
         << " bp checked, among these: "       << endl;
    cout << "         "                        << (*hap_win_marker).size()  
         << " contigs with a total size of "   << total_ctg_size_hap
         << " bp selected as over "            << min_hapctg_size
         << " bp and with a hap win ratio >= " << haplotig_win_ratio*100
         << "%. "                              << endl;
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
   	        }
   	Note, here only contigs in given list of target_contig_size will be collected
    */
    ifstream ifp;
    ifp.open(win_marker_file.c_str(), ios::in);
    if(!ifp.good())
    {
        return false;
    }
    int num_win       = 0; // given group of contigs
    int num_win_total = 0; // full assembly
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
        num_win_total ++;
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
            tmp_node.grp_link_val.clear();
            tmp_node.grp_link_ctg.clear();
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
    cout << "   Info: " << target_contig_size.size() << " target ctgs given, "          << endl;
    cout << "       : " << num_win_total             << " windows from full assembly, " << endl;
    cout << "       : " << (*win_marker).size()      << " target ctgs collected with "  << num_win 
         << " windows (note, windows not belonging to given group were skipped). "      << endl;
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
//
bool read_hic_binning_lg(string                       hic_group_file, 
                         map<string, int>*            ctg2_hic_group,
                         map<int, map<string, int> >* hic_group_2ctg)
{
    /*  function: read cluster of haplotigs from a dot file
            ctg2_hic_group........: <ctg_id, clusterid_old>
            hic_group_2ctg........: <clusterid_old, <ctg_id, clusterid_new> >                
    */
    ifstream ifp;
    ifp.open(hic_group_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open gb group file " << hic_group_file << endl;
        return false;
    }
    int hic_grpid_new = 0;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        // "\tutg001319l_pilon -- utg000269l_pilon [color=gold, fontcolor=gold, penwidth=1, label=146.1]; /* cluster 6 */"
        // cout << "   check: " << line << endl;
        if(line.find("subgraph cluster") != std::string::npos)
        {
            hic_grpid_new ++;
        }
        if(line.find(" -- ") != std::string::npos)
        {
            vector<string> lineinfo2 = split_string(line, '\t');
            vector<string> lineinfo  = split_string(lineinfo2[0], ' ');            
            string vet1  = lineinfo[0];
            string vet2  = lineinfo[2];            
            int hic_grpid_old = atoi(lineinfo[9].c_str()); // clusterid_old   
            //cout << "   check: vet1=" << vet1 << " -- vet2=" << vet2 << ", gb group = " << hic_grpid_old << endl;
            if(!add_vet_to_group(vet1, 
                                 hic_grpid_old, 
                                 hic_grpid_new, 
                                 hic_group_2ctg, 
                                 ctg2_hic_group) )
            {
                return false;
            }
            if(!add_vet_to_group(vet2, 
                                 hic_grpid_old, 
                                 hic_grpid_new, 
                                 hic_group_2ctg, 
                                 ctg2_hic_group) )
            {
                return false;
            }            
        }        
    }
    ifp.close();
    //
    return true;
}
//
bool add_vet_to_group(string                       vet, 
                      int                          hic_grpid_old,
                      int                          hic_grpid_new,
                      map<int, map<string, int> >* hic_group_2ctg,
                      map<string, int>*            ctg2_hic_group)
{
    /*
        hic_grpid_old.........: old group id: any number
        hic_grpid_new.........: new group id: 1-4
        hic_group_2ctg........: <clusterid_old, <ctg_id, clusterid_new> >
    */    
    // update map from linkage group to list of contigs 
    if( (*hic_group_2ctg).find(hic_grpid_old) == (*hic_group_2ctg).end() )
    {
        map<string, int> tmpgroup;
        tmpgroup.insert(std::pair<string, int>(vet, hic_grpid_new));
        (*hic_group_2ctg).insert(std::pair<int, map<string, int> >(hic_grpid_old, tmpgroup));
    }else
    {
        if( (*hic_group_2ctg)[hic_grpid_old].find(vet) == (*hic_group_2ctg)[hic_grpid_old].end() )
        {        
            (*hic_group_2ctg)[hic_grpid_old].insert(std::pair<string, int>(vet, hic_grpid_new));
        }else
        {
            ; // redundant 
        }
    }
    // update map from a contig to its linkage group
    if((*ctg2_hic_group).find(vet) == (*ctg2_hic_group).end())
    {
        (*ctg2_hic_group).insert(std::pair<string, int>(vet, hic_grpid_old));
    }else
    {
        // should never happen
        if( (*ctg2_hic_group)[vet] != hic_grpid_old )
        {
            cout << "   Warning: contig found in group " << (*ctg2_hic_group)[vet] << " and " << hic_grpid_old << endl;
            cout << "            if it is non-haplotig markers, it is fine; otherwise, something wrong!"     << endl;
        }
    }
    //      
    return true;
}
bool get_options(int           argc, 
                char*          argv[], 
                double*        min_contact, 
                double*        min_contact_i,
                double*        haplotig_win_ratio,         
                unsigned long* min_hapctg_size,  
                double*        max_allelic_ratio,                     
                string*        ctg_file,
                string*        win_marker_file,
                string*        allelic_file,
                string*        align_file,
                string*        out_prefix)
{
    // log cmdline
    cout << "   CMDline: ";
    int ic = 0;
    while(ic < argc)
    {
        cout << argv[ic] << " ";
        ic ++;
    }
    cout << endl << endl;
    // get values: 
    ic = 1;    
    while (ic < argc) 
    {
        string optstr = (string)argv[ic];
        if(optstr.compare("--min-hic") == 0)
        {
            ic ++;
            *min_contact = atof(argv[ic]);
            cout << "   Info: contig contact graph will be built with minimum contact cutoff of: " 
                 << *min_contact 
                 << endl;
        }else
        if(optstr.compare("--min-hic-i") == 0)
        {
            ic ++;
            *min_contact_i = atof(argv[ic]);
            cout << "   Info: minimum contig contact to insert a contig into an existing group: " 
                 << *min_contact_i 
                 << endl;
        }else        
        if(optstr.compare("--hap-win-r") == 0)
        {
            ic ++;
            *haplotig_win_ratio = (double)atof(argv[ic]);
            cout << "   Info: ratio of haplotig windows along a unitig to create initial clusters: "
                 << *haplotig_win_ratio 
                 << endl;
        }else
        if(optstr.compare("--min-hap-size") == 0)
        {
            ic ++;
            *min_hapctg_size = strtoul(argv[ic], NULL, 0);
            cout << "   Info: min size (bp) of haplotigs to build backbone of linkage groups: "
                 << *min_hapctg_size 
                 << endl;
        }else
        if(optstr.compare("--max-allelic-r") == 0)
        {
            ic ++;
            *max_allelic_ratio = (double)atof(argv[ic]);
            cout << "   Info: ratio of a ctg overlapping a group larger this would be considered as allelic: "
                 << *max_allelic_ratio
                 << endl;
        }else
        if(optstr.compare("--ctg") == 0)
        {
            ic ++;
            *ctg_file = (string)argv[ic];
            ifstream fp;
            fp.open((*ctg_file).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: file of list of target contigs...........: " 
                     << *ctg_file 
                     << endl;
            }
            else
            {
                cout << "   Error: cannot open file of list of target contigs" 
                     << *ctg_file 
                     << endl;
                return false;
            }
        }else
        if(optstr.compare("--win-marker") == 0)
        {
            ic ++;
            *win_marker_file = (string)argv[ic];
            ifstream fp;
            fp.open((*win_marker_file).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: file of coverage-defined window markers..: " 
                     << *win_marker_file 
                     << endl;
            }
            else
            {
                cout << "   Error: cannot open file of coverage-defined window markers: " 
                     << *win_marker_file 
                     << endl;
                return false;
            }
        }else
        if(optstr.compare("--gmap-allelic") == 0)
        {
            ic ++;
            *allelic_file = (string)argv[ic];
            ifstream fp;
            fp.open((*allelic_file).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: file of allelic contigs defined with gmap: " 
                     << *allelic_file 
                     << endl;
            }
            else
            {
                cout << "   Error: cannot open file of allelic contigs defined with gmap: " 
                     << *allelic_file 
                     << endl;
                return false;
            }
        }else
        if(optstr.compare("--align-paf") == 0)
        {
            ic ++;
            *align_file = (string)argv[ic];
            ifstream fp;
            fp.open((*align_file).c_str(), ios::in);
            if(fp.good())
            {
                fp.close();
                cout << "   Info: file of alignment of contigs to reference: " 
                     << *align_file 
                     << endl;
            }
            else
            {
                cout << "   Error: cannot open file of alignment of contigs to reference: " 
                     << *align_file 
                     << endl;
                return false;
            }
        }else
        if(optstr.compare("-o") == 0)
        {
            ic ++;
            *out_prefix = (string)argv[ic];
            cout << "   Info: outputs will be labeled with \""    << *out_prefix << "\"" << endl;
        }else
        if(optstr.compare("-v") == 0)
        {
            verbose = true;
            cout << "   Info: be verbose - output some intermdiate info. " << endl;
        }        
        else
        {
            cout << "   Warning: option " << argv[ic] << " was not recognized and ignored."  << endl;
        }
        // next option
        ic ++;
    }
    
    return true;
}
