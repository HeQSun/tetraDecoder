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
//
double min_contact            = 25;    // hic contact of two markers
int    N_cluster              = 16;    // number of segmental regions - if 80 Mb, only cluster regions within 5 Mb - strong signal!
int    N_cluster_final        = 4;     // number of expected linkage groups: 1 * iploidy
int    iploidy                = 4;     // ploidy level; default 4 for tetraploid
double haplotig_win_ratio     = 0.999; // haplotig win ratio along contigs; used to select "pure"-haplotig markers
unsigned long min_hapctg_size = 50000; // min sizes of contigs to build backbone of chr-groups
map<string, CONTACT> ctg_cross_link_long_hap2; // main: < ctg1#ctg2, {LL, LR, RL, RR} >
map<string, CONTACT> ctg_cross_link_other2;    // main: < ctg1#ctg2, {LL, LR, RL, RR} >      
map<string, map<unsigned long, NODE> > hap_win_marker;
map<string, map<unsigned long, NODE> > other_win_marker;
map<string, string> allelic_map;       // allelic relationship between contigs <ctg1#ctg2, "0"> - gmap/gene based inferrence
map<string, string> allelic_map_seq;   // allelic relationship between contigs <ctg1#ctg2, "0"> - ref-align based inferrence
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
                  map<string, map<unsigned long, NODE> >*  other_win_marker);
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
                  map<string, ENDREGION>                   ctg_end_reg);          
bool get_clusters(map<string, int>                         vertex,
                  map<string, vector<string> >             edge,
                  map<string, double>                      cor,
                  map<int, map<string, double> >*          cluster_edge);
int find_smallest_cluster(map<int, unsigned long>          cluster_vertex_size,
                  map<int, int>                            cluster_cannot_be_merged);
bool merge_clusters(map<string, double>                    cor_any,
                  map<int, map<string, double> >           cluster_edge,
                  map<int, map<string, int> >              cluster_vertex,
                  map<int, unsigned long>                  cluster_vertex_size,
                  map<string, unsigned long>               target_contig_size,
                  string                                   tmpfolder);
int allelic_check_gene_based(map<string, int>              smallest_cluster_vertex,
                  map<string, int>                         larger_cluster_vertex);
//
int main(int argc, char* argv[])
{
    if(argc < 9)
    {
        cout << "\nFunction: analyze Hi-C contacts of tig markers, cluster markers and bin long reads. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl
             << "\nUsage: samtools view -F 3840 x_merged_chrxx_extract.bam"
             << " |"
             << " hic_binning - contact_cutoff target_contig_size.txt window_markers.txt min_hapctg_size allelic.table.txt align_file out_prefix_str"
             << endl 
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
             << "                             note, contact at ends of paired-contigs also calculated at this resolution. "                         
             << endl
             << "      *allelic.table.txt  -- allelic conditions among contigs. "
             << endl
             << "      *align_file         -- alignment of contigs to DM reference . "
             << endl
             << "      *output_prefix      -- flag of output file." << endl 
             << endl;
        cout << "   Note-S3: to filter    dot: dot_filter_v2 xxx.dot 0.75"               << endl;
        cout << "   Note-S3: to visualize dot: circo -Tpdf xxx.dot > xxx_circo.pdf, or " << endl;
        cout << "                                fdp -Tpdf xxx.dot > xxx_fdp.pdf . "     << endl;
        cout << endl;
        //             -f 0x2              -- only properly paired read alignments. 
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
    min_contact            =    atof(argv[2]);
    cout << "   Info: contig contact graph will be built with minimum contact cutoff of: " << min_contact     << endl;
    string ctg_file        = (string)argv[3];  
    cout << "   Info: target contig list file: "                                           << ctg_file        << endl;    
    string win_marker_file =  string(argv[4]);
    cout << "   Info: coverage-defined window marker file: "                               << win_marker_file << endl;
    min_hapctg_size        = strtoul(argv[5], NULL, 0);
    cout << "   Info: min size in bp of haplotigs to build backbone of linkage groups "    << min_hapctg_size << endl; 
    string allelic_file    = (string)argv[6];
    cout << "   Info: allelic information provided in "                                    << allelic_file    << endl;
    string align_file      = (string)argv[7];
    cout << "   Info: alignment of contigs to (DM) reference provided in "                 << align_file      << endl;   
    string out_prefix      = (string)argv[8];    
    cout << "   Info: output will be labeled with "                                        << out_prefix      << endl;
    cout << endl;
    // step s1: read target contig and size info
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
    map<string, map<unsigned long, NODE> > win_marker;
    if(!read_win_marker_ctg(win_marker_file, target_contig_size, &win_marker) )
    {
        cout << "   Error: failed in reading coverage-defined window markers. " << endl;
        return 1;
    }
    cout << "   Info: step 2 - read info on window markers done. " << endl;        
    cout << endl;        
    // step s3: select long haplotigs: set up gloabl variables: hap_win_marker and other_win_marker
    cout << "   Info: step 3 - select contigs over " << min_hapctg_size 
         << " bp with "                              << haplotig_win_ratio  
         << " hap-windows. "                         << endl;        
    if(!categorize_ctg(win_marker, haplotig_win_ratio, min_hapctg_size, &hap_win_marker, &other_win_marker))
    {        
        cout << "   Error: failed in categorizing contigs. " << endl;
        return 1;
    }
    cout << "   Info: step 3 - select long haplotigs done. " << endl;
    cout << endl;            
    // step s4. read allelic contigs 
    cout << "   Info: step 4 - reading allelic relationship of contigs..." << endl; // => allelic_map=<"ctg1#ctg2", "0">
    if(!get_allelic_ctg(allelic_file))
    {        
        cout << "   Error: failed in reading allelic map. "  << endl;
        return 1;
    }
    if(!get_allelic_ctg_align_based(align_file, target_contig_size))
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
        string key1 = this_a_ctg + "#" + this_b_ctg;
        string key2 = this_b_ctg + "#" + this_a_ctg;          
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
        // left/right region of contig a
        unsigned long left_sta_a_ctg = 1;
        unsigned long left_end_a_ctg = min_hapctg_size;
        if(left_end_a_ctg > length_a_ctg)
        {
            left_end_a_ctg = length_a_ctg;
        }
        unsigned long right_sta_a_ctg = 1;
        unsigned long right_end_a_ctg = length_a_ctg;        
        if(length_a_ctg > min_hapctg_size)
        {
            right_sta_a_ctg = length_a_ctg - min_hapctg_size + 1;
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
        unsigned long left_end_b_ctg = min_hapctg_size;
        if(left_end_b_ctg > length_b_ctg)
        {
            left_end_b_ctg = length_b_ctg;
        }
        unsigned long right_sta_b_ctg = 1;
        unsigned long right_end_b_ctg = length_b_ctg;        
        if(length_b_ctg > min_hapctg_size)
        {
            right_sta_b_ctg = length_b_ctg - min_hapctg_size + 1;
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
         << " collected into given group: "                      << endl
         << "         a4="                                       << num_pair_inter_long_hap  
         << " R1-R2 cross-aligned to two long haplotigs, over "  << min_hapctg_size 
         << " bp with "                                          << haplotig_win_ratio
         << " as haplotig windows; "                             << endl
         << "         a5="                                       << num_pair_other 
         << " R1-R2 aligned as other cases. "                    << endl;
    cout << "         Note, a0=a1+a2+a3, where a3=a4+a5."        << endl; 
    cout << "   Info: step 5 - traverse bam file to get Hi-C link  done. " << endl;
    cout << endl;  
    // step s6. output raw hic links among haplotigs and other cases   
    cout << "   Info: step 6 - output raw Hi-C link info (for checking purpose)... " << endl;    
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
        cout << "   Error: cannot output raw Hi-C link info for other cases. " << endl;
        return 1;
    }
    if(!output_raw_links(ctg_cross_link_long_hap2,
                         target_contig_size,
                         tmpfolder_s6, 
                         "hap",
                         ctg_end_region) )
    {
        cout << "   Error: cannot output raw Hi-C link info for long haplotigs. " << endl;
        return 1;
    }
    cout << "   Info: step 6 - output raw Hi-C link info (for checking purpose) done. " << endl;        
    cout << endl;
    // step s7. cluster haplotigs
    cout << "   Info: step 7 - initial clustering of haplotigs using Hi-C link... " << endl;    
    if(!cluster_hap_links(ctg_cross_link_long_hap2, 
                         target_contig_size,
                         tmpfolder_s6,
                         "hap",
                         ctg_end_region) )
    {
        cout << "   Error: failed in haplotig clustering. " << endl;
        return 1;
    }
    cout << "   Info: step 7 - initial clustering of haplotigs using Hi-C link done. " << endl;        
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    //
    return 0;
}
//
bool cluster_hap_links(map<string, CONTACT>       hic_cross_link, 
                       map<string, unsigned long> target_contig_size,
                       string                     tmpfolder,
                       string                     this_label,
                       map<string, ENDREGION>     ctg_end_region)
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
    map<int, map<string, double> > cluster_edge;        // <clusterid, <ctg1-ctg2, hic.contact.value > >
    map<int, map<string, int> >    cluster_vertex;      // <clusterid, <ctg, 1 > >    
    map<int, unsigned long>        cluster_vertex_size; // <clusterid, total_ctg_length >   
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
            vector<string> keyinfo = split_string(this_key, '#');
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
            double normalized_reads = this_cnt*1.0 / (ctg1_end_size+ctg2_end_size) * 100000;  // reads per 100 kb      
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
        cout << "   Info: " << unique_contigs.size() << " raw unique contigs checked. " << endl;    
        // ss1: initial clustering of haplotigs        
        if(this_label.compare("hap")==0)
        {
            cout << "   Info: " << vertex_all.size()     << " vertices selected with "
                 <<                cor_high.size()        << " edges. " 
                 << endl;  
            // get sub-clusters = backbone of potential linkage groups
            if( !get_clusters(vertex_all, edge_all, cor_high, &cluster_edge) )
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
                    subcluster_ofp << "\t" 
                                   << (*vvitr).first 
                                   << " [color=orangered, style=filled, fillcolor=orangered, fontcolor=white];"
                                   << endl;
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
    cout << "   Info: " << cor_any.size()  << " edges collected, among these, " << endl
         << "         " << cor_high.size() << " edges showing contact > "       << min_contact << endl;        
    if(0)
    {
        //cout << "   Info: more clusters than expected number of " << N_cluster << endl
        //     << "         merge of clusters would be performed. " << endl;
        if(!merge_clusters(cor_any,
                           cluster_edge,
                           cluster_vertex,
                           cluster_vertex_size,
                           target_contig_size,
                           tmpfolder))
        {
            cout << "   Error: failed in merging clusters. " << endl;
            return false;
        }
    }else
    {
        cout << "   Info: number of clusters = expected number of "     << N_cluster << endl;        
    }      
    //  
    return true;
}
bool merge_clusters(map<string, double>            cor_any,
                    map<int, map<string, double> > cluster_edge,
                    map<int, map<string, int> >    cluster_vertex,   
                    map<int, unsigned long>        cluster_vertex_size, 
                    map<string, unsigned long>     target_contig_size,
                    string                         tmpfolder)
{
    /*
        Inputs:
            map<string, double>            cor_any;             // <ctg1-ctg2, hic.contact.value >    
            map<int, map<string, double> > cluster_edge;        // <clusterid, <ctg1-ctg2, hic.contact.value > >
            map<int, map<string, int> >    cluster_vertex;      // <clusterid, <ctg, 1> >    
            map<int, unsigned long>        cluster_vertex_size; // <clusterid, total_ctg_length >               
    */
    map<int, vector<int> > merging_record; // <clusterid, <merged-subcluster-ids> >         
    //
    int cluster_after_merge = cluster_edge.size();
    map<int, int> cluster_cannot_be_merged; // some smaller cluster possibly cannot be merged..    
    while(cluster_after_merge > N_cluster)
    {
        cout << "   Info: ---- number of current clusters " << cluster_edge.size() 
             << " > expected cluster number "               << N_cluster 
             << " (with "                                   << cluster_cannot_be_merged.size()
             << " cannot be merged). "                      << endl;
        cout << "         ---- smallest clusters will be merged into larger ones. "   
             << endl;
        cout << "   Info: ---- merging clusters..."    << endl;   
        //
        int smallest_cluster_id = find_smallest_cluster(cluster_vertex_size, cluster_cannot_be_merged);                    
        assert( cluster_vertex.find(smallest_cluster_id) != cluster_vertex.end() );
        map<string, int> smallest_cluster_vertex = cluster_vertex[smallest_cluster_id];
        // initialize largest correlation cluster
        int    best_cluster_id  = smallest_cluster_id;
        double best_cluster_cor = 0;
        string best_connect_edge= "";
        // record contigs in the smaller cluster connect to larger clusters
        map<int, vector<string> > best_larger_cluster_id;    // <best_larger_id, <smaller_ctgs> >
        map<int, vector<double> > best_larger_cluster_cor;   // <best_larger_id, <smaller_larger_edge_score> >
        map<int, vector<string> > best_larger_cluster_edge;  // <best_larger_id, <smaller_larger_edge> >          
        // check best correlation by traversing all vertices
        map<string, int>::iterator svitr;
        map<string, int>::iterator svitr_end;
        svitr     = smallest_cluster_vertex.begin();
        svitr_end = smallest_cluster_vertex.end();
        map<int, int> allelic_group;
        bool can_be_merged = false;                    
        while(svitr != svitr_end)
        {
            // get a contig id in the smallest cluster
            string sctg = (*svitr).first;
            // clusters of vertices
            map<int, map<string, int> >::iterator vclitr;
            map<int, map<string, int> >::iterator vclitr_end;
            vclitr     = cluster_vertex.begin();
            vclitr_end = cluster_vertex.end();
            while(vclitr != vclitr_end)
            {
                int              larger_cluster_id     = (*vclitr).first;
                map<string, int> larger_cluster_vertex = (*vclitr).second;    
                if(allelic_group.find( larger_cluster_id ) != allelic_group.end() )
                {
                    // next larger cluster 
                    vclitr ++;       
                    continue;                    
                }            
                // check with allelic table, if the smaller and larger clusters should avoid being merged
                int allelic_cnt = allelic_check_gene_based(smallest_cluster_vertex, larger_cluster_vertex); 
                if(allelic_cnt >= 1)
                {
                    cout << "   Info: found allelic contigs between groups " << smallest_cluster_id
                         << " and "                                          << larger_cluster_id
                         << endl;
                    allelic_group.insert(std::pair<int, int>(larger_cluster_id, smallest_cluster_id));
                    // next larger cluster 
                    vclitr ++;       
                    continue;             
                }else
                if( larger_cluster_id != smallest_cluster_id )
                {
                    map<string, int>::iterator lvitr;
                    map<string, int>::iterator lvitr_end;
                    lvitr     = larger_cluster_vertex.begin();
                    lvitr_end = larger_cluster_vertex.end();                    
                    while(lvitr != lvitr_end)
                    {
                        // get a contig id in current larger cluster
                        string lctg = (*lvitr).first;
                        // check hic contacts
                        string edge1 = sctg + "-" + lctg;
                        string edge2 = lctg + "-" + sctg;
                        if( cor_any.find(edge1) != cor_any.end() )
                        {
                            if(cor_any[ edge1 ] > best_cluster_cor)
                            {
                                best_cluster_cor  = cor_any[ edge1 ];
                                best_cluster_id   = larger_cluster_id;
                                best_connect_edge = edge1;
                                can_be_merged     = true;
                            } 
                        }else
                        if( cor_any.find(edge2) != cor_any.end() )
                        {
                            if(cor_any[ edge2 ] > best_cluster_cor)
                            {
                                best_cluster_cor  = cor_any[ edge2 ];
                                best_cluster_id   = larger_cluster_id;
                                best_connect_edge = edge2;
                                can_be_merged     = true;
                            } 
                        }else ;                  
                        // next contig in current larger cluster
                        lvitr ++;                    
                    }
                }else ;
               // next larger cluster 
                vclitr ++;
            }
            // next contig in the smallest cluster 
            svitr ++;        
        }
        if(can_be_merged == false)
        {
            cout << "   warning: small group "             << smallest_cluster_id 
                 << " cannot be merged into larger ones. " << endl;
            cluster_after_merge --;
            cluster_cannot_be_merged.insert(std::pair<int, int>(smallest_cluster_id, 1));
        }        
        // merge the smaller cluster into a larger group according to best correlation value
        if(smallest_cluster_id != best_cluster_id)
        {
            cout << "   Info: ---- cluster "              << smallest_cluster_id 
                 << " will be merged into cluster "  << best_cluster_id << endl;
            // update edge info 
            cout << "   Info: ---- cluster " 
                 << best_cluster_id 
                 << " with " 
                 << cluster_edge[best_cluster_id].size() 
                 << " edges updated as ";
            cluster_edge[best_cluster_id].insert(cluster_edge[smallest_cluster_id].begin(), 
                                                 cluster_edge[smallest_cluster_id].end());
            // add the "bridging" edge
            cluster_edge[best_cluster_id].insert(std::pair<string, double>(best_connect_edge, best_cluster_cor));
            cout << cluster_edge[best_cluster_id].size() 
                 << " edges. " 
                 << endl;
            // remove edges of smallest group
            cluster_edge.erase(smallest_cluster_id);                 
            // update vertex info 
            cout << "   Info: ---- cluster " 
                 << best_cluster_id 
                 << " with " 
                 << cluster_vertex[best_cluster_id].size() 
                 << " vertices updated as ";            
            cluster_vertex[best_cluster_id].insert(cluster_vertex[smallest_cluster_id].begin(),
                                                   cluster_vertex[smallest_cluster_id].end());
            cout << cluster_vertex[best_cluster_id].size() 
                 << " vertices. " 
                 << endl;    
            // remove vertices of smallest group 
            cluster_vertex.erase(smallest_cluster_id);                         
            // update cluster size info 
            cout << "   Info: ---- cluster " 
                 << best_cluster_id 
                 << " with size of " 
                 << cluster_vertex_size[best_cluster_id] 
                 << " bp updated as ";              
            cluster_vertex_size[best_cluster_id] += cluster_vertex_size[smallest_cluster_id];
            cout << cluster_vertex_size[best_cluster_id] 
                 << " bp. " 
                 << endl;   
            // remove sizes of smallest group
            cluster_vertex_size.erase(smallest_cluster_id);  
            // record how merging happened
            if(merging_record.find(smallest_cluster_id) != merging_record.end() )
            {
                if(merging_record.find(best_cluster_id) == merging_record.end() )
                {
                    merging_record.insert(std::pair<int, vector<int> >(best_cluster_id, merging_record[smallest_cluster_id] ));
                    merging_record[best_cluster_id].push_back(smallest_cluster_id);
                }else
                {
                    merging_record[best_cluster_id].insert(merging_record[best_cluster_id].end(),
                                                           merging_record[smallest_cluster_id].begin(),
                                                           merging_record[smallest_cluster_id].end());
                    merging_record[best_cluster_id].push_back(smallest_cluster_id);
                }               
            }else
            {
                if(merging_record.find(best_cluster_id) == merging_record.end() )
                {
                    vector<int> tmpvec;
                    tmpvec.push_back(smallest_cluster_id);
                    merging_record.insert(std::pair<int, vector<int> >(best_cluster_id, tmpvec));
                }else
                {
                    merging_record[best_cluster_id].push_back(smallest_cluster_id);
                }
            } 
            //
            cluster_after_merge --;
        }else
        {
            cout << "   warning: cannot find a proper target cluster to merge current cluster " 
                 << smallest_cluster_id
                 << ", no action." 
                 << endl;
        }
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
            subcluster_ofp << "\t" 
                           << (*vvitr).first 
                           << " [color=orangered, style=filled, fillcolor=orangered, fontcolor=white];"
                           << endl;
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
    //
    return true;
}                
//
int allelic_check_gene_based(map<string, int> smallest_cluster_vertex, 
                  map<string, int> larger_cluster_vertex)
{
    /* function: check if contigs in small list having allelic relationship with those in large list, with 
       map<string, string> allelic_map --- <ctg1#ctg2, "0">
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
            string key1  = ctg_s + "#" + ctg_l;
            string key2  = ctg_l + "#" + ctg_s;
            if(allelic_map.find(key1)!=allelic_map.end() || allelic_map.find(key2)!=allelic_map.end() )
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
int find_smallest_cluster(map<int, unsigned long> cluster_vertex_size,
                          map<int, int> cluster_cannot_be_merged)
{
    /*
        map<int, int> cluster_cannot_be_merged; // some smaller cluster possibly cannot be merged..    
    */
    int           smallest_cluster_id   = 0;
    unsigned long smallest_cluster_size = 999999999;
    map<int, unsigned long>::iterator citr;
    map<int, unsigned long>::iterator citr_end;
    citr     = cluster_vertex_size.begin();
    citr_end = cluster_vertex_size.end();
    while(citr != citr_end)
    {
        if(cluster_cannot_be_merged.find((*citr).first) != cluster_cannot_be_merged.end() )
        {
            citr ++;
            continue;
        }
        if( (*citr).second < smallest_cluster_size )
        {
            smallest_cluster_size = (*citr).second;
            smallest_cluster_id   = (*citr).first;
        }
        citr ++;
    }
    return smallest_cluster_id;
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
        vector<string> keyinfo = split_string(this_key, '#');
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
        double normalized_reads = this_cnt*1.0 / (ctg1_size+ctg2_size) * 100000;  // reads per 100 kb  
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
    cout << "   Info: " << unique_contigs.size() << " unique contigs in " << siminfo << endl;
    //  
    return true;
}
//
bool get_clusters(map<string, int> vertex,
                  map<string, vector<string> > edge,
                  map<string, double> cor,
                  map<int, map<string, double> >* cluster_edge
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
    cout << "   Info: " << (*cluster_edge).size() << " clusters found. " << endl;
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
	map<string, CONTACT> ctg_cross_link_long_hap2;           // global: < ctg1#ctg2, {LL, LR, RL, RR} >
	map<string, CONTACT> ctg_cross_link_other2;              // global: < ctg1#ctg2, {LL, LR, RL, RR} >      
	map<string, map<unsigned long, NODE> > hap_win_marker;   // global
	map<string, map<unsigned long, NODE> > other_win_marker; // global
    */
    string key1 = this_a_ctg + "#" + this_b_ctg;
    string key2 = this_b_ctg + "#" + this_a_ctg;
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
    cout << "   check: cleaning alignments... " << endl;    
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
        //
        citr ++;
    }
    cout << "   check: cleaning alignments done. " << endl;            
    // s3. analyze allelic relationship among contigs 
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
            if(overlap_size*1.0 > target_contig_size[ctg_1]*0.3 || 
               overlap_size*1.0 > target_contig_size[ctg_2]*0.3 )
            {
                cout << "   check: "      << ctg_1 
                     << " and "           << ctg_2 
                     << " overlap by "    << overlap_size 
                     << " bp, with size " << target_contig_size[ctg_1]
                     << " vs "            << target_contig_size[ctg_2]
                     << endl;
                string key1  = ctg_1 + "#" + ctg_2;
                string key2  = ctg_2 + "#" + ctg_1;
                if( allelic_map_seq.find(key1) == allelic_map_seq.end() )
                {
                    allelic_map_seq.insert(std::pair<string, string>(key1, "0"));
                }
                if( allelic_map_seq.find(key2) == allelic_map_seq.end() )
                {
                    allelic_map_seq.insert(std::pair<string, string>(key2, "0"));
                }
            }
            //
            citr_2nd ++;
        }
        //
        citr_1st ++;
    }
    //
    return true;
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
         map<string, string> allelic_map ---- gloabl <"ctg1#ctg2", "0">
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
                string key1  = ctg1 + "#" + ctg2;
                string key2  = ctg2 + "#" + ctg1;
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
