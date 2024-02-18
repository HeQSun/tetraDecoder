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
double min_contact            =  0.50;  // correlation of genotype sequences at two markers defined by sequencing depth 
int    N_cluster              =  4;     // number of expected linkage groups: 1 * iploidy
int    iploidy                =  4;     // ploidy level; default 4 for tetraploid
double haplotig_win_ratio     =  0.95;  // haplotig win ratio along contigs; used to select "pure"-haplotig markers
unsigned long min_hapctg_size = 50000;  // min sizes of contigs to build backbone of chr-groups
map<string, CONTACT> ctg_cross_link_long_hap2;      // main: < ctg1#ctg2, {LL, LR, RL, RR} >
map<string, CONTACT> ctg_cross_link_other2;         // main: < ctg1#ctg2, {LL, LR, RL, RR} >      
map<string, map<unsigned long, NODE> > hap_win_marker;
map<string, map<unsigned long, NODE> > other_win_marker;
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
void update_links(string this_a_ctg, 
                  string this_b_ctg,
                  string contact_type,                                 // LL, LR, RL, RR
                  unsigned long* num_pair_inter_long_hap,
                  unsigned long* num_pair_other
                  );
bool output_raw_links(map<string, CONTACT> hic_cross_link, 
                      map<string, unsigned long> target_contig_size,
                      string                     tmpfolder,
                      string                     this_label,
                      map<string, ENDREGION> ctg_end_region);                   
//
int main(int argc, char* argv[])
{
    if(argc < 7)
    {
        cout << "\nFunction: analyze Hi-C contacts of tig markers, cluster markers and bin long reads. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl
             << "\nUsage: samtools view -F 3840 x_merged_chrxx_extract.bam"
             << " |"
             << " hic_binning - contact_cutoff target_contig_size.txt window_markers.txt min_hapctg_size out_prefix_str"
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
    map<string, ENDREGION> ctg_end_region;              // main: < ctg, {ctg_left_sta/end, ctg_right_sta/end} >
    unsigned long numraw                   = 0;
    unsigned long num_pair_skip            = 0;
    unsigned long num_pair_inter_long_hap  = 0;
    unsigned long num_pair_other           = 0;
    unsigned long num_pair_intra           = 0;
    unsigned long num_pair_not_at_ctg_ends = 0;
    unsigned long lenNoCG  = 0; // no cigar string unmapped - length
    unsigned long numNoCG  = 0; // no cigar string unmapped - number
    string exampleNoCG     = "";     
    string last_read_name  = "";// avoid repeating counting of the same read pair on the same pair of contigs.        
    string last_ctg_pair   = "";
    std::string read_a; // first  row
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
        // skip read pairs mapped within a contig.
        if(this_b_ctg.compare("=") == 0)
        {
            num_pair_intra ++;
            num_pair_skip ++;            
            this_b_ctg = this_a_ctg;
            continue;
        }
        // only when "this_a_ctg!=this_b_ctg", the followings needed!
        assert( this_a_ctg.compare(this_b_ctg) != 0 )
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
        string key1 = this_a_ctg + "#" + this_b_ctg;
        string key2 = this_b_ctg + "#" + this_a_ctg;        
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
        int inter_ll = 0;
        int inter_lr = 0;
        int inter_rl = 0;
        int inter_rr = 0;
        int other_ll = 0;
        int other_lr = 0;
        int other_rl = 0;
        int other_rr = 0;
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
         << " reads skipped (including a1.1="                    << num_pair_intra
         << " intra-contig links and a1.2="                      << num_pair_not_at_ctg_ends
         << " not in end-regions of contigs), "                  << endl
         << "         a2="                                       << num_pair_kept            
         << " collected into given group: "                      << endl
         << "         a3="                                       << num_pair_inter_long_hap  
         << " R1-R2 cross-aligned to two long haplotigs, over "  << min_hapctg_size 
         << " bp with "                                          << haplotig_win_ratio
         << " as haplotig windows; "                             << endl
         << "         a4="                                       << num_pair_other 
         << " R1-R2 aligned to other cases. "                    << endl;
    cout << "         Note, a0=a1+a2, where a2=a3+a4."           << endl; 
    cout << "   Info: step 4 - traverse bam file to get Hi-C link  done. " << endl;
    //
    cout << "   Info: step 5 - output raw Hi-C link info (for checking purpose)... " << endl;    
    string tmpfolder_s5 = "s5_"+out_prefix+"_raw_tig_marker_cross_link_count";
    if(!create_folder(tmpfolder_s5))
    {
        return 1;
    } 
    if(!output_raw_links(ctg_cross_link_other2,
                         target_contig_size,
                         tmpfolder_s5, 
                         "other",
                         ctg_end_region) )
    {
        cout << "   Error: cannot output raw Hi-C link info for other cases. " << endl;
        return 1;
    }
    if(!output_raw_links(ctg_cross_link_long_hap2,
                         target_contig_size,
                         tmpfolder_s5, 
                         "hap",
                         ctg_end_region) )
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
void update_links(string this_a_ctg, 
                  string this_b_ctg,
                  string contact_type,                                 // LL, LR, RL, RR
                  unsigned long* n_pair_inter_long_hap,
                  unsigned long* n_pair_other
                  )
{
    /*
	map<string, CONTACT> ctg_cross_link_long_hap2;      // main: < ctg1#ctg2, {LL, LR, RL, RR} >
	map<string, CONTACT> ctg_cross_link_other2;         // main: < ctg1#ctg2, {LL, LR, RL, RR} >      
	map<string, map<unsigned long, NODE> > hap_win_marker;
	map<string, map<unsigned long, NODE> > other_win_marker;    
    */
    //
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
bool output_raw_links(map<string, CONTACT> hic_cross_link, 
                      map<string, unsigned long> target_contig_size,
                      string                     tmpfolder,
                      string                     this_label,
                      map<string, ENDREGION> ctg_end_region)
{
    /* functon: output raw Hi-C links for checking purpose */
    if(hic_cross_link.size() == 0)
    {
        cout << "   Error: size of hic_cross_link is 0 - no cross-linking found between contigs? " << endl;
        return false;
    }
    //
    map<string, int> unique_contigs;
    //
    string siminfo = tmpfolder + "/s5_cross_link_" + this_label + "_raw_full.dot\0"; // 
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
        int normalized_reads = this_cnt*1.0 / (ctg1_size+ctg2_size) * 100000;  // reads per 100 kb  
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
    //
    cout << "   Info: " << unique_contigs.size() << " unique contigs in " << siminfo << endl;
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
