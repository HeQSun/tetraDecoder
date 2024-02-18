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
   s4:
   marker with limited pollen depth......: s4p1_all_raw_unlinked_markers_with_limited_pollen_support.bed
   raw tetraplotig markers not groupped..: s4p2_all_raw_unlinked_tetraplotig_markers.bed
   raw grouped markers...................: s4p3_all_raw_linked_hap_dip_triplotig_markers.bed
   raw replotig markers not groupped.....: s4p4_all_raw_rep_markers.bed
   tmp on re-linking non-grouped-markers.: s4p5_refine_grouping_tmp.txt
   final grouping of all possible markers: s4p6_refine_grouping_final_window_markers.txt
#
Written by Hequan Sun, MPIPZ (Germany)
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
struct NODE 
{
    unsigned long           depth; // raw depth at this marker
    double              depth_nor; // normalized depth, rpkm    
    string                   type; // hap/dip/trip/tetrap/rep  
    unsigned long            size; // size of this marker
    unsigned long        ctg_size; // size of the contig having this marker
    unsigned long        genotype; // genotype-like value determined by depth
};
struct NODE2
{
    vector<int>          depth_gt; // at a single window marker
    string                   type; // hap/dip/trip/tetrap/rep
    map<long, string> highest_cor; // several markers with the highest correlation to the current marker.
};
struct BIPAT
{
    string                leftPat;
    string               rightPat;
    int                   win_num; // window markers along this contig leading to these genotype patterns - not used!
};
struct BIPAT_INT
{
    vector<int>       leftPat_int;
    vector<int>      rightPat_int;
    int                   win_num;// not used yet!
};
struct lgGROUP
{
    vector<string>           type; // type of current window marker: hap/dip/trip/tetrap...
    vector<int>               qlg; // LG of query contig = 1:48
    vector<int>              qchr; // Chr of query contig = 1:12, i.e., homocluster_id
    vector<string>           sctg; // id of linked subject contig
    vector<int>               slg; // LG of linked subject contig = 1:48
    vector<double>          qscor; // correlation value linked query and subject contigs
    vector<string>          qscas; // case info: how this linking is found
    vector<unsigned long>   msize; // size of current window marker -- not used yet!    
};
double minCOscore         =  0.64; // crossover score 0.7: 0.83666 * 0.83666; 0.5:0.7071068*0.7071068
double minCorscore        =  0.50; // correlation of genotype sequences at two markers defined by sequencing depth 
int    N_cluster          =    48; // number of expected linkage groups: 12 * iploidy
int    iploidy            =     4; // ploidy level; default 4 for tetraploid
double haplotig_win_ratio =   0.9; // haplotig win ratio along contigs; used to select pure-haplotig markers-TODO option
double hom_cor_cutoff     = -0.25; // maximum correlation for finding "homologous" haplotype-specific linkage groups.
                                   // note: "homologous" haplotype-specific LGs expected to be with negative correlation.                                   
int    min_hapctg_size    = 15000; // minimum size of haplotigs to be selected as initial pure haplotig markers for 
                                   // for building up the backbone of linkage groups
bool   lg_interctg_cor    = false; // turn on to recalculate intra-lg inter-contig correlation
// 
bool create_folder(string folder);
bool get_individual_bed_file(string                  list_bed_files, 
                  vector<string>*                    bedfiles);   
bool read_current_bed(string                         this_bed, 
                  map<string, NODE>*                 this_pollen_depth,
                  map<string, unsigned long>*        ctg_size_all,
                  vector<string>*                    marker_order);   
bool normalize_read_count(map<string, NODE>*         this_pollen_depth);  
double find_haploid_mean_depth(map<string, NODE>     this_pollen_depth);   
string get_break_pos(string                          depth_gt, 
                  string                             pmflag,
                  int                                pollenid,
                  unsigned long*                     leftp, 
                  double*                            score,
                  std::stringstream*                 ss,                     
                  string                             contigid,
                  bool                               output_ss);
string smooth_depth_gt(string                        depth_gt); // note: naively smooth only those >=10 bits     
double calculate_correlation(vector<int>*            marker_x, 
                  vector<int>*                       marker_y);  
bool greedy_cluster_haplotigs(map<string, BIPAT_INT> contigPollenGTSeq_int, 
                  map<string, unsigned long>         contigsize,
                  string                             tmpfolder,
                  map<int, int>*                     gb_group_id,
                  map<int, map<string, int> >*       gb_group);
bool get_clusters(map<string, int>                   vertex,
                  map<string, vector<string> >       edge,
                  map<string, double>                cor,
                  map<int, map<string, double> >*    cluster_edge);
bool collect_contig_size(string                      ctgsize_file, 
                  map<string, unsigned long>*        contigsize);   
bool break_clusters(map<int, map<string, int> >      cluster_vertex,  
                  map<string, BIPAT_INT>             contigPollenGTSeq_int,
                  string                             tmpfolder,
                  map<string, unsigned long>         contigsize,
                  map<int, map<string, double> >*    cluster_edge_updated,
                  map<int, unsigned long>*           cluster_vertex_size_updated,
                  map<int, map<string, int> >*       cluster_vertex_updated,
                  map<int, int>*                     mutually_exclusive_cluster);
bool calc_inter_clusters_cor(map<int, map<string, int> > cluster_vertex,  
                  map<int, unsigned long>            cluster_vertex_size_updated,
                  map<string, BIPAT_INT>             contigPollenGTSeq_int,
                  string                             tmpfolder,
                  map<string, unsigned long>         contigsize,
                  map<int, map<string, double> >*    cluster_edge_updated2,
                  map<int, unsigned long>*           cluster_vertex_size_updated2,
                  map<int, map<string, int> >*       cluster_vertex_updated2,
                  map<int, int>*                     mutually_exclusive_cluster);                   
int find_contig_set_overlap(map<string, double>      map1, 
                  map<string, double>                map2);                          
bool merge_clusters(map<int, map<string, double> >   cluster_edge, 
                  map<int, map<string, int> >        cluster_vertex,  
                  map<int, unsigned long>            cluster_vertex_size,
                  map<int, int>                      mutually_exclusive_cluster,
                  map<string, BIPAT_INT>             contigPollenGTSeq_int,  
                  map<string, unsigned long>         contigsize,
                  int                                N_cluster,
                  string                             tmpfolder,
                  map<int, int>*                     gb_group_id,
                  map<int, map<string, int> >*       gb_group);                  
int find_smallest_cluster(map<int, unsigned long>    cluster_vertex_size);     
bool get_pure_haplotig_markers(string                bed_file, 
                  map<string, unsigned long>         contigsize,
                  map<string, double>*               hap_win_ratio);   
// below is to find "homologous haplotype-specific LGs" (4 x 12 in potato Otava)                                   
bool find_hom_group(map<string, BIPAT_INT>           contigPollenGTSeq_int, 
                  map<int, map<string, int> >        gb_group,
                  map<int, int>                      gb_group_id,
                  double                             hom_cor_cutoff,
                  string                             tmpfolder);   
bool get_lg_clusters(map<string, vector<string> >    edge,
                  map<string, double>                hom_pairs_simple,
                  map<int, map<string, double> >*    cluster_edge,
                  map<int, map<string, int> >*       cluster_vertex);                  
bool read_gamete_binning_lg(string                   gb_group_file, 
                  map<string, int>*                  ctg2_gb_group,
                  map<int, map<string, int> >*       gb_group_2ctg);
bool add_vet_to_group(string                         vet, 
                  int                                gb_grpid,
                  int                                gb_grpid_new,
                  map<int, map<string, int> >*       gb_group_2ctg,
                  map<string, int>*                  ctg2_gb_group); 
bool read_gb_homo_lg(string                          homo_lg_file,
                  map<int, map<string, int> >*       cluster_id_2lg,
                  map<string, int>*                  lg2_cluster_id); 
bool get_depth_AND(vector<int>*                      marker_x, 
                  vector<int>*                       marker_y, 
                  vector<int>*                       dand);         
double calculate_correlation_diplotig(vector<int>*   dip_marker_x, 
                  vector<int>*                       hap_marker_y);
double calculate_correlation_triplotig(vector<int>*  trip_marker_x, 
                  vector<int>*                       hap_marker_y);
bool integrate_dtt_markers_to_lg(map<string, NODE2>* marker_depth_all,
                  vector<string>                     marker_order_keep,
                  map<string, double>                hap_win_ratio,
                  map<int, int>                      gb_group_id,
                  map<int, int>                      gb_group_id_rev,
                  map<int, map<string, int> >        gb_group_2ctg,
                  map<string, int>                   ctg2_gb_group,
                  map<int, map<string, int> >        homocluster_id_2lg,
                  map<string, int>                   lg2_homocluster_id,
                  map<string, BIPAT_INT>             contigPollenGTSeq_int_pure_hap,
                  map<string, unsigned long>         contigsize,
                  string                             tmpfolder_s4);  
int count_effective_gt_signal(vector<int>            gts); 
bool update_hom_into(string                          first_lgid, 
                  string                             second_lgid, 
                  map<string, double>*               hom_pairs_simple,
                  map<string, double>*               hom_pairs,
                  map<string, vector<string> >*      edge_all);    
bool update_grouping(string                          ctg_id, 
                  string                             win_sta_end,
                  string                             grouping, 
                  lgGROUP                            groupinfo,
                  map<string, map<string, map<string, lgGROUP> > >* raw_grouping);    
bool refine_grouping(map<string, map<string, map<string, lgGROUP> > >* raw_grouping,
                  map<string, NODE2>*                marker_depth_all,
                  map<int, map<string, int> >*       gb_group_2ctg,
                  map<string, int>*                  ctg2_gb_group,
                  map<int, int>*                     gb_group_id,
	          map<int, int>*                     gb_group_id_rev,
	          map<int, map<string, int> >*       homocluster_id_2lg,
     	          map<string, int>*                  lg2_homocluster_id,
     	          map<string, BIPAT_INT>*            contigPollenGTSeq_int_pure_hap,
                  string                             tmpfolder_s4);
bool recalculate_lg_contig_cor(map<string, NODE2>*   marker_depth_all,
                  string                             tmpfolder_s4);                  
//
int main(int argc, char* argv[])
{
    if(argc < 8)
    {
        // g++ gamete_binning.cpp split_string.cpp -O3 -o gamete_binning
        cout << "\nFunction: analyze sequencing depth of tig markers, cluster markers and bin long reads. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl
             << "\nUsage: gamete_binning"
                 << " list_bed_files.txt"
                 << " cor_cutoff"
                 << " contig_size.txt"
                 << " cluster_num_exp"
                 << " hom_cor_cutoff"
                 << " min_hapctg_size"
                 << " output_prefix" << endl 
             << endl
             << "      *list_bed_files.txt -- list of \"lib/barcode_win_marker_read_count_MQ5_updated.bed\". "
             << endl
             << "      *cor_cutoff         -- minimum correlation to build the contact graph of contigs [0.5]. "  
             << endl
             << "      *contig_size.txt    -- list of contig sizes (with two tab-separated columns). " 
             << endl
             << "      *cluster_num_exp    -- expected number of linkage groups. " 
             << endl         
             << "      *hom_cor_cutoff     -- maximum correlation to find \"homologous\" linkage groups [-0.25]."                 
             << endl
             << "      *min_hapctg_size    -- minimum contig size to select haplotig markers for linkage clusterings [50000]."                              
             << "      *output_prefix      -- flag of output file." << endl 
             << endl;
        cout << "   Note-S3: to filter    dot: dot_filter_v2 xxx.dot 0.75"               << endl;
        cout << "   Note-S3: to visualize dot: circo -Tpdf xxx.dot > xxx_circo.pdf, or " << endl;
        cout << "                                fdp -Tpdf xxx.dot > xxx_fdp.pdf . "     << endl;
        return 1;
    }    
    double startT= clock();
    if(argc > 8)
    {
        cout << "   Info: cor-recalc: lg_interctg_cor is asked; running time will increase by half an hour!" << endl;
        lg_interctg_cor = true;
    }
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
    minCorscore = atof(argv[2]);
    cout << "   Info: contig contact graph will be built with minimum correlation of " << minCorscore << endl;
    // read contig size info
    string ctgsizefile = (string)argv[3];  
    cout << "   Info: read contig size info... " << endl;
    map<string, unsigned long> contigsize;
    if(!collect_contig_size(ctgsizefile, &contigsize))
    {
        cout << "   Error: failed in collecting contig size info." << endl;
        return 1;
    }    
    //
    N_cluster = atoi(argv[4]);  
    cout << "   Info: expected number of linkage groups given as " << N_cluster << "." << endl; 
    // get cutoff to find "homologous haplotype-specific linkage groups. "
    hom_cor_cutoff = atof(argv[5]);
    //
    min_hapctg_size = atoi(argv[6]);
    cout << "   Info: ctg-size cutoff for selecting haplotig markers: " << min_hapctg_size << " bp."    << endl;
    cout << "         note that smaller contigs will be inserted to initialized linkage groups later. " << endl;
    //
    string outprefix = (string)argv[7];
    string tmpfolder_s = "s1_"+outprefix+"_tmp_tig_marker_depth_genotype";
    if(!create_folder(tmpfolder_s))
    {
        return 1;
    }
    cout << endl;
    cout << "   Info: starting depth analysis..." << endl;
    // get bed file
    string list_bed_files = (string)argv[1];   
    vector<string> bedfiles;
    if(!(get_individual_bed_file(list_bed_files, &bedfiles)))
    {
        cout << "   Error: cannot get bed files from " << list_bed_files << endl;
        return 1;
    }
    // prepare pure haplotig markers (at least 80% are hap windows: TODO 20200714 - set up option)
    string bed_file = *( bedfiles.begin() );
    map<string, double> hap_win_ratio;
    cout << endl;
    cout << "   Info: finding \"pure\"-haplotig markers..." << endl;
    if(!get_pure_haplotig_markers(bed_file, contigsize, &hap_win_ratio))
    {
        cout << "   Error: failed to calculate ratio of haplo-win markers along contigs." << endl;
        return 1;
    }
    cout << "   Info: " << hap_win_ratio.size() << " (haplotig & non-haplotig) markers collected. " << endl;
    cout << "   Info: finding \"pure\"-haplotig markers done. " << endl;    
    //
    cout << endl;
    cout << "   Info: analyzing depth_genotype (window_marker+pollen => contig_end_marker)..." << endl;
    vector<string>::iterator bitr;
    vector<string>::iterator bitr_end;
    bitr     = bedfiles.begin();
    bitr_end = bedfiles.end();
    int    i = 0;
    // marker_depth_all collects all depth-genotypes at all (hap-/dip-/trip-/tetraplotig) window markers.
    map<string, NODE2>             marker_depth_all; // <ctg_id\tsta\tend, depth_genotype_across_all_pollens>    
    map<int, map<string, NODE> >   pollen_depth_all; // <pollen_id, map<ctg_id\tsta\tend, {depth, marker-type}> >
    map<int, map<string, string> > pollen_ctg_gt_all;// <pollen_id, map<ctg_id, genotype_sequence_along_ctg> >    
    map<string, unsigned long>     ctg_size_all;     // <ctg_id, ctg_length> ~=~ contigsize <= redundant
    vector<string>                 marker_order_keep;// input order of all win_marker_id = ctg_id\tsta\tend  
    map<string, string>            marker_order_keep_hap_only;    
    bool                           first_file = true;
    vector<string>                 s1_tmp_files;     // collect files generated at step 1-3.
    while(bitr != bitr_end)
    {
        string this_bed = *bitr;
        // plollen int id
        i ++; 
        cout << "\n   Info: reading "         << i        << "/" << bedfiles.size() 
             << " depth data from bed file: " << this_bed << endl;
        // step 1. read this pollen depth
        map<string, NODE> this_pollen_depth;  // <ctg_id\tsta\tend, {depth, marker-type}>
        vector<string>    marker_order;            
        if(!read_current_bed(this_bed, 
                             &this_pollen_depth,
                             &ctg_size_all,
                             &marker_order)
          )
        {
            cout << "   Warning: reading failed at " << this_bed << endl;
            continue; 
        }
        assert(this_pollen_depth.find("#total_reads") != this_pollen_depth.end());
        cout << "   Info:   total read alignment observed = " 
             << this_pollen_depth["#total_reads"].depth << endl;
        // step 2. analyze current bed
        if(!normalize_read_count(&this_pollen_depth))
        {
            cout << "   Error: cannot nomalize read counts. " << endl;
            return 1;
        }
        //
        double this_hap_mean = find_haploid_mean_depth(this_pollen_depth);
        cout << "   Info: found normalized read count at haplotig marker: " << this_hap_mean 
             << "." << endl;
        if(this_hap_mean == 0)
        {
            // avoid 0
            this_hap_mean = 0.00001;
        }
        // step 3. output tmp info and update genotype-like value
        //
        vector<string> this_bed_info = split_string(this_bed, '/');
        std::stringstream ss;
        ss.str("");
        ss << tmpfolder_s << "/" << i << "th" << "_pollen_" << this_bed_info[this_bed_info.size()-1];
        ofstream ofp;
        ofp.open(ss.str().c_str(), ios::out);
        if(!ofp.good())
        {
            cout << "   Error: cannot open output file. " << endl;
            return 1;
        }
        s1_tmp_files.push_back(ss.str());
        //
        map<string, string> this_pollen_ctg_gt; // <ctg_id, gt_string>
        string last_ctg  = "";
        string last_type = "";
        //        
        vector<string>::iterator mitr;
        vector<string>::iterator mitr_end;
        mitr     = marker_order.begin();
        mitr_end = marker_order.end();    
        while(mitr != mitr_end)
        {
            string         this_key     = *mitr; // ctg\tstart\tend
            vector<string> this_keyinfo = split_string(this_key, '\t');
            assert(this_pollen_depth.find(this_key) != this_pollen_depth.end());            
            NODE tmpnode = this_pollen_depth[this_key];   
            if(this_key.compare("#total_reads") != 0)
            {
                int this_gt = round(tmpnode.depth_nor/this_hap_mean); // GT=0,1,2,3,...
                if(tmpnode.type.compare("hap") == 0 && this_gt > 1)
                {
                    this_gt = 1;
                }else
                if(tmpnode.type.compare("rep") != 0 && this_gt > 2)
                {
                    this_gt = 2;
                }else ;
                this_pollen_depth[this_key].genotype = this_gt; // update genotype value
                //
                if(tmpnode.type.compare("hap") == 0)
                {
                    if(this_pollen_ctg_gt.find(this_keyinfo[0]) != this_pollen_ctg_gt.end())
                    {
                        this_pollen_ctg_gt[ this_keyinfo[0] ] += to_string(this_gt);
                    }else
                    {
                        // collect                    
                        this_pollen_ctg_gt.insert(std::pair<string, string>(this_keyinfo[0], to_string(this_gt)));
                        last_ctg = this_keyinfo[0];
                    }
                    last_type = "hap";
                    if(i==1) // bed files have same marker info; recording only from pollen i=1 is enough.
                    {
                        map<string, string>::iterator ckitr = marker_order_keep_hap_only.find(this_keyinfo[0]);
                        if(ckitr == marker_order_keep_hap_only.end())
                        {
                            marker_order_keep_hap_only.insert(std::pair<string, string>(this_keyinfo[0], "0"));
                        }
                    }
                }
                //
                /* header: 
                           1.ctg_id                     2.win_sta           3.win_end 
                           4.genotypevalue={0,1,2,x}    5.normalized_depth  6.real_depth  
                           7.type_of_marker={hap,dip,trip,tetrap,rep}       8.ctg_size  
                */ 
                ofp << this_key          << "\t"
                    << this_gt           << "\t"
                    << tmpnode.depth_nor << "\t"
                    << tmpnode.depth     << "\t"
                    << tmpnode.type      << "\t"
                    << ctg_size_all[this_keyinfo[0] ] << endl;
                // collect markers at this pollen 
                if(first_file)
                {
                    // initialize gt string at first pollen gt at this marker
                    NODE2 tmpnode2;
                    (tmpnode2.depth_gt).push_back(this_gt);
                    tmpnode2.type     = tmpnode.type;
                    marker_depth_all.insert(std::pair<string, NODE2>(this_key, tmpnode2));
                }else
                {   
                    // update gt string with additional pollen gt at this marker
                    assert(marker_depth_all.find(this_key) != marker_depth_all.end());
                    (marker_depth_all[this_key].depth_gt).push_back(this_gt);
                }
            }
            mitr ++;  
        } 
        // check
        if(0)
        {
            map<string, string>::iterator hapitr;
            map<string, string>::iterator hapitr_end;
            hapitr     = this_pollen_ctg_gt.begin();
            hapitr_end = this_pollen_ctg_gt.end();
            while(hapitr != hapitr_end)
            {
                cout << "   check-gt-str " << (*hapitr).first << " gt = " << (*hapitr).second << endl;
                hapitr ++;
            }
        }
        // close output file 
        ofp.close();  
        // step 4. collect current bed
        pollen_depth_all.insert(std::pair<int, map<string, NODE> >(i, this_pollen_depth));   
        this_pollen_depth.clear();     
        pollen_ctg_gt_all.insert(std::pair<int, map<string, string> >(i, this_pollen_ctg_gt));
        this_pollen_ctg_gt.clear();
        if(marker_order_keep.size()==0)
        {
            marker_order_keep = marker_order;
        }
        marker_order.clear();// if not clear, it's accumualating causing redundant info.
        //
        if(first_file) first_file = false;
        // next bed file
        bitr ++;
    }
    //
    // output genotype matrix for all pollen at all markers - dot format
    std::stringstream gtfile;
    gtfile.str("");
    gtfile << tmpfolder_s << "/final_depth_genotype_matrix.txt";
    ofstream ofp2;
    ofp2.open(gtfile.str().c_str(), ios::out);
    if(!ofp2.good())
    {
        cout << "   Error: cannot open output file: " << gtfile.str() << endl;
        return 1;
    }  
    //      
    // rows...: contig window marekers  
    // columns: pollen depth genotypes
    vector<string>::iterator mitr;
    vector<string>::iterator mitr_end;
    mitr     = marker_order_keep.begin();
    mitr_end = marker_order_keep.end();
    while(mitr != mitr_end)
    {
        string this_key  = *mitr;
        assert(marker_depth_all.find(this_key) != marker_depth_all.end());
        // contig win_sta win_end        
        NODE2 tmpnode2                  = marker_depth_all[this_key];  
        std::replace( this_key.begin(), this_key.end(), '\t', ' ');
        ofp2 << this_key; 
        // hap/dip/trip/tetrap                 
        string this_type                = tmpnode2.type;
        ofp2 << " " << this_type;
        // depth_genotype at all pollen
        vector<int> this_gts  = tmpnode2.depth_gt;        
        vector<int>::iterator gtitr;
        vector<int>::iterator gtitr_end;
        gtitr     = this_gts.begin();
        gtitr_end = this_gts.end();
        while(gtitr != gtitr_end)
        {
            ofp2 << " " << *gtitr;
            gtitr ++;
        }
        ofp2 << endl;
        // next marker
        mitr ++;
    }
    //
    ofp2.close();    
    //
    cout << "\n   Info: " << pollen_depth_all.size() << " pollen bed files analyzed. " << endl;
    // step 5. create haplotig end markers indicating potential crossovers
    string tmpfolder_s2 = "s2_"+outprefix+"_tmp_haplotig_end_markers";
    if(!create_folder(tmpfolder_s2))
    {
        return 1;
    }
    std::stringstream gtinfo;
    gtinfo.str("");
    gtinfo << tmpfolder_s2
           << "/s2_genotype_contig_seq.txt\0"; //
    ofstream ctgPMofp;
    ctgPMofp.open((gtinfo.str()).c_str(), ios::out);
    if(!ctgPMofp.good())
    {
        cout   << "   Error: cannot open file " << gtinfo.str() << endl;
        return 1;
    }
    // to record <ctg_id, {leftPat, rightPat}>
    map<string, BIPAT> contigPollenGTSeq;
    map<string, BIPAT_INT> contigPollenGTSeq_int;          // all hap-win markers
    map<string, BIPAT_INT> contigPollenGTSeq_int_pure_hap; // where all contig-wins are hap
    //      
    map<string, string>::iterator moitr;
    map<string, string>::iterator moitr_end;
    moitr     = marker_order_keep_hap_only.begin();
    moitr_end = marker_order_keep_hap_only.end();
    bool output_ss = false;
    while(moitr != moitr_end)
    {
        string this_key  = (*moitr).first;  // contig id             
        cout     << "#genotype at "        << this_key << " (hap-marker, raw size: " 
                 << ctg_size_all[this_key] << " bp)"   << endl << endl;
        ofstream ctgGTofp;                 
        if(output_ss)
        {
            std::stringstream ctg_gt_file;
            ctg_gt_file.str("");
            ctg_gt_file << tmpfolder_s2 << "/" 
                        << this_key     << "_genotype_seq_details_all_pollen.txt\0"; //
                ctgGTofp.open((ctg_gt_file.str()).c_str(), ios::out);
            if(!ctgGTofp.good())
            {
                cout   << "   Error: cannot open file " << ctg_gt_file.str() << endl;
                return 1;
            }         
            ctgGTofp << "#genotype at "        << this_key << " (hap-marker, raw size: " 
                     << ctg_size_all[this_key] << " bp)"   << endl << endl;             
        }
        //
        string leftGateSeq           = "";
        string rightGateSeq          = ""; 
        string recombined            = "";    
        //
        vector<int> leftGateSeq_int;
        vector<int> rightGateSeq_int;  
        //
        map<int, map<string, string> >::iterator pitr;
        map<int, map<string, string> >::iterator pitr_end;
        pitr     = pollen_ctg_gt_all.begin();
        pitr_end = pollen_ctg_gt_all.end();
        while(pitr != pitr_end)
        {
            if((*pitr).second.find(this_key) != (*pitr).second.end())
            {
                unsigned long leftp    = 0;
                double        score    = 0;
                std::stringstream outss;
                outss.str("");                                     
                string        depth_gt = (*pitr).second[this_key];
                int           pollenid = (*pitr).first;
                string pmstring        = get_break_pos(depth_gt, 
                                                       "U", 
                                                       pollenid, 
                                                       &leftp, 
                                                       &score, 
                                                       &outss, 
                                                       this_key, 
                                                       output_ss);
                //
                if(output_ss)
                {
                    ctgGTofp << outss.str();
                }        
                //
                assert(pmstring.size()==2);
                leftGateSeq    += pmstring.substr(0, 1);
                rightGateSeq   += pmstring.substr(1, 1);         
                if(pmstring.substr(0, 1).compare( pmstring.substr(1, 1) ) != 0)
                {
                    recombined     += "R";
                }
                else
                {
                    recombined     += "|";
                } 
                //
                leftGateSeq_int.push_back(atoi(pmstring.substr(0, 1).c_str()));
                rightGateSeq_int.push_back(atoi(pmstring.substr(1, 1).c_str()));
                // cout << (*pitr).second[this_key] << " p_" << (*pitr).first << "_x" << endl;
            }
            pitr ++;
        }
        cout     << endl;
        if(output_ss)
        {
            ctgGTofp << endl;
            ctgGTofp.close();
        }
        //
        ctgPMofp << leftGateSeq  << "\t" << this_key << "\tleft-GT"  << endl;  
        ctgPMofp << recombined   << "\t" << this_key << "\tCO-info:" << ctg_size_all[this_key] << "bp" << endl;          
        ctgPMofp << rightGateSeq << "\t" << this_key << "\tright-GT" << endl;     
        // collect
        BIPAT tmpbipat;
        tmpbipat.leftPat  = leftGateSeq;
        tmpbipat.rightPat = rightGateSeq;
        contigPollenGTSeq.insert(std::pair<string, BIPAT>(this_key, tmpbipat));
        // collect int
        BIPAT_INT tmpbipat_int;
        tmpbipat_int.leftPat_int  = leftGateSeq_int;
        tmpbipat_int.rightPat_int = rightGateSeq_int;
        contigPollenGTSeq_int.insert(std::pair<string, BIPAT_INT>(this_key, tmpbipat_int));   
        // here we have one level of haplotig marker selection: at least min_hapctg_size bp haplotig markers used to build a backbone!
        assert( hap_win_ratio.find(this_key) != hap_win_ratio.end() );
        if(hap_win_ratio[this_key] >= haplotig_win_ratio && ctg_size_all[this_key]>=min_hapctg_size)
        {
            contigPollenGTSeq_int_pure_hap.insert(std::pair<string, BIPAT_INT>(this_key, tmpbipat_int));               
        }      
        // next marker
        moitr ++;
    }
    // close file 
    ctgPMofp.close();
    // step 6: LG analysis: build contig contact graph with correlation of genotype sequences
    string tmpfolder_s3 = "s3_"+outprefix+"_tmp_haplotig_end_markers_LG_analysis";
    if(!create_folder(tmpfolder_s3))
    {
        return 1;
    }  
    // contigPollenGTSeq_int including many with only smaller ratio of haplotigs causing mis-joining of different groups
    // thus here we use only "pure"-haplotigs, where hap-wins take up >= 80% contig size.
    map<int, int>               gb_group_id; // <old_group_id:to_dot, new_group_id>
    map<int, map<string, int> > gb_group;    // <old_group_id, <ctg_id, ref_grp_id > >
    if(!greedy_cluster_haplotigs(contigPollenGTSeq_int_pure_hap, 
                                 contigsize, 
                                 tmpfolder_s3,
                                 &gb_group_id,
                                 &gb_group))
    {
        cout << "    Error: greedy clustering of haplotigs was not finished as expected. " << endl;
        return 1;
    }
    // create a reverse searching map: new_group_id as key, old_group_id as value.
    map<int, int> gb_group_id_rev; // <new_group_id, old_group_id>
    map<int, int>::iterator iditr;
    map<int, int>::iterator iditr_end;
    iditr     = gb_group_id.begin();
    iditr_end = gb_group_id.end();
    while(iditr != iditr_end)
    {
        gb_group_id_rev.insert(std::pair<int, int>((*iditr).second, (*iditr).first));
        iditr ++;
    }
    // read linkage group information
    cout << "   Info: re-collect output of function: merge_clusters ..." << endl;
    string gb_group_file = tmpfolder_s3 + "/s3_genotype_haplotig_GT_similarity_matrix_subclusters_final.dot\0"; // ;
    map<string, int>            ctg2_gb_group; // <ctg_id, old_group_id>
    map<int, map<string, int> > gb_group_2ctg; // <old_group_id, <ctg_id, new_group_id> >
    if(!read_gamete_binning_lg(gb_group_file, 
                               &ctg2_gb_group,
                               &gb_group_2ctg) )
    {
        return 1;
    }
    cout << "   Info: re-collect output of function: merge_clusters done. " << endl;    
    // read homologous clusters of linkage groups
    cout << "   Info: re-collect output of function: find_hom_group.." << endl;
    string ohomfilename = tmpfolder_s3 + "/s3_res_homologous_linkgage_groups.dot";    
    // homocluster_id=1:12 links to new_group_id: new_group_id="1:48"
    map<int, map<string, int> > homocluster_id_2lg; // <homocluster_id=1:12, <new_group_id="1:48", 1> > 
    map<string, int>            lg2_homocluster_id; // <new_group_id="1:48", homocluster_id=1:12>......
    if(!read_gb_homo_lg(ohomfilename,
                        &homocluster_id_2lg,
                        &lg2_homocluster_id))
    {
        return 1;
    }
    cout << "   Info: re-collect output of function: find_hom_group done." << endl;    
    // step 7: LG analysis: integrate all dip/trip/tetraplotig window markers to linkage groups
    //
    string tmpfolder_s4 = "s4_"+outprefix+"_tmp_integrating_all_ctg_markers_to_LGs";
    if(!create_folder(tmpfolder_s4))
    {
        return 1;
    }   
    if(!integrate_dtt_markers_to_lg(&marker_depth_all,
                                    marker_order_keep,
                                    hap_win_ratio,
                                    gb_group_id,
                                    gb_group_id_rev,
                                    gb_group_2ctg,
                                    ctg2_gb_group,
                                    homocluster_id_2lg,
                                    lg2_homocluster_id,
                                    contigPollenGTSeq_int_pure_hap,
                                    contigsize,
                                    tmpfolder_s4))
    {
        cout << "   Error: cannot integrate dip/trip/tetraplotig window markers. " << endl;
        return 1;
    } 
    // default: lg_interctg_cor==false! This is for checking purpose!
    if(lg_interctg_cor && !recalculate_lg_contig_cor(&marker_depth_all, tmpfolder_s4) )
    {
        cout << "   Error: failed in re-calculating intra-LG inter-CTG correlation. " << endl;
        return 1;
    }
    //
    double finishT= clock();
    cout << "   Time consumed: " << (finishT-startT)/1000000 << " seconds" << endl << endl;
    //
    return 0;
}
/*
ctg			start	end	#reads	ctg		 start	end	#raw.avg.reads	#windows ctg_len type
utg000001l_pilon	150001	200000	0	utg000001l_pilon 150001	200000	115.379		227	 2264094 hap
*/
bool recalculate_lg_contig_cor(map<string, NODE2>*   marker_depth_all,
                               string                tmpfolder_s4)
{
    /*
        Function: re-calculate intra-LG inter-CTG correlation: 
            if a marker shows negative correlation to ~50% of the LG ctgs, it needs to be re-assigned to another LG.
        
        marker_depth_all.....: <ctg_id\twin_sta\twin_end, {<depth_gt_int>, hap/dip/trip/tetra, ..} >....................
    */
    // step 1. read refined window markers
    string ifilenamefinal = tmpfolder_s4 + "/s4p6_refine_grouping_final_window_markers.txt";
    ifstream ifp_win_final;
    ifp_win_final.open(ifilenamefinal.c_str(), ios::in);
    if(!ifp_win_final.good())
    {
        cout << "   Error: cannot open file " << ifilenamefinal << endl;
        return false;
    }     
    map<string, map<string, string> > refined_grouping; // <LGID-string, <ctg_id\tsta\tend, group_line > >
    while(ifp_win_final.good())
    {
        string line("");
        getline(ifp_win_final, line);
        if(line.size()==0 || line[0]=='#') continue;
        /* example lines:
		utg000001l_pilon    1   50000   hap 1   11  utg000001l_pilon    1.1 1   0   linked  case_1.1
		utg000001l_pilon    50001   100000  hap 1   11  utg000001l_pilon    1.1 1   0   linked  case_1.1
		utg000001l_pilon    100001  150000  hap 1   11  utg000001l_pilon    1.1 1   0   linked  case_1.1        
		...
        */
        vector<string> lineinfo = split_string(line, '\t');
        string this_LG  = lineinfo[4];
        string this_key = lineinfo[0] + "\t" + lineinfo[1] + "\t" + lineinfo[2]; // ctg_id\tsta\end
        string this_val = line;
        if(refined_grouping.find(this_LG) == refined_grouping.end() )
        {
            map<string, string> tmp_win_marker;
            tmp_win_marker.insert(std::pair<string, string>(this_key, this_val));
            refined_grouping.insert(std::pair<string, map<string, string> >(this_LG, tmp_win_marker));
        }else
        {
            assert(refined_grouping[this_LG].find(this_key) == refined_grouping[this_LG].end() );
            refined_grouping[this_LG].insert(std::pair<string, string>(this_key, this_val));
        }
    }
    ifp_win_final.close();
    // step 2. check correlation within each LG
    map<string, map<string, string> >::iterator lgitr;
    map<string, map<string, string> >::iterator lgitr_end;
    lgitr     = refined_grouping.begin();
    lgitr_end = refined_grouping.end();
    while(lgitr != lgitr_end)
    {
        cout << "    recalc-check: linkage group id " << (*lgitr).first 
             << " with "                              << ((*lgitr).second).size() 
             << " window markers." << endl;
        lgitr ++;
    }    
    lgitr     = refined_grouping.begin();
    lgitr_end = refined_grouping.end();
    while(lgitr != lgitr_end)
    {
        string this_LG;
        map<string, string> this_LG_win_marker; // <ctg_id\tsta\tend, markerline>
        this_LG            = (*lgitr).first;
        this_LG_win_marker = (*lgitr).second;
        if(this_LG.compare("-1") == 0)
        {
            lgitr ++;
            continue;
        }        
        map<string, string>::iterator witr1;
        map<string, string>::iterator witr1_end;
        witr1     = this_LG_win_marker.begin();
        witr1_end = this_LG_win_marker.end();
        while(witr1 != witr1_end)
        {
            string wmarker1 = (*witr1).first; // ctg_id\tsta\tend
            vector<string> wm1info = split_string(wmarker1, '\t');
            string ctg1_id = wm1info[0];
            cout << "   recalc-check: LG " << this_LG << ":" << (*lgitr).first << " at line " << (*witr1).second << endl;            
            //
            assert((*marker_depth_all).find(wmarker1) != (*marker_depth_all).end());
            NODE2 tmpnode1  = (*marker_depth_all)[wmarker1];  
            // get type........: hap/dip/trip/tetrap              
            string this_type1     = tmpnode1.type;
            // get genotype seq: depth_genotype at all pollen
            vector<int> this_gts1 = tmpnode1.depth_gt;            
            //
            map<string, string>::iterator witr2;
            map<string, string>::iterator witr2_end;
            witr2     = this_LG_win_marker.begin();
            witr2_end = this_LG_win_marker.end();
            while(witr2 != witr2_end)
            {
                string wmarker2 = (*witr2).first; // ctg_id\tsta\tend                
                if( wmarker1.compare(wmarker2) == 0 )
                {
                    // next window marker 2
                    witr2 ++;  
                    continue;                  
                }
                vector<string> wm2info = split_string(wmarker2, '\t');
                string ctg2_id = wm2info[0];
                if(ctg1_id.compare(ctg2_id) == 0)
                {
                    witr2 ++;
                    continue;
                }
                // 
                assert((*marker_depth_all).find(wmarker2) != (*marker_depth_all).end());
                NODE2 tmpnode2  = (*marker_depth_all)[wmarker2];  
                // get type........: hap/dip/trip/tetrap              
                string this_type2     = tmpnode2.type;
                // get genotype seq: depth_genotype at all pollen
                vector<int> this_gts2 = tmpnode2.depth_gt;     
                // calculate corralation
                double cor = -200.0;    
                bool recalc=true;            
                if(this_type1.compare("hap") == 0 &&
                   this_type2.compare("hap") == 0 )
                {
                    // hap1-hap2
                    cor = calculate_correlation(&this_gts1,  &this_gts2);                     
                }else
                if(this_type1.compare("hap") == 0 &&
                   this_type2.compare("dip") == 0 )
                {
                    // hap1-dip2
                    cor = calculate_correlation_diplotig(&this_gts2,  &this_gts1);                     
                }else                
                if(this_type1.compare("hap") == 0 &&
                   this_type2.compare("trip") == 0 )
                {
                    // hap1-trip2
                    cor = calculate_correlation_triplotig(&this_gts2,  &this_gts1);                      
                    cor = cor * -1.0;
                }else                
                if(this_type1.compare("dip") == 0 &&
                   this_type2.compare("hap") == 0 )
                {
                    // dip1-hap2
                    cor = calculate_correlation_diplotig(&this_gts1,  &this_gts2);                     
                }else                 
                if(this_type1.compare("trip") == 0 &&
                   this_type2.compare("hap") == 0 )
                {
                    // trip1-hap2
                    cor = calculate_correlation_triplotig(&this_gts1,  &this_gts2);   
                    cor = cor * -1.0;                                      
                }else 
                {
                    //                                         hap-tetrap;
                    //                dip-dip;    dip-trip;    dip-tetrap; 
                    //               trip-dip;   trip-trip;   trip-tetrap;
                    // tetrap-hap; tetrap-dip; tetrap-trip; tetrap-tetrap.     
                    recalc = false;    
                }   
                //
                if(recalc)
                cout << "   recalc-check: LG " << this_LG 
                     << " "                    << wmarker1 << " " << this_type1 << " "
                     << " -- "                 << wmarker2 << " " << this_type2 << " "
                     << " cor.val "            << cor 
                     << endl;
                // next window marker 2
                witr2 ++;
            }
            // next window marker 1
            witr1 ++;
        }
        // next LG
        lgitr ++;
    }
    //
    return true;
}                  
// check how many pollens showing positive signal of "1"s: haplotig marker only
int count_effective_gt_signal(vector<int> gts)
{
    int pcount = 0; // number of positive bits    
    vector<int>::iterator cditr;
    vector<int>::iterator cditr_end;
    cditr      = gts.begin();
    cditr_end  = gts.end();                
    while(cditr != cditr_end)
    {   
        if(*cditr > 0)
        { 
            pcount ++;
        }
        cditr ++;
    }  
    return pcount;     
}
//
bool update_grouping(string  ctg_id, 
                     string  win_sta_end,
                     string  grouping, 
                     lgGROUP groupinfo,
                     map<string, map<string, map<string, lgGROUP> > >* raw_grouping)
{
    /*
        ctg_id......: contig id 
        grouping....: linked/unlinked_tetra/unlinked_rep/unlinked_limited_gamete
        groupinfo...: detailes on this window marker, see lgGROUP
        raw_grouping: map of grouping info :: <ctg_id, <"linked/...", <win_sta_end, lgGROUP > > > 	
    */
    if( (*raw_grouping).find(ctg_id) == (*raw_grouping).end() )
    {
        // contig not found; initialize for contig 
        map<string, lgGROUP> tmp_win; // <win_sta_end, groupinfo >
        tmp_win.insert(std::pair<string, lgGROUP>(win_sta_end, groupinfo));
        //
        map<string, map<string, lgGROUP> > tmp_grouping; // <grouping, tmp_win>
        tmp_grouping.insert(std::pair<string, map<string, lgGROUP> >(grouping, tmp_win));    
        //            
        (*raw_grouping).insert(std::pair<string, map<string, map<string, lgGROUP> > >(ctg_id, tmp_grouping));
    }else
    {
        // contig found: work with (*raw_grouping)[ctg_id]
        if( (*raw_grouping)[ctg_id].find(grouping) == (*raw_grouping)[ctg_id].end() )
        {
            // grouping not found; initialize for grouping 
            map<string, lgGROUP> tmp_win; // <win_sta_end, groupinfo >
            tmp_win.insert(std::pair<string, lgGROUP>(win_sta_end, groupinfo));
            //
            (*raw_grouping)[ctg_id].insert(std::pair<string, map<string, lgGROUP> >(grouping, tmp_win));                      
        }else
        {
            // grouping found: work with (*raw_grouping)[ctg_id][grouping]
            if( (*raw_grouping)[ctg_id][grouping].find(win_sta_end) == (*raw_grouping)[ctg_id][grouping].end() )
            {
                // win_sta_end not found; initialize for win_sta_end
                (*raw_grouping)[ctg_id][grouping].insert(std::pair<string, lgGROUP>(win_sta_end, groupinfo));
            }else
            {            
                // win_sta_end found: work with (*raw_grouping)[ctg_id][grouping][win_sta_end]
                (*raw_grouping)[ctg_id][grouping][win_sta_end].type.push_back(  groupinfo.type[0] );
                (*raw_grouping)[ctg_id][grouping][win_sta_end].qlg.push_back(   groupinfo.qlg[0] );     
                (*raw_grouping)[ctg_id][grouping][win_sta_end].qchr.push_back(  groupinfo.qchr[0] );            
                (*raw_grouping)[ctg_id][grouping][win_sta_end].sctg.push_back(  groupinfo.sctg[0] );            
                (*raw_grouping)[ctg_id][grouping][win_sta_end].slg.push_back(   groupinfo.slg[0] );            
                (*raw_grouping)[ctg_id][grouping][win_sta_end].qscor.push_back( groupinfo.qscor[0] );            
                (*raw_grouping)[ctg_id][grouping][win_sta_end].qscas.push_back( groupinfo.qscas[0] );               
            }
        }                
    }
    //    
    return true;
}                     
// function: integrate dip-/trip-/tetraplotig window markers to linkage groups
bool integrate_dtt_markers_to_lg(map<string, NODE2>*         marker_depth_all,
                                 vector<string>              marker_order_keep,
                                 map<string, double>         hap_win_ratio,
                                 map<int, int>               gb_group_id,
                                 map<int, int>               gb_group_id_rev,
                                 map<int, map<string, int> > gb_group_2ctg,
                                 map<string, int>            ctg2_gb_group,
                                 map<int, map<string, int> > homocluster_id_2lg,
                                 map<string, int>            lg2_homocluster_id,
                                 map<string, BIPAT_INT>      contigPollenGTSeq_int_pure_hap,
                                 map<string, unsigned long>  contigsize,
                                 string                      tmpfolder_s4)
{
    /*  
        dtt..................: diplotig, triplotig, tetraplotig.......; these are the markers to integrate..............
        marker_depth_all.....: <ctg_id\twin_sta\twin_end, {<depth_gt_int>, hap/dip/trip/tetrap, ..} >...................
        marker_order_keep....; <ctg_id\tsta\tend>.....................; input order of key = ctg_id\tsta\tend...........
        hap_win_ratio........: <ctg_id, hap-win-ratio>................; how much of a contig is haplotig windows........
        gb_group_id..........: <old_group_id, new_group_id=1:48>......; ................................................      
        gb_group_id_rev......: <new_group_id=1:48, old_group_id>......; ................................................
        gb_group_2ctg........: <old_group_id, <ctg_id, new_group_id> >; list of contigs belong to one LG................    
        ctg2_gb_group........: <ctg_id, old_group_id>.................; using contig to find old LG id..................     
        homocluster_id_2lg...: <homocluster_id=1:12, <new_group_id="1:48", 1> >.; list of homologous LGs................
        lg2_homocluster_id...: <new_group_id="1:48", homocluster_id=1:12>.......; new_group_id to find homo-cluster id.. 
        contigPollenGTSeq_int_pure_hap: <ctg_id, {leftPat_int, rightPat_int}>...; GT. across all pollen at contig ends..
        contigsize...........: <ctg_id, ctg_length>
        
        Note: contigs of ctg2_gb_group is a subset of contigPollenGTSeq_int_pure_hap       
    */
    // contig-wise variable for collecting outputs: <contig, <"linked/unlinked_tetra/...", <sta-end, {type, qlg, qchr, ...} > > >
    map<string, map<string, map<string, lgGROUP> > > raw_grouping;
    unsigned long case_collect = 0;
    // prepare files for collecting outputs: <tmpfolder_s4+"problematic/unlinked/linked_.bed", *ofp >    
    // markers where signal as '1/2's in depth genotype sequence is supported by less than 1% pollen
    string ofilename_limited = tmpfolder_s4 + "/s4p1_all_raw_unlinked_markers_with_limited_pollen_support.bed";
    ofstream ofp_limited = ofstream(ofilename_limited.c_str(), ios::out);
    if(!ofp_limited.good())
    {
        cout << "   Error: cannot set up output file " << ofilename_limited << endl;
        return false;
    }
    ofp_limited << "#header: ctg\tstart\tend\tmarker_type\tpollen_support" << endl;
    // tetraplotig markers that can not be linked to existing linkage groups
    string ofilename_tetrap = tmpfolder_s4 + "/s4p2_all_raw_unlinked_tetraplotig_markers.bed";
    ofstream ofp_tetra = ofstream(ofilename_tetrap.c_str(), ios::out);
    if(!ofp_tetra.good())
    {
        cout << "   Error: cannot set up output file " << ofilename_tetrap << endl;
        return false;
    }
    //  
    string ofilename_linked = tmpfolder_s4 + "/s4p3_all_raw_linked_hap_dip_triplotig_markers.bed";
    ofstream ofp_linked = ofstream(ofilename_linked.c_str(), ios::out);
    if(!ofp_linked.good())
    {
        cout << "   Error: cannot set up output file " << ofilename_linked << endl;
        return false;
    }
    //  
    string ofilename_rep = tmpfolder_s4 + "/s4p4_all_raw_rep_markers.bed";
    ofstream ofp_rep = ofstream(ofilename_rep.c_str(), ios::out);
    if(!ofp_rep.good())
    {
        cout << "   Error: cannot set up output file " << ofilename_rep << endl;
        return false;
    } 
    // filter subject genotype sequences: 
    //        1. remove those with very low signals with limited number of pollen showing "1"s
    //           filtering because those with limited number of '1's will misguiding identification of matching.
    //        2. remove those not in existing LGs.
    //        3. remove those less than 50 kb
    cout << "   Info: " << contigPollenGTSeq_int_pure_hap.size() << " initial haplotig-end-genotypes provided. "<< endl;
    cout << "   Info: " << ctg2_gb_group.size()                  << " initial haplotigs in existing LGs"        << endl;
    cout << "         the latter subset will be used as subject to anchor contigs not in existing LGs. "        << endl;
    map<string, BIPAT_INT>::iterator epitr;
    map<string, BIPAT_INT>::iterator epitr_end;
    epitr     = contigPollenGTSeq_int_pure_hap.begin();
    epitr_end = contigPollenGTSeq_int_pure_hap.end();
    int cerased1 = 0;
    int cerased2 = 0;
    int cerased3 = 0;
    int ckept    = 0;
    while(epitr != epitr_end)
    {
        BIPAT_INT pat;
        pat.leftPat_int   = (*epitr).second.leftPat_int;
        int exp_size      = pat.leftPat_int.size();
        int left_cnt      = count_effective_gt_signal(pat.leftPat_int);
        pat.rightPat_int  = (*epitr).second.rightPat_int; 
        int right_cnt     = count_effective_gt_signal(pat.rightPat_int);   
        if(ctg2_gb_group.find( (*epitr).first ) == ctg2_gb_group.end() )
        {
            cerased1 ++;                 
            contigPollenGTSeq_int_pure_hap.erase(epitr ++);   
        }else
        if( right_cnt*1.0/exp_size<0.135 || left_cnt*1.0/exp_size<0.135) // 100/720
        {
            cerased2 ++;        
            contigPollenGTSeq_int_pure_hap.erase(epitr ++);                       
        }else
        if(contigsize[ (*epitr).first ] < 40000) // marker too short might be noises
        {
            cerased3 ++;
            contigPollenGTSeq_int_pure_hap.erase(epitr ++);            
        }
        else
        {
            ckept ++;
            epitr ++; 
        }
        //
    }
    cout << "   Info: " << cerased1 << " contigs removed from subject list as not in existing LGs. " << endl
         << "         " << cerased2 << " contigs removed from subject list due to limited \"1\"s. "  << endl
         << "         " << cerased3 << " contigs removed from subject list due to <40 kb "           << endl         
         << "         " << ckept    << " contigs kept in the subject list. "                         << endl;
    cout << "   Info: " << contigPollenGTSeq_int_pure_hap.size() 
         << " good haplotig-end-genotypes remained for matching dip/trip/tetraplotigs. " << endl;  
    //
    vector<string>::iterator mitr;
    vector<string>::iterator mitr_end;
    mitr     = marker_order_keep.begin();
    mitr_end = marker_order_keep.end();
    while(mitr != mitr_end)
    {
        // get key.........: ctg_id\tsta\tend
        string this_key = *mitr;
        vector<string> keyinfo = split_string(this_key, '\t');
        string ctg_id = keyinfo[0];
        // get node........: {<depth_gt_int>, hap/dip/trip/tetrap, ..}
        assert((*marker_depth_all).find(this_key) != (*marker_depth_all).end());
        NODE2 tmpnode2  = (*marker_depth_all)[this_key];  
        // get type........: hap/dip/trip/tetrap              
        string this_type     = tmpnode2.type;
        // get genotype seq: depth_genotype at all pollen
        vector<int> this_gts = tmpnode2.depth_gt;
        // check signal strength at this marker. At some markers, interval is too small and depth gt are all '0'.
        // e.g., utg000005l_pilon	1	10000	111	1	10518574	hap
        vector<int>::iterator cditr;
        vector<int>::iterator cditr_end;
        cditr      = this_gts.begin();
        cditr_end  = this_gts.end();                
        int pcount = 0; // number of positive bits
        while(cditr != cditr_end)
        {   
            if(this_type.compare("trip") == 0)
            {
                if(*cditr>0 && 2-*cditr > 0)
                {
                    pcount ++;
                }
            }else
            if(*cditr > 0)
            { 
                pcount ++;
            }
            cditr ++;
        }        
        if( pcount*1.0/this_gts.size() < 0.01 )
        {
            // this needs to go to a no_signal_cluster
            cout << "   output0: problematic win-marker="    << this_key 
                 << ", type="                                << this_type  
                 << ", with limited "                        << pcount
                 << " pollens with depth signal of \"1\" "
                 << endl; 
            ofp_limited << this_key  << "\t"
                        << this_type << "\t"
                        << pcount    << "\t"
                        << "case_0"  << endl;            
            // collect grouping info 
            // string ctg_id      = keyinfo[0];
            string win_sta_end = keyinfo[1]+"\t"+keyinfo[2];
            string grouping    = "unlinked_limited_gamete";            
            lgGROUP groupinfo;
            groupinfo.type.push_back( this_type ); // type of current window marker: hap/dip/trip/tetrap...
            groupinfo.qlg.push_back(  -1 );        // LG of query contig  = 1:48
            groupinfo.qchr.push_back( -1 );        // Chr of query contig = 1:12, i.e., homocluster_id
            groupinfo.sctg.push_back("NA.ctg" );   // id of linked subject contig
            groupinfo.slg.push_back(  -1 );        // LG of linked subject contig = 1:48
            groupinfo.qscor.push_back( -100.0 );   // correlation value linked query and subject contigs
            groupinfo.qscas.push_back( "case_0" ); // case info: how this linking is found                                            
            update_grouping(ctg_id, 
                            win_sta_end,
                            grouping, 
                            groupinfo,
                            &raw_grouping);   
            case_collect ++;                                                 
            // next marker please!
            mitr ++;
            continue;
            //           
        }             
        // two major cases: 1. current marker partially in an existing LG; 2. not in existing LGs at all.
        map<string, int>::iterator clgitr = ctg2_gb_group.find(ctg_id);
        if(clgitr != ctg2_gb_group.end())
        {
            // found ctg in existing an linkage group (old_group_id)
            // get old group id
            int this_old_group_id = (*clgitr).second;
            std::stringstream this_old_group_id_ss;
            this_old_group_id_ss.str("");
            this_old_group_id_ss << this_old_group_id;
            // get new group id
            int this_new_group_id = gb_group_id[ this_old_group_id ] ;
            std::stringstream this_new_group_id_ss;
            this_new_group_id_ss.str("");
            this_new_group_id_ss << this_new_group_id;            
            // classify the marker to one pre-define linkage group
            cout << "   output1.1: win-marker=" << this_key 
                 << ", type="                   << this_type  
                 << " will go to LG="           << this_new_group_id
                 << " (old_group_id:"           << this_old_group_id << ")"                  
                 << endl;
            ofp_linked << this_key             << "\t"
                       << this_type            << "\t"
                       << this_new_group_id    << "\t"
                       << lg2_homocluster_id[ this_new_group_id_ss.str() ] << "\t"
                     //<< this_old_group_id    << "\t"
                       << "Exist_in_LG"        << "\t"
                       << "case_1.1"           << endl;
            // collect grouping info 
            // string ctg_id      = keyinfo[0];
            string win_sta_end = keyinfo[1]+"\t"+keyinfo[2];
            string grouping    = "linked";            
            lgGROUP groupinfo;
            groupinfo.type.push_back( this_type );         // type of current window marker: hap/dip/trip/tetrap...
            groupinfo.qlg.push_back(  this_new_group_id ); // LG of query contig  = 1:48
            groupinfo.qchr.push_back( lg2_homocluster_id[ this_new_group_id_ss.str() ] ); // Chr of query contig = 1:12, i.e., homocluster_id
            groupinfo.sctg.push_back( ctg_id );            // id of linked subject contig
            groupinfo.slg.push_back(  this_new_group_id ); // LG of linked subject contig = 1:48
            groupinfo.qscor.push_back( 1.1 );              // correlation value linked query and subject contigs
            groupinfo.qscas.push_back( "case_1.1" );       // case info: how this linking is found                                            
            update_grouping(ctg_id, 
                            win_sta_end,
                            grouping, 
                            groupinfo,
                            &raw_grouping); 
            case_collect ++;                                                                              
            /* TODO: send one copy to LG */            
            if(this_type.compare("dip")==0)
            {
                /* send another copy to one of the homologous LGs of this_new_group_id */                            
                // get homocluster_id
                assert( lg2_homocluster_id.find( this_new_group_id_ss.str() )  != lg2_homocluster_id.end() );
                int homocluster_id = lg2_homocluster_id[ this_new_group_id_ss.str() ];                
                // get list of homologous linkage groups
                assert( homocluster_id_2lg.find( homocluster_id ) != homocluster_id_2lg.end() );
                map<string, int>  homo_lgs = homocluster_id_2lg[ homocluster_id ]; // <"1:48", 1>
                // find best correlation between target widnow marker and contigs in these homolougous linkage groups
                double best_ctg_cor_this_homocluster = 0; 
                string best_ctg_id_this_homocluster  = "";
                int    best_linkage_group_id         = 0;
                map<string, int>::iterator hcitr;
                map<string, int>::iterator hcitr_end;
                hcitr     = homo_lgs.begin();
                hcitr_end = homo_lgs.end();
                while(hcitr != hcitr_end)
                {
                    //
                    if( ( (*hcitr).first).compare( this_new_group_id_ss.str() ) ==0 )
                    {
                        hcitr ++;
                        continue;
                    }
                    // get new group id
                    int tmpLG_id_new = atoi( (*hcitr).first.c_str() ); // 1:48
                    assert( gb_group_id_rev.find( tmpLG_id_new ) != gb_group_id_rev.end() );
                    // get old group id
                    int tmpLG_id_old = gb_group_id_rev[ tmpLG_id_new ];
                    // get list of contigs in the linkage group
                    assert( gb_group_2ctg.find( tmpLG_id_old ) != gb_group_2ctg.end() );                   
                    map<string, int> this_LG_ctgs = gb_group_2ctg[ tmpLG_id_old ];
                    // find best correlation between target window marker and contigs in this linkage group
                    double best_ctg_cor_this_lg = 0; 
                    string best_ctg_id_this_lg  = "";                                       
                    map<string, int>::iterator citr;
                    map<string, int>::iterator citr_end;
                    citr     = this_LG_ctgs.begin();
                    citr_end = this_LG_ctgs.end();
                    while(citr != citr_end)
                    {
                        // get contig id 
                        string this_ctg_id = (*citr).first;
                        // get contig end marker
                        map<string, BIPAT_INT>::iterator pitr = contigPollenGTSeq_int_pure_hap.find(this_ctg_id);
                        if( pitr == contigPollenGTSeq_int_pure_hap.end() )
                        {
                            citr ++;
                            continue;
                        }
                        BIPAT_INT pat;
                        pat.leftPat_int  = (*pitr).second.leftPat_int;
                        pat.rightPat_int = (*pitr).second.rightPat_int;  
                        // calculate correlation
                        double cor[2];
                        cor[0] = calculate_correlation_diplotig(&this_gts, &pat.leftPat_int);      
                        cor[1] = calculate_correlation_diplotig(&this_gts, &pat.rightPat_int);      
                        for(int ci=0; ci<2; ci++)
                        {
                            if(cor[ci] > best_ctg_cor_this_lg)
                            {
                                best_ctg_cor_this_lg = cor[ci];
                                best_ctg_id_this_lg  = this_ctg_id;
                            }
                        }
                        // next contig
                        citr ++;
                    }
                    //
                    if(best_ctg_cor_this_lg > best_ctg_cor_this_homocluster)
                    {
                        best_ctg_cor_this_homocluster = best_ctg_cor_this_lg;
                        best_ctg_id_this_homocluster  = best_ctg_id_this_lg;
                        best_linkage_group_id         = tmpLG_id_old;
                    }
                    // next homologous cluster
                    hcitr ++;
                }
                //
                cout << "   output1.2: win-marker="      << this_key 
                     << ", type="                        << this_type  
                     << " will also go to LG="           << gb_group_id[ best_linkage_group_id ]
                     << " (old_group_id:"                << best_linkage_group_id 
                     << ") with best matching haplotig=" << best_ctg_id_this_homocluster
                     << ", correlation.value="           << best_ctg_cor_this_homocluster
                     << endl;
                ofp_linked << this_key                             << "\t"
                           << this_type                            << "\t"
                           << gb_group_id[ best_linkage_group_id ] << "\t"
                           << lg2_homocluster_id[ this_new_group_id_ss.str() ] << "\t"
                       //  << best_linkage_group_id                << "\t"
                           << best_ctg_id_this_homocluster         << "\t"
                           << best_ctg_cor_this_homocluster        << "\t"
                           << "Corr_to_LG"                         << "\t"
                           << "case_1.2"                           << endl; 
                // collect grouping info                            
                // string ctg_id      = keyinfo[0];
                string win_sta_end = keyinfo[1]+"\t"+keyinfo[2];
                string grouping    = "linked";            
                lgGROUP groupinfo;
                groupinfo.type.push_back( this_type );      // type of current window marker: hap/dip/trip/tetrap...
                groupinfo.qlg.push_back(  gb_group_id[ best_linkage_group_id ] ); // LG of query contig  = 1:48
                groupinfo.qchr.push_back( lg2_homocluster_id[ this_new_group_id_ss.str() ] ); // Chr of query contig = 1:12, i.e., homocluster_id
                groupinfo.sctg.push_back( ctg_id );         // id of linked subject contig
                groupinfo.slg.push_back(  gb_group_id[ best_linkage_group_id ] ); // LG of linked subject contig = 1:48
                groupinfo.qscor.push_back(1.2);             // correlation value linked query and subject contigs
                groupinfo.qscas.push_back("case_1.2");      // case info: how this linking is found                                            
                update_grouping(ctg_id, 
                                win_sta_end,
                                grouping, 
                                groupinfo,
                                &raw_grouping);   
                case_collect ++;                                                                  
            }else
            if(this_type.compare("trip")==0)
            {
                // if trip marker 
                /* send another copy to two of the homologous LGs of this_new_group_id */
                // get homocluster_id
                assert( lg2_homocluster_id.find( this_new_group_id_ss.str() )  != lg2_homocluster_id.end() );
                int homocluster_id = lg2_homocluster_id[ this_new_group_id_ss.str() ];                
                // get list of homologous linkage groups
                assert( homocluster_id_2lg.find( homocluster_id ) != homocluster_id_2lg.end() );
                map<string, int>  homo_lgs = homocluster_id_2lg[ homocluster_id ]; // <"1:48", 1>
                // find best correlation between target widnow marker and contigs in these homolougous linkage groups
                double best_ctg_cor_this_homocluster = 0; 
                string best_ctg_id_this_homocluster  = "";
                int    best_linkage_group_id         = 0;
                map<string, int>::iterator hcitr;
                map<string, int>::iterator hcitr_end;
                hcitr     = homo_lgs.begin();
                hcitr_end = homo_lgs.end();
                while(hcitr != hcitr_end)
                {
                    //
                    if( ( (*hcitr).first).compare( this_new_group_id_ss.str() ) ==0 )
                    {
                        hcitr ++;
                        continue;
                    }
                    // get new group id
                    int tmpLG_id_new = atoi( (*hcitr).first.c_str() ); // 1:48
                    assert( gb_group_id_rev.find( tmpLG_id_new ) != gb_group_id_rev.end() );
                    // get old group id
                    int tmpLG_id_old = gb_group_id_rev[ tmpLG_id_new ];
                    // get list of contigs in the linkage group
                    assert( gb_group_2ctg.find( tmpLG_id_old ) != gb_group_2ctg.end() );                   
                    map<string, int> this_LG_ctgs = gb_group_2ctg[ tmpLG_id_old ];
                    // find best correlation between target window marker and contigs in this linkage group
                    double best_ctg_cor_this_lg = 0; 
                    string best_ctg_id_this_lg  = "";                                       
                    map<string, int>::iterator citr;
                    map<string, int>::iterator citr_end;
                    citr     = this_LG_ctgs.begin();
                    citr_end = this_LG_ctgs.end();
                    while(citr != citr_end)
                    {
                        // get contig id 
                        string this_ctg_id = (*citr).first;
                        // get contig end marker
                        map<string, BIPAT_INT>::iterator pitr = contigPollenGTSeq_int_pure_hap.find(this_ctg_id);
                        if( pitr == contigPollenGTSeq_int_pure_hap.end() )
                        {
                            citr ++;
                            continue;
                        }
                        BIPAT_INT pat;
                        pat.leftPat_int  = (*pitr).second.leftPat_int;
                        pat.rightPat_int = (*pitr).second.rightPat_int;  
                        // calculate correlation
                        double cor[2];
                        cor[0] = calculate_correlation_triplotig(&this_gts, &pat.leftPat_int);      
                        cor[1] = calculate_correlation_triplotig(&this_gts, &pat.rightPat_int);      
                        for(int ci=0; ci<2; ci++)
                        {
                            if(cor[ci] > best_ctg_cor_this_lg)
                            {
                                best_ctg_cor_this_lg = cor[ci];
                                best_ctg_id_this_lg  = this_ctg_id;
                            }
                        }
                        // next contig
                        citr ++;
                    }
                    //
                    if(best_ctg_cor_this_lg > best_ctg_cor_this_homocluster)
                    {
                        best_ctg_cor_this_homocluster = best_ctg_cor_this_lg;
                        best_ctg_id_this_homocluster  = best_ctg_id_this_lg;
                        best_linkage_group_id         = tmpLG_id_old;
                    }
                    // next homologous cluster
                    hcitr ++;
                }
                // 
                cout << "   output1.3: haplotig complementing the tripoltig is in " 
                                                          << gb_group_id[ best_linkage_group_id ]
                     << " (old_group_id:"                 << best_linkage_group_id 
                     << "), with best matching haplotig=" << best_ctg_id_this_homocluster
                     << ", correlation.value="            << best_ctg_cor_this_homocluster       
                     << ", therefore, "              
                     << endl;
                // The best match would not get reads; as it is the haplotig complement to this triplotig.
                hcitr     = homo_lgs.begin();
                hcitr_end = homo_lgs.end();
                while(hcitr != hcitr_end)
                {
                    if( ( (*hcitr).first).compare( this_new_group_id_ss.str() ) ==0 )
                    {
                        hcitr ++;
                        continue;
                    }    
                    std::stringstream best_linkage_new_group_id_ss;
                    best_linkage_new_group_id_ss.str("");
                    best_linkage_new_group_id_ss << gb_group_id[ best_linkage_group_id ];                      
                    if( ( (*hcitr).first).compare( best_linkage_new_group_id_ss.str() ) ==0 )
                    {
                        hcitr ++;
                        continue;
                    }                                    
                    cout << "   output1.4: win-marker="        << this_key 
                         << ", type="                          << this_type  
                         << " will also go to LG="             << (*hcitr).first
                         << " (old_group_id:"                  << gb_group_id_rev[ atoi( (*hcitr).first.c_str() ) ]
                         << ") due to best matching haplotig=" << best_ctg_id_this_homocluster
                         << ", correlation.value="             << best_ctg_cor_this_homocluster
                         << ", in linkage group "              << gb_group_id[ best_linkage_group_id ]
                         << " (old_group_id:"                  << best_linkage_group_id << ")"
                         << endl; 
                    ofp_linked << this_key                                          << "\t"
                               << this_type                                         << "\t"
                               << (*hcitr).first                                    << "\t"
                               << lg2_homocluster_id[ this_new_group_id_ss.str() ] << "\t"                               
                           //  << gb_group_id_rev[ atoi( (*hcitr).first.c_str() ) ] << "\t"
                               << best_ctg_id_this_homocluster                      << "\t"
                               << best_ctg_cor_this_homocluster                     << "\t"
                               << "Corr_to_LG"                                      << "\t"
                               << "case_1.4"                                        << endl; 
                    // collect grouping info                                
                    // string ctg_id      = keyinfo[0];
                    string win_sta_end = keyinfo[1]+"\t"+keyinfo[2];
                    string grouping    = "linked";            
                    lgGROUP groupinfo;
                    groupinfo.type.push_back( this_type );      // type of current window marker: hap/dip/trip/tetrap...
                    groupinfo.qlg.push_back( atoi ( (*hcitr).first.c_str() ) );  // LG of query contig  = 1:48
                    groupinfo.qchr.push_back( lg2_homocluster_id[ this_new_group_id_ss.str() ] ); // Chr of query contig = 1:12, i.e., homocluster_id
                    groupinfo.sctg.push_back( ctg_id );         // id of linked subject contig
                    groupinfo.slg.push_back( atoi ( (*hcitr).first.c_str() ) );  // LG of linked subject contig = 1:48
                    groupinfo.qscor.push_back(1.4);             // correlation value linked query and subject contigs
                    groupinfo.qscas.push_back("case_1.4");      // case info: how this linking is found                                            
                    update_grouping(ctg_id, 
                                    win_sta_end,
                                    grouping, 
                                    groupinfo,
                                    &raw_grouping);  
                    case_collect ++;                                                                                            
                    //
                    hcitr ++;    
                }                                                                   
            }else            
            if(this_type.compare("tetrap") == 0)
            {
                //  if tetrap marker and linked with homologous cluster, reads go to each homologous LG
                /* send another copy to three of the homologous LGs of this_new_group_id */
                // get homocluster_id
                assert( lg2_homocluster_id.find( this_new_group_id_ss.str() )  != lg2_homocluster_id.end() );
                int homocluster_id = lg2_homocluster_id[ this_new_group_id_ss.str() ];                
                // get list of homologous linkage groups
                assert( homocluster_id_2lg.find( homocluster_id ) != homocluster_id_2lg.end() );
                map<string, int>  homo_lgs = homocluster_id_2lg[ homocluster_id ]; // <"1:48", 1>                
                //
                map<string, int>::iterator hcitr;
                map<string, int>::iterator hcitr_end;                
                hcitr     = homo_lgs.begin();
                hcitr_end = homo_lgs.end();
                while(hcitr != hcitr_end)
                {
                    if( ( (*hcitr).first).compare( this_new_group_id_ss.str() ) ==0 )
                    {
                        hcitr ++;
                        continue;
                    }                                   
                    cout << "   output1.5: win-marker="  << this_key 
                         << ", type="                    << this_type  
                         << " will also go to LG="       << (*hcitr).first
                         << " (old_group_id:"            << gb_group_id_rev[ atoi( (*hcitr).first.c_str() ) ]
                         << ") due to tetraplotig. "     << endl;
                    ofp_linked << this_key                                          << "\t"
                               << this_type                                         << "\t"
                               << (*hcitr).first                                    << "\t"
                               << lg2_homocluster_id[ this_new_group_id_ss.str() ]  << "\t"                               
                            // << gb_group_id_rev[ atoi( (*hcitr).first.c_str() ) ] << "\t"
                               << "NA.ctg"                                          << "\t"
                               << "NA.cor"                                          << "\t"
                               << "Corr_to_LG"                                      << "\t"
                               << "case_1.5"                                        << endl;          
                    // collect grouping info
                    // string ctg_id      = keyinfo[0];
                    string win_sta_end = keyinfo[1]+"\t"+keyinfo[2];
                    string grouping    = "linked";
                    lgGROUP groupinfo;
                    groupinfo.type.push_back( this_type );      // type of current window marker: hap/dip/trip/tetrap...
                    groupinfo.qlg.push_back( atoi ( (*hcitr).first.c_str() ) );  // LG of query contig  = 1:48
                    groupinfo.qchr.push_back( lg2_homocluster_id[ this_new_group_id_ss.str() ] ); // Chr of query contig = 1:12, i.e., homocluster_id
                    groupinfo.sctg.push_back( ctg_id );         // id of linked subject contig
                    groupinfo.slg.push_back( atoi ( (*hcitr).first.c_str() ) );  // LG of linked subject contig = 1:48
                    groupinfo.qscor.push_back( 1.5 );           // correlation value linked query and subject contigs
                    groupinfo.qscas.push_back("case_1.5");      // case info: how this linking is found                                            
                    update_grouping(ctg_id, 
                                    win_sta_end,
                                    grouping, 
                                    groupinfo,
                                    &raw_grouping); 
                    case_collect ++;                                                                                            
                    //
                    hcitr ++;    
                }                
            }else
            if(this_type.compare("rep") == 0)
            {
                // collect to non-classified rep groups (to give to all linkage groups?)
                cout << "   output1.6: win-marker="                << this_key 
                     << ", type="                                  << this_type  
                     << " will go to all (which) four LGs - TODO " << endl;
                ofp_rep << this_key     << "\t"
                        << this_type    << "\t"
                        << -1           << "\t"
                        << -1           << "\t"
                        << "NA.ctg"     << "\t"
                        << "NA.cor"     << "\t"
                        << "NA.newLG"   << "\t"
                        << "NA.oldLG"   << "\t"
                        << "TODO_to_LG" << "\t"
                        << "case_2.6"   << endl;  
                    // collect grouping info                                
                    // string ctg_id      = keyinfo[0];
                    string win_sta_end = keyinfo[1]+"\t"+keyinfo[2];
                    string grouping    = "unlinked_rep";            
                    lgGROUP groupinfo;
                    groupinfo.type.push_back( this_type ); // type of current window marker: hap/dip/trip/tetrap...
                    groupinfo.qlg.push_back( -1 );         // LG of query contig  = 1:48
                    groupinfo.qchr.push_back( -1 );        // Chr of query contig = 1:12, i.e., homocluster_id
                    groupinfo.sctg.push_back("NA.ctg");    // id of linked subject contig
                    groupinfo.slg.push_back(-1);           // LG of linked subject contig = 1:48
                    groupinfo.qscor.push_back(-1.6);       // correlation value linked query and subject contigs
                    groupinfo.qscas.push_back("case_1.6"); // case info: how this linking is found                                            
                    update_grouping(ctg_id, 
                                    win_sta_end,
                                    grouping, 
                                    groupinfo,
                                    &raw_grouping);  
                    case_collect ++;                                                                                 
            }
            else 
            {
                // nothing to be further done for haplotig markers
                ;
            }
        }else
        {
            //if hap/dip/trip/tetrap/rep not found in existing linkage group: 
            if(this_type.compare("hap") == 0)
            {
                // "hap" not in existing LGs
                double best_ctg_cor = 0; 
                string best_ctg_id  = "";
                int    best_this_old_group_id = 0;                
                int    best_this_new_group_id = 1;                
                // find one best ctg
                map<string, BIPAT_INT>::iterator citr;
                map<string, BIPAT_INT>::iterator citr_end;
                citr     = contigPollenGTSeq_int_pure_hap.begin();
                citr_end = contigPollenGTSeq_int_pure_hap.end();
                while(citr != citr_end)
                {
                    BIPAT_INT pat;
                    pat.leftPat_int  = (*citr).second.leftPat_int;
                    pat.rightPat_int = (*citr).second.rightPat_int;     
                    double cor[2];
                    cor[0] = calculate_correlation_diplotig(&this_gts, &pat.leftPat_int);      
                    cor[1] = calculate_correlation_diplotig(&this_gts, &pat.rightPat_int);                          
                    //
                    double max_cor = cor[0];
                    double max_i   = 0;
                    for(int i = 0; i < 2; i ++)
                    {
                        if(cor[i]>best_ctg_cor)
                        {
                            best_ctg_cor = cor[i];
                            best_ctg_id  = (*citr).first;
                        }
                    }
                    //
                    citr ++;
                }
                // find best old_group_id
                assert(best_ctg_id.size() > 0);
                assert( ctg2_gb_group.find(best_ctg_id) != ctg2_gb_group.end() );
                best_this_old_group_id = ctg2_gb_group[ best_ctg_id ];
                assert( gb_group_id.find(best_this_old_group_id) != gb_group_id.end() );
                best_this_new_group_id = gb_group_id[ best_this_old_group_id ];
                // classify the marker to one pre-define linkage group
                cout << "   output2.1: win-marker="     << this_key 
                     << ", type="                       << this_type  
                     << " will go to LG="               << best_this_new_group_id
                     << " (old_group_id:"               << best_this_old_group_id 
                     << "), with best matching contig " << best_ctg_id
                     << " showing corralation of "      << best_ctg_cor                 
                     << endl;    
                std::stringstream best_this_new_group_id_ss2;
                best_this_new_group_id_ss2.str();
                best_this_new_group_id_ss2 << best_this_new_group_id;
                assert( lg2_homocluster_id.find( best_this_new_group_id_ss2.str() ) != lg2_homocluster_id.end() );                
                ofp_linked << this_key               << "\t"
                           << this_type              << "\t"
                           << best_this_new_group_id << "\t"
                           << lg2_homocluster_id[ best_this_new_group_id_ss2.str() ] << "\t"                           
                         //<< best_this_old_group_id << "\t"
                           << best_ctg_id            << "\t"
                           << best_ctg_cor           << "\t"
                           << "Corr_to_LG"           << "\t"
                           << "case_2.1"             << endl;
                // collect grouping info                                
                // string ctg_id      = keyinfo[0];
                string win_sta_end = keyinfo[1]+"\t"+keyinfo[2];
                string grouping    = "linked";            
                lgGROUP groupinfo;
                groupinfo.type.push_back(  this_type );              // type of current window marker: hap/dip/trip/tetrap...
                groupinfo.qlg.push_back(   best_this_new_group_id ); // LG of query contig  = 1:48
                groupinfo.qchr.push_back(  lg2_homocluster_id[ best_this_new_group_id_ss2.str() ] ); // Chr of query contig = 1:12, i.e., homocluster_id
                groupinfo.sctg.push_back(  best_ctg_id );            // id of linked subject contig
                groupinfo.slg.push_back(   best_this_new_group_id ); // LG of linked subject contig = 1:48
                groupinfo.qscor.push_back( best_ctg_cor );           // correlation value linked query and subject contigs
                groupinfo.qscas.push_back("case_2.1");               // case info: how this linking is found                                            
                update_grouping(ctg_id, 
                                win_sta_end,
                                grouping, 
                                groupinfo,
                                &raw_grouping);  
                case_collect ++;                                                           
            }else
            if(this_type.compare("dip") == 0)
            {                     
                // "dip" not in existing LGs
                // find two best linkage groups;
                map<int, double> best_ctg_cor; // <old_group_id, cor.value> 
                map<int, string> best_ctg_id;  // <old_group_id, ctg_id>
                // initialize correlation score to each linkage group (old_group_id)
                map<int, int>::iterator ogitr;
                map<int, int>::iterator ogitr_end;
                ogitr     = gb_group_id.begin();
                ogitr_end = gb_group_id.end();
                while(ogitr != ogitr_end)
                {
                    int this_old_group_id = (*ogitr).first;
                    best_ctg_cor.insert(std::pair<int, double>(this_old_group_id, 0.0));
                    best_ctg_id.insert(std::pair<int, string>(this_old_group_id, ""));
                    //
                    ogitr ++;
                }           
                // find one best ctg
                map<string, BIPAT_INT>::iterator citr;
                map<string, BIPAT_INT>::iterator citr_end;
                citr     = contigPollenGTSeq_int_pure_hap.begin();
                citr_end = contigPollenGTSeq_int_pure_hap.end();
                while(citr != citr_end)
                {
                    // tmp ctg id to compare
                    string ctg_id_tmp = (*citr).first;                   
                    // tmp linkage group 
                    assert(ctg2_gb_group.find(ctg_id_tmp) != ctg2_gb_group.end() );                    
                    int old_group_id_tmp = ctg2_gb_group[ctg_id_tmp];                    
                    // tmp ctg id genotype sequence
                    BIPAT_INT pat;
                    pat.leftPat_int  = (*citr).second.leftPat_int;
                    pat.rightPat_int = (*citr).second.rightPat_int; 
                    // calculate correlation of the query window marker with tmp ctg      
                    double cor[2];
                    cor[0] = calculate_correlation_diplotig(&this_gts, &pat.leftPat_int);      
                    cor[1] = calculate_correlation_diplotig(&this_gts, &pat.rightPat_int);                          
                    //
                    double max_cor = cor[0];
                    double max_i   = 0;
                    for(int i = 0; i < 2; i ++)
                    {
                        if(cor[i]>best_ctg_cor[old_group_id_tmp])
                        {
                            best_ctg_cor[old_group_id_tmp] = cor[i];
                            best_ctg_id[old_group_id_tmp]  = (*citr).first;
                        }
                    }
                    // next tmp ctg
                    citr ++;
                }
                // top 1: get the linkage group with the highest correlation value
                double best_ctg_cor_1st = 0;                 
                int    best_this_old_group_id = 0;                
                int    best_this_new_group_id = 1;                  
                ogitr     = gb_group_id.begin();
                ogitr_end = gb_group_id.end();
                while(ogitr != ogitr_end)
                {
                    int this_old_group_id = (*ogitr).first;
                    if( best_ctg_cor[this_old_group_id]> best_ctg_cor_1st)
                    {
                        best_ctg_cor_1st       = best_ctg_cor[this_old_group_id];
                        best_this_old_group_id = this_old_group_id;
                    }
                    //
                    ogitr ++;
                }   
                assert( gb_group_id.find(best_this_old_group_id) != gb_group_id.end() );
                best_this_new_group_id = gb_group_id[ best_this_old_group_id ];
                // classify the marker to top 1 matching linkage group                
                cout << "   output2.2: win-marker="       << this_key 
                     << ", type="                         << this_type  
                     << " will go to LG="                 << best_this_new_group_id
                     << " (old_group_id:"                 << best_this_old_group_id 
                     << "), with best matching haplotig " << best_ctg_id[best_this_old_group_id]
                     << " showing correlation value of "  << best_ctg_cor_1st                
                     << endl;    
                std::stringstream best_this_new_group_id_ss2;
                best_this_new_group_id_ss2.str("");
                best_this_new_group_id_ss2 << best_this_new_group_id;  
                assert( lg2_homocluster_id.find( best_this_new_group_id_ss2.str() ) != lg2_homocluster_id.end() );                                                      
                ofp_linked << this_key                            << "\t"
                           << this_type                           << "\t"
                           << best_this_new_group_id              << "\t"
                           << lg2_homocluster_id[ best_this_new_group_id_ss2.str() ] << "\t"                           
                         //<< best_this_old_group_id              << "\t"
                           << best_ctg_id[best_this_old_group_id] << "\t"
                           << best_ctg_cor_1st                    << "\t"
                           << "Corr_to_LG"                        << "\t"
                           << "case_2.2"                          << endl;  
                // collect grouping info                                
                // string ctg_id      = keyinfo[0];
                string win_sta_end = keyinfo[1]+"\t"+keyinfo[2];
                string grouping    = "linked";            
                lgGROUP groupinfo;
                groupinfo.type.push_back(  this_type );              // type of current window marker: hap/dip/trip/tetrap...
                groupinfo.qlg.push_back(   best_this_new_group_id ); // LG of query contig  = 1:48
                groupinfo.qchr.push_back(  lg2_homocluster_id[ best_this_new_group_id_ss2.str() ] ); // Chr of query contig = 1:12, i.e., homocluster_id
                groupinfo.sctg.push_back(  best_ctg_id[best_this_old_group_id] );                    // id of linked subject contig
                groupinfo.slg.push_back(   best_this_new_group_id ); // LG of linked subject contig = 1:48
                groupinfo.qscor.push_back( best_ctg_cor_1st );       // correlation value linked query and subject contigs
                groupinfo.qscas.push_back("case_2.2");               // case info: how this linking is found                                            
                update_grouping(ctg_id, 
                                win_sta_end,
                                grouping, 
                                groupinfo,
                                &raw_grouping);  
                case_collect ++;                                                                             
                // top 1: find homologous cluster: best linkage group related homocluster_id
                int homocluster_id_top1                 = lg2_homocluster_id[ best_this_new_group_id_ss2.str() ];
                assert( homocluster_id_2lg.find(homocluster_id_top1)    != homocluster_id_2lg.end() );
                map<string, int> homocluster_id_all4LGs = homocluster_id_2lg[homocluster_id_top1]; // <"subset(1:48)",1>
                // top 2: get the linkage group with the 2nd-highest correlation value
                double best_ctg_cor_2nd = 0;                 
                int    best_this_old_group_id_2nd = 0;
                int    best_this_new_group_id_2nd = 1;
                map<string, int>::iterator hcitr_tmp;
                map<string, int>::iterator hcitr_tmp_end;
                hcitr_tmp     = homocluster_id_all4LGs.begin();
                hcitr_tmp_end = homocluster_id_all4LGs.end();
                while(hcitr_tmp != hcitr_tmp_end)
                {
                    int new_group_id_tmp = atoi( ((*hcitr_tmp).first).c_str() );
                    int old_group_id_tmp = gb_group_id_rev[ new_group_id_tmp ];
                    if(new_group_id_tmp == best_this_new_group_id)
                    { 
                        hcitr_tmp ++;
                        continue;
                    }
                    //
                    if( best_ctg_cor[old_group_id_tmp]> best_ctg_cor_2nd)
                    {
                        best_ctg_cor_2nd           = best_ctg_cor[old_group_id_tmp];
                        best_this_old_group_id_2nd = old_group_id_tmp;
                    } 
                    //
                    hcitr_tmp ++;
                }
                // classify the marker to top-2 matching linkage group
                assert( gb_group_id.find(best_this_old_group_id_2nd) != gb_group_id.end() );
                best_this_new_group_id_2nd = gb_group_id[ best_this_old_group_id_2nd ];                
                cout << "   output2.3: win-marker="       << this_key 
                     << ", type="                         << this_type  
                     << " will also go to LG="            << best_this_new_group_id_2nd
                     << " (old_group_id:"                 << best_this_old_group_id_2nd 
                     << "), with best matching haplotig " << best_ctg_id[best_this_old_group_id_2nd]
                     << " showing correlation value of "  << best_ctg_cor_2nd                
                     << endl;
                ofp_linked << this_key                                << "\t"
                           << this_type                               << "\t"
                           << best_this_new_group_id_2nd              << "\t"
                           << lg2_homocluster_id[ best_this_new_group_id_ss2.str() ] << "\t"                           
                        // << best_this_old_group_id_2nd              << "\t"
                           << best_ctg_id[best_this_old_group_id_2nd] << "\t"
                           << best_ctg_cor_2nd                        << "\t"
                           << "Corr_to_LG"                            << "\t"
                           << "case_2.3"                              << endl;     
                // collect grouping info                                
                // string ctg_id      = keyinfo[0];
                string win_sta_end2 = keyinfo[1]+"\t"+keyinfo[2];
                string grouping2    = "linked";            
                lgGROUP groupinfo2;
                groupinfo2.type.push_back(  this_type );                  // type of current window marker: hap/dip/trip/tetrap...
                groupinfo2.qlg.push_back(   best_this_new_group_id_2nd ); // LG of query contig  = 1:48
                groupinfo2.qchr.push_back(  lg2_homocluster_id[ best_this_new_group_id_ss2.str() ] ); // Chr of query contig = 1:12, i.e., homocluster_id
                groupinfo2.sctg.push_back(  best_ctg_id[best_this_old_group_id_2nd] );                // id of linked subject contig
                groupinfo2.slg.push_back(   best_this_new_group_id_2nd ); // LG of linked subject contig = 1:48
                groupinfo2.qscor.push_back( best_ctg_cor_2nd );           // correlation value linked query and subject contigs
                groupinfo2.qscas.push_back("case_2.3");                   // case info: how this linking is found                                            
                update_grouping(ctg_id, 
                                win_sta_end2,
                                grouping2, 
                                groupinfo2,
                                &raw_grouping);       
                case_collect ++;
            }else
            if(this_type.compare("trip") == 0)
            {
                // find three best linkage groups
                double best_ctg_cor = 0; 
                string best_ctg_id  = "";
                int    best_this_old_group_id = 0;                
                int    best_this_new_group_id = 1;                
                // find one best ctg
                map<string, BIPAT_INT>::iterator citr;
                map<string, BIPAT_INT>::iterator citr_end;
                citr     = contigPollenGTSeq_int_pure_hap.begin();
                citr_end = contigPollenGTSeq_int_pure_hap.end();
                while(citr != citr_end)
                {                   
                    BIPAT_INT pat;
                    pat.leftPat_int  = (*citr).second.leftPat_int;
                    pat.rightPat_int = (*citr).second.rightPat_int;       
                    double cor[2];
                    cor[0] = calculate_correlation_triplotig(&this_gts, &pat.leftPat_int);      
                    cor[1] = calculate_correlation_triplotig(&this_gts, &pat.rightPat_int);                          
                    //
                    double max_cor = cor[0];
                    double max_i   = 0;
                    for(int i = 0; i < 2; i ++)
                    {
                        if(cor[i]>best_ctg_cor)
                        {
                            best_ctg_cor = cor[i];
                            best_ctg_id  = (*citr).first;
                        }
                    }
                    //
                    citr ++;
                }
                // find best old_group_id
                assert(best_ctg_id.size() > 0);
                assert( ctg2_gb_group.find(best_ctg_id) != ctg2_gb_group.end() );
                best_this_old_group_id = ctg2_gb_group[ best_ctg_id ];
                assert( gb_group_id.find(best_this_old_group_id) != gb_group_id.end() );
                best_this_new_group_id = gb_group_id[ best_this_old_group_id ];                
                // find homologous cluster: best linkage group related homocluster_id
                std::stringstream best_this_new_group_id_ss2;
                best_this_new_group_id_ss2.str("");
                best_this_new_group_id_ss2 << best_this_new_group_id;
                assert( lg2_homocluster_id.find( best_this_new_group_id_ss2.str() ) != lg2_homocluster_id.end() );
                int homocluster_id_top1                 = lg2_homocluster_id[ best_this_new_group_id_ss2.str() ];
                assert( homocluster_id_2lg.find(homocluster_id_top1)    != homocluster_id_2lg.end() );
                map<string, int> homocluster_id_all4LGs = homocluster_id_2lg[homocluster_id_top1]; // <"subset(1:48)",1>                
                // classify the marker to 3 complementing linkage groups (to the one with the highest correlation score)
                map<string, int>::iterator hcitr_tmp;
                map<string, int>::iterator hcitr_tmp_end;
                hcitr_tmp     = homocluster_id_all4LGs.begin();
                hcitr_tmp_end = homocluster_id_all4LGs.end();
                while(hcitr_tmp != hcitr_tmp_end)
                {
                    int new_group_id_tmp = atoi( ((*hcitr_tmp).first).c_str() );
                    int old_group_id_tmp = gb_group_id_rev[ new_group_id_tmp ];
                    if(new_group_id_tmp == best_this_new_group_id)
                    { 
                        // current group id is the one complementing the haplotig
                        hcitr_tmp ++;
                        continue;
                    }
                    //
                    cout << "   output2.4: win-marker="              << this_key 
                         << ", type="                                << this_type  
                         << " will also go to LG="                   << new_group_id_tmp
                         << " (old_group_id:"                        << old_group_id_tmp 
                         << "), due to best complementing haplotig " << best_ctg_id
                         << " showing correlation value of "         << best_ctg_cor 
                         << " in LG="                                << best_this_new_group_id
                         << " (old_group_id:"                        << gb_group_id_rev[ best_this_new_group_id ]
                         << ")"
                         << endl; 
                   std::stringstream ss3;
                   ss3.str("");
                   ss3 << new_group_id_tmp;
                   ofp_linked << this_key         << "\t"
                              << this_type        << "\t"
                              << new_group_id_tmp << "\t"
                              << lg2_homocluster_id[ ss3.str() ]           << "\t"
                          //  << old_group_id_tmp << "\t"
                              << best_ctg_id      << "\t"
                              << best_ctg_cor     << "\t"
                              << best_this_new_group_id                    << "\t"
                              << gb_group_id_rev[ best_this_new_group_id ] << "\t"
                              << "Corr_to_LG"     << "\t"
                              << "case_2.4"       << endl;    
                    // collect grouping info                                
                    // string ctg_id      = keyinfo[0];
                    string win_sta_end = keyinfo[1]+"\t"+keyinfo[2];
                    string grouping    = "linked";            
                    lgGROUP groupinfo;
                    groupinfo.type.push_back(  this_type );                       // type of current window marker: hap/dip/trip/tetrap...
                    groupinfo.qlg.push_back(   new_group_id_tmp );                // LG of query contig  = 1:48
                    groupinfo.qchr.push_back(  lg2_homocluster_id[ ss3.str() ] ); // Chr of query contig = 1:12, i.e., homocluster_id
                    groupinfo.sctg.push_back(  best_ctg_id );                     // id of linked subject contig
                    groupinfo.slg.push_back(   best_this_new_group_id );          // LG of linked subject contig = 1:48
                    groupinfo.qscor.push_back( best_ctg_cor );                    // correlation value linked query and subject contigs
                    groupinfo.qscas.push_back("case_2.4");                        // case info: how this linking is found                                            
                    update_grouping(ctg_id, 
                                    win_sta_end,
                                    grouping, 
                                    groupinfo,
                                    &raw_grouping);   
                    case_collect ++;                                                                                                                                                                  
                    //
                    hcitr_tmp ++;
                }
            }else
            if(this_type.compare("tetrap") == 0)
            {
                // collect to non-classified tetrap groups (to give to all linkage groups?)
                cout << "   output2.5: win-marker="                << this_key 
                     << ", type="                                  << this_type  
                     << " will go to all (which) four LGs - TODO " << endl;
                ofp_tetra << this_key     << "\t"
                          << this_type    << "\t"
                          << -1           << "\t"
                          << -1           << "\t"
                          << "NA.ctg"     << "\t"
                          << "NA.cor"     << "\t"
                          << "NA.newLG"   << "\t"
                          << "NA.oldLG"   << "\t"
                          << "TODO_to_LG" << "\t"
                          << "case_2.5"   << endl;   
                // collect grouping info                                
                // string ctg_id      = keyinfo[0];
                string win_sta_end = keyinfo[1]+"\t"+keyinfo[2];
                string grouping    = "unlinked_tetrap";            
                lgGROUP groupinfo;
                groupinfo.type.push_back(  this_type );// type of current window marker: hap/dip/trip/tetrap...
                groupinfo.qlg.push_back(   -1 );       // LG of query contig  = 1:48
                groupinfo.qchr.push_back(  -1 );       // Chr of query contig = 1:12, i.e., homocluster_id
                groupinfo.sctg.push_back(  "NA.ctg" ); // id of linked subject contig
                groupinfo.slg.push_back(   -1 );       // LG of linked subject contig = 1:48
                groupinfo.qscor.push_back( -2.5 );     // correlation value linked query and subject contigs
                groupinfo.qscas.push_back("case_2.5"); // case info: how this linking is found                                            
                update_grouping(ctg_id, 
                                win_sta_end,
                                grouping, 
                                groupinfo,
                                &raw_grouping);
                case_collect ++;                                                                                                           
            }else
            if(this_type.compare("rep") == 0)
            {
                // collect to non-classified rep groups (to give to all linkage groups?)
                cout << "   output2.6: win-marker="                << this_key 
                     << ", type="                                  << this_type  
                     << " will go to all (which) four LGs - TODO " << endl;
                ofp_rep << this_key     << "\t"
                        << this_type    << "\t"
                        << -1           << "\t"
                        << -1           << "\t"
                        << "NA.ctg"     << "\t"
                        << "NA.cor"     << "\t"
                        << "NA.newLG"   << "\t"
                        << "NA.oldLG"   << "\t"
                        << "TODO_to_LG" << "\t"
                        << "case_2.6"   << endl;
                // collect grouping info   
                // string ctg_id      = keyinfo[0];
                string win_sta_end = keyinfo[1]+"\t"+keyinfo[2];
                string grouping    = "unlinked_rep";
                lgGROUP groupinfo;
                groupinfo.type.push_back(  this_type );// type of current window marker: hap/dip/trip/tetrap...
                groupinfo.qlg.push_back(   -1 );       // LG of query contig  = 1:48
                groupinfo.qchr.push_back(  -1 );       // Chr of query contig = 1:12, i.e., homocluster_id
                groupinfo.sctg.push_back(  "NA.ctg" ); // id of linked subject contig
                groupinfo.slg.push_back(   -1 );       // LG of linked subject contig = 1:48
                groupinfo.qscor.push_back( -2.6 );     // correlation value linked query and subject contigs
                groupinfo.qscas.push_back("case_2.6"); // case info: how this linking is found                                            
                update_grouping(ctg_id, 
                                win_sta_end,
                                grouping, 
                                groupinfo,
                                &raw_grouping); 
                case_collect ++;                                                                           
            }            
            else
            {
		/* for checking purpose:
		map<string, NODE2>*         marker_depth_all,
	        vector<string>              marker_order_keep,
	        map<string, double>         hap_win_ratio,
	        map<int, int>               gb_group_id,
	        map<int, int>               gb_group_id_rev,
	        map<int, map<string, int> > gb_group_2ctg,
	        map<string, int>            ctg2_gb_group,
	        map<int, map<string, int> > homocluster_id_2lg,
	        map<string, int>            lg2_homocluster_id,
	        map<string, BIPAT_INT>      contigPollenGTSeq_int_pure_hap
	        
	        dtt..................: diplotig, triplotig, tetraplotig.......; these are the markers to integrate..............
	        marker_depth_all.....: <ctg_id\twin_sta\twin_end, {<depth_gt_int>, hap/dip/trip/tetrap, ..} >...................
	        marker_order_keep....; <ctg_id\tsta\tend>.....................; input order of key = ctg_id\tsta\tend...........
	        hap_win_ratio........: <ctg_id, hap-win-ratio>................; how much of a contig is haplotig windows........
	        gb_group_id..........: <old_group_id, new_group_id=1:48>......; ................................................      
	        gb_group_id_rev......: <new_group_id=1:48, old_group_id>......; ................................................
	        gb_group_2ctg........: <old_group_id, <ctg_id, new_group_id> >; list of contigs belong to one LG................    
	        ctg2_gb_group........: <ctg_id, old_group_id>.................; using contig to find old LG id..................     
	        homocluster_id_2lg...: <homocluster_id=1:12, <new_group_id="1:48", 1> >.; list of homologous LGs................
        	lg2_homocluster_id...: <new_group_id="1:48", homocluster_id=1:12>.......; new_group_id to find homo-cluster id.. 
	        contigPollenGTSeq_int_pure_hap: <ctg_id, {leftPat_int, rightPat_int}>...; GT. across all pollen at contig ends..
		*/                  
                cout << "   Warning: this should never happen!" << endl;
            }
        }        
        // next marker please!
        mitr ++;
    } 
    //
    cout << "   Info: " << case_collect << " initial window marker linking information collected. " << endl;
    // refine grouping          
    if( !refine_grouping(&raw_grouping,
                         marker_depth_all,
                         &gb_group_2ctg,
                         &ctg2_gb_group,
                         &gb_group_id,
                         &gb_group_id_rev,
                         &homocluster_id_2lg,
                         &lg2_homocluster_id,
                         &contigPollenGTSeq_int_pure_hap,
                         tmpfolder_s4) )
    {
        cout << "   Error: refining grouping of contig-win markers failed. " << endl;
        return false;
    }                 
    // close files
    ofp_limited.close();
    ofp_tetra.close();
    ofp_linked.close();   
    ofp_rep.close(); // how to classify repetitive regions? TODO-20200720
    //
    return true;
}
// refine raw grouping of ctg-window markers based on a complete view of grouping of a contig
bool refine_grouping(map<string, map<string, map<string, lgGROUP> > >* raw_grouping,
                     map<string, NODE2>*                               marker_depth_all,
                     map<int, map<string, int> >*                      gb_group_2ctg,
                     map<string, int>*                                 ctg2_gb_group,
                     map<int, int>*                                    gb_group_id,
	             map<int, int>*                                    gb_group_id_rev,
	             map<int, map<string, int> >*                      homocluster_id_2lg,
     	             map<string, int>*                                 lg2_homocluster_id,
     	             map<string, BIPAT_INT>*                           contigPollenGTSeq_int_pure_hap,
                     string                                            tmpfolder_s4)
{
    /*
        raw_grouping: map of grouping info :: <ctg_id, <"linked/...", <win_sta_end, lgGROUP > > > 
    */
    string ofilename = tmpfolder_s4 + "/s4p5_refine_grouping_tmp.txt";
    ofstream ofp_refine;
    ofp_refine.open(ofilename.c_str(), ios::out);
    if(!ofp_refine.good())
    {
        cout << "   Error: cannot open file " << ofilename << endl;
        return false;
    }
    map<string, map<string, map<string, lgGROUP> > >::iterator citr;
    map<string, map<string, map<string, lgGROUP> > >::iterator citr_end;
    citr     = (*raw_grouping).begin();
    citr_end = (*raw_grouping).end();
    while(citr != citr_end)
    {
        string this_ctg = (*citr).first;
        map<string, map<string, lgGROUP> > typeinfo = (*citr).second;
        ofp_refine << endl;
        ofp_refine << this_ctg << endl; // ctg     
        // linking info and LGs
        map<int, double> homID_best_cor;            // <ID=1:12, highest_correlation>
        map<int, string> homID_best_mkr;            // <ID=1:12, highest_cor_marker="linking/unlinked_*\twin_sta_end">
        map<int, unsigned long> homID_best_mkr_len; // <ID=1:12, marker_size_giving_highest_cor>
        // raw grouping types
        map<string, map<string, lgGROUP> >::iterator gtitr;
        map<string, map<string, lgGROUP> >::iterator gtitr_end;
        gtitr = typeinfo.begin();
        gtitr_end = typeinfo.end();
        while(gtitr != gtitr_end)
        {
            // linked, unlinked_limited_gamete, unlinked_tetra, unlinked_rep
            string this_grouping         = (*gtitr).first; 
            map<string, lgGROUP> wininfo = (*gtitr).second;
            ofp_refine << "\t";
            ofp_refine << this_grouping << endl;  // linked/unlinked_*     
            //
            map<string, lgGROUP>::iterator witr;
            map<string, lgGROUP>::iterator witr_end;
            witr     = wininfo.begin();
            witr_end = wininfo.end();
            while(witr != witr_end)
            {
                string win_sta_end = (*witr).first;
                lgGROUP groupinfo  = (*witr).second;
                vector<string> win_info = split_string(win_sta_end, '\t');
                unsigned long win_sta = strtoul(win_info[0].c_str(), NULL, 0);
                unsigned long win_end = strtoul(win_info[1].c_str(), NULL, 0);
                // when selecting the best homID, only consider larger windows!
                //
                ofp_refine << "\t\t";
                ofp_refine << win_sta_end << endl;// win_sta_end                
                //                
                vector<string>  type = groupinfo.type;  // type of current window marker: hap/dip/trip/tetrap...
	        vector<int>      qlg = groupinfo.qlg;   // LG of query contig = 1:48
	        vector<int>     qchr = groupinfo.qchr;  // Chr of query contig = 1:12, i.e., homocluster_id
	        vector<string>  sctg = groupinfo.sctg;  // id of linked subject contig
	        vector<int>      slg = groupinfo.slg;   // LG of linked subject contig = 1:48
     	        vector<double> qscor = groupinfo.qscor; // correlation value linked query and subject contigs
                vector<string> qscas = groupinfo.qscas; // case info: how this linking is found                                                
                for(int wi=0; wi<type.size(); wi++)
                {
                    ofp_refine << "\t\t\t";
                    ofp_refine << type[wi]  << "\t"
                               << qlg[wi]   << "\t"
                               << qchr[wi]  << "\t"
                               << sctg[wi]  << "\t"
                               << slg[wi]   << "\t"
                               << qscor[wi] << "\t"  
                               << qscas[wi] << endl;
                    //
                    if(homID_best_cor.find( qchr[wi] ) == homID_best_cor.end() )
                    {
                        homID_best_cor.insert(std::pair<int, double>(qchr[wi], qscor[wi] ));
                        homID_best_mkr.insert(std::pair<int, string>(qchr[wi], this_grouping + "\t" + win_sta_end) );
                        homID_best_mkr_len.insert(std::pair<int, unsigned long>(qchr[wi], win_end - win_sta + 1));
                    }else
                    {
                        // using the longer marker to compare the best!
                        if(homID_best_mkr_len[ qchr[wi] ] < win_end-win_sta+1)
                        {
                            homID_best_cor[ qchr[wi] ]     = qscor[wi];
                            homID_best_mkr[ qchr[wi] ]     = this_grouping + "\t" + win_sta_end;
                            homID_best_mkr_len[ qchr[wi] ] = win_end - win_sta + 1;
                        }else
                        if(homID_best_mkr_len[ qchr[wi] ] == win_end-win_sta+1)
                        {
                            if(homID_best_cor[ qchr[wi] ] < qscor[wi])
                            {
                                homID_best_cor[ qchr[wi] ]     = qscor[wi];
                                homID_best_mkr[ qchr[wi] ]     = this_grouping + "\t" + win_sta_end;
                                homID_best_mkr_len[ qchr[wi] ] = win_end - win_sta + 1;                            
                            }
                        }else ;
                    }
                }
                // next window marker 
                witr ++;
            }
            // next raw group type 
            gtitr ++;
        }
        // check 
        // >1: more than only "-1" found so there is a way to relink
        if(homID_best_cor.find(-1)!=homID_best_cor.end() && homID_best_cor.size() == 1)
        {
            ofp_refine << "\tCase 1: re-linking unlinked impossible. " << endl;                
        }else
        if(homID_best_cor.size() > 1)
        {
            ofp_refine << "\tCase 2: re-linking un-/linked due to " << homID_best_cor.size() 
                       << " hom-clusters found. "                   << endl;     
            // find the ID with the highest cor score; and match the others in the respecitve 4 LGs linking this ID
            // map<int, double> homID_best_cor;   // <ID=1:12, highest_correlation>
            // map<int, string> homID_best_mkr;   // <ID=1:12, highest_cor_marker="linking/unlinked_*\twin_sta_end">            
            int    best_homID = -1;
            double best_cor   = -1.0;
            long target_win_size = 50000;// NOTE and TODO: we need an input option here!                  
            while(best_homID == -1 && target_win_size>=0)
            {
                ofp_refine << "      check: looking for best homID under target_win_size = " << target_win_size << endl;
                map<int, double>::iterator hitr;
                map<int, double>::iterator hitr_end;
                hitr     = homID_best_cor.begin();
                hitr_end = homID_best_cor.end();            
                while(hitr != hitr_end)
                {
                    assert(homID_best_mkr_len.find( (*hitr).first ) != homID_best_mkr_len.end() );
                    if((*hitr).second > best_cor && homID_best_mkr_len[ (*hitr).first ] >= target_win_size)
                    {
                        ofp_refine << "       check: new cor="     << (*hitr).second 
                                   << " at homologous cluster="    << (*hitr).first
                                   << " versus current best "      << best_cor 
                                   << " at homologous cluster="    << best_homID
                                   << endl;                    
                        best_homID = (*hitr).first;
                        best_cor   = (*hitr).second;
                        ofp_refine << "       updated: current best cor=" << best_cor 
                                   << " at homologous cluster="           << best_homID
                                   << endl;                    
                    }
                    hitr ++;
                }
                //
                target_win_size -= 10000;
                if(best_homID != -1) break;
            }                         
            assert(best_homID != -1);         
            if(1)
            {       
                ofp_refine << "\t\tBest hom id is " 
                           << best_homID 
                           << "; markers not assigned here need to be re-linked. " 
                           << endl;  
                // get LGs related to best_homID: homocluster_id_2lg...: <homocluster_id=1:12, <new_group_id="1:48", 1> >.; 
                //     list of homologous LGs
                assert( (*homocluster_id_2lg).find( best_homID ) != (*homocluster_id_2lg).end() );
                map<string, int> hom_LGs = (*homocluster_id_2lg)[best_homID];
                // all below need to work with iterator to update raw_grouping; 
                map<string, map<string, lgGROUP> >::iterator gtitr; // <"linked/unlinked_*", <win_sta_end, lgGROUP > >
                map<string, map<string, lgGROUP> >::iterator gtitr_end;
                gtitr     = (*citr).second.begin();
                gtitr_end = (*citr).second.end();
                while(gtitr != gtitr_end)
                {
                    if( ((*gtitr).first).compare("unlinked_limited_gamete") == 0 || 
                        ((*gtitr).first).compare("unlinked_rep") == 0 )
                    {
                        ofp_refine << "\t\tsuch a case refers to organelle genomes / repeats - not updated!" 
                                   << endl;                  
                    }else
                    if( ((*gtitr).first).compare("unlinked_tetrap") == 0 )
                    {
                        map<string, lgGROUP>::iterator witr;
                        map<string, lgGROUP>::iterator witr_end;
                        witr     = (*gtitr).second.begin() ;
                        witr_end = (*gtitr).second.end();
                        while(witr != witr_end)
                        {
                            string win_sta_end = (*witr).first;	                
	                    (*witr).second.type.clear();  // type of current window marker: hap/dip/trip/tetrap...
    	                    (*witr).second.qlg.clear();   // LG of query contig = 1:48
	                    (*witr).second.qchr.clear();  // Chr of query contig = 1:12, i.e., homocluster_id
	                    (*witr).second.sctg.clear();  // id of linked subject contig
	                    (*witr).second.slg.clear();   // LG of linked subject contig = 1:48
 	                    (*witr).second.qscor.clear(); // correlation value linked query and subject contigs
    	                    (*witr).second.qscas.clear(); // case info: how this linking is found    
    	                    // assigned this window marker to each of the 4 hom LGs    	                
    	                    map<string, int>::iterator update_lg_itr;
    	                    map<string, int>::iterator update_lg_itr_end;
    	                    update_lg_itr     = hom_LGs.begin();
    	                    update_lg_itr_end = hom_LGs.end();
    	                    while(update_lg_itr != update_lg_itr_end)
    	                    {
                                string grouping    = "unlinked_tetrap";  // already linked now!!
                                lgGROUP groupinfo;
                                groupinfo.type.push_back(  "tetrap" );   // type of current window marker: hap/dip/trip/tetrap...
                                groupinfo.qlg.push_back(   atoi((*update_lg_itr).first.c_str()) );// LG of query contig  = 1:48
                                groupinfo.qchr.push_back(  best_homID ); // Chr of query contig = 1:12, i.e., homocluster_id
                                groupinfo.sctg.push_back(  "NA.ctg" );   // id of linked subject contig
                                groupinfo.slg.push_back(   atoi((*update_lg_itr).first.c_str()) );// LG of linked subject contig = 1:48
                                groupinfo.qscor.push_back( 3.0 );        // correlation value linked query and subject contigs
                                groupinfo.qscas.push_back("reln_3.0");   // case info: how this linking is found
                                update_grouping(this_ctg, 
                                                win_sta_end,
                                                grouping, 
                                                groupinfo,
                                                raw_grouping);   
                                ofp_refine << "\t\tthis tetrap marker " << win_sta_end
                                           << " re-linked to "          << atoi((*update_lg_itr).first.c_str())
                                           << " in hom-cluster "        << best_homID 
                                           << endl;
    	                        //
    	                        update_lg_itr ++;
    	                    }
    	                    //                
                            witr ++;
                        }
                    }else
                    {
                        // correct existing linkage group assignment: (*witr).second = lgGROUP
                        // if it is not consistent with best_homID, re-link the window markers
                        map<string, lgGROUP>::iterator witr;
                        map<string, lgGROUP>::iterator witr_end;
                        witr     = (*gtitr).second.begin() ;
                        witr_end = (*gtitr).second.end();
                        while(witr != witr_end)
                        {
                            string win_sta_end = (*witr).first;	                
                            // not in the correct homologous cluster of linkage groups
                            //                
	                    vector<string>  type = ((*witr).second).type;  // type of current window marker: hap/dip/trip/tetrap...
	                    vector<int>      qlg = ((*witr).second).qlg;   // LG of query contig = 1:48
	                    vector<int>     qchr = ((*witr).second).qchr;  // Chr of query contig = 1:12, i.e., homocluster_id
	                    vector<string>  sctg = ((*witr).second).sctg;  // id of linked subject contig
	                    vector<int>      slg = ((*witr).second).slg;   // LG of linked subject contig = 1:48
 	                    vector<double> qscor = ((*witr).second).qscor; // correlation value linked query and subject contigs
                            vector<string> qscas = ((*witr).second).qscas; // case info: how this linking is found                                                
                            for(int wi=0; wi<1; wi++) // no need to repeat with wi<type.size()
                            {
                                if( qchr[wi] != best_homID)
                                {
                                    ofp_refine << "\t\tthis "                  << type[wi] 
                                               << " marker "                   << win_sta_end 
                                               << " currently in hom-cluster " << qchr[wi] 
                                               << " need to re-link to "       << best_homID
                                               << endl;
                                    /* relinking current win marker */
                                    // relinking begins
                                    // get its depth genotype sequence 
                                    // find best correlation to related linkage groups: need to check hap/dip/trip/tetrap
                                    // get key.........: ctg_id\tsta\tend, where win_sta_end="sta\tend"
                                    string this_key = this_ctg + "\t" + win_sta_end;     
                                    // get node........: {<depth_gt_int>, hap/dip/trip/tetrap, ..}                                                                   
                                    assert((*marker_depth_all).find(this_key) != (*marker_depth_all).end());
                                    NODE2 tmpnode2  = (*marker_depth_all)[this_key];  
                                    // get type........: hap/dip/trip/tetrap              
                                    string this_type     = tmpnode2.type;
                                    // get genotype seq: depth_genotype at all pollen
                                    vector<int> this_gts = tmpnode2.depth_gt;  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                                    
//////////////////////////////////////////////////// re-linking begins /////////////////////////////////////////////////
                                    //if hap/dip/trip/tetrap/rep not found in existing linkage group: 
                                    if(this_type.compare("hap") == 0)
                                    {
                                        // "hap" not in existing LGs
                                        double best_ctg_cor = 0; 
                                        string best_ctg_id  = "";
                                        int    best_this_old_group_id = 0;                
                                        int    best_this_new_group_id = 1;                
                                        // find one best ctg
                                        map<string, BIPAT_INT>::iterator re_citr;
                                        map<string, BIPAT_INT>::iterator re_citr_end;
                                        re_citr     = (*contigPollenGTSeq_int_pure_hap).begin();
                                        re_citr_end = (*contigPollenGTSeq_int_pure_hap).end();
                                        while(re_citr != re_citr_end)
                                        {
                                            // check if ctg is in the best homologous cluster: hom_LGs <"1:48", 1>
                                            string subj_ctg_id = (*re_citr).first;
                                            assert((*ctg2_gb_group).find( subj_ctg_id ) != (*ctg2_gb_group).end());
                                            int subj_old_group_id = (*ctg2_gb_group)[subj_ctg_id];
                                            assert((*gb_group_id).find( subj_old_group_id ) != (*gb_group_id).end() );
                                            int subj_new_group_id = (*gb_group_id)[ subj_old_group_id ];
                                            std::stringstream subj_new_group_id_ss;
                                            subj_new_group_id_ss << subj_new_group_id;
                                            // hom_LGs = (*homocluster_id_2lg)[best_homID]
                                            if(hom_LGs.find( subj_new_group_id_ss.str() ) == hom_LGs.end() )
                                            {
                                                re_citr ++;
                                                continue;
                                            }
                                            //
                                            BIPAT_INT pat;
                                            pat.leftPat_int  = (*re_citr).second.leftPat_int;
                                            pat.rightPat_int = (*re_citr).second.rightPat_int;     
                                            double cor[2];
                                            cor[0] = calculate_correlation_diplotig(&this_gts, &pat.leftPat_int);      
                                            cor[1] = calculate_correlation_diplotig(&this_gts, &pat.rightPat_int);                          
                                            //
                                            double max_cor = cor[0];
                                            double max_i   = 0;
                                            for(int i = 0; i < 2; i ++)
                                            {
                                                if(cor[i]>best_ctg_cor)
                                                {
                                                    best_ctg_cor = cor[i];
                                                    best_ctg_id  = (*re_citr).first;
                                                }
                                            }
                                            //
                                            re_citr ++;
                                        }
                                        // find best old_group_id
                                        assert(best_ctg_id.size() > 0);
                                        assert( (*ctg2_gb_group).find(best_ctg_id) != (*ctg2_gb_group).end() );
                                        best_this_old_group_id = (*ctg2_gb_group)[ best_ctg_id ];
                                        assert( (*gb_group_id).find(best_this_old_group_id) != (*gb_group_id).end() );
                                        best_this_new_group_id = (*gb_group_id)[ best_this_old_group_id ];
                                        // classify the marker to one pre-define linkage group
                                        cout << "   output4.1: win-marker="     << this_key 
                                             << ", type="                       << this_type  
                                             << " will go to LG="               << best_this_new_group_id
                                             << " (old_group_id:"               << best_this_old_group_id 
                                             << "), with best matching contig " << best_ctg_id
                                             << " showing corralation of "      << best_ctg_cor                 
                                             << endl;    
                                        std::stringstream best_this_new_group_id_ss2;
                                        best_this_new_group_id_ss2.str();
                                        best_this_new_group_id_ss2 << best_this_new_group_id;
                                        assert( (*lg2_homocluster_id).find( best_this_new_group_id_ss2.str() ) != (*lg2_homocluster_id).end() );                
                                        ofp_refine << "\t\t\t"
                                                   << this_key               << "\t"
                                                   << this_type              << "\t"
                                                   << best_this_new_group_id << "\t"
                                                   << (*lg2_homocluster_id)[ best_this_new_group_id_ss2.str() ] << "\t"                           
                                                 //<< best_this_old_group_id << "\t"
                                                   << best_ctg_id            << "\t"
                                                   << best_ctg_cor           << "\t"
                                                   << "Corr_to_LG"           << "\t"
                                                   << "case_4.1"             << endl;
                                        ofp_refine << "\t\t\t updating wi="    << wi << endl;
                                        // update grouping info
	                                ((*witr).second).type[wi]  = this_type;    // type of current window marker: hap/dip/trip/tetrap...
	                                ((*witr).second).qlg[wi]   = best_this_new_group_id; // LG of query contig = 1:48
	                                ((*witr).second).qchr[wi]  = (*lg2_homocluster_id)[ best_this_new_group_id_ss2.str() ]; // Chr of query contig = 1:12, i.e., homocluster_id
	                                ((*witr).second).sctg[wi]  = best_ctg_id;  // id of linked subject contig
	                                ((*witr).second).slg[wi]   = best_this_new_group_id; // LG of linked subject contig = 1:48
 	                                ((*witr).second).qscor[wi] = best_ctg_cor; // correlation value linked query and subject contigs
                                        ((*witr).second).qscas[wi] = "case_4.1";   // case info: how this linking is found                                         
                                    }else
                                    if(this_type.compare("dip") == 0)
                                    {                     
                                        // "dip" not in existing LGs
                                        // find two best linkage groups;
                                        map<int, double> best_ctg_cor; // <old_group_id, cor.value> 
                                        map<int, string> best_ctg_id;  // <old_group_id, ctg_id>
                                        // initialize correlation score to each linkage group (old_group_id)
                                        map<int, int>::iterator ogitr;
                                        map<int, int>::iterator ogitr_end;
                                        ogitr     = (*gb_group_id).begin();
                                        ogitr_end = (*gb_group_id).end();
                                        while(ogitr != ogitr_end)
                                        {
                                            int this_old_group_id = (*ogitr).first;
                                            best_ctg_cor.insert(std::pair<int, double>(this_old_group_id, 0.0));
                                            best_ctg_id.insert(std::pair<int, string>(this_old_group_id, ""));
                                            //
                                            ogitr ++;
                                        }           
                                        // find one best ctg
                                        map<string, BIPAT_INT>::iterator re_citr;
                                        map<string, BIPAT_INT>::iterator re_citr_end;
                                        re_citr     = (*contigPollenGTSeq_int_pure_hap).begin();
                                        re_citr_end = (*contigPollenGTSeq_int_pure_hap).end();
                                        while(re_citr != re_citr_end)
                                        {
                                            // check if ctg is in the best homologous cluster: hom_LGs <"1:48", 1>
                                            string subj_ctg_id = (*re_citr).first;
                                            assert((*ctg2_gb_group).find( subj_ctg_id ) != (*ctg2_gb_group).end());
                                            int subj_old_group_id = (*ctg2_gb_group)[subj_ctg_id];
                                            assert((*gb_group_id).find( subj_old_group_id ) != (*gb_group_id).end() );
                                            int subj_new_group_id = (*gb_group_id)[ subj_old_group_id ];
                                            std::stringstream subj_new_group_id_ss;
                                            subj_new_group_id_ss << subj_new_group_id;
                                            // hom_LGs = (*homocluster_id_2lg)[best_homID]
                                            if(hom_LGs.find( subj_new_group_id_ss.str() ) == hom_LGs.end() )
                                            {
                                                re_citr ++;
                                                continue;
                                            }                                                                                                          
                                            // tmp linkage group 
                                            int old_group_id_tmp = subj_old_group_id;                    
                                            // tmp ctg id genotype sequence
                                            BIPAT_INT pat;
                                            pat.leftPat_int  = (*re_citr).second.leftPat_int;
                                            pat.rightPat_int = (*re_citr).second.rightPat_int; 
                                            // calculate correlation of the query window marker with tmp ctg      
                                            double cor[2];
                                            cor[0] = calculate_correlation_diplotig(&this_gts, &pat.leftPat_int);      
                                            cor[1] = calculate_correlation_diplotig(&this_gts, &pat.rightPat_int);                          
                                            //
                                            double max_cor = cor[0];
                                            double max_i   = 0;
                                            for(int i = 0; i < 2; i ++)
                                            {
                                                if(cor[i]>best_ctg_cor[old_group_id_tmp])
                                                {
                                                    best_ctg_cor[old_group_id_tmp] = cor[i];
                                                    best_ctg_id[old_group_id_tmp]  = (*re_citr).first;
                                                }
                                            }
                                            // next tmp ctg
                                            re_citr ++;
                                        }
                                        // top 1: get the linkage group with the highest correlation value
                                        double best_ctg_cor_1st = 0;                 
                                        int    best_this_old_group_id = 0;                
                                        int    best_this_new_group_id = 1;                  
                                        ogitr     = (*gb_group_id).begin();
                                        ogitr_end = (*gb_group_id).end();
                                        while(ogitr != ogitr_end)
                                        {
                                            int this_old_group_id = (*ogitr).first;
                                            if( best_ctg_cor[this_old_group_id]> best_ctg_cor_1st)
                                            {
                                                best_ctg_cor_1st       = best_ctg_cor[this_old_group_id];
                                                best_this_old_group_id = this_old_group_id;
                                            }
                                            //
                                            ogitr ++;
                                        }   
                                        assert( (*gb_group_id).find(best_this_old_group_id) != (*gb_group_id).end() );
                                        best_this_new_group_id = (*gb_group_id)[ best_this_old_group_id ];
                                        // classify the marker to top 1 matching linkage group                
                                        cout << "   output4.2: win-marker="       << this_key 
                                             << ", type="                         << this_type  
                                             << " will go to LG="                 << best_this_new_group_id
                                             << " (old_group_id:"                 << best_this_old_group_id 
                                             << "), with best matching haplotig " << best_ctg_id[best_this_old_group_id]
                                             << " showing correlation value of "  << best_ctg_cor_1st                
                                             << endl;    
                                        std::stringstream best_this_new_group_id_ss2;
                                        best_this_new_group_id_ss2.str("");
                                        best_this_new_group_id_ss2 << best_this_new_group_id;  
                                        assert( (*lg2_homocluster_id).find( best_this_new_group_id_ss2.str() ) != (*lg2_homocluster_id).end() );                                                      
                                        ofp_refine << "\t\t\t"
                                                   << this_key                            << "\t"
                                                   << this_type                           << "\t"
                                                   << best_this_new_group_id              << "\t"
                                                   << (*lg2_homocluster_id)[ best_this_new_group_id_ss2.str() ] << "\t"                           
                                                 //<< best_this_old_group_id              << "\t"
                                                   << best_ctg_id[best_this_old_group_id] << "\t"
                                                   << best_ctg_cor_1st                    << "\t"
                                                   << "Corr_to_LG"                        << "\t"
                                                   << "case_4.2"                          << endl;  
                                        ofp_refine << "\t\t\t updating wi="         << wi << endl;                                                   
                                        // update grouping info
	                                ((*witr).second).type[0]  = this_type;    // type of current window marker: hap/dip/trip/tetrap...
	                                ((*witr).second).qlg[0]   = best_this_new_group_id; // LG of query contig = 1:48
	                                ((*witr).second).qchr[0]  = (*lg2_homocluster_id)[ best_this_new_group_id_ss2.str() ]; // Chr of query contig = 1:12, i.e., homocluster_id
	                                ((*witr).second).sctg[0]  = best_ctg_id[best_this_old_group_id];  // id of linked subject contig
	                                ((*witr).second).slg[0]   = best_this_new_group_id; // LG of linked subject contig = 1:48
 	                                ((*witr).second).qscor[0] = best_ctg_cor_1st; // correlation value linked query and subject contigs
                                        ((*witr).second).qscas[0] = "case_4.2";   // case info: how this linking is found 
                                        //                                                                            
                                        // top 1: find homologous cluster: best linkage group related homocluster_id
                                        int homocluster_id_top1    = (*lg2_homocluster_id)[ best_this_new_group_id_ss2.str() ];
                                        assert( (*homocluster_id_2lg).find(homocluster_id_top1) != (*homocluster_id_2lg).end() );
                                        map<string, int> homocluster_id_all4LGs = (*homocluster_id_2lg)[homocluster_id_top1]; // <"subset(1:48)",1>
                                        // top 2: get the linkage group with the 2nd-highest correlation value
                                        double best_ctg_cor_2nd = 0;                 
                                        int    best_this_old_group_id_2nd = 0;
                                        int    best_this_new_group_id_2nd = 1;
                                        map<string, int>::iterator hcitr_tmp;
                                        map<string, int>::iterator hcitr_tmp_end;
                                        hcitr_tmp     = homocluster_id_all4LGs.begin();
                                        hcitr_tmp_end = homocluster_id_all4LGs.end();
                                        while(hcitr_tmp != hcitr_tmp_end)
                                        {
                                            int new_group_id_tmp = atoi( ((*hcitr_tmp).first).c_str() );
                                            int old_group_id_tmp = (*gb_group_id_rev)[ new_group_id_tmp ];
                                            if(new_group_id_tmp == best_this_new_group_id)
                                            { 
                                                hcitr_tmp ++;
                                                continue;
                                            }
                                            //
                                            if( best_ctg_cor[old_group_id_tmp]> best_ctg_cor_2nd)
                                            {
                                                best_ctg_cor_2nd           = best_ctg_cor[old_group_id_tmp];
                                                best_this_old_group_id_2nd = old_group_id_tmp;
                                            } 
                                            //
                                            hcitr_tmp ++;
                                        }
                                        // classify the marker to top-2 matching linkage group
                                        assert( (*gb_group_id).find(best_this_old_group_id_2nd) != (*gb_group_id).end() );
                                        best_this_new_group_id_2nd = (*gb_group_id)[ best_this_old_group_id_2nd ];                
                                        cout << "   output4.3: win-marker="       << this_key 
                                             << ", type="                         << this_type  
                                             << " will also go to LG="            << best_this_new_group_id_2nd
                                             << " (old_group_id:"                 << best_this_old_group_id_2nd 
                                             << "), with best matching haplotig " << best_ctg_id[best_this_old_group_id_2nd]
                                             << " showing correlation value of "  << best_ctg_cor_2nd                
                                             << endl;
                                        ofp_refine << this_key                                << "\t"
                                                   << this_type                               << "\t"
                                                   << best_this_new_group_id_2nd              << "\t"
                                                   << (*lg2_homocluster_id)[ best_this_new_group_id_ss2.str() ] << "\t"                           
                                                // << best_this_old_group_id_2nd              << "\t"
                                                   << best_ctg_id[best_this_old_group_id_2nd] << "\t"
                                                   << best_ctg_cor_2nd                        << "\t"
                                                   << "Corr_to_LG"                            << "\t"
                                                   << "case_4.3"                              << endl;     
                                        ofp_refine << "\t\t\t updating wi="             << wi << endl;                                                    
                                        // update grouping info
	                                ((*witr).second).type[1]  = this_type;    // type of current window marker: hap/dip/trip/tetrap...
	                                ((*witr).second).qlg[1]   = best_this_new_group_id_2nd; // LG of query contig = 1:48
	                                ((*witr).second).qchr[1]  = (*lg2_homocluster_id)[ best_this_new_group_id_ss2.str() ]; // Chr of query contig = 1:12, i.e., homocluster_id
	                                ((*witr).second).sctg[1]  = best_ctg_id[best_this_old_group_id_2nd];  // id of linked subject contig
	                                ((*witr).second).slg[1]   = best_this_new_group_id_2nd; // LG of linked subject contig = 1:48
 	                                ((*witr).second).qscor[1] = best_ctg_cor_2nd; // correlation value linked query and subject contigs
                                        ((*witr).second).qscas[1] = "case_4.3";   // case info: how this linking is found 
                                    }else
                                    if(this_type.compare("trip") == 0)
                                    {
                                        // find three best linkage groups
                                        double best_ctg_cor = 0; 
                                        string best_ctg_id  = "";
                                        int    best_this_old_group_id = 0;                
                                        int    best_this_new_group_id = 1;                
                                        // find one best ctg
                                        map<string, BIPAT_INT>::iterator re_citr;
                                        map<string, BIPAT_INT>::iterator re_citr_end;
                                        re_citr     = (*contigPollenGTSeq_int_pure_hap).begin();
                                        re_citr_end = (*contigPollenGTSeq_int_pure_hap).end();
                                        while(re_citr != re_citr_end)
                                        {                   
                                            // check if ctg is in the best homologous cluster: hom_LGs <"1:48", 1>
                                            string subj_ctg_id = (*re_citr).first;
                                            assert((*ctg2_gb_group).find( subj_ctg_id ) != (*ctg2_gb_group).end());
                                            int subj_old_group_id = (*ctg2_gb_group)[subj_ctg_id];
                                            assert((*gb_group_id).find( subj_old_group_id ) != (*gb_group_id).end() );
                                            int subj_new_group_id = (*gb_group_id)[ subj_old_group_id ];
                                            std::stringstream subj_new_group_id_ss;
                                            subj_new_group_id_ss << subj_new_group_id;
                                            // hom_LGs = (*homocluster_id_2lg)[best_homID]
                                            if(hom_LGs.find( subj_new_group_id_ss.str() ) == hom_LGs.end() )
                                            {
                                                re_citr ++;
                                                continue;
                                            }
                                            BIPAT_INT pat;
                                            pat.leftPat_int  = (*re_citr).second.leftPat_int;
                                            pat.rightPat_int = (*re_citr).second.rightPat_int;       
                                            double cor[2];
                                            cor[0] = calculate_correlation_triplotig(&this_gts, &pat.leftPat_int);      
                                            cor[1] = calculate_correlation_triplotig(&this_gts, &pat.rightPat_int);                          
                                            //
                                            double max_cor = cor[0];
                                            double max_i   = 0;
                                            for(int i = 0; i < 2; i ++)
                                            {
                                                if(cor[i]>best_ctg_cor)
                                                {
                                                    best_ctg_cor = cor[i];
                                                    best_ctg_id  = (*re_citr).first;
                                                }
                                            }
                                            //
                                            re_citr ++;
                                        }
                                        // find best old_group_id
                                        assert(best_ctg_id.size() > 0);
                                        assert( (*ctg2_gb_group).find(best_ctg_id) != (*ctg2_gb_group).end() );
                                        best_this_old_group_id = (*ctg2_gb_group)[ best_ctg_id ];
                                        assert( (*gb_group_id).find(best_this_old_group_id) != (*gb_group_id).end() );
                                        best_this_new_group_id = (*gb_group_id)[ best_this_old_group_id ];                
                                        // find homologous cluster: best linkage group related homocluster_id
                                        std::stringstream best_this_new_group_id_ss2;
                                        best_this_new_group_id_ss2.str("");
                                        best_this_new_group_id_ss2 << best_this_new_group_id;
                                        assert( (*lg2_homocluster_id).find( best_this_new_group_id_ss2.str() ) != (*lg2_homocluster_id).end() );
                                        int homocluster_id_top1                 = (*lg2_homocluster_id)[ best_this_new_group_id_ss2.str() ];
                                        assert( (*homocluster_id_2lg).find(homocluster_id_top1)    != (*homocluster_id_2lg).end() );
                                        map<string, int> homocluster_id_all4LGs = (*homocluster_id_2lg)[homocluster_id_top1]; // <"subset(1:48)",1>                
                                        // classify the marker to 3 complementing linkage groups (to the one with the highest correlation score)
                                        map<string, int>::iterator hcitr_tmp;
                                        map<string, int>::iterator hcitr_tmp_end;
                                        hcitr_tmp     = homocluster_id_all4LGs.begin();
                                        hcitr_tmp_end = homocluster_id_all4LGs.end();
                                        int upi = 0;
                                        while(hcitr_tmp != hcitr_tmp_end)
                                        {
                                            int new_group_id_tmp = atoi( ((*hcitr_tmp).first).c_str() );
                                            int old_group_id_tmp = (*gb_group_id_rev)[ new_group_id_tmp ];
                                            if(new_group_id_tmp == best_this_new_group_id)
                                            { 
                                                // current group id is the one complementing the haplotig
                                                hcitr_tmp ++;
                                                continue;
                                            }
                                            //
                                            cout << "   output4.4: win-marker="              << this_key 
                                                 << ", type="                                << this_type  
                                                 << " will also go to LG="                   << new_group_id_tmp
                                                 << " (old_group_id:"                        << old_group_id_tmp 
                                                 << "), due to best complementing haplotig " << best_ctg_id
                                                 << " showing correlation value of "         << best_ctg_cor 
                                                 << " in LG="                                << best_this_new_group_id
                                                 << " (old_group_id:"                        << (*gb_group_id_rev)[ best_this_new_group_id ]
                                                 << ")"
                                                 << endl; 
                                           std::stringstream ss3;
                                           ss3.str("");
                                           ss3 << new_group_id_tmp;
                                           ofp_refine << "\t\t\t"
                                                      << this_key         << "\t"
                                                      << this_type        << "\t"
                                                      << new_group_id_tmp << "\t"
                                                      << (*lg2_homocluster_id)[ ss3.str() ]           << "\t"
                                                  //  << old_group_id_tmp << "\t"
                                                      << best_ctg_id      << "\t"
                                                      << best_ctg_cor     << "\t"
                                                      << best_this_new_group_id                       << "\t"
                                                      << (*gb_group_id_rev)[ best_this_new_group_id ] << "\t"
                                                      << "Corr_to_LG"     << "\t"
                                                      << "case_4.4"       << endl;  
                                            ofp_refine << "\t\t\t updating wi="                 << wi << endl;   
                                            // update grouping info
	                                    ((*witr).second).type[upi]  = this_type;    // type of current window marker: hap/dip/trip/tetrap...
	                                    ((*witr).second).qlg[upi]   = new_group_id_tmp; // LG of query contig = 1:48
	                                    ((*witr).second).qchr[upi]  = (*lg2_homocluster_id)[ ss3.str() ]; // Chr of query contig = 1:12, i.e., homocluster_id
	                                    ((*witr).second).sctg[upi]  = best_ctg_id;  // id of linked subject contig
  	                                    ((*witr).second).slg[upi]   = best_this_new_group_id; // LG of linked subject contig = 1:48
    	                                    ((*witr).second).qscor[upi] = best_ctg_cor; // correlation value linked query and subject contigs
                                            ((*witr).second).qscas[upi] = "case_4.4";   // case info: how this linking is found 
                                            // the 0/1/2nd linkage group info
                                            upi ++; 
                                            hcitr_tmp ++;
                                        }
                                    }else
                                    {
                                        ; // should never happen!              
                                    }
///////////////////////////////////////////////////// re-linking ends //////////////////////////////////////////////////                                    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                                                          
                                    // relinking ends
                                }
                            }
        	            // next window marker    
                            witr ++;
                        }                
                    }
                    // next grouping type 
                    gtitr ++;
                // end of while(gtitr != gtitr_end)
                }
            //  end of if(1)
            }                                
        }else
        {
            ofp_refine << "\tCase 3: re-linking not needed. " << endl;
        }
        // next contig 
        citr ++;
    }
    ofp_refine.close();
    // output finally refined window markers
    string ofilenamefinal = tmpfolder_s4 + "/s4p6_refine_grouping_final_window_markers.txt";
    ofstream ofp_win_final;
    ofp_win_final.open(ofilenamefinal.c_str(), ios::out);
    if(!ofp_win_final.good())
    {
        cout << "   Error: cannot open file " << ofilenamefinal << endl;
        return false;
    }   
    // map<string, map<string, map<string, lgGROUP> > >::iterator citr;
    // map<string, map<string, map<string, lgGROUP> > >::iterator citr_end;
    citr     = (*raw_grouping).begin();
    citr_end = (*raw_grouping).end();
    while(citr != citr_end)
    {
        string this_ctg_id = (*citr).first;
        //
        map<string, map<string, lgGROUP> > tmp_grouping;    // <grouping, tmp_win>
        map<string, map<string, lgGROUP> >::iterator gtitr; // <"linked/unlinked_*", <win_sta_end, lgGROUP > >
        map<string, map<string, lgGROUP> >::iterator gtitr_end;
        gtitr     = (*citr).second.begin();
        gtitr_end = (*citr).second.end();
        while(gtitr != gtitr_end)
        {
            string link_unlinked = (*gtitr).first;
            //
            map<string, lgGROUP>::iterator witr;
            map<string, lgGROUP>::iterator witr_end;
            witr     = (*gtitr).second.begin() ;
            witr_end = (*gtitr).second.end();
            while(witr != witr_end)
            {
                string win_sta_end = (*witr).first;
                string this_key    = this_ctg_id + "\t" + win_sta_end;	                
	        vector<string>  type = ((*witr).second).type;  // type of current window marker: hap/dip/trip/tetrap...
	        vector<int>      qlg = ((*witr).second).qlg;   // LG of query contig = 1:48
	        vector<int>     qchr = ((*witr).second).qchr;  // Chr of query contig = 1:12, i.e., homocluster_id
	        vector<string>  sctg = ((*witr).second).sctg;  // id of linked subject contig
	        vector<int>      slg = ((*witr).second).slg;   // LG of linked subject contig = 1:48
 	        vector<double> qscor = ((*witr).second).qscor; // correlation value linked query and subject contigs
                vector<string> qscas = ((*witr).second).qscas; // case info: how this linking is found                                                
                for(int wi=0; wi<type.size(); wi++)
                {                    
                    ofp_win_final << this_key        << "\t"
                                  << type[wi]        << "\t"
                                  << qlg[wi]         << "\t"
                                  << qchr[wi]        << "\t"
                                  << sctg[wi]        << "\t"
                                  << qscor[wi]       << "\t"
                                  << slg[wi]         << "\t"
                                  << (*gb_group_id_rev)[ qlg[wi] ] << "\t"
                                  << link_unlinked   << "\t"
                                  << qscas[wi]       << endl;
                }
                // next window:sta\tend
                witr ++;            
            }
            // next link/unlinked_*
            gtitr ++;        
        }
        // next contig
        citr ++;
    }    
    //
    ofp_win_final.close(); 
    //
    return true;
}                
//
double calculate_correlation_triplotig(vector<int>* trip_marker_x, vector<int>* hap_marker_y)
{
    /*
        Theoretically, depth values at a triplotig is a sum of its three neighbouring homologous haplotig markers; 
              [there can be left/right multiple solutions]
              haplotig ............a......................: 111000111
              haplotig ............b......................: 000111010
              haplotig ............c......................: 001011101
              triplotig collapsing ab neighbouring regions: 112122222
              Then, as all homolgous hapltigs make........: 222222222, 222222222 - 112122222 => the complement
              haplotig ............d......................: 110100000
           Thus, cor(2...2 - triplotig_ab, haplotig_d) should reach 1.0 in an ideal case!
        trip_marker_x: made up by {1,2}; 
        hap_marker_y: made up by {0,1}
        this calculation is based on double calculate_correlation(vector<int>* marker_x, vector<int>* marker_y)
    */
    double trip_cor_xy = 0.0;    
    vector<int> dand;
    //
    for(int d=0; d<(*trip_marker_x).size(); d++)
    {
        int dtmp = 0;
        if((*trip_marker_x)[d]>0)
        {
            dtmp = 2 - (*trip_marker_x)[d];
            assert(dtmp>=0);
        }else
        {
            // if '0' occurs on triplotig, meaning missing values -- halfly ignore such cases = no change!
            dtmp = 0; 
        }
        dand.push_back(dtmp);
    }
    //        
    trip_cor_xy = calculate_correlation(&dand, hap_marker_y);
    //
    return trip_cor_xy;
}
//
double calculate_correlation_diplotig(vector<int>* dip_marker_x, vector<int>* hap_marker_y)
{
    /*
        Theoretically, depth values at a diplotig is a sum of its two neighbouring homologous haplotig markers; 
              [there can be left/right multiple solutions]
              haplotig ............a.....................: 100110010
              haplotig ............b.....................: 010101010
              diplotig collapsing ab neighbouring regions: 110211020
           Thus, cor(diplotig_ab & haplotig_a/b, haplotig_a/b) should reach 1.0 in an ideal case!
        dip_marker_x: made up by {0,1,2}
        hap_marker_y: made up by {0,1}
        this calculation is based on double calculate_correlation(vector<int>* marker_x, vector<int>* marker_y)
    */
    double dip_cor_xy = 0.0;    
    vector<int> dand;
    if(!get_depth_AND(dip_marker_x, 
                      hap_marker_y, 
                      &dand))
    {
        // dand <- dip_marker_x;
        dand.assign( (*dip_marker_x).begin(), (*dip_marker_x).end() );
        cout << "   Warning: cannot find AND values of two depth genotypes; diplotig values used! " << endl;
    }        
    dip_cor_xy = calculate_correlation(&dand, hap_marker_y);
    //
    return dip_cor_xy;
}
//
bool get_depth_AND(vector<int>* marker_x, vector<int>* marker_y, vector<int>* dand)
{
    /*
       dand: when both depth > 0, it is 1; otherwise it is 0
    */
    if( (*marker_x).size() != (*marker_y).size() )
    {
        cout << "   Error: unequal length of depth genotypes. " << endl;
        return false;
    }
    //
    for(int d=0; d<(*marker_x).size(); d++)
    {
        int dtmp = 0;
        if((*marker_x)[d]>0 && (*marker_y)[d]>0)
        {
            dtmp = 1;
        }
        (*dand).push_back(dtmp);
    }
    //
    return true;
}
//
bool read_gamete_binning_lg(string             gb_group_file, 
                  map<string, int>*            ctg2_gb_group,
                  map<int, map<string, int> >* gb_group_2ctg)
{
    /*
        ctg2_gb_group........: <ctg_id, old_group_id>
        gb_group_2ctg........: <old_group_id, <ctg_id, new_group_id> >
    */
    ifstream ifp;
    ifp.open(gb_group_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open gb group file " << gb_group_file << endl;
        return false;
    }
    int gb_grpid_new = 0;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        // "\tutg000001l_pilon -- utg000070l_pilon [color=gold, penwidth=1, arrowsize=1, label=0.969478]; /* cluster 0 */"
        // cout << "   check: " << line << endl;
        if(line.find(" -- ") != std::string::npos)
        {
            vector<string> lineinfo2 = split_string(line, '\t');
            vector<string> lineinfo  = split_string(lineinfo2[0], ' ');            
            string vet1  = lineinfo[0];
            string vet2  = lineinfo[2];            
            int gb_grpid_old = atoi(lineinfo[9].c_str()); // old_group_id   
            //cout << "   check: vet1=" << vet1 << " -- vet2=" << vet2 << ", gb group = " << gb_grpid_old << endl;
            if(!add_vet_to_group(vet1, 
                                 gb_grpid_old, 
                                 gb_grpid_new, 
                                 gb_group_2ctg, 
                                 ctg2_gb_group) )
            {
                return false;
            }
            if(!add_vet_to_group(vet2, 
                                 gb_grpid_old, 
                                 gb_grpid_new, 
                                 gb_group_2ctg, 
                                 ctg2_gb_group) )
            {
                return false;
            }            
        }        
    }
    ifp.close();
    //
    return true;
}
// function: read homologous "linkage groups"; in tetraploid potato, every four make one cluster!
bool read_gb_homo_lg(string homo_lg_file,
                     map<int, map<string, int> >* cluster_id_2lg,
                     map<string, int>*            lg2_cluster_id)
{
    /*
        cluster_id_2lg: <cluster_id, <lgx, cluster_id> >
        lg2_cluster_id: <lgx, cluster_id>
    */
    ifstream ifp;
    ifp.open(homo_lg_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "    Error: cannot open file " << homo_lg_file << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);	
        if(line.size()==0 || line.find("//")!=std::string::npos || line[0]=='#') continue;
        if(line.find("--") == std::string::npos) continue;
	//\t38 -- 1 [color=gold, fontcolor=gold, penwidth=1, label="-"]; /* cluster 1 */
 	//\t11 -- 4 [color=gold, fontcolor=gold, penwidth=1, label="-"]; /* cluster 2 */
	vector<string> lineinfo2 = split_string(line, '\t');
	vector<string> lineinfo  = split_string(lineinfo2[0], ' ');
	
	int cluster_id = atoi(lineinfo[9].c_str() );	
	if( (*cluster_id_2lg).find(cluster_id) == (*cluster_id_2lg).end() )
	{
	    map<string, int> tmpcluster;
	    tmpcluster.insert(std::pair<string, int>(lineinfo[0], 1));
	    tmpcluster.insert(std::pair<string, int>(lineinfo[2], 1));
	    (*cluster_id_2lg).insert(std::pair<int, map<string, int> >(cluster_id, tmpcluster));
	}else
	{
	    if( (*cluster_id_2lg)[cluster_id].find(lineinfo[0]) == (*cluster_id_2lg)[cluster_id].end() )
	    {
	        (*cluster_id_2lg)[cluster_id].insert(std::pair<string, int>(lineinfo[0], 1));
	    }
	    if( (*cluster_id_2lg)[cluster_id].find(lineinfo[2]) == (*cluster_id_2lg)[cluster_id].end() )
	    {
	        (*cluster_id_2lg)[cluster_id].insert(std::pair<string, int>(lineinfo[2], 1));
	    }	    
	}
	//
	if( (*lg2_cluster_id).find(lineinfo[0]) == (*lg2_cluster_id).end() )
	{
	    (*lg2_cluster_id).insert(std::pair<string, int>(lineinfo[0], cluster_id));
	}else
	{
	    if( (*lg2_cluster_id)[lineinfo[0] ] != cluster_id)
	    {
	        cout << "   Error: linkage group " << lineinfo[0] << "  found in different homologous clusters. " << endl;
	        return false;
	    }
	}
	if( (*lg2_cluster_id).find(lineinfo[2]) == (*lg2_cluster_id).end() )
	{
	    (*lg2_cluster_id).insert(std::pair<string, int>(lineinfo[2], cluster_id));
	}else
	{
	    if( (*lg2_cluster_id)[lineinfo[2] ] != cluster_id)
	    {
	        cout << "   Error: linkage group " << lineinfo[2] << "  found in different homologous clusters. " << endl;
	        return false;
	    }	
	}	
    }
    ifp.close();
    //
    bool check_output = true;
    if(check_output)
    {
        map<int, map<string, int> >::iterator hmitr;
        map<int, map<string, int> >::iterator hmitr_end;
        hmitr     = (*cluster_id_2lg).begin();
        hmitr_end = (*cluster_id_2lg).end();
        while(hmitr != hmitr_end)
        {
            cout << "   check: homologous cluster " << (*hmitr).first << " including linkage groups: " << endl;
            map<string, int> hcluster = (*hmitr).second;
            map<string, int>::iterator lgitr;
            map<string, int>::iterator lgitr_end;
            lgitr     = hcluster.begin();
            lgitr_end = hcluster.end();
            while(lgitr != lgitr_end)
            {
                cout << "       " << (*lgitr).first << endl;
                lgitr ++;
            }
            //
            hmitr ++;
        }
        //
        map<string, int>::iterator lgitr;
        map<string, int>::iterator lgitr_end;
        lgitr     = (*lg2_cluster_id).begin();
        lgitr_end = (*lg2_cluster_id).end();
        while(lgitr != lgitr_end)
        {
            cout << "    check: linkage group " << (*lgitr).first 
                 << " in homologous cluster "   << (*lgitr).second << endl;
            lgitr ++;
        }
    }
    //
    return true;
}
//
bool add_vet_to_group(string                       vet, 
                      int                          gb_grpid_old,
                      int                          gb_grpid_new,
                      map<int, map<string, int> >* gb_group_2ctg,
                      map<string, int>*            ctg2_gb_group)
{
    /*
        gb_grpid_old.........: old group id: any number
        gb_grpid_new.........: new group id: 1-48
        gb_group_2ctg........: <old_group_id, <ctg_id, new_group_id> >
    */    
    // update map from linkage group to list of contigs 
    if( (*gb_group_2ctg).find(gb_grpid_old) == (*gb_group_2ctg).end() )
    {
        map<string, int> tmpgroup;
        tmpgroup.insert(std::pair<string, int>(vet, gb_grpid_new));
        (*gb_group_2ctg).insert(std::pair<int, map<string, int> >(gb_grpid_old, tmpgroup));
    }else
    {
        if( (*gb_group_2ctg)[gb_grpid_old].find(vet) == (*gb_group_2ctg)[gb_grpid_old].end() )
        {        
            (*gb_group_2ctg)[gb_grpid_old].insert(std::pair<string, int>(vet, gb_grpid_new));
        }else
        {
            ; // redundant 
        }
    }
    // update map from a contig to its linkage group
    if((*ctg2_gb_group).find(vet) == (*ctg2_gb_group).end())
    {
        (*ctg2_gb_group).insert(std::pair<string, int>(vet, gb_grpid_old));
    }else
    {
        // should never happen
        if( (*ctg2_gb_group)[vet] != gb_grpid_old )
        {
            cout << "   Warning: contig found in group " << (*ctg2_gb_group)[vet] << " and " << gb_grpid_old << endl;
            cout << "            if it is non-haplotig markers, it is fine; otherwise, something wrong!"     << endl;
        }
    }
    //      
    return true;
}
//
bool greedy_cluster_haplotigs(map<string, BIPAT_INT>       contigPollenGTSeq_int, 
                              map<string, unsigned long>   contigsize, 
                              string                       tmpfolder,
                              map<int, int>*               gb_group_id,
                              map<int, map<string, int> >* gb_group)
{
    /*        
        cluster the haplotig markers into linkage groups using a greedy way:
            here only correlation > 0.6 will generate an edge!
        contigPollenGTSeq_int: <ctg_id, {leftPat, rightPat}> .... a list of pure-haplotig markers: hap_win_ratio>0.8           
        //
        // later visualization:
        circo -Tpdf -s108 xxx.dot > xxx.circo.pdf
    */
    map<string, int> vertex_all;           // <contig, 1>
    map<string, vector<string> > edge_all; // <contig1, <contigs-connecting-contig1> >
    map<string, double> cor_all;           // <contig1-contig2, correlation.value >
    // prepare dot output file
    string siminfo = tmpfolder + "/s3_genotype_haplotig_GT_similarity_matrix.dot\0"; // 
    ofstream sofp;
    sofp.open(siminfo.c_str(), ios::out);
    if(!sofp.good())
    {
        cout   << "   Error: cannot open file " << siminfo << endl;
        return false;
    }
    // un-directed graph
    sofp << "/* Contact graph by correlation of genotypes defined by sequencing depth */" << endl;
    sofp << "graph\tGraph_1 {" << endl;
    //
    map<string, BIPAT_INT>::iterator citr1;
    map<string, BIPAT_INT>::iterator citr1_end;
    citr1     = contigPollenGTSeq_int.begin();
    citr1_end = contigPollenGTSeq_int.end();
    while(citr1 != citr1_end)
    {
        BIPAT_INT pat1;
        pat1.leftPat_int  = (*citr1).second.leftPat_int;
        pat1.rightPat_int = (*citr1).second.rightPat_int;        
        map<string, BIPAT_INT>::iterator citr2;
        map<string, BIPAT_INT>::iterator citr2_end;
        citr2     = citr1;
        citr2     ++;
        citr2_end = contigPollenGTSeq_int.end();
        while(citr2 != citr2_end)
        {
            BIPAT_INT pat2;
            pat2.leftPat_int  = (*citr2).second.leftPat_int;
            pat2.rightPat_int = (*citr2).second.rightPat_int;            
            double cor[4];
            cor[0] = calculate_correlation(&pat1.leftPat_int, &pat2.leftPat_int);      
            cor[1] = calculate_correlation(&pat1.leftPat_int, &pat2.rightPat_int);          
            cor[2] = calculate_correlation(&pat1.rightPat_int, &pat2.leftPat_int);      
            cor[3] = calculate_correlation(&pat1.rightPat_int, &pat2.rightPat_int);            
            //
            double max_cor = cor[0];
            double max_i   = 0;
            for(int i = 1; i < 4; i ++)
            {
                if(cor[i]>max_cor)
                {
                    max_cor = cor[i];
                    max_i   = i;
                }
            }
            //
            if(max_cor > minCorscore) // to set up: default 0.6: here we have another level of haplotig marker selection
            {
                string vet1  = (*citr1).first;
                string vet2  = (*citr2).first;
                string edge1 = vet1 + "-" + vet2;
                string edge2 = vet2 + "-" + vet1;
                if(max_i == 0)
                {
                    sofp << "\t"
                         << (*citr1).first << " -- " << (*citr2).first << " "            
                         << "[color=gold, penwidth=1, arrowsize=1, label=" << max_cor << "];" 
                         << " /* RLLR */"  << endl;                       
                }else
                if(max_i == 1)
                {
                    sofp << "\t"
                         << (*citr1).first << " -- " << (*citr2).first << " "            
                         << "[color=gold, penwidth=1, arrowsize=1, label=" << max_cor << "];" 
                         << " /* RLRL */"  << endl;                             
                }else 
                if(max_i == 2)
                {
                    sofp << "\t"
                         << (*citr1).first << " -- " << (*citr2).first << " "            
                         << "[color=gold, penwidth=1, arrowsize=1, label=" << max_cor << "];" 
                         << " /* LRLR */"  << endl;                             
                }else        
                if(max_i == 3)
                {
                    sofp << "\t"
                         << (*citr1).first << " -- " << (*citr2).first << " "            
                         << "[color=gold, penwidth=1, arrowsize=1, label=" << max_cor << "];" 
                         << " /* LRRL */"  << endl;                             
                }else  ; 
                // collect vertex <contig, 1>
                if( vertex_all.find(vet1) == vertex_all.end() )
                {
                    vertex_all.insert(std::pair<string, int>(vet1, 1));
                }
                if( vertex_all.find(vet2) == vertex_all.end() )
                {
                    vertex_all.insert(std::pair<string, int>(vet2, 1));
                }  
                // collect edge <contig, <contig> >
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
                // collect correlation: <contig1-contig2, cor>             
                if(cor_all.find(edge1) == cor_all.end() )
                {
                    cor_all.insert(std::pair<string, double>(edge1, max_cor) );
                }else
                {
                    // should never happen
                    cout << "   Warning: edge " << edge1 << " seen more than 1 time. " << endl;
                }
            }
            //
            citr2 ++;                
        }                
        citr1 ++;
    }
    sofp << "}" << endl; 
    //
    sofp.close();       
    //
    cout << "   Info: " << vertex_all.size() << " vertices collected with "
         <<                cor_all.size()    << " edges. " 
         << endl;
    // get sub-clusters = backbone of potential linkage groups 
    map<int, map<string, double> > cluster_edge; 
    if( !get_clusters(vertex_all, edge_all, cor_all, &cluster_edge) )
    {
        return false;
    }
    // output subclusters: this is raw, can be more than expected number of linkage groups thus needs further merging.
    string subclusterinfo = tmpfolder + "/s3_genotype_haplotig_GT_similarity_matrix_subclusters_raw.dot\0"; // 
    ofstream subcluster_ofp;
    subcluster_ofp.open(subclusterinfo.c_str(), ios::out);
    if(!subcluster_ofp.good())
    {
        cout   << "   Error: cannot open file " << subclusterinfo << endl;
        return false;
    } 
    subcluster_ofp << "/* Here are the subclusters of contigs */" << endl;
    subcluster_ofp << "graph\tGraph_1 {"                           << endl;
    // to collect all clusters-specific vertex
    map<int, map<string, int> > cluster_vertex;      // <clusterid, <contig, 1> >    
    map<int, unsigned long>     cluster_vertex_size; // <clusterid, total_ctg_length >        
    map<int, map<string, double> >::iterator clitr;
    map<int, map<string, double> >::iterator clitr_end;
    clitr     = cluster_edge.begin();
    clitr_end = cluster_edge.end();
    unsigned long total_contig_size = 0;
    unsigned long total_contig_numb = 0;
    while(clitr != clitr_end)
    {
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
            << endl;  
            //
            if(subcluster_vertex.find(edgeinfo[0]) == subcluster_vertex.end() )
            {
                subcluster_vertex.insert(std::pair<string, int>(edgeinfo[0], 1));
                assert(contigsize.find(edgeinfo[0]) != contigsize.end());
                cluster_ctg_size += contigsize[edgeinfo[0]];
            }
            if(subcluster_vertex.find(edgeinfo[1]) == subcluster_vertex.end() )
            {
                subcluster_vertex.insert(std::pair<string, int>(edgeinfo[1], 1));
                assert(contigsize.find(edgeinfo[1]) != contigsize.end());
                cluster_ctg_size += contigsize[edgeinfo[1]];
            }            
            //           
            eitr ++;
        }
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
    cout << "   Info: "         << total_contig_numb << " contigs of " << total_contig_size << " bp clustered into "
         << cluster_edge.size() << " clusters. "     << endl;
    // analyze cluster; 
    //    all members of a good cluster show ">0" correlation values
    //    many "cor<-0.1" means mis-joining, thus need breaking
    map<int, map<string, double> > cluster_edge_updated;
    map<int, unsigned long>        cluster_vertex_size_updated;
    map<int, map<string, int> >    cluster_vertex_updated;   
    map<int, int>                  mutually_exclusive_cluster;
    if(!break_clusters(cluster_vertex,
                       contigPollenGTSeq_int,
                       tmpfolder,
                       contigsize,
                       &cluster_edge_updated,
                       &cluster_vertex_size_updated,
                       &cluster_vertex_updated,
                       &mutually_exclusive_cluster))
    {
        cout << "   Error: failed in breaking clusters. " << endl;
        return false;
    } 
    //
    map<int, map<string, double> >    cluster_edge_updated2;
    map<int, unsigned long>           cluster_vertex_size_updated2;
    map<int, map<string, int> >       cluster_vertex_updated2;
    //map<int, int>                   mutually_exclusive_cluster;
    if(!calc_inter_clusters_cor(cluster_vertex_updated,
                                cluster_vertex_size_updated,
                                contigPollenGTSeq_int,
                                tmpfolder,
                                contigsize,
                                &cluster_edge_updated2,
                                &cluster_vertex_size_updated2,
                                &cluster_vertex_updated2,
                                &mutually_exclusive_cluster))
    {
        cout << "   Error: failed in calculating inter-cluster correlations. " << endl;
        return false;
    }
    // find LGs by merging raw-advanced clusters according to expected number: 
    //    chromosome-number X ploidy, in potato 12 x 4 = 48. 
    // build map<int, int>*               gb_group_id; // <old_group_id, new_group_id>
    // build map<int, map<string, int> >* gb_group;    // <old_group_id, <ctg_id, 1> >
    if( !merge_clusters(cluster_edge_updated2, 
                        cluster_vertex_updated2, 
                        cluster_vertex_size_updated2,
                        mutually_exclusive_cluster,
                        contigPollenGTSeq_int, 
                        contigsize,
                        N_cluster,
                        tmpfolder,
                        gb_group_id,
                        gb_group) )
    {
        cout << "   Error: merging cluser was not finished as expected. " << endl;
        return false;
    }
    // find "homologous haplotype-specific linkage groups", i.e., "ploidy" LGs of the same chromosome.
    // TODO: hom_cor_cutoff to setup input option
    if(!find_hom_group(contigPollenGTSeq_int, 
                       *gb_group, 
                       *gb_group_id, 
                       hom_cor_cutoff,
                       tmpfolder) )
    {
        return false;
    }    
    //
    return true;
}
//
//
bool get_lg_clusters(map<string, vector<string> > edge,
                  map<string, double>             hom_pairs_simple,
                  map<int, map<string, double> >* cluster_edge,
                  map<int, map<string, int> >*    cluster_vertex)
{
    // function: traverse LG1-LG2 graph, connected LGs will be collected as one cluster, 
    //           i.e., homologous LGs = LGs refer to the same chromosome.
    //
    // vertex, edge and correlation
    // map<string, vector<string> > edge;           // <group_id, <connected-groups> >
    // map<string, double> hom_pairs_simple;        // <group1-group2, correlation.value >
    // map<int, map<string, double> > cluster_edge; // <cluster_id, <group1-group2, correlation.value> >
    // map<int, map<string, int> > cluster_vertex;  // <cluster_id, <group_id, 1> >
    int cluster_id = -1;
    bool check_on  = false;
    map<string, int> visited_vet;
    while(hom_pairs_simple.size() > 0)
    {
        // get an <group1-group2, cor>
        map<string, double>::iterator citr = hom_pairs_simple.begin();
        string         this_edge    = (*citr).first; // group1-group2
        vector<string> this_vetinfo = split_string(this_edge, '-');
        // initialize a cluster id
        cluster_id ++;
        if(check_on) cout << "   check: building a cluster of homolougous linkage groups " << cluster_id << endl;
        // initialize edge of this cluster
        map<string, double> tmpcor; // <group1-group2, cor>
        tmpcor.insert(std::pair<string, double>( (*citr).first, (*citr).second) );
        (*cluster_edge).insert(std::pair<int, map<string, double> >(cluster_id, tmpcor) );
        // initialize vertex of this cluster
        map<string, int> tmpvet;
        tmpvet.insert(std::pair<string, int>(this_vetinfo[0], 1));
        tmpvet.insert(std::pair<string, int>(this_vetinfo[1], 1));
        (*cluster_vertex).insert(std::pair<int, map<string, int> >(cluster_id, tmpvet));
        // collect visited vertex
        visited_vet.insert(std::pair<string, int>(this_vetinfo[0], 1));
        visited_vet.insert(std::pair<string, int>(this_vetinfo[1], 1));
        // remove this <group1-group2, cor>
        hom_pairs_simple.erase(citr);
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
                    if((*cluster_vertex)[cluster_id].find(*vitr) == (*cluster_vertex)[cluster_id].end())
                    {
                        (*cluster_vertex)[cluster_id].insert(std::pair<string, int>(*vitr, 1));
                    }
                    map<string, int>::iterator existence_itr = visited_vet.find(*vitr);
                    if(existence_itr == visited_vet.end())
                    {
                        vet_queue.push_back(*vitr);
                    }
                    // update edge with cor
                    string tmp_edge1 = this_vet + "-" + *vitr; 
                    string tmp_edge2 = *vitr    + "-" + this_vet;                                         
                    citr = hom_pairs_simple.find(tmp_edge1);
                    if(citr == hom_pairs_simple.end())
                    {
                        citr = hom_pairs_simple.find(tmp_edge2);
                    }
                    if(citr != hom_pairs_simple.end() )
                    {
                        string         tmp_edge    = (*citr).first; // contig1-contig2
                        vector<string> tmp_vetinfo = split_string(tmp_edge, '-');
                        // update edge
                        (*cluster_edge)[cluster_id].insert(std::pair<string, double>( (*citr).first, (*citr).second) );
                        // update visited vertex 
                        visited_vet.insert(std::pair<string, int>(tmp_vetinfo[0], 1));
                        visited_vet.insert(std::pair<string, int>(tmp_vetinfo[1], 1));                        
                        // remove this <contig1-contig2, cor>
                        hom_pairs_simple.erase(citr);
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
bool find_hom_group(map<string, BIPAT_INT>     contigPollenGTSeq_int, 
                  map<int, map<string, int> >  gb_group,
                  map<int, int>                gb_group_id,
                  double                       hom_cor_cutoff,
                  string                       tmpfolder)
{
    /*
        contigPollenGTSeq_int: <ctg_id, {leftPat, rightPat}>
        gb_group.............: <old_group_id, <ctg_id, 1> >
        gb_group_id..........: <old_group_id:to_dot, new_group_id>
    */
    map<string, double> hom_pairs;                // <group1-group2, correlation> that are defined as homologous
    map<string, double> hom_pairs_simple;         // <group1-group2, correlation> that are defined as homologous
    map<string, vector<string> > hom_edge_all;    // <group_id, <list-groups> >       
    map<string, map<string, double> > out_degree; // <group1, <group2, score> >    
    //
    map<int, map<string, int> >::iterator gbitr1;
    map<int, map<string, int> >::iterator gbitr1_end;
    gbitr1     = gb_group.begin();
    gbitr1_end = gb_group.end();
    while(gbitr1 != gbitr1_end)
    {
        int              grp1_id  = (*gbitr1).first;
        map<string, int> grp1_ctg = (*gbitr1).second;
        //
        map<int, map<string, int> >::iterator gbitr2;
        map<int, map<string, int> >::iterator gbitr2_end;        
        //gbitr2   = gbitr1;
        //gbitr2   ++;
        gbitr2     = gb_group.begin();
        gbitr2_end = gb_group.end();
        while(gbitr2 != gbitr2_end && gbitr1 != gbitr2)
        {
            int              grp2_id  = (*gbitr2).first;
            map<string, int> grp2_ctg = (*gbitr2).second;
            //
            map<string, int>::iterator grp1_citr;
            map<string, int>::iterator grp1_citr_end;
            grp1_citr     = grp1_ctg.begin();
            grp1_citr_end = grp1_ctg.end();
            while(grp1_citr != grp1_citr_end)
            {
                string this_grp1_ctg = (*grp1_citr).first;
                assert(contigPollenGTSeq_int.find(this_grp1_ctg) != contigPollenGTSeq_int.end() );
                BIPAT_INT pat1;
                pat1.leftPat_int  = contigPollenGTSeq_int[this_grp1_ctg].leftPat_int;
                pat1.rightPat_int = contigPollenGTSeq_int[this_grp1_ctg].rightPat_int;                
                //
                map<string, int>::iterator grp2_citr;
                map<string, int>::iterator grp2_citr_end;
                grp2_citr     = grp2_ctg.begin();
                grp2_citr_end = grp2_ctg.end();    
                while(grp2_citr != grp2_citr_end)
                {
                    string this_grp2_ctg = (*grp2_citr).first;   
                    assert(contigPollenGTSeq_int.find(this_grp2_ctg) != contigPollenGTSeq_int.end() );
                    BIPAT_INT pat2;
                    pat2.leftPat_int  = contigPollenGTSeq_int[this_grp2_ctg].leftPat_int;
                    pat2.rightPat_int = contigPollenGTSeq_int[this_grp2_ctg].rightPat_int;
                    // calculate correlation: negative expected
                    double cor[4];
                    cor[0] = calculate_correlation(&pat1.leftPat_int, &pat2.leftPat_int);      
                    cor[1] = calculate_correlation(&pat1.leftPat_int, &pat2.rightPat_int);          
                    cor[2] = calculate_correlation(&pat1.rightPat_int, &pat2.leftPat_int);      
                    cor[3] = calculate_correlation(&pat1.rightPat_int, &pat2.rightPat_int);  
                    double min_cor = cor[0]; // [-1, 0)
                    double min_i   = 0;
                    for(int i = 1; i < 4; i ++)
                    {
                        if(cor[i]<min_cor && cor[i]>-100)
                        {
                            min_cor = cor[i];
                            min_i   = i;
                        }
                    }
                    // selection 
                    if(min_cor < hom_cor_cutoff && min_cor>-100) // negative expected between homologuous haplotigs
                    {
                        string vet1  = this_grp1_ctg;
                        string vet2  = this_grp2_ctg;
                        string edge1 = vet1 + "-" + vet2;
                        cout << "   Info: " << gb_group_id[grp1_id] << " hom " << gb_group_id[grp2_id] << endl; 
                        std::stringstream ss;
                        ss.str("");
                        ss << gb_group_id[grp1_id] << " -- " << gb_group_id[grp2_id];
                        if(hom_pairs.find(ss.str()) == hom_pairs.end() )
                        {
                            hom_pairs.insert(std::pair<string, double>(ss.str(), min_cor));
                        }else
                        {
                            hom_pairs[ss.str()] += min_cor;
                        }
                        //
                        ss.str(""); 
                        ss << gb_group_id[grp1_id] << "-" << gb_group_id[grp2_id];                        
                        if(1)
                        {
                            if(hom_pairs_simple.find(ss.str()) == hom_pairs_simple.end() )
                            {
                                hom_pairs_simple.insert(std::pair<string, double>(ss.str(), min_cor));
                            }else
                            {
                                hom_pairs_simple[ss.str()] += min_cor;
                            }
                            // new: group-wise checking 20200731
                            std::stringstream ss1;
                            ss1.str("");
                            ss1 << gb_group_id[grp1_id];
                            std::stringstream ss2;
                            ss2.str("");
                            ss2 << gb_group_id[grp2_id];                         
                            if(out_degree.find( ss1.str() ) == out_degree.end() )
                            {
                                map<string, double> tmp_outd;
                                tmp_outd.insert(std::pair<string, double>( ss2.str(), min_cor));
                                out_degree.insert(std::pair<string, map<string, double> >( ss1.str(), tmp_outd));
                            }else
                            {
                                if( out_degree[ ss1.str() ].find( ss2.str() ) == out_degree[ ss1.str() ].end() )
                                {
                                    out_degree[ ss1.str() ].insert(std::pair<string, double>( ss2.str(), min_cor));
                                }else
                                {
                                    out_degree[ ss1.str() ][ ss2.str() ] += min_cor;
                                }
                            }
                            if(out_degree.find( ss2.str() ) == out_degree.end() )
                            {
                                map<string, double> tmp_outd;
                                tmp_outd.insert(std::pair<string, double>( ss1.str(), min_cor));
                                out_degree.insert(std::pair<string, map<string, double> >( ss2.str(), tmp_outd));
                            }else
                            {
                                if( out_degree[ ss2.str() ].find( ss1.str() ) == out_degree[ ss2.str() ].end() )
                                {
                                    out_degree[ ss2.str() ].insert(std::pair<string, double>( ss1.str(), min_cor));
                                }else
                                {
                                    out_degree[ ss2.str() ][ ss1.str() ] += min_cor;
                                }
                            } // new end 20200731                                                                                       
                            // collect edge <group, <list-groups> >
                            /*
                            std::stringstream ss1;
                            ss1.str("");
                            ss1 << gb_group_id[grp1_id];
                            std::stringstream ss2;
                            ss2.str("");
                            ss2 << gb_group_id[grp2_id];   
                            */                     
                            if( hom_edge_all.find(ss1.str()) == hom_edge_all.end() )
                            {
                                vector<string> tmpedge;
                                tmpedge.push_back(ss2.str());
                                hom_edge_all.insert(std::pair<string, vector<string> >(ss1.str(), tmpedge) );
                            }else
                            {
                                if(std::find(hom_edge_all[ ss1.str() ].begin(), 
                                             hom_edge_all[ ss1.str() ].end(),
                                             ss2.str() ) == hom_edge_all[ ss1.str() ].end() )
                                {             
                                    hom_edge_all[ ss1.str() ].push_back( ss2.str() );
                                }
                            }
                            if(hom_edge_all.find( ss2.str() ) == hom_edge_all.end() )
                            {
                                vector<string> tmpedge;
                                tmpedge.push_back( ss1.str() );
                                hom_edge_all.insert(std::pair<string, vector<string> >(ss2.str(), tmpedge) );
                            }else
                            {
                                if(std::find(hom_edge_all[ ss2.str() ].begin(), 
                                             hom_edge_all[ ss2.str() ].end(),
                                             ss1.str() ) == hom_edge_all[ ss2.str() ].end() )
                                {                               
                                    hom_edge_all[ ss2.str() ].push_back( ss1.str() );
                                }
                            }                              
                        }else
                        {
                            cout << "   Warning: need a cluster-broken-merge algorithm here!!!!!!!!!!!!!!!!!!!!!!! " 
                                 << endl;
                            cout << "            TODO-20200714: using top 25\% long contigs to find homologous lgs. " 
                                 << endl;         
                            cout << "    Caution: still need a tiny method to break connections between 12 and 33/47!"
                                 << "             check src code. " << endl;                                                
                        }                                             
                        //
                        if(min_i == 0)
                        {
                            cout << "\t"
                                 << vet1           << " -- "  << vet2    << " "            
                                 << "[color=gold, penwidth=1, arrowsize=1, label=" << min_cor << "];" 
                                 << " /* RLLR */"  << endl;                       
                        }else
                        if(min_i == 1)
                        {
                            cout << "\t"
                                 << vet1 << " -- " << vet2 << " "            
                                 << "[color=gold, penwidth=1, arrowsize=1, label=" << min_cor << "];" 
                                 << " /* RLRL */"  << endl;                             
                        }else 
                        if(min_i == 2)
                        {
                            cout << "\t"
                                 << vet1 << " -- " << vet2 << " "            
                                 << "[color=gold, penwidth=1, arrowsize=1, label=" << min_cor << "];" 
                                 << " /* LRLR */"  << endl;                             
                        }else        
                        if(min_i == 3)
                        {
                            cout << "\t"
                                 << vet1 << " -- " << vet2 << " "            
                                 << "[color=gold, penwidth=1, arrowsize=1, label=" << min_cor << "];" 
                                 << " /* LRRL */"  << endl;                             
                        }else  ; 
                    }                                                                                                  
                    //
                    grp2_citr ++;
                }                            
                //
                grp1_citr ++;
            }
            //
            gbitr2 ++;
        }
        //
        gbitr1 ++;
    }
    // begin new 20200731
    // post-processing: 1. break clusters with more than iploidy-1 linkage groups; 
    //                  2. TODO: merge clusters with less than iploidy linkage groups?
    // map<string, double> hom_pairs_simple;         // <group1-group2, correlation> that are defined as homologous
    // map<string, vector<string> > hom_edge_all;    // <group_id, <list-groups> >
    // map<string, map<string, double> > out_degree; // <group1, <group2, score> >
    map<string, map<string, double> >::iterator oditr;
    map<string, map<string, double> >::iterator oditr_end;
    oditr     = out_degree.begin();
    oditr_end = out_degree.end();
    while(oditr != oditr_end)
    {
        map<string, double> tmp_outd = (*oditr).second;
        string first_lgid = (*oditr).first;
        if(tmp_outd.size() > iploidy-1)
        {
            cout << "    Info: "                                 << first_lgid
                 << " linked with "                              << tmp_outd.size() 
                 << " LGs, thus prunning and updating hom-pairs" << endl;
            // find best three = with minimum negative cor 
            double lowest_cor[3] = {0, 0, 0};
            vector<string> best_hom_lg;
            best_hom_lg.push_back("n");
            best_hom_lg.push_back("n");
            best_hom_lg.push_back("n"); 
            map<string, int> edge_to_be_removed;
            //                       
            map<string, double>::iterator titr;
            map<string, double>::iterator titr_end;
            titr     = tmp_outd.begin();
            titr_end = tmp_outd.end();
            while(titr != titr_end)
            {
                string second_lgid = (*titr).first;
                //
                if((*titr).second < lowest_cor[0])
                {
                    // if there is something in 3rd; taking it into removed class
                    if(best_hom_lg[2].compare("n") != 0)
                    update_hom_into(first_lgid, 
                                    best_hom_lg[2], 
                                    &hom_pairs_simple,
                                    &hom_pairs,                                    
                                    &hom_edge_all);                    
                    // move 2nd to 3rd
                    lowest_cor[2]  = lowest_cor[1];
                    best_hom_lg[2] = best_hom_lg[1];                        
                    // move 1st to 2nd
                    lowest_cor[1]  = lowest_cor[0];
                    best_hom_lg[1] = best_hom_lg[0];                          
                    // update 1st
                    lowest_cor[0]  = (*titr).second;
                    best_hom_lg[0] = second_lgid;
                }else
                if((*titr).second < lowest_cor[1])
                {
                    // if there is something in 3rd; taking it into removed class
                    if(best_hom_lg[2].compare("n") != 0)                    
                    update_hom_into(first_lgid, 
                                    best_hom_lg[2], 
                                    &hom_pairs_simple,
                                    &hom_pairs,                                    
                                    &hom_edge_all);                                         
                    // move 2nd to 3rd
                    lowest_cor[2]  = lowest_cor[1];
                    best_hom_lg[2] = best_hom_lg[1];                    
                    // update 2nd
                    lowest_cor[1]  = (*titr).second;
                    best_hom_lg[1] = second_lgid;
                }else
                if((*titr).second < lowest_cor[2])
                {
                    // if there is something in 3rd; taking it into removed class
                    if(best_hom_lg[2].compare("n") != 0)                    
                    update_hom_into(first_lgid, 
                                    best_hom_lg[2], 
                                    &hom_pairs_simple,
                                    &hom_pairs,                                    
                                    &hom_edge_all);                                        
                    // update 3rd
                    lowest_cor[2]  = (*titr).second;
                    best_hom_lg[2] = second_lgid;
                }else 
                {   
                    // collect group1-group2 edge to be removed.
                    edge_to_be_removed;
                    string ssid = first_lgid + "-" + second_lgid;
                    if(hom_pairs_simple.find(ssid) != hom_pairs_simple.begin() )
                    {
                        update_hom_into(first_lgid, 
                                        second_lgid, 
                                        &hom_pairs_simple,
                                        &hom_pairs,                                        
                                        &hom_edge_all);                                                                                          
                    }
                    ssid = second_lgid + "-" + first_lgid; 
                    if(hom_pairs_simple.find(ssid) != hom_pairs_simple.begin() )
                    {
                        update_hom_into(first_lgid, 
                                        second_lgid, 
                                        &hom_pairs_simple,
                                        &hom_pairs,
                                        &hom_edge_all);
                    }
                }
                //
                titr ++;
            }
        }
        //
        oditr ++;
    }    
    // end   new 20200731
    // get cluster of linkage groups
    map<int, map<string, double> > cluster_edge;   // <cluster_id, <group1-group2, cor> >
    map<int, map<string, int> >    cluster_vertex; // <cluster_id, <groupx, 1> >
    if(!get_lg_clusters(hom_edge_all, hom_pairs_simple, &cluster_edge, &cluster_vertex))
    {
        cout << "   Error: cannot find cluster of homologous linkage groups. " << endl;
        return false;
    }    
    // output homologuous groups as a graph: this has been integrated into gamete binning
    string ohomfilename = tmpfolder + "/s3_res_homologous_linkgage_groups.dot";
    ofstream ofp;
    ofp.open(ohomfilename.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open output file " << ohomfilename << endl;
        return false;
    }
    ofp << "/* Here are the homolougous linkage groups */" << endl;
    ofp << "graph\tGraph_1 {"                           << endl;    
    //
    map<int, map<string, double> >::iterator hcitr;
    map<int, map<string, double> >::iterator hcitr_end;
    hcitr     = cluster_edge.begin();
    hcitr_end = cluster_edge.end();
    while(hcitr != hcitr_end)
    {
        map<string, double>::iterator homitr;
        map<string, double>::iterator homitr_end;
        homitr     = ((*hcitr).second).begin();
        homitr_end = ((*hcitr).second).end();
        ofp << "\tsubgraph cluster_" << (*hcitr).first+1 << " {" << endl;        
        while(homitr != homitr_end)
        {
            vector<string> edgeinfo = split_string((*homitr).first, '-');
            
            ofp << "\t" << edgeinfo[0] << " -- " << edgeinfo[1] << " " 
                    << "[color=gold, fontcolor=gold, penwidth=1, label=\"-\"];" 
                    << " /* cluster " << (*hcitr).first+1 << " */" << endl;
            //
            homitr ++;
        }
        // sub cluster id: homologous linkage groups = hLG
        ofp << "\tlabel=\"hLGs_" << (*hcitr).first+1 << "\";" << endl;
        ofp << "\tfontsize=90;"            << endl; // Make title stand out by giving a large font size
        ofp << "\tfontcolor=orangered;"    << endl;
        ofp << "\tcolor=gray;"             << endl;        
        ofp << "\t}" << endl;
        //
        hcitr ++;
    }
    // output node/vertex as group id
    for(int ni=1; ni<=gb_group.size(); ni ++)
    {
        ofp << "\t" 
            << ni
            << " [color=orangered, style=filled, fillcolor=orangered, fontcolor=white];"
            << endl;    
    }
    //
    ofp << "}" << endl;
    ofp.close();
    cout << "   Info: to visualize the graph\n      circo -Tpdf -s108 xxx.dot > xxx.circo.pdf" << endl;
    //
    return true;
}
// new 20200731
bool update_hom_into(string first_lgid, 
                     string second_lgid, 
                     map<string, double>* hom_pairs_simple,
                     map<string, double>* hom_pairs,                     
                     map<string, vector<string> >* edge_all)
{
    // first_lgid/second_lgid                  // id of linkage groups showing initial hom-relationships
    // map<string, double>* hom_pairs_simple;  // <group1-group2, correlation> that are defined as homologous
    // map<string, vector<string> >* edge_all; // <group_id, <list-groups> >   
    //   
    string ssid;
    //
    ssid = first_lgid + "-" + second_lgid;    
    if((*hom_pairs_simple).find(ssid) != (*hom_pairs_simple).end() )
    {
        (*hom_pairs_simple).erase(ssid);
        cout << "   info: " << ssid << " removed from hom_pairs_simple. " << endl;         
    }
    ssid = first_lgid + " -- " + second_lgid;  
    if( (*hom_pairs).find(ssid) != (*hom_pairs).end() )
    {
        (*hom_pairs).erase(ssid);
        cout << "   info: " << ssid << " removed from hom_pairs. "        << endl;   
    }  
    //
    ssid = second_lgid + "-" + first_lgid;    
    if((*hom_pairs_simple).find(ssid) != (*hom_pairs_simple).end() )
    {
        (*hom_pairs_simple).erase(ssid);
        cout << "   info: " << ssid << " removed from hom_pairs_simple. " << endl;         
    }
    ssid = second_lgid + " -- " + first_lgid;
    if( (*hom_pairs).find(ssid) != (*hom_pairs).end() )
    {
        (*hom_pairs).erase(ssid);
        cout << "   info: " << ssid << " removed from hom_pairs. "        << endl;   
    }      
    //
    vector<string>::iterator tvitr;
    tvitr = std::find((*edge_all)[first_lgid].begin(), 
                      (*edge_all)[first_lgid].end(),
                      second_lgid);
    if(tvitr != (*edge_all)[second_lgid].end() )                                                                                                    
    (*edge_all)[first_lgid].erase(tvitr);   
    //   
    tvitr = std::find((*edge_all)[second_lgid].begin(), 
                      (*edge_all)[second_lgid].end(),
                      first_lgid);     
    if(tvitr != (*edge_all)[second_lgid].end() )                                                          
    (*edge_all)[second_lgid].erase(tvitr);
    //
    return true;
}
//
bool calc_inter_clusters_cor(map<int, map<string, int> > cluster_vertex,  
                  map<int, unsigned long>            cluster_vertex_size_updated,
                  map<string, BIPAT_INT>             contigPollenGTSeq_int,
                  string                             tmpfolder,
                  map<string, unsigned long>         contigsize,
                  map<int, map<string, double> >*    cluster_edge_updated2,
                  map<int, unsigned long>*           cluster_vertex_size_updated2,
                  map<int, map<string, int> >*       cluster_vertex_updated2,
                  map<int, int>*                     mutually_exclusive_cluster)
{
    /*
        this function should be applied after break_clusters(): 
              merge small clusters of [2 Mb, 10 Mb] into larger ones of (10 Mb, +++)!
        idea: calculate inter-cluster correlation values;
              case 1. most paired-ctgs of two clusters showing correlation of 0 meaning the two clusters are non-relevant
              case 2. most paired-ctgs of two clusters showing correlation of >cutoff meaning the two clusters should be merged
              case 3. most paired-ctgs of two clusters showing correlation of <0 meaning the two clusters are homologous linkage groups
        //
        cluster_vertex..............: <clusterid, <ctg_id, 1> >
        cluster_vertex_size_updated.: <clusterid, cluster-contig-length>; raw-advanced
        contigPollenGTSeq_int.......: <ctg_id, {leftPat, rightPat}>; here is indeed selected pure haplotig markers
        tmpfolder...................: output directory: s3_        
        contigsize..................: <ctg_id, ctg_size>
        cluster_edge_updated2.......: <clusterid, <contig1-contig2, correlation value> >
        cluster_vertex_size_updated2: <clusterid, cluster-contig-length>; updated in this process
        cluster_vertex_updated2.....: <clusterid, <ctg_id, 1> >;  updated       
        mutually_exclusive_cluster..: <cluster_id1, cluster_id2>; cluster_id1/2 broken from the same old cluster id
    */
    string interLG_cor_file = tmpfolder + "/s3_genotype_haplotig_GT_similarity_matrix_subclusters_corr_check_torm.txt\0";
    bool LG_Cor_verbose  = true;
    ofstream ofp_LGcor;
    if(LG_Cor_verbose)
    {
        ofp_LGcor.open(interLG_cor_file.c_str(), ios::out);
        if(!ofp_LGcor.good())
        {
            cout   << "   Error: cannot open file " << interLG_cor_file << endl;
            return false;
        }
        ofp_LGcor  << "   Note: merge small clusters of [1 Mb, 10 Mb] into larger ones of [10 Mb, +++)! "
                   << endl
                   << "         Not all paired inter-clusters are checked."
                   << endl
                   << "         Further merging will happen in merge_clusters(...)"
                   << endl
                   << endl;
    } 
    // initialize variable collecting pairs of clusters to be merged: one larger might be merged with more smallers
    map<int, map<int, int> > clusters_to_be_merged; // <cluster_id2_larger, <cluster_id1_smallers, >1 >
    map<int, int> smaller_to_be_excluded;           // <cluster_id1_smaller, 1>; excluded in cluster_vertex_updated2
    // calculate inter-LG correlations    
    map<int, map<string, int> >::iterator clitr1;
    map<int, map<string, int> >::iterator clitr1_end;
    clitr1     = cluster_vertex.begin();
    clitr1_end = cluster_vertex.end();    
    while(clitr1 != clitr1_end)
    {
        int cluster1_id = (*clitr1).first;
        map<string, int> cluster1_contigs = (*clitr1).second; // <ctg_id, 1>
        // note only targeting clusters of [2, 10] Mb
        assert( cluster_vertex_size_updated.find(cluster1_id) != cluster_vertex_size_updated.end() );
        unsigned long cluster1_total_size = cluster_vertex_size_updated[cluster1_id];
        if(cluster1_total_size<2000000 || cluster1_total_size>10000000)
        {
            clitr1 ++; // skipping
            continue;
        }
        //
        int high_positive_cnt        =  0; // number of paired-ctgs show high positive cor>=0.25 (i.e., count_plus_3)
        double high_positive_ratio   =  0; // count_plus_3 / count_total
        int high_positive_cluster_id = -1; // the respective cluster 
        // 
        map<int, map<string, int> >::iterator clitr2;
        map<int, map<string, int> >::iterator clitr2_end;
        clitr2     = cluster_vertex.begin();
        clitr2_end = cluster_vertex.end(); 
        while(clitr2 != clitr2_end)
        {
            int cluster2_id = (*clitr2).first;
            if(cluster1_id == cluster2_id)
            {
                clitr2 ++;
                continue;
            }
            // check if cluster1 and cluster2 should avoid being merged 
            if((*mutually_exclusive_cluster).find(cluster1_id) != (*mutually_exclusive_cluster).end() && 
               (*mutually_exclusive_cluster)[cluster1_id] == cluster2_id )
            {
                clitr2 ++;
                continue;
            } 
            // note only targeting clusters of [10, +++) Mb
            assert( cluster_vertex_size_updated.find(cluster2_id) != cluster_vertex_size_updated.end() );
            unsigned long cluster2_total_size = cluster_vertex_size_updated[cluster2_id];
            if(cluster2_total_size<=10000000)
            {
                clitr2 ++; // skipping
                continue;
            }            
            //
            map<string, int> cluster2_contigs = (*clitr2).second; // <ctg_id, 1>   
            //
            if(LG_Cor_verbose)               
            {
                ofp_LGcor << "   inter-LG-cor-to-check: cluster " << cluster1_id 
                          << " vs. cluster "                      << cluster2_id 
                          << endl;
            }
            //
            int count_minus_3 = 0; // how many paired-ctgs of two clusters showed negative     correlation (-,      -0.25]
            int count_minus_2 = 0; // how many paired-ctgs of two clusters showed negative     correlation (-0.25,  -0.10]
            int count_minus_1 = 0; // how many paired-ctgs of two clusters showed negative     correlation (-0.10,      0]
            int count_plus_1  = 0; // how many paired-ctgs of two clusters showed negative     correlation (0,      +0.10]
            int count_plus_2  = 0; // how many paired-ctgs of two clusters showed negative     correlation (+0.10,  +0.25]
            int count_plus_3  = 0; // how many paired-ctgs of two clusters showed negative     correlation (+0.25,      +)            
            int count_total   = 0; // how many paired-ctgs of two clusters are calculated
            // calculate correlation
            map<string, int>::iterator c1itr;
            map<string, int>::iterator c1itr_end;
            c1itr     = cluster1_contigs.begin();
            c1itr_end = cluster1_contigs.end();
            while(c1itr != c1itr_end)
            {
                // get ctg1 genotype information in cluster1 
                string this_ctg1 = (*c1itr).first;
                map<string, BIPAT_INT>::iterator gt_citr1 = contigPollenGTSeq_int.find(this_ctg1);
                assert(gt_citr1 != contigPollenGTSeq_int.end() );
                BIPAT_INT pctg1;
                pctg1.leftPat_int  = (*gt_citr1).second.leftPat_int;
                pctg1.rightPat_int = (*gt_citr1).second.rightPat_int;                  
                //
                map<string, int>::iterator c2itr;
                map<string, int>::iterator c2itr_end;
                c2itr     = cluster2_contigs.begin();
                c2itr_end = cluster2_contigs.end();                
                while(c2itr != c2itr_end)
                {
                    // get ctg2 genotype information            
                    string this_ctg2 = (*c2itr).first;
                    if(this_ctg1.compare(this_ctg2) == 0)
                    {
                        cout << "   Error: same haplotig " << this_ctg1 << " found in two different clusters. " << endl;
                        return false;
                    }
                    map<string, BIPAT_INT>::iterator gt_citr2 = contigPollenGTSeq_int.find(this_ctg2);
                    assert(gt_citr2 != contigPollenGTSeq_int.end() );
                    BIPAT_INT pctg2;
                    pctg2.leftPat_int  = (*gt_citr2).second.leftPat_int;
                    pctg2.rightPat_int = (*gt_citr2).second.rightPat_int;  
                    // calculate correlation
                    double cor[4];
                    cor[0] = calculate_correlation(&pctg1.leftPat_int,  &pctg2.leftPat_int);      
                    cor[1] = calculate_correlation(&pctg1.leftPat_int,  &pctg2.rightPat_int);          
                    cor[2] = calculate_correlation(&pctg1.rightPat_int, &pctg2.leftPat_int);      
                    cor[3] = calculate_correlation(&pctg1.rightPat_int, &pctg2.rightPat_int);            
                    // find highest correlation: can be minus
                    int    best_i   = 0;
                    double best_cor = 0; 
                    for(int i = 0; i < 4; i ++)
                    {                     
                        if(fabs(cor[i])>best_cor)
                        {
                            best_cor = fabs(cor[i]);
                            best_i   = i;                                                         
                        }
                    }
                    best_cor = cor[best_i]; // re-collect minus/plus
                    //
                    if(best_cor <= -0.25)
                    {
                        count_minus_3 ++;
                    }else                    
                    if(best_cor <= -0.10)
                    {
                        count_minus_2 ++;
                    }else
                    if(best_cor <=  0.00)
                    {
                        count_minus_1 ++;
                    }
                    else
                    if(best_cor <= 0.10)
                    {
                        count_plus_1 ++;
                    }else
                    if(best_cor <= 0.25)
                    {
                        count_plus_2 ++;
                    }else                    
                    {
                        count_plus_3 ++;
                    }
                    count_total ++;
                    if(LG_Cor_verbose)               
                    {
                        ofp_LGcor << "       "        << this_ctg1 
                                  << " -- "           << this_ctg2
                                  << " : best_cor = " << best_cor
                                  << endl;
                    }                                               
                    // next ctg in cluster 2
                    c2itr ++;
                }
                // next ctg in cluster 1
                c1itr ++;
            }
            //
            if(count_plus_3 > high_positive_cnt)
            {
                high_positive_cnt        = count_plus_3;
                high_positive_ratio      = count_plus_3*1.0/count_total;
                high_positive_cluster_id = cluster2_id;
            }
            if(LG_Cor_verbose)               
            {
                ofp_LGcor << "   summary: cluster " << cluster1_id 
                          << " vs. cluster "        << cluster2_id 
                          << ": "
                          << "      total correlations checked: "  << count_total 
                          << "; among these, "                     << endl
                          << "          cor in (-----,   -0.25]: " << count_minus_3 << endl                          
                          << "          cor in (-0.25,   -0.10]: " << count_minus_2 << endl
                          << "          cor in (-0.10,   00.00]: " << count_minus_1 << endl                          
                          << "          cor in (00.00,   +0.10]: " << count_plus_1  << endl
                          << "          cor in (+0.10,   +0.25]: " << count_plus_2  << endl
                          << "          cor in (+0.25,   +++++): " << count_plus_3  << endl;
            }            
            // next cluster2
            clitr2 ++;
        }
        // high_positive_cluster_id will collect contigs in smaller cluster1_id
        if(high_positive_cnt>=50 && high_positive_ratio>=0.25)
        {            
            if(clusters_to_be_merged.find(high_positive_cluster_id) != clusters_to_be_merged.end() )
            {
                clusters_to_be_merged[ high_positive_cluster_id ].insert(std::pair<int, int>(cluster1_id, 1));
            }else
            {
                map<int, int> smaller_set;
                smaller_set.insert(std::pair<int, int>(cluster1_id, 1));
                clusters_to_be_merged.insert(std::pair<int, map<int, int> >(high_positive_cluster_id, smaller_set));
            }
            //
            smaller_to_be_excluded.insert(std::pair<int, int>(cluster1_id, 1));
            //
            if(LG_Cor_verbose)               
            {
                ofp_LGcor << "   conclusion: cluster "         << cluster1_id 
                          << " having "                        << cluster_vertex[cluster1_id].size()
                          << " contigs of "                    << cluster_vertex_size_updated[cluster1_id]
                          << " bp will be merged into cluster "<< high_positive_cluster_id 
                          << " having "                        << cluster_vertex[high_positive_cluster_id].size()
                          << " contigs of "                    << cluster_vertex_size_updated[high_positive_cluster_id]
                          << " bp with inter-cluster-paired-contgis of high-positive-correlation: " << high_positive_cnt
                          << " (ratio of count_plus_3 / count_total = " << high_positive_ratio
                          << ")"
                          << endl << endl;                
            }            
        }   
        // next cluster1
        clitr1 ++;
    }
    //
    if(LG_Cor_verbose)
    {
        ofp_LGcor.close();
    } 
    // update cluster with merging of size-selected smaller clusters [2, 10 Mb] into size-selected larger ones (over 10 Mb). 
    map<int, map<string, int> >::iterator clitr;
    map<int, map<string, int> >::iterator clitr_end;    
    clitr     = cluster_vertex.begin();
    clitr_end = cluster_vertex.end(); 
    while(clitr != clitr_end)
    {
        int cluster_id = (*clitr).first;
        if( smaller_to_be_excluded.find( cluster_id ) != smaller_to_be_excluded.end() )
        {
            clitr ++;
            if(LG_Cor_verbose)
            {
                ofp_LGcor << "   Update: this cluster " << cluster_id << " would be merged into others. " << endl;
            }            
            continue;
        }
        // collect vertex set 
        map<string, int> this_cluster_vertex = (*clitr).second;        
        (*cluster_vertex_updated2).insert(std::pair<int, map<string, int> >(cluster_id, this_cluster_vertex));
        //
        if(clusters_to_be_merged.find(cluster_id) != clusters_to_be_merged.end() )
        {
            map<int, int>::iterator mergeitr;
            map<int, int>::iterator mergeitr_end;
            mergeitr     = clusters_to_be_merged[cluster_id].begin();
            mergeitr_end = clusters_to_be_merged[cluster_id].end();
            while(mergeitr != mergeitr_end)
            {
                int smaller_cluster_id = (*mergeitr).first;
                (*cluster_vertex_updated2)[cluster_id].insert(cluster_vertex[ smaller_cluster_id ].begin(), 
                                                           cluster_vertex[ smaller_cluster_id ].end());
                if(LG_Cor_verbose)
                {
                    ofp_LGcor << "   Update: small cluster "      << smaller_cluster_id 
                              << " has been merged into cluster " << cluster_id
                              << endl;
                }                                               
                //
                mergeitr ++;
            }
        }
        //
        clitr ++;
    }    
    // re-build sets of edges and so on after merging some targeted clusters
    //  contigsize..................: <ctg_id, ctg_size>
    //  cluster_edge_updated2.......: <clusterid, <contig1-contig2, correlation value> >
    //  cluster_vertex_size_updated2: <clusterid, cluster-contig-length>  
    //  cluster_vertex_updated2.....: <clusterid, <ctg_id, 1> >;  updated       
    // build cluster_edge_updated2
    // map<int, map<string, int> >::iterator clitr;
    // map<int, map<string, int> >::iterator clitr_end;
    clitr     = (*cluster_vertex_updated2).begin();
    clitr_end = (*cluster_vertex_updated2).end();            
    while(clitr != clitr_end)
    {
        int cluster_id = (*clitr).first;
        map<string, int> this_cluster_vertex = (*clitr).second;
        // initiaze edge variable
        assert( (*cluster_edge_updated2).find(cluster_id) == (*cluster_edge_updated2).end() );
        map<string, double> this_cluster_edge;
        (*cluster_edge_updated2).insert(std::pair<int, map<string, double> >(cluster_id, this_cluster_edge));
        // initialize total vertex size 
        assert( (*cluster_vertex_size_updated2).find(cluster_id) == (*cluster_vertex_size_updated2).end() );
        (*cluster_vertex_size_updated2).insert(std::pair<int, unsigned long>(cluster_id, 0));
        //
        map<string, int>::iterator ctgitr1;
        map<string, int>::iterator ctgitr1_end;
        ctgitr1     = this_cluster_vertex.begin();
        ctgitr1_end = this_cluster_vertex.end();
        while(ctgitr1 != ctgitr1_end)
        {
            string this_ctg1 = (*ctgitr1).first;
            // update total vertex size 
            assert( contigsize.find(this_ctg1) != contigsize.end() );
            (*cluster_vertex_size_updated2)[ cluster_id ] += contigsize[ this_ctg1 ];
            //
            map<string, BIPAT_INT>::iterator gt_citr1 = contigPollenGTSeq_int.find(this_ctg1);
            assert(gt_citr1 != contigPollenGTSeq_int.end() );
            BIPAT_INT pctg1;
            pctg1.leftPat_int  = (*gt_citr1).second.leftPat_int;
            pctg1.rightPat_int = (*gt_citr1).second.rightPat_int;
            //
            map<string, int>::iterator ctgitr2;
            map<string, int>::iterator ctgitr2_end;
            ctgitr2     = this_cluster_vertex.begin();
            ctgitr2_end = this_cluster_vertex.end();    
            while(ctgitr2 != ctgitr2_end)
            {
                string this_ctg2 = (*ctgitr2).first;   
                if(this_ctg1.compare(this_ctg2) == 0)
                {
                    ctgitr2 ++;
                    continue;
                }
                //
                string this_edge_tmp1 = this_ctg1 + "-" + this_ctg2;                
                string this_edge_tmp2 = this_ctg2 + "-" + this_ctg1;
                if((*cluster_edge_updated2)[cluster_id].find(this_edge_tmp1) != (*cluster_edge_updated2)[cluster_id].end() || 
                   (*cluster_edge_updated2)[cluster_id].find(this_edge_tmp2) != (*cluster_edge_updated2)[cluster_id].end() )
                {
                    ctgitr2 ++;
                    continue;
                }
                //                
                map<string, BIPAT_INT>::iterator gt_citr2 = contigPollenGTSeq_int.find(this_ctg2);
                assert(gt_citr2 != contigPollenGTSeq_int.end() );
                BIPAT_INT pctg2;
                pctg2.leftPat_int  = (*gt_citr2).second.leftPat_int;
                pctg2.rightPat_int = (*gt_citr2).second.rightPat_int;                   
                // calculate correlation
                double cor[4];
                cor[0] = calculate_correlation(&pctg1.leftPat_int,  &pctg2.leftPat_int);      
                cor[1] = calculate_correlation(&pctg1.leftPat_int,  &pctg2.rightPat_int);          
                cor[2] = calculate_correlation(&pctg1.rightPat_int, &pctg2.leftPat_int);      
                cor[3] = calculate_correlation(&pctg1.rightPat_int, &pctg2.rightPat_int);            
                // update larger cluster with better correlation
                int best_i              = 0;
                double best_cluster_cor = 0;
                for(int i = 0; i < 4; i ++)
                {                     
                    if(cor[i]>best_cluster_cor)
                    {
                        best_cluster_cor = cor[i];
                        best_i           = i;                                                         
                    }
                }
                // update edge 
                if(best_cluster_cor > minCorscore)
                {
                    (*cluster_edge_updated2)[cluster_id].insert(std::pair<string, double>(this_edge_tmp1, best_cluster_cor));
                }
                // next ctg2
                ctgitr2 ++;
            }
            // next ctg1
            ctgitr1 ++;
        }
        // next cluster 
        clitr ++;
    }
    // output updated raw-advanced-merged clusters: 
    // this is processed x 2 raw, can be more than expected number of linkage groups thus needs further merging.
    string subclusterinfo = tmpfolder + "/s3_genotype_haplotig_GT_similarity_matrix_subclusters_raw_advanced_merge.dot\0"; // 
    ofstream subcluster_ofp;
    subcluster_ofp.open(subclusterinfo.c_str(), ios::out);
    if(!subcluster_ofp.good())
    {
        cout   << "   Error: cannot open file " << subclusterinfo << endl;
        return false;
    } 
    subcluster_ofp << "/* Here are the subclusters of contigs with correlation based cluster break&merge */" << endl;
    subcluster_ofp << "graph\tGraph_1 {"                           << endl;
    // to collect all clusters-specific vertex
    map<int, map<string, double> >::iterator cleitr;
    map<int, map<string, double> >::iterator cleitr_end;
    cleitr     = (*cluster_edge_updated2).begin();
    cleitr_end = (*cluster_edge_updated2).end();
    unsigned long total_contig_size = 0;
    unsigned long total_contig_numb = 0;
    while(cleitr != cleitr_end)
    {
        subcluster_ofp << "\tsubgraph cluster_" << (*cleitr).first << " {" << endl;
        unsigned long cluster_ctg_size = 0;
        map<string, int> subcluster_vertex;
        //
        map<string, double>::iterator eitr;
        map<string, double>::iterator eitr_end;
        eitr     = (*cleitr).second.begin();
        eitr_end = (*cleitr).second.end();
        while(eitr != eitr_end)
        {
            vector<string> edgeinfo = split_string((*eitr).first, '-');
            subcluster_ofp << "\t"
            << edgeinfo[0] << " -- " << edgeinfo[1]
            << " "            << "[color=gold, penwidth=1, arrowsize=1, label=" << (*eitr).second << "];" 
            << endl;  
            //
            if(subcluster_vertex.find(edgeinfo[0]) == subcluster_vertex.end() )
            {
                subcluster_vertex.insert(std::pair<string, int>(edgeinfo[0], 1));
                assert(contigsize.find(edgeinfo[0]) != contigsize.end());
                cluster_ctg_size += contigsize[edgeinfo[0]];
            }
            if(subcluster_vertex.find(edgeinfo[1]) == subcluster_vertex.end() )
            {
                subcluster_vertex.insert(std::pair<string, int>(edgeinfo[1], 1));
                assert(contigsize.find(edgeinfo[1]) != contigsize.end());
                cluster_ctg_size += contigsize[edgeinfo[1]];
            }            
            //           
            eitr ++;
        }
        //
        subcluster_ofp << "\t/* "          << subcluster_vertex.size() << " contigs with total size of " 
                       << cluster_ctg_size << " bp */"              << endl;
        if(clusters_to_be_merged.find( (*cleitr).first ) != clusters_to_be_merged.end() )
        {
            subcluster_ofp << "\t/* this cluster has been merged with smaller clusters of "; 
            map<int, int>::iterator mergeitr;
            map<int, int>::iterator mergeitr_end;
            mergeitr     = clusters_to_be_merged[ (*cleitr).first ].begin();
            mergeitr_end = clusters_to_be_merged[ (*cleitr).first ].end();
            while(mergeitr != mergeitr_end)
            {
                int smaller_cluster_id = (*mergeitr).first;    
                if(mergeitr != clusters_to_be_merged[ (*cleitr).first ].begin() )
                {
                    subcluster_ofp << "," << smaller_cluster_id;                    
                }else
                {
                    subcluster_ofp << smaller_cluster_id;
                }                
                //
                mergeitr ++;
            }
            subcluster_ofp << " */"              << endl;                                    
        }
        subcluster_ofp << "\t}" << endl;     
        //
        total_contig_numb += subcluster_vertex.size();
        total_contig_size += cluster_ctg_size;  
        // next cluster of edges
        cleitr ++;
    }
    //
    subcluster_ofp << "}" << endl;
    subcluster_ofp.close();    
    //
    cout << "   Info: after seletively merging clusters: " << total_contig_numb 
         << " contigs of "                                 << total_contig_size 
         << " bp clustered into "                          << (*cluster_edge_updated2).size() 
         << " clusters. "                                  << endl;     
    // 
    return true;
}                  
//
bool break_clusters(map<int, map<string, int> >   cluster_vertex,  
                  map<string, BIPAT_INT>          contigPollenGTSeq_int,
                  string                          tmpfolder,
                  map<string, unsigned long>      contigsize,
                  map<int, map<string, double> >* cluster_edge_updated,
                  map<int, unsigned long>*        cluster_vertex_size_updated,
                  map<int, map<string, int> >*    cluster_vertex_updated,
                  map<int, int>*                  mutually_exclusive_cluster)
{
    /*
        this function should be applied before merge_clusters()!
        idea: contigs in the same Lg should all show positive correlation to each other;
              once a cluster is broken, find a most closely correlation one directly
        //
        cluster_vertex.............: <clusterid, <ctg_id, 1> >
        contigPollenGTSeq_int......: <ctg_id, {leftPat, rightPat}>; here is indeed selected pure haplotig markers
        tmpfolder..................: output directory: s3_        
        contigsize.................: <ctg_id, ctg_size>
        cluster_edge_updated.......: <clusterid, <contig1-contig2, correlation value> >
        cluster_vertex_size_updated: <clusterid, cluster-contig-length>  
        cluster_vertex_updated.....: <clusterid, <ctg_id, 1> >;  updated       
        mutually_exclusive_cluster.: <cluster_id1, cluster_id2>; cluster_id1/2 broken from the same old cluster id
    */ 
    string breakLGfile = tmpfolder + "/s3_genotype_haplotig_GT_similarity_matrix_subclusters_breakLG_check_torm.txt\0";
    bool break_verbose  = false;
    ofstream ofp_break;
    if(break_verbose)
    {
        ofp_break.open(breakLGfile.c_str(), ios::out);
        if(!ofp_break.good())
        {
            cout   << "   Error: cannot open file " << breakLGfile << endl;
            return false;
        }
    }
    // traverse each cluster, and try to break a cluster into two if many vertices showing negative correlation
    map<int, map<string, int> >::iterator clitr;
    map<int, map<string, int> >::iterator clitr_end;
    clitr     = cluster_vertex.begin();
    clitr_end = cluster_vertex.end();    
    while(clitr != clitr_end)
    {
        int cluster_id = (*clitr).first;        
        map<string, int> this_cluster_vertex = (*clitr).second;
        // if cluster size is already small, skip breaking it 
        if(this_cluster_vertex.size() < 5)
        {
            if(break_verbose)            
            ofp_break << "   recalc: skipped inter-contig cor within cluster "  << cluster_id
                      << " with "                                               << this_cluster_vertex.size() 
                      << " vertices; collected as is. "                         << endl; 
            // collect this cluster as is
            (*cluster_vertex_updated).insert(std::pair<int, map<string, int> >(cluster_id, this_cluster_vertex));
            // next cluster 
            clitr ++;
            continue;            
        }
        //
        if(break_verbose)        
        ofp_break << "   recalc: checking inter-contig cor within cluster "  << cluster_id
                  << " with "      << this_cluster_vertex.size() 
                  << " vertices. " << endl;        
        //
        map<string, map<string, double> > positive_inter_cor; // <ctg1, <ctg2, cor(>=0) > > 
        map<string, map<string, double> > negative_inter_cor; // <ctg1, <ctg2, cor(< 0) > > 
        //
        map<string, int>::iterator ctgitr1;
        map<string, int>::iterator ctgitr1_end;
        ctgitr1     = this_cluster_vertex.begin();
        ctgitr1_end = this_cluster_vertex.end();
        while(ctgitr1 != ctgitr1_end)
        {
            // get ctg1 genotype information
            string this_ctg1 = (*ctgitr1).first;
            map<string, BIPAT_INT>::iterator gt_citr1 = contigPollenGTSeq_int.find(this_ctg1);
            assert(gt_citr1 != contigPollenGTSeq_int.end() );
            BIPAT_INT pctg1;
            pctg1.leftPat_int  = (*gt_citr1).second.leftPat_int;
            pctg1.rightPat_int = (*gt_citr1).second.rightPat_int;             
            //
            map<string, int>::iterator ctgitr2;
            map<string, int>::iterator ctgitr2_end;
            ctgitr2     = this_cluster_vertex.begin();
            ctgitr2_end = this_cluster_vertex.end();  
            while(ctgitr2 != ctgitr2_end)
            {
                // get ctg2 genotype information            
                string this_ctg2 = (*ctgitr2).first;
                if(this_ctg1.compare(this_ctg2) == 0)
                {
                    ctgitr2 ++;
                    continue;
                }
                map<string, BIPAT_INT>::iterator gt_citr2 = contigPollenGTSeq_int.find(this_ctg2);
                assert(gt_citr2 != contigPollenGTSeq_int.end() );
                BIPAT_INT pctg2;
                pctg2.leftPat_int  = (*gt_citr2).second.leftPat_int;
                pctg2.rightPat_int = (*gt_citr2).second.rightPat_int;                   
                // calculate correlation
                double cor[4];
                cor[0] = calculate_correlation(&pctg1.leftPat_int,  &pctg2.leftPat_int);      
                cor[1] = calculate_correlation(&pctg1.leftPat_int,  &pctg2.rightPat_int);          
                cor[2] = calculate_correlation(&pctg1.rightPat_int, &pctg2.leftPat_int);      
                cor[3] = calculate_correlation(&pctg1.rightPat_int, &pctg2.rightPat_int);            
                // update larger cluster with better correlation
                int best_i              = 0;
                double best_cluster_cor = 0;
                for(int i = 0; i < 4; i ++)
                {                     
                    if(fabs(cor[i])>best_cluster_cor)
                    {
                        best_cluster_cor = fabs(cor[i]);
                        best_i           = i;                                                         
                    }
                }       
                best_cluster_cor = cor[best_i]; // note it might become negative if two contigs belong to different LGs
                if(break_verbose)                
                ofp_break << "   check: best cor of " << this_ctg1 
                          << " with "                 << this_ctg2
                          << " updated as "           << best_cluster_cor 
                          << endl;                
                if(best_cluster_cor>=0)
                {
                    if(positive_inter_cor.find(this_ctg1) == positive_inter_cor.end() )
                    {
                        map<string, double> tmpcor;
                        tmpcor.insert(std::pair<string, double>(this_ctg2, best_cluster_cor));
                        positive_inter_cor.insert(std::pair<string, map<string, double> >(this_ctg1, tmpcor));
                    }else
                    {
                        assert(positive_inter_cor[this_ctg1].find(this_ctg2) == positive_inter_cor[this_ctg1].end());
                        positive_inter_cor[this_ctg1].insert(std::pair<string, double>(this_ctg2, best_cluster_cor));
                    }
                }else
                {
                    if(negative_inter_cor.find(this_ctg1) == negative_inter_cor.end() )
                    {
                        map<string, double> tmpcor;
                        tmpcor.insert(std::pair<string, double>(this_ctg2, best_cluster_cor));
                        negative_inter_cor.insert(std::pair<string, map<string, double> >(this_ctg1, tmpcor));
                    }else
                    {
                        assert(negative_inter_cor[this_ctg1].find(this_ctg2) == negative_inter_cor[this_ctg1].end());
                        negative_inter_cor[this_ctg1].insert(std::pair<string, double>(this_ctg2, best_cluster_cor));
                    }               
                }                                                
                // next ctg2
                ctgitr2 ++;
            }      
            // check 
            if(break_verbose)            
            if(negative_inter_cor.find(this_ctg1) != negative_inter_cor.end())
            {
                ofp_break << "   check: "                     << negative_inter_cor[ this_ctg1 ].size() 
                          << "/"                              << this_cluster_vertex.size()
                          << " ctgs showing negative cor to " << this_ctg1 << endl;                 
            }
            if(break_verbose)            
            if(positive_inter_cor.find(this_ctg1) != positive_inter_cor.end())
            {
                ofp_break << "   check: "                     << positive_inter_cor[ this_ctg1 ].size()
                          << "/"                              << this_cluster_vertex.size()                
                          << " ctgs showing positive cor to " << this_ctg1 << endl;                
            }            
            // next ctg1
            ctgitr1 ++;
        }
        // check if mutually exclusive contigs exist 
        map<string, int> this_cluster_positive_set; // in case separation occurs, this is one     cluster 
        map<string, int> this_cluster_negative_set; // in case separation occurs, this is another cluster 
        map<string, int>::iterator meitr1;
        map<string, int>::iterator meitr1_end;
        meitr1     = this_cluster_vertex.begin();
        meitr1_end = this_cluster_vertex.end();
        while(meitr1 != meitr1_end)
        {
            string this_ctg1 = (*meitr1).first;
            if(this_cluster_positive_set.size() == 0)
            {
                // initialize a positive cluster
                this_cluster_positive_set.insert(std::pair<string, int>(this_ctg1, 1));
            }else
            {
                int which_to_check = 0;
                int this_ctg1_positive_support = 0;
                if(positive_inter_cor.find(this_ctg1) != positive_inter_cor.end() )
                {
                    this_ctg1_positive_support = positive_inter_cor[this_ctg1].size();
                }
                int this_ctg1_negative_support = 0;                
                if(negative_inter_cor.find(this_ctg1) != negative_inter_cor.end() )
                {
                    this_ctg1_negative_support = negative_inter_cor[this_ctg1].size();
                }
                if(this_ctg1_positive_support > this_ctg1_negative_support)
                {
                    which_to_check = 1;  // check positive set of this_ctg1                    
                }else
                if(this_ctg1_positive_support < this_ctg1_negative_support)                
                {
                    which_to_check = -1; // check negative set of this_ctg1                                        
                }else
                if(this_ctg1_positive_support > 0)
                {
                    // equal size: select any one would be okay!
                    which_to_check = 1;
                }else 
                {
                    cout << "   Error: no correlation contigs found? " << endl;
                    return false;
                }                
                // case 1: comparing using its positive set
                bool goto_positive_case1 = false;
                bool goto_negative_case1 = false;                
                map<string, double> positive_ctg1;
                if(positive_inter_cor.find(this_ctg1) != positive_inter_cor.end() && which_to_check==1)
                {
                    // check positive contig list of this ctg1
                    positive_ctg1 = positive_inter_cor[this_ctg1]; // <ctg1's ctg-list, ">=0">
                    //
                    map<string, int>::iterator meitr_tmp;
                    map<string, int>::iterator meitr_tmp_end;                    
                    // compare with existing positive set 
                    meitr_tmp     = this_cluster_positive_set.begin();
                    meitr_tmp_end = this_cluster_positive_set.end();
                    int goto_positive_set = 0; // 
                    int goto_negative_set = 0; // times supporting this positive_ctg1 goes to positive/negative set
                    int goto_positive_cnt = 0; // 
                    int goto_negative_cnt = 0; // ctgs  supporting this negative_ctg1 goes to positive/negative set                         
                    while(meitr_tmp != meitr_tmp_end)
                    {
                        string tmp_ctg = (*meitr_tmp).first;
                        int pos_common = 0;
                        int neg_common = 0;
                        //
                        if(positive_inter_cor.find(tmp_ctg) != positive_inter_cor.end() )
                        {
                            // compare with positive contig list of this tmp contig: positive_ctg1 vs positive_tmp_ctg
                            map<string, double> positive_tmp_ctg;
                            positive_tmp_ctg = positive_inter_cor[tmp_ctg]; // <tmp_ctg's ctg-list, ">=0">
                            pos_common = find_contig_set_overlap(positive_ctg1, positive_tmp_ctg);                            
                        }
                        if(negative_inter_cor.find(tmp_ctg) != negative_inter_cor.end() )
                        {
                            // compare with negative contig list of this tmp contig: positive_ctg1 vs negative_tmp_ctg 
                            map<string, double> negative_tmp_ctg;
                            negative_tmp_ctg = negative_inter_cor[tmp_ctg]; // <tmp_ctg's ctg-list, "<0">
                            neg_common = find_contig_set_overlap(positive_ctg1, negative_tmp_ctg);                                                        
                        }
                        //
                        if(break_verbose)                        
                        ofp_break << "   check: positive set of ctg "    << this_ctg1 
                                  << " common with positive set of ctg " << tmp_ctg 
                                  << ": "                                << pos_common
                                  << endl;
                        if(break_verbose)                                  
                        ofp_break << "   check: positive set of ctg "    << this_ctg1 
                                  << " common with negative set of ctg " << tmp_ctg 
                                  << ": "                                << neg_common
                                  << endl;
                        // checking positive set of this_ctg1 at the moment, if overlapping more + of subject: ++ => +                                  
                        if(pos_common > neg_common)
                        {
                            goto_positive_set ++;
                            goto_positive_cnt += pos_common;   
                            if(break_verbose)                                                     
                            ofp_break << "   check: this comparion supports "                 << this_ctg1 
                                      << " going to positive set: this_cluster_positive_set." << endl;                            
                        }else
                        if(pos_common < neg_common)
                        {
                            goto_negative_set ++;
                            goto_negative_cnt += neg_common;          
                            if(break_verbose)                                              
                            ofp_break << "   check: this comparion supports " << this_ctg1 
                                      << " going to negative set: this_cluster_negative_set." << endl; 
                        }else ;                                                   
                        //
                        meitr_tmp ++;
                    }
                    // collect this_ctg1
                    if(goto_positive_set > goto_negative_set)
                    {
                        goto_positive_case1 = true;
                    }else
                    if(goto_positive_set < goto_negative_set)                    
                    {
                        goto_negative_case1 = true;                    
                    }else
                    {
                        if(goto_positive_cnt==goto_negative_cnt) // compare absolute number of ctgs supporting moving!
                        {
                            cout      << "   Error: equality of going to both sets should never happen. " << endl;
                            if(break_verbose)                            
                            ofp_break << "   Error: equality of going to both sets should never happen. " << endl;                        
                            return false;
                        }else
                        if(goto_positive_cnt > goto_negative_cnt)
                        {
                            goto_positive_case1 = true;
                        }else
                        {
                            goto_negative_case1 = true;                    
                        }
                    }
                    // compare with existing negative set 
                }
                // 
                // case 2: comparing using its negative set  
                bool goto_positive_case2 = false;
                bool goto_negative_case2 = false;                               
                map<string, double> negative_ctg1;
                if(negative_inter_cor.find(this_ctg1) != negative_inter_cor.end() && which_to_check==-1)
                {
                    // check negative contig list of this ctg1
                    negative_ctg1 = negative_inter_cor[this_ctg1]; // <ctg1's ctg-list, "<0">
                    //
                    map<string, int>::iterator meitr_tmp;
                    map<string, int>::iterator meitr_tmp_end;                    
                    // compare with existing positive set 
                    meitr_tmp     = this_cluster_positive_set.begin();
                    meitr_tmp_end = this_cluster_positive_set.end();
                    int goto_positive_set = 0; // 
                    int goto_negative_set = 0; // times supporting this negative_ctg1 goes to positive/negative set
                    int goto_positive_cnt = 0; // 
                    int goto_negative_cnt = 0; // ctgs  supporting this negative_ctg1 goes to positive/negative set                    
                    while(meitr_tmp != meitr_tmp_end)
                    {
                        string tmp_ctg = (*meitr_tmp).first;
                        int pos_common = 0;
                        int neg_common = 0;
                        //
                        if(positive_inter_cor.find(tmp_ctg) != positive_inter_cor.end() )
                        {
                            // compare with positive contig list of this tmp contig: negative_ctg1 vs positive_tmp_ctg
                            map<string, double> positive_tmp_ctg;
                            positive_tmp_ctg = positive_inter_cor[tmp_ctg]; // <tmp_ctg's ctg-list, ">=0">
                            pos_common = find_contig_set_overlap(negative_ctg1, positive_tmp_ctg);                            
                        }
                        if(negative_inter_cor.find(tmp_ctg) != negative_inter_cor.end() )
                        {
                            // compare with negative contig list of this tmp contig: negative_ctg1 vs negative_tmp_ctg 
                            map<string, double> negative_tmp_ctg;
                            negative_tmp_ctg = negative_inter_cor[tmp_ctg]; // <tmp_ctg's ctg-list, "<0">
                            neg_common = find_contig_set_overlap(negative_ctg1, negative_tmp_ctg);                                                        
                        }
                        //
                        if(break_verbose)                        
                        ofp_break << "   check: negative set of ctg "    << this_ctg1 
                                  << " common with positive set of ctg " << tmp_ctg 
                                  << ": "                                << pos_common
                                  << endl;
                        if(break_verbose)                                  
                        ofp_break << "   check: negative set of ctg "    << this_ctg1 
                                  << " common with negative set of ctg " << tmp_ctg 
                                  << ": "                                << neg_common
                                  << endl;
                        // checking negative set of this_ctg1 at the moment, if overlapping more + of subject: -+ => -  
                        // note below needs to be reserved for selecting positive/negative set!                                
                        if(pos_common > neg_common)
                        {
                            goto_negative_set ++;
                            goto_negative_cnt += pos_common;
                            if(break_verbose)                            
                            ofp_break << "   check: this comparion supports "                  << this_ctg1 
                                      << " going to negative set: this_cluster_negative_set. " << endl;                            
                        }else   
                        if(pos_common < neg_common)
                        {
                            goto_positive_set ++;
                            goto_positive_cnt += neg_common;                            
                            if(break_verbose)                            
                            ofp_break << "   check: this comparion supports " << this_ctg1 
                                      << " going to positive set: this_cluster_positive_set "  << endl;                             
                        }else ;
                        //
                        meitr_tmp ++;
                    }
                    // collect this_ctg1
                    if(goto_positive_set > goto_negative_set)
                    {
                        goto_positive_case2 = true;
                    }else
                    if(goto_positive_set < goto_negative_set)                    
                    {
                        goto_negative_case2 = true;                    
                    }else
                    {
                        if(goto_positive_cnt == goto_negative_cnt) // compare absolute number of ctgs supporting moving!
                        {
                            cout      << "   Error: equality of going to both sets should never happen. " << endl;
                            if(break_verbose)                            
                            ofp_break << "   Error: equality of going to both sets should never happen. " << endl;                                                
                            return false;
                        }else
                        if(goto_positive_cnt > goto_negative_cnt)
                        {
                            goto_positive_case2 = true;                            
                        }else
                        {
                            goto_negative_case2 = true;
                        }
                    }
                    // compare with existing negative set 
                }
                // merge checking result:
                if( positive_inter_cor.find(this_ctg1) != positive_inter_cor.end() && which_to_check==1)
                {
                    if(goto_positive_case1==true)
                    {
                        if(break_verbose)                    
                        ofp_break << "   check: case 2.1 - final decision of "                        << this_ctg1 
                                  << " going to positive set: this_cluster_positive_set. " << endl;
                        this_cluster_positive_set.insert(std::pair<string, int>(this_ctg1, 1));                        
                    }else
                    {
                        if(break_verbose)                    
                        ofp_break << "   check: case 2.2 - final decision of "                        << this_ctg1 
                                  << " going to negative set: this_cluster_negative_set. " << endl; 
                        this_cluster_negative_set.insert(std::pair<string, int>(this_ctg1, 1));                        
                    }
                }else
                if( negative_inter_cor.find(this_ctg1) != negative_inter_cor.end() && which_to_check==-1)
                {
                    if(goto_positive_case2==true)
                    {
                        if(break_verbose)                    
                        ofp_break << "   check: case 3.1 - final decision of "                        << this_ctg1 
                                  << " going to positive set: this_cluster_positive_set. " << endl;
                        this_cluster_positive_set.insert(std::pair<string, int>(this_ctg1, 1));
                    }else
                    {
                        if(break_verbose)                    
                        ofp_break << "   check: case 3.2 - final decision of "                        << this_ctg1 
                                  << " going to negative set: this_cluster_negative_set. " << endl; 
                        this_cluster_negative_set.insert(std::pair<string, int>(this_ctg1, 1)); 
                    }
                }else 
                { 
                    cout << "   Error: no valid cluster breaking. " << endl;
                    return false;
                }
            }
            //
            meitr1 ++;
        }
        //
        if(break_verbose)        
        ofp_break << "   check: for this cluster "   << (*clitr).first                   << ":"
                  << " vertices in positive set = "  << this_cluster_positive_set.size() << " vs. "
                  << " vertices in negative set = "  << this_cluster_negative_set.size()
                  << endl;
        // update clusters of vertices: 
        if(this_cluster_positive_set.size()>0 && this_cluster_negative_set.size()>0)
        {
            // current raw cluster has been broken into two 
            // collect positive set using original cluster id 
            (*cluster_vertex_updated).insert(std::pair<int, map<string, int> >(cluster_id, this_cluster_positive_set));
            // create another cluster id for negative set 
            int new_cluster_id = 0;
            while(cluster_vertex.find(new_cluster_id) != cluster_vertex.end() )
            {
                new_cluster_id ++;
            }
            (*cluster_vertex_updated).insert(std::pair<int, map<string, int> >(new_cluster_id, this_cluster_negative_set));  
            if(break_verbose)                 
            ofp_break << "   check: this cluster "            << cluster_id
                      << " has been broken into old cluster " << cluster_id << " and "
                      << " new cluster "                      << new_cluster_id
                      << endl;   
            // create mutually exclusive pairs of cluster ids; they should aovid being merged later! 
            (*mutually_exclusive_cluster).insert(std::pair<int, int>(cluster_id, new_cluster_id));
            (*mutually_exclusive_cluster).insert(std::pair<int, int>(new_cluster_id, cluster_id));             
        }else
        if(this_cluster_positive_set.size()>0)
        {
            (*cluster_vertex_updated).insert(std::pair<int, map<string, int> >(cluster_id, this_cluster_positive_set)); 
            if(break_verbose)            
            ofp_break << "   check: this cluster "               << cluster_id
                      << " has been recollected as old cluster " << cluster_id
                      << endl;                       
        }else
        if(this_cluster_negative_set.size()>0)
        {
            (*cluster_vertex_updated).insert(std::pair<int, map<string, int> >(cluster_id, this_cluster_negative_set));  
            if(break_verbose)            
            ofp_break << "   check: this cluster "               << cluster_id
                      << " has been recollected as old cluster " << cluster_id
                      << endl;                             
        }else
        {
            ;
        }
        if(break_verbose)
        {
            // check purpose
            map<string, int>::iterator checkitr;
            map<string, int>::iterator checkitr_end;
            ofp_break << "      positive set: " << endl;            
            checkitr     = this_cluster_positive_set.begin();
            checkitr_end = this_cluster_positive_set.end();
            while(checkitr != checkitr_end)
            {
                ofp_break << "         " << (*checkitr).first << endl;
                checkitr ++;
            }
            ofp_break << "      negative set: " << endl;            
            checkitr     = this_cluster_negative_set.begin();
            checkitr_end = this_cluster_negative_set.end();
            while(checkitr != checkitr_end)
            {
                ofp_break << "         " << (*checkitr).first << endl;
                checkitr ++;
            }            
        }        
        // next cluster please
        clitr ++;
    }
    //
    if(break_verbose)   
    ofp_break.close();
    // find proper clusters    
    // re-build sets of edges and so on
    //  contigsize.................: <ctg_id, ctg_size>
    //  cluster_edge_updated.......: <clusterid, <contig1-contig2, correlation value> >
    //  cluster_vertex_size_updated: <clusterid, cluster-contig-length>  
    //  cluster_vertex_updated.....: <clusterid, <ctg_id, 1> >;  updated       
    //  mutually_exclusive_cluster.: <cluster_id1, cluster_id2>; cluster_id1/2 broken from the same old cluster id 
    // build cluster_edge_updated
    // map<int, map<string, int> >::iterator clitr;
    // map<int, map<string, int> >::iterator clitr_end;
    clitr     = (*cluster_vertex_updated).begin();
    clitr_end = (*cluster_vertex_updated).end();            
    while(clitr != clitr_end)
    {
        int cluster_id = (*clitr).first;
        map<string, int> this_cluster_vertex = (*clitr).second;
        // initiaze edge variable
        assert( (*cluster_edge_updated).find(cluster_id) == (*cluster_edge_updated).end() );
        map<string, double> this_cluster_edge;
        (*cluster_edge_updated).insert(std::pair<int, map<string, double> >(cluster_id, this_cluster_edge));
        // initialize total vertex size 
        assert( (*cluster_vertex_size_updated).find(cluster_id) == (*cluster_vertex_size_updated).end() );
        (*cluster_vertex_size_updated).insert(std::pair<int, unsigned long>(cluster_id, 0));
        //
        map<string, int>::iterator ctgitr1;
        map<string, int>::iterator ctgitr1_end;
        ctgitr1     = this_cluster_vertex.begin();
        ctgitr1_end = this_cluster_vertex.end();
        while(ctgitr1 != ctgitr1_end)
        {
            string this_ctg1 = (*ctgitr1).first;
            // update total vertex size 
            assert( contigsize.find(this_ctg1) != contigsize.end() );
            (*cluster_vertex_size_updated)[ cluster_id ] += contigsize[ this_ctg1 ];
            //
            map<string, BIPAT_INT>::iterator gt_citr1 = contigPollenGTSeq_int.find(this_ctg1);
            assert(gt_citr1 != contigPollenGTSeq_int.end() );
            BIPAT_INT pctg1;
            pctg1.leftPat_int  = (*gt_citr1).second.leftPat_int;
            pctg1.rightPat_int = (*gt_citr1).second.rightPat_int;             
            //
            map<string, int>::iterator ctgitr2;
            map<string, int>::iterator ctgitr2_end;
            ctgitr2     = this_cluster_vertex.begin();
            ctgitr2_end = this_cluster_vertex.end();    
            while(ctgitr2 != ctgitr2_end)
            {
                string this_ctg2 = (*ctgitr2).first;   
                if(this_ctg1.compare(this_ctg2) == 0)
                {
                    ctgitr2 ++;
                    continue;
                }
                //
                string this_edge_tmp1 = this_ctg1 + "-" + this_ctg2;                
                string this_edge_tmp2 = this_ctg2 + "-" + this_ctg1;
                if((*cluster_edge_updated)[cluster_id].find(this_edge_tmp1) != (*cluster_edge_updated)[cluster_id].end() || 
                   (*cluster_edge_updated)[cluster_id].find(this_edge_tmp2) != (*cluster_edge_updated)[cluster_id].end() )
                {
                    ctgitr2 ++;
                    continue;
                }
                //                
                map<string, BIPAT_INT>::iterator gt_citr2 = contigPollenGTSeq_int.find(this_ctg2);
                assert(gt_citr2 != contigPollenGTSeq_int.end() );
                BIPAT_INT pctg2;
                pctg2.leftPat_int  = (*gt_citr2).second.leftPat_int;
                pctg2.rightPat_int = (*gt_citr2).second.rightPat_int;                   
                // calculate correlation
                double cor[4];
                cor[0] = calculate_correlation(&pctg1.leftPat_int,  &pctg2.leftPat_int);      
                cor[1] = calculate_correlation(&pctg1.leftPat_int,  &pctg2.rightPat_int);          
                cor[2] = calculate_correlation(&pctg1.rightPat_int, &pctg2.leftPat_int);      
                cor[3] = calculate_correlation(&pctg1.rightPat_int, &pctg2.rightPat_int);            
                // update larger cluster with better correlation
                int best_i              = 0;
                double best_cluster_cor = 0;
                for(int i = 0; i < 4; i ++)
                {                     
                    if(cor[i]>best_cluster_cor)
                    {
                        best_cluster_cor = cor[i];
                        best_i           = i;                                                         
                    }
                }
                // update edge 
                if(best_cluster_cor > minCorscore)
                {
                    (*cluster_edge_updated)[cluster_id].insert(std::pair<string, double>(this_edge_tmp1, best_cluster_cor));
                }
                // next ctg2
                ctgitr2 ++;
            }
            // next ctg1
            ctgitr1 ++;
        }
        // next cluster 
        clitr ++;
    }
    // output updated raw clusters: 
    // this is processed raw, can be more than expected number of linkage groups thus needs further merging.
    string subclusterinfo = tmpfolder + "/s3_genotype_haplotig_GT_similarity_matrix_subclusters_raw_advanced.dot\0"; // 
    ofstream subcluster_ofp;
    subcluster_ofp.open(subclusterinfo.c_str(), ios::out);
    if(!subcluster_ofp.good())
    {
        cout   << "   Error: cannot open file " << subclusterinfo << endl;
        return false;
    } 
    subcluster_ofp << "/* Here are the subclusters of contigs with correlation based cluster broken */" << endl;
    subcluster_ofp << "graph\tGraph_1 {"                           << endl;
    // to collect all clusters-specific vertex
    map<int, map<string, double> >::iterator cleitr;
    map<int, map<string, double> >::iterator cleitr_end;
    cleitr     = (*cluster_edge_updated).begin();
    cleitr_end = (*cluster_edge_updated).end();
    unsigned long total_contig_size = 0;
    unsigned long total_contig_numb = 0;
    while(cleitr != cleitr_end)
    {
        subcluster_ofp << "\tsubgraph cluster_" << (*cleitr).first << " {" << endl;
        unsigned long cluster_ctg_size = 0;
        map<string, int> subcluster_vertex;
        //
        map<string, double>::iterator eitr;
        map<string, double>::iterator eitr_end;
        eitr     = (*cleitr).second.begin();
        eitr_end = (*cleitr).second.end();
        while(eitr != eitr_end)
        {
            vector<string> edgeinfo = split_string((*eitr).first, '-');
            subcluster_ofp << "\t"
            << edgeinfo[0] << " -- " << edgeinfo[1]
            << " "            << "[color=gold, penwidth=1, arrowsize=1, label=" << (*eitr).second << "];" 
            << endl;  
            //
            if(subcluster_vertex.find(edgeinfo[0]) == subcluster_vertex.end() )
            {
                subcluster_vertex.insert(std::pair<string, int>(edgeinfo[0], 1));
                assert(contigsize.find(edgeinfo[0]) != contigsize.end());
                cluster_ctg_size += contigsize[edgeinfo[0]];
            }
            if(subcluster_vertex.find(edgeinfo[1]) == subcluster_vertex.end() )
            {
                subcluster_vertex.insert(std::pair<string, int>(edgeinfo[1], 1));
                assert(contigsize.find(edgeinfo[1]) != contigsize.end());
                cluster_ctg_size += contigsize[edgeinfo[1]];
            }            
            //           
            eitr ++;
        }
        //
        subcluster_ofp << "\t/* "          << subcluster_vertex.size() << " contigs with total size of " 
                       << cluster_ctg_size << " bp */"              << endl;
        subcluster_ofp << "\t}" << endl;     
        //
        total_contig_numb += subcluster_vertex.size();
        total_contig_size += cluster_ctg_size;  
        // next cluster of edges
        cleitr ++;
    }
    //
    subcluster_ofp << "}" << endl;
    subcluster_ofp.close();    
    //
    cout << "   Info: after breaking clusters: " << total_contig_numb 
         << " contigs of "                       << total_contig_size 
         << " bp clustered into "                << (*cluster_edge_updated).size() 
         << " clusters. "                        << endl;          
    //       
    return true;
}  
// 
int find_contig_set_overlap(map<string, double> map1, map<string, double> map2)
{
    /*
        function: find how many contigs are common in the two maps
        map1/2: <contig-id, correlation>; map1/2 has the set of contigs show correlations to ctg1/2 respectively
                                        ; if map1 and map2 have many common ctgs, ctg1 and ctg2 belong to same cluster
    */
    int common_ctg_num = 0;
    map<string, double>::iterator map1_itr;
    map<string, double>::iterator map1_itr_end;
    map1_itr     = map1.begin();
    map1_itr_end = map1.end();
    while(map1_itr != map1_itr_end)
    {
        map<string, double>::iterator map2_itr = map2.find((*map1_itr).first);
        if(map2_itr != map2.end() )
        {
            common_ctg_num ++;
        }
        //        
        map1_itr ++;
    }
    //
    return common_ctg_num;
}
// 
bool merge_clusters(map<int, map<string, double> > cluster_edge, 
                    map<int, map<string, int> >    cluster_vertex,  
                    map<int, unsigned long>        cluster_vertex_size,
                    map<int, int>                  mutually_exclusive_cluster,
                    map<string, BIPAT_INT>         contigPollenGTSeq_int,  
                    map<string, unsigned long>     contigsize,
                    int                            N_cluster,
                    string                         tmpfolder,
                    map<int, int>*                 gb_group_id,
                    map<int, map<string, int> >*   gb_group)
{
    /*
        this function works on the output of function, or raw-but-advanced clusters: break_clusters(...) 
        cluster_edge..............: <clusterid, <contig1-contig2, correlation value> >
        cluster_vertex............: <clusterid, <ctg_id, 1> >
        cluster_vertex_size.......: <clusterid, cluster-contig-length> 
        mutually_exclusive_cluster: <clusterid1, clusterid2>; clusters to avoid being merged; see break_clusters(...)
        contigPollenGTSeq_int.....: <ctg_id, {leftPat, rightPat}>  
        contigsize................: <ctg_id, individual-contig-length>  
        N_cluster.................: expected number of linkage groups/clusters
        tmpfolder.................: output directory
        gb_group..................: <old_group_id, <ctg_id, 1> >
        gb_group_id...............: <old_group_id, new_group_id>
    */
    map<int, vector<int> > merging_record; // <clusterid, <merged-subcluster-ids> >    
    if(cluster_edge.size() < N_cluster)
    {
        cout << "   Error: number of current clusters "    << cluster_edge.size() 
             << " < expected cluster number " << N_cluster << endl;
        cout << "          Hint: you can slowly increase value of option cor_cutoff to get more clusters. "   
             << endl;
        return false;
    }else
    if(cluster_edge.size() == N_cluster)
    {
        cout << "   Info: number of current clusters = expected cluster number = " << N_cluster << endl;
        cout << "         No further merging needed. "     << endl;
        return true;
    }else
    {
        cout << "   Info: number of current clusters "     << cluster_edge.size() 
             << " > expected cluster number " << N_cluster << endl;
        cout << "         Smallest clusters will be merged into larger ones. "   
             << endl;
        cout << "   Info: merging clusters..."     << endl;
        while(cluster_edge.size() > N_cluster)
        {
            int smallest_cluster_id = find_smallest_cluster(cluster_vertex_size);            
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
            while(svitr != svitr_end)
            {
                // get a contig id in the smallest cluster
                string sctg = (*svitr).first;
                map<string, BIPAT_INT>::iterator sctg_citr = contigPollenGTSeq_int.find(sctg);
                assert(sctg_citr != contigPollenGTSeq_int.end() );
                BIPAT_INT psmall;
                psmall.leftPat_int  = (*sctg_citr).second.leftPat_int;
                psmall.rightPat_int = (*sctg_citr).second.rightPat_int;  
                // clusters of vertices
                map<int, map<string, int> >::iterator vclitr;
                map<int, map<string, int> >::iterator vclitr_end;
                vclitr     = cluster_vertex.begin();
                vclitr_end = cluster_vertex.end();
                while(vclitr != vclitr_end)
                {
                    int larger_cluster_id  = (*vclitr).first;
                    // check if the smaller and larger clusters should avoid being merged 
                    if(mutually_exclusive_cluster.find(smallest_cluster_id) != mutually_exclusive_cluster.end() && 
                       mutually_exclusive_cluster[smallest_cluster_id] == larger_cluster_id )
                    {
                        vclitr ++;
                        continue;
                    }
                    //
                    if( larger_cluster_id != smallest_cluster_id )
                    {
                        map<string, int> larger_cluster_vertex = (*vclitr).second;
                        map<string, int>::iterator lvitr;
                        map<string, int>::iterator lvitr_end;
                        lvitr     = larger_cluster_vertex.begin();
                        lvitr_end = larger_cluster_vertex.end();
                        while(lvitr != lvitr_end)
                        {
                            // get a contig id in current larger cluster                             
                            string lctg = (*lvitr).first;      
                            map<string, BIPAT_INT>::iterator lctg_citr = contigPollenGTSeq_int.find(lctg);
                            assert(lctg_citr != contigPollenGTSeq_int.end() );
                            BIPAT_INT plarger;
                            plarger.leftPat_int  = (*lctg_citr).second.leftPat_int;
                            plarger.rightPat_int = (*lctg_citr).second.rightPat_int;  
                            // calculate correlation 
                            double cor[4];
                            cor[0] = calculate_correlation(&psmall.leftPat_int,  &plarger.leftPat_int);      
                            cor[1] = calculate_correlation(&psmall.leftPat_int,  &plarger.rightPat_int);          
                            cor[2] = calculate_correlation(&psmall.rightPat_int, &plarger.leftPat_int);      
                            cor[3] = calculate_correlation(&psmall.rightPat_int, &plarger.rightPat_int);            
                            // update larger cluster with better correlation
                            for(int i = 0; i < 4; i ++)
                            {
                                if(cor[i]>best_cluster_cor)
                                {
                                    best_cluster_cor = cor[i];
                                    best_cluster_id  = larger_cluster_id;
                                    best_connect_edge= sctg + "-" + lctg;
                                }
                            }
                            // next contig id in current larger cluster 
                            lvitr ++;
                        }
                    }
                    // next larger cluster 
                    vclitr ++;
                }
                //
                if(0)
                {
                    // this does not work as expected so not used 20200714.
                    if(best_larger_cluster_id.find(best_cluster_id) == best_larger_cluster_id.end())
                    {
                        // contig in smaller cluster showing best cor to larger cluster
                        vector<string> sctg_id;
                        sctg_id.push_back(sctg);
                        best_larger_cluster_id.insert(std::pair<int, vector<string> >(best_cluster_id, sctg_id));
                        // best correlation of contigs in smaller and larger clusters
                        vector<double> slcluster_cor;
                        slcluster_cor.push_back(best_cluster_cor);
                        best_larger_cluster_cor.insert(std::pair<int, vector<double> >(best_cluster_id, slcluster_cor));
                        // edge with best correlation connecting smaller-larger clusters
                        vector<string> slcluster_edge;
                        slcluster_edge.push_back(best_connect_edge);
                        best_larger_cluster_edge.insert(std::pair<int, vector<string> >(best_cluster_id, slcluster_edge));
                    }else
                    {
                        // update 
                        best_larger_cluster_id[best_cluster_id].push_back(sctg);
                        best_larger_cluster_cor[best_cluster_id].push_back(best_cluster_cor);
                        best_larger_cluster_edge[best_cluster_id].push_back(best_connect_edge);
                    }
                }
                // next contig id in the smallest cluster 
                svitr ++;
            }
            // 
            if(0)
            {   
                // this does not work as expected so not used 20200714.
                map<int, vector<string> >::iterator bestmitr;
                map<int, vector<string> >::iterator bestmitr_end;
                bestmitr     = best_larger_cluster_id.begin();
                bestmitr_end = best_larger_cluster_id.end();
                int major_support = 0;; // how many contigs from the smaller cluster connecting to this larger cluster
                while(bestmitr != bestmitr_end)
                {
                    if( ( (*bestmitr).second ).size() > major_support )
                    {
                        best_cluster_id = (*bestmitr).first;
                        // find best edge with cor
                        assert( best_larger_cluster_cor.find(best_cluster_id)  != best_larger_cluster_cor.end() );
                        assert( best_larger_cluster_edge.find(best_cluster_id) != best_larger_cluster_edge.end() );
                        best_cluster_cor= 0;
                        vector<double> slcluster_cor  = best_larger_cluster_cor[best_cluster_id];
                        vector<string> slcluster_edge = best_larger_cluster_edge[best_cluster_id];
                        for(int coi = 0; coi < slcluster_cor.size(); coi++)
                        {
                            if(slcluster_cor[coi] > best_cluster_cor)
                            {
                                best_cluster_cor  = slcluster_cor[coi];
                                best_connect_edge = slcluster_edge[coi];
                            }else
                            if(slcluster_cor[coi] == best_cluster_cor)
                            {
                                cout << "   Warning: equal cor score found, but not considered. " << endl;
                            }else ;
                        }
                        //
                        major_support   = ( (*bestmitr).second ).size();
                    }
                    bestmitr ++;
                }
            }
            // merge the smaller cluster into a larger group according to best correlation value
            cout << "   Info: cluster "              << smallest_cluster_id 
                 << " will be merged into cluster "  << best_cluster_id << endl;
            // update edge info 
            cout << "   Info: cluster " 
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
            cout << "   Info: cluster " 
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
            cout << "   Info: cluster " 
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
        }
        cout << "   Info: merging clusters done. " << endl;
    } 
    //    
    // output subclusters: this is final, merged from raw subclusters thus equals expected number of numbers.
    string subclusterinfo = tmpfolder + "/s3_genotype_haplotig_GT_similarity_matrix_subclusters_final.dot\0"; // 
    ofstream subcluster_ofp;
    subcluster_ofp.open(subclusterinfo.c_str(), ios::out);
    if(!subcluster_ofp.good())
    {
        cout   << "   Error: cannot open file " << subclusterinfo << endl;
        return false;
    } 
    subcluster_ofp << "/* Here are the final subclusters of contigs */" << endl;
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
                assert(contigsize.find(edgeinfo[0]) != contigsize.end());
                cluster_ctg_size += contigsize[edgeinfo[0]];
            }
            if(subcluster_vertex.find(edgeinfo[1]) == subcluster_vertex.end() )
            {
                subcluster_vertex.insert(std::pair<string, int>(edgeinfo[1], 1));
                assert(contigsize.find(edgeinfo[1]) != contigsize.end());
                cluster_ctg_size += contigsize[edgeinfo[1]];
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
        int gb_grpid = (*clitr).first;
        assert((*gb_group).find(gb_grpid) == (*gb_group).end() );
        (*gb_group).insert(std::pair<int, map<string, int> >(gb_grpid, subcluster_vertex));   
        (*gb_group_id).insert( std::pair<int, int>( (*clitr).first, count_ci) ); // <old_group_id:to_dot, new_group_id>    
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
int find_smallest_cluster(map<int, unsigned long> cluster_vertex_size)
{
    int           smallest_cluster_id   = 0;
    unsigned long smallest_cluster_size = 999999999;
    map<int, unsigned long>::iterator citr;
    map<int, unsigned long>::iterator citr_end;
    citr     = cluster_vertex_size.begin();
    citr_end = cluster_vertex_size.end();
    while(citr != citr_end)
    {
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
bool get_clusters(map<string, int> vertex,
                  map<string, vector<string> > edge,
                  map<string, double> cor,
                  map<int, map<string, double> >* cluster_edge
                 )
{
    // function: traverse contig1-contig2 graph, connected vertices will be collected as one group = potential LGs.
    //
    // vertex, edge and correlation
    // map<string, int> vertex;           // <contig, 1>
    // map<string, vector<string> > edge; // <contig, <list-contig> >
    // map<string, double> cor;           // <contig1-contig2, correlation.value >
    // map<int, map<string, double> > cluster_edge;   // <clusterid, <contig1-contig2, correlation.value> >
    map<int, map<string, int> > cluster_vertex; // <clusterid, <contig, 1> >
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
string get_break_pos(string             depth_gt, 
                     string             pmflag,
                     int                pollenid,
                     unsigned long*     leftp, 
                     double*            score,
                     std::stringstream* ss,
                     string             contigid,
                     bool               output_ss)
{
    // used after 20200604
    // depth_gt: is the pollen depth-defined genotype sequence at haplotig markers
    // pmflag  : reserved, not really used yet.
    // leftp   : would be the i-th   marker position along the haplotig
    // rightp  : would be the i+1-th marker position along the haplotig
    // score   : would be score for the reported co
    // minimum number of markers hard-coded: at least 3
    //
    string pmstring    = "";
    string poPatbk     = depth_gt;
    // note the smoothing is only for genotype sequences >=10 bits
    string depth_gt_sm = smooth_depth_gt(depth_gt);
    if(depth_gt_sm.size()>=10 && false)
    cout << "   Check2rm: " << depth_gt    << " raw " << endl
         << "   Check2rm: " << depth_gt_sm << " smt " << endl;
    // bool output_ss = false;
    // remove non-effective marekers with value: 'u': not existing with current genotype definition
    poPatbk.erase(std::remove(poPatbk.begin(), poPatbk.end(), 'u'), poPatbk.end());
    if(poPatbk.size()==0) 
    {
        pmstring = "uu";// not sure on co; this should not never happen with current genotype definition
        *leftp   = 0;
        *score   = 0;
        cout << "   Check2rm: " << depth_gt       << " ---- determined as " 
             << pmstring        << " for pollen " << pollenid 
             << "_x\t0_"        << "U" 
             << " (no effective markers) "        << endl; 
        if(output_ss)             
        *ss  << "   Check2rm: " << depth_gt       << " ---- determined as " 
             << pmstring        << " for pollen " << pollenid 
             << "_x\t0_"        << "U" 
             << " (no effective markers) "        << endl;
        return pmstring;
    }
    if(poPatbk.size()==1)
    {
        if(poPatbk.compare("0")==0)
        {
            pmstring = "00";
            pmflag   = "M";
        }else
        if(poPatbk.compare("1")==0)
        {
            pmstring = "11";
            pmflag   = "P";            
        }else ;
        *leftp   = 0; // no co
        *score   = 0; // no co     
        cout << "   Check2rm: " << depth_gt       << " ---- determined as "  
             << pmstring        << " for pollen " << pollenid 
             << "_x\t0_"        << pmflag 
             << " ( 1 effective marker) "         << endl; 
        if(output_ss)                          
        *ss  << "   Check2rm: " << depth_gt       << " ---- determined as "  
             << pmstring        << " for pollen " << pollenid 
             << "_x\t0_"        << pmflag 
             << " ( 1 effective marker) "         << endl;              
        return pmstring;  
    }
    //
    *score        = 0.0;
    size_t c1Lef  = 0;   // Left  P cnt for highest score
    size_t c2Lef  = 0;   // Left  M cnt for highest score
    size_t c1Rig  = 0;   // right P cnt for highest score
    size_t c2Rig  = 0;   // right M cnt for highest score
    // later control on number of markers, and no need checking all pos.
    for(size_t mi = 1; mi <= depth_gt_sm.size()-1; mi ++) 
    {
        string strLef = depth_gt_sm.substr(0, mi); // val @ pos of 0,...,mi-1
        string strRig = depth_gt_sm.substr(mi);    // val @ pos of mi,mi+1,...
        // current break point: note: 'u' excluded from consideration
        size_t c1Leftmp  = std::count(strLef.begin(), strLef.end(), '1'); // Left  P:1
        size_t c2Leftmp  = std::count(strLef.begin(), strLef.end(), '0'); // Left  M:0: note: 0 represent for 3 haps.
        size_t c1Rigtmp  = std::count(strRig.begin(), strRig.end(), '1'); // Right P:1
        size_t c2Rigtmp  = std::count(strRig.begin(), strRig.end(), '0'); // Right M:0
        double af1Leftmp = (double)c1Leftmp*1.0/(c1Leftmp+c2Leftmp);
        double af2Leftmp = (double)c2Leftmp*1.0/(c1Leftmp+c2Leftmp);        
        double af1Rigtmp = (double)c1Rigtmp*1.0/(c1Rigtmp+c2Rigtmp);
        double af2Rigtmp = (double)c2Rigtmp*1.0/(c1Rigtmp+c2Rigtmp);
        // left as P right as M
        double score1 = af1Leftmp * af2Rigtmp; 
        // left as M right as P
        double score2 = af2Leftmp * af1Rigtmp; 
        double scotmp = score1>score2?score1:score2;   
        if(scotmp >= *score)
        {
            *score  = scotmp;
            // (*leftp)-th marker is the left break marker, (*leftp+1) is the right break marker; 1-based positioning
            *leftp  = mi;  
            if(score1>score2)
            { 
                pmstring  = "10";
                pmflag    = "P";                
            }else
            {
                pmstring  = "01";
                pmflag    = "P";                
            }  
            //
            c1Lef  = c1Leftmp;
            c2Lef  = c2Leftmp;
            c1Rig  = c1Rigtmp;
            c2Rig  = c2Rigtmp;                     
        }
    }
    // 
    string checkfrom = "";
    // control on marker number: MuMMMMuuuuuuPuPPPu
    if(*score > 0 && poPatbk.size()<10)
    {
            // found "1" on left end;            
            if(pmstring.compare("10") == 0)
            { 
               // when considering "0", its number is required larger; as "0" might be from miss-sequencing.
               // at least 1/3 "1"; first position was "1" + 
               // at least 1/3 "0"; last  position was "0"
               if(depth_gt_sm.compare("010") == 0)
               {
                    pmstring  = "11"; 
                    pmflag    = "P";                      
                    checkfrom = "1.1";                   
               }else
               if( (poPatbk.substr(0, 1)).compare("1") == 0 && 
                   c1Lef >= poPatbk.size()/3.0              &&
                   (poPatbk.substr(poPatbk.size()-1, 1)).compare("0") == 0 &&
                   c2Rig >  poPatbk.size()/3.0 && c2Rig>=3                 
                 )
                {
                    checkfrom = "1.2";
                }else
                {
                    if(c1Lef+c1Rig >= c2Lef+c2Rig)
                    {
                        // as long as there is "1", assume this region is fully covered; "0": missed-sequencing
                        pmstring  = "11"; 
                        pmflag    = "P";  
                    }else
                    {
                        pmstring  = "10"; 
                        pmflag    = "P";                          
                    }
                    *leftp    = 0; // no co    
                    //*score  = 0; // no co
                    checkfrom = "2";
                }
            }
            // found "1" on right end
            if(pmstring.compare("01") == 0)
            {
               // at least 1/3 "0"; first position was "0" + at least three "0"
               // at least 1/3 "1"; last  position was "1"
               if( (poPatbk.substr(0, 1)).compare("0") == 0 && 
                   c2Lef >  poPatbk.size()/3.0 && c2Lef>=3  &&
                   (poPatbk.substr(poPatbk.size()-1, 1)).compare("1") == 0 &&
                   c1Rig >= poPatbk.size()/3.0
                 )
                {
                    checkfrom = "3";
                }else
                {
                    // as long as there is "1", assume this region is fully covered; "0": missed-sequencing
                    if(c1Lef+c1Rig >= c2Lef+c2Rig)
                    {
                        // as long as there is "1", assume this region is fully covered; "0": missed-sequencing
                        pmstring  = "11"; 
                        pmflag    = "P";  
                    }else
                    {
                        pmstring  = "01"; 
                        pmflag    = "P";                          
                    }
                    *leftp    = 0; // no co
                    //*score  = 0; // no co 
                    checkfrom = "4";
                }
            }
    }else
    if(*score>=minCOscore && poPatbk.size()>10)
    {
        // found "1" on left end
        if(pmstring.compare("10") == 0)
        {            
            // at least 3 effective "1" plus at least 6 effective "0" (as "0" can be due to missed-sequencing)            
            if( c1Lef >= 3 && c2Rig >= 6 )
            {
                checkfrom = "9";
            }else
            {
                if(c2Rig > c1Lef)
                {
                    pmstring = "00";
                    pmflag   = "M";                    
                }else
                {
                    pmstring = "11";
                    pmflag   = "P";                    
                }
                ////*leftp   = 0;  // no co
                //*score  = 0; // no co
                checkfrom = "10";
            }          
        }
        // found "1" on right end
        if(pmstring.compare("01") == 0)
        {
            // at least 6 effective "0" plus at least 3 effective "1" 
            // (poPatbk.substr(0, 1)).compare("0") == 0 &&   
            // (poPatbk.substr(poPatbk.size()-1, 1)).compare("1") == 0     
            if( c2Lef >= 6 && c1Rig >= 3 )
            {
                checkfrom = "11";
            }else
            {
                if(c2Lef > c1Rig)
                {
                    pmstring = "00";
                    pmflag   = "M";
                }else
                {
                    pmstring = "11";
                    pmflag   = "P";
                }
                ////*leftp   = 0;  // no co
                //*score  = 0; // no co
                checkfrom = "10";
            } 
        }
    }else
    {
        if(poPatbk.size()<=12 && c1Rig + c1Lef >= poPatbk.size()/4.0)
        {
            pmstring = "11";
            pmflag   = "P";            
        }else
        if(c1Rig + c1Lef >= poPatbk.size()/3.0)
        {
            pmstring = "11";
            pmflag   = "P";            
        }else
        if(c2Rig + c2Lef > c1Rig + c1Lef)
        {
            pmstring = "00";
            pmflag   = "M";
        }else
        {
            pmstring = "11";
            pmflag   = "P";            
        }
        checkfrom = "16";
        *leftp   = 0; // no co
        //*score = 0; // no co  
    }
    //
    cout << "   Check2rm: " << depth_gt       << " ---- determined as "  
         << pmstring        << " for pollen " << pollenid 
         << "_x\t0_"        << pmflag 
         << " (CO after "   << *leftp         << "th marker: score "     
         << *score          << " - case "     << checkfrom << " ) "; // *** 
    if(output_ss)                      
    *ss  << "   Check2rm: " << depth_gt       << " ---- determined as "  
         << pmstring        << " for pollen " << pollenid 
         << "_x\t0_"        << pmflag 
         << " (CO after "   << *leftp         << "th marker: score "     
         << *score          << " - case "     << checkfrom << " ) "; // ***          
    // output co interval
    if(*leftp != 0)
    {
        int leftindex  = *leftp-1;
        while(leftindex>0)
        {
            if(depth_gt.substr(leftindex, 1).compare("u") == 0)
            {
                leftindex --;
            }else
            {
                break;
            }
        }
        int rightindex = *leftp;
        while(rightindex<depth_gt.size()-1)
        {
            if(depth_gt.substr(rightindex, 1).compare("u") == 0)
            {
                rightindex ++;
            }else
            {
                break;
            }
        }
        // output CO real positions of markers: increase to 1-based positioning <= leftindex+1; real position in depth_gt
        unsigned long leftmarker  = 0; // not really used
        unsigned long rightmarker = 0; // not really used
        cout << " (CO interval " 
             << contigid 
             << " "
             << leftindex+1  << "th.mkr=" << leftmarker 
             << "-"
             << rightindex+1 << "th.mkr=" << rightmarker
             << ", size="
             << rightmarker - leftmarker + 50000 << ")"
             << endl;
        if(output_ss)                          
        *ss  << " (CO interval " 
             << contigid 
             << " "
             << leftindex+1  << "th.mkr=" << leftmarker 
             << "-"
             << rightindex+1 << "th.mkr=" << rightmarker
             << ", size="
             << rightmarker - leftmarker + 50000 << ")"
             << endl;             
    }else
    {
        cout << endl; // to end ***
        if(output_ss)        
        *ss  << endl;
    }
    return pmstring;
}
//
string smooth_depth_gt(string depth_gt)
{
    string depth_gt_smoothed = depth_gt;
    if(depth_gt.size()<10)
    {
        return depth_gt_smoothed;
    }
    size_t pos = depth_gt.find("101", 0);    
    while(pos != std::string::npos)
    {
        depth_gt_smoothed.replace(pos+1, 1, "1");
        pos = depth_gt.find("101", pos+1);
    }
    return depth_gt_smoothed;
}
//
double calculate_correlation(vector<int>* marker_x, vector<int>* marker_y) 
{
    /*
        cor_xy <- sum[ (xi-mean_x) * (yi-mean_y) ] / [ sqrt( sum(xi-mean_x)^2 ) * sqrt( sum(yi-mean_y)^2 ) ] 
    */
    double cor_xy = 0.0;
    double mean_x = 0.0;
    double mean_y = 0.0;
    double len_x  = (*marker_x).size();
    double len_y  = (*marker_y).size();
    assert(len_x == len_y && len_x > 0);
    for(vector<int>::iterator it = (*marker_x).begin(); it != (*marker_x).end(); ++it)
    {
        mean_x += *it;
    }
    if(mean_x < 3)
    {
        cor_xy = -100.0; // indicating inaccuracy     
        return cor_xy;   
    }
    mean_x = mean_x/len_x;
    for(vector<int>::iterator it = (*marker_y).begin(); it != (*marker_y).end(); ++it)
    {
        mean_y += *it;
    }
    if(mean_y < 3)
    {
        cor_xy = -100.0; // indicating inaccuracy 
        return cor_xy;               
    }    
    mean_y = mean_y/len_y;
    //
    double numerator    = 0;
    double denominatorx = 0;
    double denominatory = 0;   
    for(int i = 0; i < len_x; i ++)
    {
        numerator    += ((*marker_x)[i]-mean_x)*( (*marker_y)[i]-mean_y);
        denominatorx += pow( (*marker_x)[i]-mean_x, 2);
        denominatory += pow( (*marker_y)[i]-mean_y, 2);
    }
    //
    if(denominatorx==0 || denominatory==0)
    {
        //cout   << "   Warning: denominator* == 0. " << endl;
        cor_xy = -100.0; // indicating inaccuracy
    }else
    {
        cor_xy = numerator / ( sqrt(denominatorx) * sqrt(denominatory) );
    }
    //
    return cor_xy;
}
//
double find_haploid_mean_depth(map<string, NODE> this_pollen_depth)
{
    double this_mean = 0.0; // most likely with our depth
    //
    double cutoffreadcnt = 4*(this_pollen_depth["#total_reads"].depth*1.0/557309) * 1.0; 
    cout << "   Info: raw but scaled cutoff to determine haplotig mean reads: [" << cutoffreadcnt << ", ..)" << endl;
    //
    map<string, NODE>::iterator mitr;
    map<string, NODE>::iterator mitr_end;
    mitr     = this_pollen_depth.begin();
    mitr_end = this_pollen_depth.end();    
    int dused= 0;
    while(mitr != mitr_end)
    {
        if(((*mitr).second).type.compare("hap") == 0)
        {
            double this_depth_nor = ((*mitr).second).depth_nor; 
            double this_depth_raw = ((*mitr).second).depth;        
            if(this_depth_raw > cutoffreadcnt)
            {
                this_mean += this_depth_nor;
                dused ++;
            }
        }
        mitr ++;  
    }
    //
    this_mean = this_mean/dused;
    //
    return this_mean;
}
//
bool normalize_read_count(map<string, NODE>* this_pollen_depth)
{
    if((*this_pollen_depth).size()==0)
    {
        return false;
    }
    //
    unsigned long total_reads = (*this_pollen_depth)["#total_reads"].depth;
    assert(total_reads > 0);
    //
    map<string, NODE>::iterator mitr;
    map<string, NODE>::iterator mitr_end;
    mitr     = (*this_pollen_depth).begin();
    mitr_end = (*this_pollen_depth).end();    
    while(mitr != mitr_end)
    {
        ((*mitr).second).depth_nor = ((*mitr).second).depth * 1.0 / 
                                     ((*mitr).second).size  * 10000 /
                                     total_reads * 1000000; // rp10km
        mitr ++;  
    }
    //
    return true;
}
//
bool read_current_bed(string                      this_bed, 
                      map<string, NODE>*          this_pollen_depth,
                      map<string, unsigned long>* ctg_size_all,
                      vector<string>*             marker_order)
{
    ifstream ifp;
    ifp.open(this_bed.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Warning: skipping this file bed file as it coud not be opended. " << endl;
        return false;
    }    
    //
    unsigned long total_reads = 0;
    unsigned long len = 0;        
    string        lstr="";             
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;            
        vector<string> lineinfo = split_string(line, '\t');
        //
        string              key = lineinfo[0] + "\t" + lineinfo[1] + "\t" + lineinfo[2]; // ctg_id\tstart\tend
        unsigned long this_size = strtoul(lineinfo[1].c_str(), NULL, 0);
        this_size = strtoul(lineinfo[2].c_str(), NULL, 0) - this_size + 1;
        //
        if((*ctg_size_all).find(lineinfo[0]) == (*ctg_size_all).end())
        {
            len  = strtoul(lineinfo[9].c_str(), NULL, 0);  
            lstr = lineinfo[9];            
            (*ctg_size_all).insert(std::pair<string, unsigned long>(lineinfo[0], len));
            // cout << "   check: ctg " << lineinfo[0] << " size  collected as " << len << endl; 
        }
        //
        NODE   val;
        val.depth    = strtoul(lineinfo[3].c_str(), NULL, 0);
        val.type     = lineinfo[10];   
        val.size     = this_size;   
        val.ctg_size = len;     
        (*this_pollen_depth).insert(std::pair<string, NODE>(key, val));
        (*marker_order).push_back(key);
        total_reads += val.depth;
    }
    // total read count in this pollen for later normalization
    NODE dtmp;
    dtmp.depth    = total_reads;
    dtmp.type     = "info";
    dtmp.size     = len;
    dtmp.ctg_size = len;
    (*this_pollen_depth).insert(std::pair<string, NODE>("#total_reads", dtmp)); 
    //
    return true;
}
//
bool get_individual_bed_file(string list_bed_files, vector<string>* bedfiles)
{
    ifstream ifp;
    ifp.open(list_bed_files.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << list_bed_files << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;        
        (*bedfiles).push_back(line);
    }   
    ifp.close();
    return true;
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
bool get_pure_haplotig_markers(string                     bed_file, 
                               map<string, unsigned long> contigsize, 
                               map<string, double>*       hap_win_ratio)
{
    /*
        calculate ratio of hap-wins along a contig 
        ...
	utg000014l_pilon    1   50000   0   utg000014l_pilon    1   50000   126.87  138 2441410 hap
	utg000014l_pilon    50001   100000  0   utg000014l_pilon    50001   100000  126.87  138 2441410 hap
	utg000014l_pilon    100001  150000  0   utg000014l_pilon    100001  150000  126.87  138 2441410 hap 
	utg000014l_pilon    1350001 1380000 0   utg000014l_pilon    1350001 1380000 126.87  138 2441410 hap
	utg000014l_pilon    1380001 1390000 0   utg000014l_pilon    1380001 1390000 179 1   2441410 dip
	utg000014l_pilon    1390001 1440000 0   utg000014l_pilon    1390001 1440000 122.415 106 2441410 hap	
	...       
    */    
    ifstream ifp;
    ifp.open(bed_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open file " << bed_file << endl;
        return false;
    }
    map<string, double> hap_marker;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 11)
        {
            cout << "   Warning: skipped insufficient line " << line << endl;
            continue;
        }
        if(lineinfo[10].find("hap") != std::string::npos)
        {
            double len = atof(lineinfo[2].c_str()) - atof(lineinfo[1].c_str()) + 1;
            if( hap_marker.find(lineinfo[0]) != hap_marker.end() )
            {
                hap_marker[lineinfo[0]] += len;
            }else
            {
                hap_marker.insert(std::pair<string, double>(lineinfo[0], len));
            }
        }else
        {
            double len = 0.0;
            if( hap_marker.find(lineinfo[0]) != hap_marker.end() )
            {
                hap_marker[ lineinfo[0] ] += len;
            }else
            {
                hap_marker.insert(std::pair<string, double>(lineinfo[0], len));
            }            
        }        
    }
    ifp.close();
    // find highly confident haplotig markers 
    map<string, double>::iterator witr;
    map<string, double>::iterator witr_end;
    witr     = hap_marker.begin();
    witr_end = hap_marker.end();
    while(witr != witr_end)
    {
        string ctg_id = (*witr).first;
        assert(contigsize.find(ctg_id) != contigsize.end());
        double ratio = (*witr).second / contigsize[ctg_id];
        (*hap_win_ratio).insert(std::pair<string, double>(ctg_id, ratio));        
        if(ratio >= haplotig_win_ratio)
        {
            cout << "   check: " << ctg_id << " selected as haplotig marker, hap-ratio = " << ratio << endl;
        }else
        {
            cout << "   check: " << ctg_id << " selected as non-haplotig marker, hap-ratio = " << ratio << endl;
        }
        //
        witr ++;
    }
    //
    return true;
}
//
bool collect_contig_size(string ctgsize_file, map<string, unsigned long>* contigsize)
{
    cout << "   Info: reading contig size info from " << ctgsize_file << endl;
    ifstream ifp;
    ifp.open(ctgsize_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open contig size file " << ctgsize_file << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size() == 0) continue;
        if(line[0]=='#')     continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size()<2) continue;
        string        ctg_id   = lineinfo[0];
        unsigned long ctgsize = strtoul(lineinfo[1].c_str(), NULL, 0);
        (*contigsize).insert(std::pair<string, unsigned long>(ctg_id, ctgsize));
        //// cout << ctg_id << "\t" << ctgsize << endl;
    }
    if((*contigsize).size() == 0)
    {
        cout << "   Erorr: no info collected from contig size file. " << endl;
        return false;
    }else
    {
        cout << "   Info: " << (*contigsize).size() << " contig size info collected. " << endl;
    }
    ifp.close();
    //
    return true;
}
