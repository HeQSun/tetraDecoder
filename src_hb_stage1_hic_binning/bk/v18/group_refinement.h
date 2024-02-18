/* this function refines initial hic-based linkage grouping among windows of each contig */
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
bool group_refinement(string hbmarker_file, 
                      string out_folder, 
                      map<string, map<unsigned long, map<unsigned long, WHIC> > > hic_matrix,
                      map<string, unsigned long> target_contig_size);
