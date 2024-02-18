/* this tool  

       given a .paf file, is aimed to 
         group contigs into reference linkage groups.
         
   Written by Hequan Sun, MPIPZ
   Email: sunhequan@gmail.com/sun@mpipz.mpg.de
   Date.: 20220826
   
   update: 20220905 unique size of query (excluding multiple alignments of the same query region) 
   
*/
#include        <stdio.h>
#include       <string.h>
#include       <stdlib.h>
#include        <fstream>
#include       <iostream>
#include        <iomanip>
#include        <sstream>
#include         <vector>
#include            <map>
#include         <time.h>
#include       <assert.h>
#include "split_string.h"
//
//
bool read_paf(string paf_file, 
              map<string, unsigned long>* ctg_size, 
              map<string, map<string, map<unsigned long, unsigned long> > >* aln,
              map<string, map<string, map<unsigned long, unsigned long> > >* aln_query);
bool find_best_lg(map<string, unsigned long> ctg_size, 
                  map<string, map<string, map<unsigned long, unsigned long> > > aln,
                  map<string, unsigned long> chrs,
                  map<string, string>* ctg_to_lg,
                  map<string, map<string, unsigned long> >* lg_with_ctg,
                  map<string, unsigned long>* lg_grouped_ctg_size );  
bool read_chrid(string chr_file, map<string, unsigned long>* chrs);                       
//
int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        printf("\nFunction: given a .paf file, group contigs into reference linkage groups.\n");
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl << endl;
        cout << "Usage: ref_linkage_grouper ref_qry.paf chr.list" << endl;        
        cout << "   chr.list gives the list of chr ids to calculate final alignment stats.\n\n";
        return 1;
    }
    time_t ftime;
    struct tm* tinfo;
    time(&ftime);
    tinfo = localtime(&ftime);
    printf("\nGrouping contigs started on %s\n", asctime(tinfo));    
    // s0. get inputs
    string paf_file = (string)argv[1];
    string chr_file = (string)argv[2];
    // s1. read chr id 
    map<string, unsigned long> chrs;
    if(!read_chrid(chr_file, &chrs))
    {
        cout << "   Error: failed in reading the chr id file. " << endl;
        return 1;
    }
    // s2. read paf
    // aln......: <qry_ctg, <ref_lg, <ref_chr_sta, ref_chr_end> > > // no overlap allowed between {ref_chr_sta, ref_chr_end}
    // aln_query: <qry_ctg, <qry_lg, <qry_chr_sta, qry_chr_end> > > // no overlap allowed between {qry_chr_sta, qry_chr_end}
    map<string, unsigned long> ctg_size;
    map<string, map<string, map<unsigned long, unsigned long> > > aln; 
    map<string, map<string, map<unsigned long, unsigned long> > > aln_query;
    if( !read_paf(paf_file, &ctg_size, &aln, &aln_query) )
    {
        cout << "   Error: failed in reading the paf file. " << endl;
        return 1;
    } 
    // s2. find best linkage group for each contig 
    map<string, string> ctg_to_lg;
    map<string, map<string, unsigned long> > lg_with_ctg;
    map<string, unsigned long> lg_grouped_ctg_size;
    if(!find_best_lg(ctg_size, aln, chrs, &ctg_to_lg, &lg_with_ctg, &lg_grouped_ctg_size) )
    {
        cout << "   Error: failed in finding best linkage groups. " << endl;
        return 1;
    }
    //
    time_t endftime;
    struct tm* endtinfo;
    time(&endftime);
    endtinfo = localtime(&endftime);
    printf("\nGrouping contigs successfully finished on %s\n", asctime(endtinfo));
    return 0;
}
//
bool find_best_lg(map<string, unsigned long> ctg_size, 
                  map<string, map<string, map<unsigned long, unsigned long> > > aln,
                  map<string, unsigned long> chrs,
                  map<string, string>* ctg_to_lg,
                  map<string, map<string, unsigned long> >* lg_with_ctg,
                  map<string, unsigned long>* lg_grouped_ctg_size)
{
    /* for each ctg, find the best matching linkage group according to alignment stats */
    // ctg_to_lg............. <ctg, ref_lg>
    // lg_with_ctg........... <ref_lg, <ctg, ctg_size> >
    // lg_grouped_ctg_size... <ref_lg, total_ctg_size >
    // chrs.................. <chrid, 1/chrsizes>
    // aln................... <qry_ctg, <ref_lg, <ref_chr_sta, ref_chr_end> > >
    //
    map<string, map<string, map<unsigned long, unsigned long> > >::iterator citr;
    map<string, map<string, map<unsigned long, unsigned long> > >::iterator citr_end;
    citr     = aln.begin();
    citr_end = aln.end();
    while(citr != citr_end)
    {
        string this_ctg_id = (*citr).first;
        map<string, map<unsigned long, unsigned long> > this_ctg_aln = (*citr).second;
        //
        map<string, unsigned long> cov_len; // <chr/scaffold, bp_covering_this_lg>
        string        best_ref_lg = "";
        unsigned long best_ref_cov = 0;
        //
        map<string, map<unsigned long, unsigned long> >::iterator litr;
        map<string, map<unsigned long, unsigned long> >::iterator litr_end;
        litr     = this_ctg_aln.begin();
        litr_end = this_ctg_aln.end();
        while(litr != litr_end)
        {
            string        this_ref_id = (*litr).first;
            unsigned long this_ref_cov = 0;            
            map<unsigned long, unsigned long> tmp_pos = (*litr).second;
            map<unsigned long, unsigned long>::iterator pitr;
            map<unsigned long, unsigned long>::iterator pitr_end;  
            pitr     = tmp_pos.begin();
            pitr_end = tmp_pos.end();
            while(pitr != pitr_end)
            {
                this_ref_cov += (*pitr).second - (*pitr).first + 1;
                //
                pitr ++;
            }
            if(this_ref_cov > best_ref_cov)
            {
                best_ref_cov = this_ref_cov;
                best_ref_lg  = this_ref_id;
            }
            //
            litr ++;
        }
        //
        assert((*ctg_to_lg).find(this_ctg_id) == (*ctg_to_lg).end());
        (*ctg_to_lg).insert(std::pair<string, string>(this_ctg_id, best_ref_lg));
        cout << "   res: "      << this_ctg_id 
             << " size "        << ctg_size[this_ctg_id]
             << " assigned to " << best_ref_lg 
             << " with cov of " << best_ref_cov 
             << endl;
        //
        if((*lg_with_ctg).find(best_ref_lg) != (*lg_with_ctg).end())
        {
            assert( (*lg_with_ctg)[best_ref_lg].find(this_ctg_id) == (*lg_with_ctg)[best_ref_lg].end() );
            (*lg_with_ctg)[best_ref_lg].insert(std::pair<string, unsigned long>(this_ctg_id, ctg_size[this_ctg_id] ));
            (*lg_grouped_ctg_size)[best_ref_lg] += ctg_size[this_ctg_id];
        }else
        {
            map<string, unsigned long> tmp_ctg_info;
            tmp_ctg_info.insert(std::pair<string, unsigned long>(this_ctg_id, ctg_size[this_ctg_id] ));
            (*lg_with_ctg).insert(std::pair<string, map<string, unsigned long> >(best_ref_lg, tmp_ctg_info));
            //
            (*lg_grouped_ctg_size).insert(std::pair<string, unsigned long>(best_ref_lg, ctg_size[this_ctg_id] ));
        }
        //
        citr ++;
    }    
    //
    if(1)
    {
        unsigned long total_ctg_size      = 0;
        unsigned long total_ctg_size_12lg = 0;
        int lg_i   = 0;
        map<string, unsigned long>::iterator sitr;
        map<string, unsigned long>::iterator sitr_end;
        sitr     = (*lg_grouped_ctg_size).begin();
        sitr_end = (*lg_grouped_ctg_size).end();
        while(sitr != sitr_end)
        {
            cout << "   final: ref " << (*sitr).first 
                 << " grouped "      << (*sitr).second 
                 << " bp. "          << endl; 
            total_ctg_size += (*sitr).second;
            //
            lg_i ++;
            if(chrs.find( (*sitr).first ) != chrs.end())
            {
                total_ctg_size_12lg += (*sitr).second;
            }
            //
            sitr ++;
        }
        cout << "   info: total grouped size " << total_ctg_size      << " bp."           << endl;
        cout << "   info: total grouped size " << total_ctg_size_12lg << " bp in 12 LGs." << endl;        
    }
    //
    return true;
}
//
bool read_paf(string paf_file, 
              map<string, unsigned long>* ctg_size, 
              map<string, map<string, map<unsigned long, unsigned long> > >* aln,
              map<string, map<string, map<unsigned long, unsigned long> > >* aln_query)
{
    // read the alignment info and collect them into the corresponding linkage group
    // aln......: <qry_ctg, <ref_lg, <ref_chr_sta, ref_chr_end> > >
    // aln_query: <qry_ctg, <qry_lg, <qry_chr_sta, qry_chr_end> > > 
    /* PAF line 
	Col	Type	Description
	0	string	Query sequence name
	1	int	Query sequence length
	2	int	Query start (0-based; BED-like; closed)
	3	int	Query end   (0-based; BED-like; open)
	4	char	Relative strand: "+" or "-"
	5	string	Target sequence name
	6	int	Target sequence length
	7	int	Target start on original strand (0-based)
	8	int	Target end   on original strand (0-based)
	9	int	Number of residue matches
	10	int	Alignment block length
	11	int	Mapping quality (0-255; 255 for missing)
    */    
    bool this_check_on = true;
    ifstream ifp;
    ifp.open(paf_file.c_str(), ios::in);
    if(!ifp.good())
    {
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#')
        {
            continue;
        }
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size() < 12) 
        {
            cout << "   Warning: skipping line with insufficient info: " << line << endl;
            continue;
        }
        string          this_ctg_id = lineinfo[0];
        unsigned long this_ctg_size = strtoul(lineinfo[1].c_str(), NULL, 0);
        string          this_ref_id = lineinfo[5];
        unsigned long  this_ref_sta = strtoul(lineinfo[7].c_str(), NULL, 0);
        unsigned long  this_ref_end = strtoul(lineinfo[8].c_str(), NULL, 0);     
        // collect size 
        if((*ctg_size).find(this_ctg_id) == (*ctg_size).end())
        {
            (*ctg_size).insert(std::pair<string, unsigned long>(this_ctg_id, this_ctg_size));
        }else
        {
            assert((*ctg_size)[this_ctg_id] == this_ctg_size);
        }
        // 1. insert new pos info to query     based coordinates
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double ovl_ratio = 0;        
        //
        string         this_qry_id  = this_ctg_id;
        unsigned long  this_qry_sta = strtoul(lineinfo[2].c_str(), NULL, 0);
        unsigned long  this_qry_end = strtoul(lineinfo[3].c_str(), NULL, 0);
        this_qry_end                = this_qry_end - 1; // paf [sta, end)
        if(this_check_on)
        {
            cout << "   check: ctg " << this_qry_id << " " << this_qry_sta << " " << this_qry_end << endl;    
        }    
        bool redundant_ctg_region   = false;
        if((*aln_query).find(this_ctg_id) == (*aln_query).end() )
        {
            // new ctg, first qry chr, first pos
            map<unsigned long, unsigned long> tmp_pos;
            tmp_pos.insert(std::pair<unsigned long, unsigned long>(this_qry_sta, this_qry_end));
            map<string, map<unsigned long, unsigned long> > tmp_chr;
            tmp_chr.insert(std::pair<string, map<unsigned long, unsigned long> >(this_qry_id, tmp_pos));
            (*aln_query).insert(std::pair<string, map<string, map<unsigned long, unsigned long> > > (this_ctg_id, tmp_chr));
        }else
        {
            // 1. update qry coordinates
            if((*aln_query)[this_ctg_id].find(this_qry_id) == (*aln_query)[this_ctg_id].end() )
            {
                // existing ctg, new qry chr, first pos
                map<unsigned long, unsigned long> tmp_pos;
                tmp_pos.insert(std::pair<unsigned long, unsigned long>(this_qry_sta, this_qry_end));
                (*aln_query)[this_ctg_id].insert(std::pair<string, map<unsigned long, unsigned long> >(this_qry_id, tmp_pos));
            }else
            {
                // existing ctg, existing qry chr, new pos: need check if pos overlapping existing ones
                map<unsigned long, unsigned long>::iterator pitr;
                map<unsigned long, unsigned long>::iterator pitr_end;
                pitr     = (*aln_query)[this_ctg_id][this_qry_id].begin();
                pitr_end = (*aln_query)[this_ctg_id][this_qry_id].end();
                bool update_record = true;  
                bool erase_ex      = false;      
                unsigned long erase_ex_sta = 0;        
                while(pitr != pitr_end)
                {
                    unsigned long ex_sta = (*pitr).first;
                    unsigned long ex_end = (*pitr).second;
                    if(ex_end < this_qry_sta)
                    {
                        pitr ++;
                        continue;
                    }
                    if(ex_sta > this_qry_end)
                    {
                        break;
                    }
                    //
                    if(ex_sta<=this_qry_sta && this_qry_sta<=ex_end && ex_end<=this_qry_end)
                    {
                        /*
                                             this_sta---------------------this_end fixed  
                                 ex_sta----------------------------ex_end                                                   
                        */
                        ovl_ratio = (ex_end - this_qry_sta + 1)*1.0 / (ex_end - ex_sta + 1);                      
                        //
                        if(this_check_on)
                        cout << "   check: qry-overlapping1 at "<< this_qry_id 
                             << " "                         << ex_sta
                             << " "                         << ex_end
                             << " vs "                      << this_qry_sta
                             << " "                         << this_qry_end
                             << " by "                      << ovl_ratio
                             << endl;
                        /* update [this_qry_sta, this_qry_end] as [this_qry_sta', this_qry_end ] */
                        this_qry_sta = ex_end+1;
                        if( ovl_ratio > 0.80 ) // hardcoded cutoff!
                        {
                            redundant_ctg_region = true;
                            break;
                        }                          
                    }else
                    if(ex_sta<=this_qry_sta && this_qry_end<=ex_end)
                    {
                        /*
                                             this_sta---------------------this_end fixed     
                                 ex_sta----------------------------------------------ex_end                                                   
                        */
                        ovl_ratio = (this_qry_end - this_qry_sta + 1)*1.0 / (ex_end - ex_sta + 1);                     
                        //                        
                        if(this_check_on)
                        cout << "   check: qry-overlapping2 at "<< this_qry_id
                             << " "                         << ex_sta
                             << " "                         << ex_end
                             << " vs "                      << this_qry_sta
                             << " "                         << this_qry_end
                             << " by "                      << ovl_ratio                             
                             << endl;
                        update_record = false;
                        /* no need to record [this_qry_sta, this_qry_end] */
                        if( ovl_ratio > 0.80 ) // hardcoded cutoff!
                        {
                            redundant_ctg_region = true;
                            break;
                        }                          
                    }else                
                    if(this_qry_sta<=ex_sta && ex_sta<=this_qry_end && this_qry_end<=ex_end)
                    {
                        /*
                                             this_sta---------------------this_end fixed     
                                                      ex_sta-------------------------ex_end
                        */   
                        ovl_ratio = (this_qry_end - ex_sta + 1)*1.0 / (ex_end - ex_sta + 1);                 
                        if(this_check_on)                 
                        cout << "   check: qry-overlapping3 at "<< this_qry_id 
                             << " "                         << ex_sta
                             << " "                         << ex_end
                             << " vs "                      << this_qry_sta
                             << " "                         << this_qry_end
                             << " by "                      << ovl_ratio                                                          
                             << endl;
                        /* update [this_qry_sta, this_qry_end] as [this_qry_sta, this_qry_end' ] */
                        this_qry_end = ex_sta - 1;   
                        if( ovl_ratio > 0.80 ) // hardcoded cutoff!
                        {
                            redundant_ctg_region = true;
                            break;
                        }                                          
                    }else
                    if(this_qry_sta<=ex_sta && ex_sta<=this_qry_end && this_qry_sta<=ex_end && ex_end<=this_qry_end)
                    {
                        /*
                                             this_sta---------------------this_end fixed
                                                      ex_sta------ex_end
                        */
                        ovl_ratio = (this_qry_end - this_qry_sta + 1)*1.0 / (ex_end - ex_sta + 1);
                        //                        
                        if(this_check_on)
                        cout << "   check: qry-overlapping4 at "<< this_qry_id 
                             << " "                         << ex_sta
                             << " "                         << ex_end
                             << " vs "                      << this_qry_sta
                             << " "                         << this_qry_end
                             << " by "                      << ovl_ratio                             
                             << endl;
                        /* erase [ex_sta, ex_end] and record [this_qry_sta, this_qry_end] */
                        erase_ex     = true;
                        erase_ex_sta = ex_sta;
                        if( ovl_ratio > 0.80 ) // hardcoded cutoff!
                        {
                            redundant_ctg_region = true;
                            break;
                        }                          
                    }else ;                       
                    // 
                    pitr ++;
                }
                //
                if(erase_ex)
                {
                    (*aln_query)[this_ctg_id][this_qry_id].erase(erase_ex_sta);
                }
                if(update_record)
                {
                    (*aln_query)[this_ctg_id][this_qry_id].insert(std::pair<unsigned long, unsigned long>(this_qry_sta, this_qry_end));
                }
            }  
        }
        //
        if(ovl_ratio > 0.80)
        {
            // no need to check again this region
            cout << "        : skipping ctg region " << this_qry_id 
                 << " "                              << this_qry_sta 
                 << " "                              << this_qry_end 
                 << endl;
            continue;
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////        
        
        // 2. insert new pos info to reference based coordinates
        if((*aln).find(this_ctg_id) == (*aln).end() )
        {
            // new ctg, first ref chr, first pos
            map<unsigned long, unsigned long> tmp_pos;
            tmp_pos.insert(std::pair<unsigned long, unsigned long>(this_ref_sta, this_ref_end));
            map<string, map<unsigned long, unsigned long> > tmp_chr;
            tmp_chr.insert(std::pair<string, map<unsigned long, unsigned long> >(this_ref_id, tmp_pos));
            (*aln).insert(std::pair<string, map<string, map<unsigned long, unsigned long> > > (this_ctg_id, tmp_chr));
        }else
        {
            // 1. update ref coordinates
            if((*aln)[this_ctg_id].find(this_ref_id) == (*aln)[this_ctg_id].end() )
            {
                // existing ctg, new ref chr, first pos
                map<unsigned long, unsigned long> tmp_pos;
                tmp_pos.insert(std::pair<unsigned long, unsigned long>(this_ref_sta, this_ref_end));
                (*aln)[this_ctg_id].insert(std::pair<string, map<unsigned long, unsigned long> >(this_ref_id, tmp_pos));
            }else
            {
                // existing ctg, existing ref chr, new pos: need check if pos overlapping existing ones
                map<unsigned long, unsigned long>::iterator pitr;
                map<unsigned long, unsigned long>::iterator pitr_end;
                pitr     = (*aln)[this_ctg_id][this_ref_id].begin();
                pitr_end = (*aln)[this_ctg_id][this_ref_id].end();
                bool update_record = true;  
                bool erase_ex      = false;      
                unsigned long erase_ex_sta = 0;        
                while(pitr != pitr_end)
                {
                    unsigned long ex_sta = (*pitr).first;
                    unsigned long ex_end = (*pitr).second;
                    if(ex_end < this_ref_sta)
                    {
                        pitr ++;
                        continue;
                    }
                    if(ex_sta > this_ref_end)
                    {
                        break;
                    }
                    //
                    if(ex_sta<=this_ref_sta && this_ref_sta<=ex_end && ex_end<=this_ref_end)
                    {
                        /*
                                             this_sta---------------------this_end fixed  
                                 ex_sta----------------------------ex_end                                                   
                        */
                        if(this_check_on)
                        cout << "   check: ref-overlapping1 at "<< this_ref_id 
                             << " "                             << ex_sta
                             << " "                             << ex_end
                             << " vs "                          << this_ref_sta
                             << " "                             << this_ref_end
                             << endl;
                        /* update [this_ref_sta, this_ref_end] as [this_ref_sta', this_ref_end ] */
                        this_ref_sta = ex_end+1;
                    }else
                    if(ex_sta<=this_ref_sta && this_ref_end<=ex_end)
                    {
                        /*
                                             this_sta---------------------this_end fixed     
                                 ex_sta----------------------------------------------ex_end                                                   
                        */
                        if(this_check_on)
                        cout << "   check: ref-overlapping2 at "<< this_ref_id
                             << " "                             << ex_sta
                             << " "                             << ex_end
                             << " vs "                          << this_ref_sta
                             << " "                             << this_ref_end
                             << endl;
                        update_record = false;
                        /* no need to record [this_ref_sta, this_ref_end] */
                    }else                
                    if(this_ref_sta<=ex_sta && ex_sta<=this_ref_end && this_ref_end<=ex_end)
                    {
                        /*
                                             this_sta---------------------this_end fixed     
                                                      ex_sta-------------------------ex_end
                        */   
                        if(this_check_on)                 
                        cout << "   check: ref-overlapping3 at "<< this_ref_id 
                             << " "                             << ex_sta
                             << " "                             << ex_end
                             << " vs "                          << this_ref_sta
                             << " "                             << this_ref_end
                             << endl;
                        /* update [this_ref_sta, this_ref_end] as [this_ref_sta, this_ref_end' ] */
                        this_ref_end = ex_sta - 1;
                    }else
                    if(this_ref_sta<=ex_sta && ex_sta<=this_ref_end && this_ref_sta<=ex_end && ex_end<=this_ref_end)
                    {
                        /*
                                             this_sta---------------------this_end fixed
                                                      ex_sta------ex_end
                        */
                        if(this_check_on)
                        cout << "   check: ref-overlapping4 at "<< this_ref_id 
                             << " "                             << ex_sta
                             << " "                             << ex_end
                             << " vs "                          << this_ref_sta
                             << " "                             << this_ref_end
                             << endl;
                        /* erase [ex_sta, ex_end] and record [this_ref_sta, this_ref_end] */
                        erase_ex     = true;
                        erase_ex_sta = ex_sta;
                    }else ;                       
                    // 
                    pitr ++;
                }
                //
                if(erase_ex)
                {
                    (*aln)[this_ctg_id][this_ref_id].erase(erase_ex_sta);
                }
                if(update_record)
                {
                    (*aln)[this_ctg_id][this_ref_id].insert(std::pair<unsigned long, unsigned long>(this_ref_sta, this_ref_end));
                }
            }  
        }
    }
    ifp.close();
    //
    return true;
}
//
bool read_chrid(string chr_file, map<string, unsigned long>* chrs)
{
    // read chr ids (not including scaffolds etc).
    ifstream ifp;
    ifp.open(chr_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: failed in reading chr list file: " << chr_file << endl;
        return false;
    }
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo.size()>=2)
        {
            string chrid = lineinfo[0];
            unsigned long chrsize = strtoul(lineinfo[1].c_str(), NULL, 0);
            (*chrs).insert(std::pair<string, unsigned long>(chrid, chrsize));
        }else
        {
            string chrid = lineinfo[0];
            (*chrs).insert(std::pair<string, unsigned long>(chrid, 1));
        }
    }
    ifp.close();
    //
    return true;
}





























