/* annotate SVs with a given gff for genes and its surrounding 3kb regions
   if a mutation overlaps with such a region, set the annotation info for it.
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
#include "split_string.h"
using namespace std;
//
bool sv_reader(string inputfile, int offset, map<string, map<unsigned long, unsigned long> >* svs, map<string, map<unsigned long, string> >* svsrawline);
//
int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        // g++ sv_annotator.cpp split_string.cpp -O3 -o sv_annotator
        cout << "\nFunction: annotate svs with a given gff. ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "Usage: gene_sv_annotator svs.txt file.gff\n"                   << endl;
        cout << "       argv[1]=svs.txt gives svs with first three columns as chr start end."  << endl;
        cout << "       argv[2]=file.gff   is the corresponding gff file. "                          << endl;
        return 1;
    }
    cout << "   pre-warning: chr-ids in snp/sv and gff file must be consistent." << endl;
    // input
    string varfile = (string)argv[1];
    string gfffile = (string)argv[2];
    // read var info
    map<string, map<unsigned long, unsigned long> > vars;
    map<string, map<unsigned long, string> >        varsrawline;
    vector<string> gene_only;
    cout << "   Reading svs, file provided: " << varfile << endl;
    int offset = 0;
    if(!sv_reader(varfile, offset, &vars, &varsrawline))
    {
      return 1;
    }
    else
    {
      cout << "   Reading svs done."          << endl;
    }
    // reading gff info
    ifstream gffp;
    gffp.open(argv[2], ios::in);
    while(gffp.good())
    {
        string line("");
        getline(gffp, line);
        if(line.size()==0) continue;
        if(line[0] == '#') continue;
        // chr1	ensembl	upstr/gene/downstr	7731	10730	.	+	.	pID=01s0011g00010_l
        vector<string> lineinfo = split_string(line, '\t');
        if(lineinfo[2].compare("gene") != 0)
        continue;
        if(vars.find(lineinfo[0]) == vars.end()) continue;
        unsigned long gene_start = strtoul(lineinfo[3].c_str(), NULL, 0);
        unsigned long gene_end   = strtoul(lineinfo[4].c_str(), NULL, 0);
        string        gene_strand=lineinfo[6]; // "+/-"
        size_t pos1              = lineinfo[8].find("ID=");
        size_t pos2              = lineinfo[8].find(";");
        string        gene_name  = lineinfo[8].substr(pos1+3, pos2-pos1-3);
        //
        long down_sta = gene_end+1;
        long down_end = gene_end+1000; // overflow at chr ends?
        long up_sta   = gene_start>1?(gene_start-1):1;
        long up_end   = gene_start>1000?(gene_start-1000):up_sta;
        if(gene_strand.compare("-")==0)
        {
            up_sta = gene_end+1;
            up_end = gene_end+1000;
            down_sta   = gene_start>1?(gene_start-1):1;
            down_end   = gene_start>1000?(gene_start-1000):down_sta;
        }
        //
        map<unsigned long, unsigned long>::iterator seg_itr;
        map<unsigned long, unsigned long>::iterator seg_itr_end;
        seg_itr     = vars[ lineinfo[0] ].begin();
        seg_itr_end = vars[ lineinfo[0] ].end();
        while(seg_itr != seg_itr_end)
        {
            unsigned long var_start = (*seg_itr).first;
            unsigned long var_end   = (*seg_itr).second;
            if(var_start > gene_end+5000) break;
            // if gene completely in sv, skip it
            if(gene_start>=var_start && gene_end<=var_end)
            {
              seg_itr ++;
              continue;
            }
            //
            if( (var_start>=gene_start && var_start<=gene_end) ||
                (var_end>=gene_start   && var_end<=gene_end) )
            {
                // gene_body
                varsrawline[ lineinfo[0] ][var_start] = varsrawline[ lineinfo[0] ][var_start] + "\t" + gene_name + ":" + lineinfo[3] + "-" + lineinfo[4] + ":" + gene_strand + ":gene_body";
                gene_only.push_back(gene_name + "\t" + lineinfo[3] + "\t" + lineinfo[4] + "\t" + gene_strand + "\tgene_body");
            }else
            if( (var_start>=up_sta && var_start<=up_end) ||
                (var_end>=up_sta   && var_end<=up_end) )
            {
                // gene upstream
                varsrawline[ lineinfo[0] ][var_start] = varsrawline[ lineinfo[0] ][var_start] + "\t" + gene_name + ":" + lineinfo[3] + "-" + lineinfo[4] + ":" + gene_strand + ":gene_up";
                gene_only.push_back(gene_name + "\t" + lineinfo[3] + "\t" + lineinfo[4] + "\t" + gene_strand + "\tgene_up");
            }else
            if( (var_start>=down_sta && var_start<=down_end) ||
                (var_end>=down_sta   && var_end<=down_end) )
            {
                // gene downstream
                varsrawline[ lineinfo[0] ][var_start] = varsrawline[ lineinfo[0] ][var_start] + "\t" + gene_name + ":" + lineinfo[3] + "-" + lineinfo[4] + ":" + gene_strand + ":gene_down";
                gene_only.push_back(gene_name + "\t" + lineinfo[3] + "\t" + lineinfo[4] + "\t" + gene_strand + "\tgene_down");
            };
            seg_itr ++;
        }
    }
    gffp.close();

    // output annotation 1
    vector<string> infileinfo = split_string(varfile, '/');
    string ofilename("");
    ofilename = "annotated_" + infileinfo[infileinfo.size()-1];
    ofstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open file " << ofilename << " to output annotation." << endl;
        return 1;
    }
    map<string, map<unsigned long, string> >::iterator chr_itr;
    map<string, map<unsigned long, string> >::iterator chr_itr_end;
    chr_itr     = varsrawline.begin();
    chr_itr_end = varsrawline.end();
    while(chr_itr != chr_itr_end)
    {
        map<unsigned long, string>::iterator seg_itr;
        map<unsigned long, string>::iterator seg_itr_end;
        seg_itr     = (*chr_itr).second.begin();
        seg_itr_end = (*chr_itr).second.end();
        while(seg_itr != seg_itr_end)
        {
            vector<string> annoinfo = split_string((*seg_itr).second, '\t');
            if(annoinfo.size()<5)
            {
              (*seg_itr).second += "\tintergenic";
            }
            ofp << (*seg_itr).second << endl;
            seg_itr ++;
        }
        chr_itr ++;
    }
    cout << "   Annotated variations have be recorded in " << ofilename  << ".\n" << endl;
    ofp.close();
    // output annotation 2
    string ofilename2("");
    ofilename2 = "annotated_sv_hit_genes_" + infileinfo[infileinfo.size()-1];
    ofstream ofp2;
    ofp2.open(ofilename2.c_str(), ios::out);
    if(!ofp2.good())
    {
        cout << "   Error: cannot open file " << ofilename2 << " to output annotation." << endl;
        return 1;
    }
    vector<string>::iterator gene_itr;
    vector<string>::iterator gene_itr_end;
    gene_itr     = gene_only.begin();
    gene_itr_end = gene_only.end();
    while(gene_itr != gene_itr_end)
    {
        ofp2 << (*gene_itr) << endl;
        gene_itr ++;
    }
    cout << "   Annotated genes have be recorded in " << ofilename2  << ".\n" << endl;
    ofp2.close();
    return 0;
}
#
bool sv_reader(string inputfile, int offset, map<string, map<unsigned long, unsigned long> >* svs, map<string, map<unsigned long, string> >* svsrawline)
{
    ifstream ifp;
    ifp.open(inputfile.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot open " << inputfile << " to read info of svs." << endl;
        return false;
    }
    unsigned long num_seg = 0;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0) continue;
        if(line[0] == '#') continue;

        // [xxx]	chr1	126618	126618 xxx [xxx]... : pos of chr can be tuned by offset value
        vector<string> lineinfo = split_string(line, '\t');
        string chr              = lineinfo[0+offset];
        unsigned long  start    = strtoul(lineinfo[1+offset].c_str(), NULL, 0);
        unsigned long  end      = strtoul(lineinfo[2+offset].c_str(), NULL, 0);

        if((*svs).find(chr) == (*svs).end())
        {
            map<unsigned long, unsigned long> tmppair;
            tmppair.insert(std::pair<unsigned long, unsigned long>(start, end));
            (*svs).insert(std::pair<string, map<unsigned long, unsigned long> >(chr, tmppair));
            // sv line info
            map<unsigned long, string> tmppair2;
            tmppair2.insert(std::pair<unsigned long, string>(start, line));
            (*svsrawline).insert(std::pair<string, map<unsigned long, string> >(chr, tmppair2));
        }
        else
        {
            map<string, map<unsigned long, unsigned long> >::iterator chritr;
            chritr = (*svs).find(chr);
            (*chritr).second.insert(std::pair<unsigned long, unsigned long>(start, end));
            // sv line info
            map<string, map<unsigned long, string> >::iterator chritr2;
            chritr2 = (*svsrawline).find(chr);
            (*chritr2).second.insert(std::pair<unsigned long, string>(start, line));
        }
        num_seg ++;
    }
    ifp.close();
    cout << "    number of svs recorded: " << num_seg << endl;
    return true;
}
