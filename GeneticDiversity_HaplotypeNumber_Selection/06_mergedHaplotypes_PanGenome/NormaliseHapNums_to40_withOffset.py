#Script to take one of Sergio's haplotypeGroup files and
import sys

inFile = sys.argv[1] #/Users/craig/Documents/Sergio_HapVis/haplotypeGroups_min10SS_chr06
outFile= inFile+"_normed40_withOffset.tsv"
out=open(outFile,"w+")
#windowSize=10000 #update this
current_start="0"
HapList=[] #List of haplotypes seen in a given window
Lines=[] #Store lines of the file for one window.
hapOffset=0.02 #An offset added to each subsequent member of the same haplotype cluster, so that the lines in R don't all draw on top of eachother
maxGap = float(sys.argv[2])

#Genomes= ["A_hap21","A_hap22","A_hap23","A_hap24","B_hap21","B_hap22","B_hap23","B_hap24","C_hap21","C_hap22","C_hap23","C_hap24","D_hap21","D_hap22","D_hap23","D_hap24","E_hap21","E_hap22","E_hap23","E_hap24","F_hap21","F_hap22","F_hap23","F_hap24","G_hap21","G_hap22","G_hap23","G_hap24","H_hap21","H_hap22","H_hap23","H_hap24","I_hap21","I_hap22","I_hap23","I_hap24","J_hap21","J_hap22","J_hap23","J_hap24"]

#Go through lines of the haplotypeGroup file
for line in open(inFile+".txt","r"):

        chr, window_start, genome, hap, missingFraction = line.split("\t")

        #If we have some new window, output the old stuff
        if window_start != current_start:
                Sorted_haps = sorted(set(HapList)) #sort the list of haplotype groups observed for this window
                seenHaps={} #Dictionary which will record the offset needed for drawing lines later, given how often we've seen a particular haplotype group in a window.

                for l in Lines: #for each line seen in this window
                        tmp_chr, tmp_window_start, tmp_genome, tmp_hap, tmp_missingFraction = l.rstrip().split("\t")

                        if float(tmp_missingFraction)>=maxGap: # assign a -1 hap to missing values
                                tmp_hap="-1"
                        else: # otherwise replace the hap number with the ordinal position of the haplotype among observed haplotypes
                                tmp_hap = str(Sorted_haps.index(tmp_hap)+1) #eg if observed haplotypes are [3,6,24], then 3->1, 6->2, 24->3 respectively

                        #Update the list of haplotypes we have seen.
                        if tmp_hap not in seenHaps:
                                seenHaps[tmp_hap]=hapOffset # If this is the first time we've seen this haplotype, add the offset ready for the next one.
                        else: #if we have seen this before, increase the offset
                                new_tmp_hap = str(float(tmp_hap) + seenHaps[tmp_hap])
                                seenHaps[tmp_hap]=seenHaps[tmp_hap]+hapOffset
                                tmp_hap=new_tmp_hap

                        #output
                        out.write("\t".join([tmp_chr, tmp_window_start, tmp_genome, tmp_hap, tmp_missingFraction])+"\n")

                        #Reset some vars
                        current_start=window_start
                        HapList=[]
                        Lines=[]
        # else, if we are still getting info for this same window
        if window_start == current_start:

                Lines.append(line) #add this line to the list of lines seen.
                if float(missingFraction)<maxGap: #only add haps that aren't missing
                        HapList.append(hap) #build the list of haplotypes



#catch the last window here
Sorted_haps = sorted(HapList)
for l in Lines:
        tmp_chr, tmp_window_start, tmp_genome, tmp_hap, tmp_missingFraction = l.rstrip().split("\t")
        if float(tmp_missingFraction)>=maxGap: # assign a -1 hap to missing values
                tmp_hap="-1"
        else: # replace the hap nurmber with the ordinal position of the haplotype among observed haplotypes
                tmp_hap = str(Sorted_haps.index(tmp_hap)+1)
        out.write("\t".join([tmp_chr, tmp_window_start, tmp_genome, tmp_hap, tmp_missingFraction])+"\n")


