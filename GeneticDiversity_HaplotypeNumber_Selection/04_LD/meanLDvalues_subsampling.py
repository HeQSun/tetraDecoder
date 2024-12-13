import os
import sys
import glob
import random
import math
import numpy as np



chrID = str("chr01")
win_size = int(10000)
chr_size = int(90000000)
output_name = str("test")

chrID = str(sys.argv[1])
win_size = int(sys.argv[2])
chr_size = int(sys.argv[3])
output_name = str(sys.argv[4])

FilenamesList = glob.glob('LD_NoSingletons_'+chrID+'_maxDis100Mb/LD_'+chrID+'_*.txt')

NWin = math.floor(int(chr_size) / win_size)+5
sumR2 = np.zeros((NWin, NWin))
NLines = np.zeros((NWin, NWin))

output_subsampling = open(f'Subsampling_python_{chrID}_maxDis100Mb_WinSize{win_size}.txt', "w")
output_meanVal = open(f'{output_name}_{chrID}_maxDis100Mb_WinSize{win_size}.txt', "w")

for file in FilenamesList:
    print(file)
    for line in open(file, "r"):
        if line[0:5] == "chr\ts":
            continue
        else:
            if random.random() < 0.0005:
                output_subsampling.write(f'{line}')
            line_values = line.strip().split("\t")
            win1 = math.floor(int(line_values[1]) / win_size)
            win2 = math.floor(int(line_values[2]) / win_size)
            sumR2[win1, win2] += float(line_values[8])
            NLines[win1, win2] += 1


for win1 in range(NWin):
    for win2 in range(NWin):
        if win2 >= win1:
            if NLines[win1, win2] > 0:
                output_meanVal.write(f'{win1*win_size}\t{win2*win_size}\t{sumR2[win1, win2]/NLines[win1, win2]}\t{NLines[win1, win2]}\n')



output_subsampling.close()
output_meanVal.close()





