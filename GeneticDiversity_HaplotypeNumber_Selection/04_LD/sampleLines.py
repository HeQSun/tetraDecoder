#!/bin/python3
import os
import sys
import random
import argparse

parser = argparse.ArgumentParser(description="Returns random lines from text files, as fast as possible.")
parser.add_argument('--file', required=True, help="the file to parse")
parser.add_argument('-n', required=True, help="amount of lines you want to return")
parser.add_argument('--min-len', dest='minLen', type=int, default=0, help="minimum line length, defaults to 0")
parser.add_argument('--max-len', dest='maxLen', type=int, default=float('inf'), help="maximum line length, defaults to ∞")
args = parser.parse_args()

# len mistake catch
if args.minLen > args.maxLen:
    sys.exit("I can't think of any strings that are longer than " + str(args.minLen) + " while being shorter than " + str(args.maxLen) + ".")

filename = args.file
filesizeBytes = os.path.getsize(filename)
#print(filesizeBytes)

# you can find this with: sudo blockdev --getbsz $(df --output=source $(pwd) | grep '/')
blockSize = 4096               # Bytes

import random
f = open(filename, "rt")
count = 0
bytesRead = 0
while (count < 10000):
  buffer = f.read(blockSize)
  if not buffer: break
  count += buffer.count('\n')
  bytesRead += blockSize
f.close()

# less then 1 blocksize file fix
if bytesRead > filesizeBytes:
    bytesRead = filesizeBytes

bytesPerLine = bytesRead/count
#print(bytesPerLine)

totalLinesEst = filesizeBytes / bytesPerLine
#print("Estimated lines in file: " + str(totalLinesEst))

# exact output count
resultlines = 0

fh = open(filename)
while resultlines < int(args.n):
    linestart = random.randint(0, int(totalLinesEst))
    readstart = (linestart * bytesPerLine) - bytesPerLine
    if readstart < bytesPerLine:
        readstart = 0
    else:
        fh.seek(readstart)
    scratch = fh.readline()
    line = fh.readline().replace('\n','')
    while not (args.minLen <= len(line) <= args.maxLen):
        line = fh.readline().replace('\n','')
    print(line)
    resultlines += 1
fh.close()


