import re
import sys

inFile = open(sys.argv[1], "r")
seq_valid = re.compile("[ACTG]+")
for line in inFile:
    line = line.strip()
    if line[0] == ">":
        header = line
    elif re.fullmatch(seq_valid, line):
        print(header)
        print(line)
inFile.close()        
            