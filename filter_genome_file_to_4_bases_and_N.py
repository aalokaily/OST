import sys
import os
import fileinput
 
fasta_file = sys.argv[1]
out_file =  open(os.path.basename(fasta_file).rsplit(".",1)[0] + ".one_line", 'wa')

# get all variants from each file without intersection
for line in fileinput.input(fasta_file):
    if line[0] == ">":
        pass
    else:
        for c in line:
            if c in ["A", "C", "G", "T", "N"]:
                out_file.write(c)
