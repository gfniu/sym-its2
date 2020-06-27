import os
import re
import argparse

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--inputfile_path', type=str, default=None,
                    help="xfasta format.")
args = parser.parse_args()
inputfile=args.inputfile_path
process = os.popen('RNAeval -i '+inputfile)
result=process.read();
lineList=result.split("\n")
print(len(lineList))
i=0
while i<len(lineList):
    line=lineList[i]
    print(str(i))
    print(line)
    index=line.find(" (")
    if index>-1:
        freeEnergy=line[0-(len(line)-index-1):]
        freeEnergy=re.sub("[ (]|[)]","",freeEnergy)
        lineList[i-2]=lineList[i-2]+" dG="+freeEnergy+"kcal/mol"
        lineList[i]=line[0:index]
    i=i+1
filename = "core set secondary structures.xfasta"
out_handle = open(filename, "w")
out_handle.write("\n".join(lineList))
out_handle.close()
process.close()
