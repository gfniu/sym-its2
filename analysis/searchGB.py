import sys
from Bio import Entrez
from Bio import SeqIO
import os
#output filename
filename = "gb_its2.gb"
Entrez.email = "niugaofeng0808@163.com"

print ("Searching...")
search_handle = Entrez.esearch(db="nucleotide",term="Symbiodinium[Organism] AND (spacer 2[All Fields])",retmax=90000)
record = Entrez.read(search_handle)

id_list =record["IdList"]
count=len(id_list)
print (str(count)+" records.downloading...")
out_handle = open(filename, "w")
n=200
num=0
for down_list in [id_list[i:i + n] for i in range(0, len(id_list), n)]:
    net_handle = Entrez.efetch(db="nucleotide", id=down_list, rettype="gb", retmode="text")
    out_handle.write(net_handle.read())
    num=num+len(down_list)
    remaining=count-num
    print (str(num)+" downloaded"+";\t"+str(remaining)+" remaining;\t waiting...")
out_handle.close()
SeqIO.convert(filename, "genbank", "gb_its2.fasta", "fasta")
print (str(count)+" Saved.")
