# Creating Fasta Gene by Gene

This file is referred to in the [molecular_dating.md](molecular_dating.md) script
This creates one file per gene, with all the exons concatenated into one sequence per species with a single N separating exons.

This script is certainly not the most straightforward way of doing it, but it is reproducible.

First, we transform the Sperm whale original annotation into a BED file with only coding sequences (CDS).

```sh
awk '$3 == "CDS" {print}' GCF_002837175.2_ASM283717v2_genomic.gff > GCF_002837175.2_ASM283717v2_genomicCDS.gff ## CDS
convert2bed --input=gff --output=bed  < GCF_002837175.2_ASM283717v2_genomicCDS.gff > GCF_002837175.2_ASM283717v2_genomicCDS.bed #convert to bed is part of the BEDOPS software.


```

Using python, I remove extra fields, keeping only coordinates of exons and gene names.


```python
 output=open("GCF_002837175.2_ASM283717v2_genomicCDSlight.bed","w")
with open("GCF_002837175.2_ASM283717v2_genomicCDS.bed") as f:
	for line in f:
   		output.write(("\t".join(line.split("\t")[0:3])+"\t"+line.split("\t")[-1].split("gene=")[1].split(";")[0]+"\n"))
output.close()
```
Then I remove all the exons that overlap with another exon in bash.

```
#bash  
cat GCF_002837175.2_ASM283717v2_genomicCDSlight.bed |     mergeBed -d 1 -i stdin  -o distinct -c 4   -delim AAAAA |     grep -v 'AAAA' | sort -k 4 > GCF_002837175.2_ASM283717v2_genomicCDSNOOVERLAPS.bed
#198666 positions
```


Removing the extra digit in the scaffold name of the annotation that is not inside our reference fasta (in Python) 

```python
output=open("GCF_002837175.2_ASM283717v2_genomicCDSNOOVERLAPSCLEAN.bed","w")
with open("GCF_002837175.2_ASM283717v2_genomicCDSNOOVERLAPS.bed") as f:
	for line in f:
		scaffold,start,end,gene= line.split()
		if scaffold.endswith(".1"):
			scaffold=scaffold.split(".1")[0]
			output.write("\t".join([scaffold,start,end,gene])+"\n")
output.close()
```
## Splitting

one fasta sequence per gene, one file per species.

```python

import os
#os.mkdir("genebygene_perspecies/")
filenames =[file for file in os.listdir(".") if  "HardMasked" in file and not file.endswith(".fai")][6:]
for file in filenames:
	print(file)
	os.system("bedtools  getfasta -fullHeader -fi "+file+" -bed GCF_002837175.2_ASM283717v2_genomicCDSNOOVERLAPSCLEAN.bed -fo "+"genebygene_perspecies/"+file+"_splitted.fa" )
```


one file per sequence, with all the species :
```python
import os
from Bio import SeqIO
#os.mkdir("genebygene_windows/")
#os.system("rm genebygene_windows/* ")
for filename in os.listdir("genebygene_perspecies"):
	print(filename)
	species=filename.split(".")[0]
	sequences = SeqIO.parse("genebygene_perspecies/"+filename, "fasta")
	nseq=0
	for seq in sequences:
		nseq+=1
		if nseq%5000==0: print(nseq)#ca 65000 for 20k
		seqname = species #for now the seqname is only the species but that might have to change once I understand the phylogenetic software
		output=open("genebygene_windows/"+seq.name.replace(":","_")+".fasta","a")
		output.write(">"+seqname+"\n"+str(seq.seq)+"\n")
		#print(">"+seqname+"\n"+str(seq.seq)+"\n")
		output.close()
```


add gene names to coordinates.

```python
#cd genebygene_windows

import os
i=0
with open("../GCF_002837175.2_ASM283717v2_genomicCDSNOOVERLAPSCLEAN.bed") as f:
	for line in f:
		i+=1
		if i%5000==0: print(i)#c
		scaffold,start,end,gene= line.split()
		oldname=scaffold+"_"+start+"-"+end+".fasta"
		newname=gene+"_"+scaffold+"_"+start+"-"+end+".fasta"
		os.system("mv "+ oldname+ " "+newname)
```

Join all exons at once per species.

```python
#!mkdir genebygene_windows_concat
#cut -f 4 GCF_002837175.2_ASM283717v2_genomicCDSNOOVERLAPSCLEAN.bed | uniq > gene_list_names.txt #20196 gene names

ipython
import os
i=0
with open("gene_list_names.txt") as f:
	for line in f:
		i+=1
		print(i)
		gene=line.strip()
		files =  ["genebygene_windows/"+file for file in os.listdir("genebygene_windows") if file.split("_")[0]==gene]
		os.system ("paste -d N "+  " ".join(files)+">  genebygene_windows_concat/"+gene+".bed")
		# rename lines to avoid messy line names
		os.system("sed -s -E 's/>[A-Za-z0-9]+_consensusN//g' genebygene_windows_concat/"+gene+".bed  -i")
```

extension bed should be fasta (fixing mistake)

```python
filenames =[ file for file in  os.listdir(".")]
for file in filenames:
	i+=1
os.system("mv "+ file+" "+file.split(".bed")[0]+".fa")
```
