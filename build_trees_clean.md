## Build trees

Modeled on [Arnason et al. 2018](https://advances.sciencemag.org/content/4/4/eaap9873):

Árnason, Úlfur, et al. "Whole-genome sequencing of the blue whale and other rorquals finds signatures for introgressive gene flow." Science advances 4.4 (2018): eaap9873.



Below is the wrapper submitting jobs for 1000 windows at a time.

```python
import os
#os.mkdir("jobs/")
#os.mkdir("trees_results/")

i=0
j=1
jobfile = open("jobs/job_njob"+str(j)+".sh","w")
jobfile.write("#!/bin/sh\n")
for file in os.listdir("20kwindow_by_window_nonmasked//"):
  i+=1
  jobfile.write(file+"\n")
  jobfile.write("perl fasta2phylip.pl 20kwindow_by_window_nonmasked/"+file+ " trees_results/"+file[:-6]+".phy\n")

  jobfile.write("raxmlHPC-PTHREADS-AVX -T2 -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s  trees_results/"+file[:-6]+".phy"+ " -o Sperm_consensus   -w /home/ludovic.dutoit/nobackup/whale_clean/ludo_calling/trees_results/  -n "+ file[:-6]+"\n")
  if i==1000:
    jobfile.close()
    print("sbatch -A uoo00116 -t 1-00:00:00 --mem=4G -c 2 jobs/job_njob"+str(j)+".sh")
    os.system("sbatch -A uoo00116 -t 1-00:00:00 --mem=4G -c 2 jobs/job_njob"+str(j)+".sh")#submitting specific to cluster
    j+=1
    i=0
    jobfile = open("jobs/job_njob"+str(j)+".sh","w")
    jobfile.write("#!/bin/sh\n")

jobfile.close()#because it is not exactly 1000 jobs
os.system("sbatch -A uoo00116 -t 10:00:00 --mem=4G -c 2 jobs/job_njob"+str(j)+".sh")

```


## Identify low quality windows

Identify low quality windows as windows where any species has less than 8000 non-N sites, remove them.

Check how many sites per windows, one file per species:

``` python
#Python
import os
from Bio import SeqIO
filenames =[file for file in os.listdir(".") if  "fastanovariantsHardMaskedMINCOV5" in file and not file.endswith(".fai")]
length=20000

for file in filenames:
	print(file)
	windows = open("windows"+str(length)+file.split(".")[0]+".bed","w")
	sequences = SeqIO.parse(file, "fasta")
	for seq in sequences:
		print(seq)
		total_valid = 0
		start = 0
		for pos,char in enumerate(str(seq.seq)):
			#print(pos,char)
			total_valid+=1
			if total_valid ==length:
				window = str(seq.name)+"\t"+str(start)+"\t"+str(pos+1)+"\t"+str(str(seq.seq)[start:pos+1].count("A")+str(seq.seq)[start:pos+1].count("T")+str(seq.seq)[start:pos+1].count("G")+str(seq.seq)[start:pos+1].count("C"))+"\n"
				windows.write(window)
				##DO A CHECK HERE

				#print(window)
				start=pos+1
				total_valid =0
	windows.close()

```

Group that into a single table with all species:

```python
import sys
import os

length=20000
output_file=open("VALIDSITES"+str(length)+".txt","w")
filenames =[file for file in os.listdir(".") if  "fastanovariantsHardMaskedMINCOV5" in file and not file.endswith(".fai")]

col = 4
res = {}
filenames = [file for file in os.listdir(".") if "windows"+str(length) in file and  not file == "windows"+str(length)+".bed" and not file == "windows"+str(length)+"0.bed"]

for file_name in filenames:
    for line_nr, line in enumerate(open(file_name)):
        res.setdefault(line_nr, []).append(line.strip().split('\t')[col-1])

output=open("temp","w")
for line_nr in sorted(res):
    output.write('\t'.join(res[line_nr])+"\n")
output.close()

os.system(' cut -f 1-3  windows'+str(length)+'.bed | paste -d "\t" - temp >temp2')
output_file.write("scaf\tstart\tend\t"+"\t".join([name.split(".bed")[0].split(str(length))[1] for name in filenames])+"\n")
with open("temp2") as f:
	for line in f:
		output_file.write(line)
output_file.close()
```

Identify poor windows and remove them. by copying only OK windows to a folder with which we'll do the analyses.

```R
#R

input_folder ="trees_results50k/"

 output_folder="clean_50Kwithoutlowdata"

	 data<-read.table("VALIDSITES50000.txt",h=T)

colnames(data)[10]
#[1] "Sperm"

dir.create(output_folder)





threshold_total=0
threshold_species=8000
lowdata<-c()
nlowdata=0
total=0
for (i in 1:dim(data)[1]){
	total= total+1
	if(sum(data[i,c(4:9,11:15)])<threshold_total | any(data[i,c(4:9,11:15)]<threshold_species)) {
		nlowdata = nlowdata +1 
		lowdata<-rbind(lowdata,data[i,])
	}else{
		print("copying!")
		system(paste("cp ", input_folder,"/RAxML_bipartitionsBranchLabels.",data[i,1],"_",data[i,2],"-",data[i,3]," ",output_folder,sep=""))
		print(c(nlowdata,"low data out of",total))

	}
}
print(nlowdata)
```


#20k 
#[1] "copying!"
[1] "6157"            "low data out of" "122377"

Concatenate them before running astral:

```
cat clean_50Kwithoutlowdata/RAxML_bipartitions* > ALLRAxML_bipartitions.tre
wc -l ALLRAxML_bipartitions.tre
```

## Running Astral



We run the commands and then describe below the output.


```
java -jar astral.5.7.5.jar -i ALLRAxML_bipartitions.tre -o ALLTHEGF_bipartitions.tre --outgroup Sperm_consensus

#annotate each branch
java -jar astral.5.7.5.jar -q ALLTHEGF_bipartitions.tre -i ALLRAxML_bipartitions.tre -o Astral_Scores.tre -t 2  --outgroup Sperm_consensus

#outputting a freqQuad.csv file
java -jar astral.5.7.5.jar -q ALLTHEGF_bipartitions.tre -i ALLRAxML_bipartitions.tre -o Astral_Scores_16.tre -t 16  --outgroup Sperm_consensus
rm Astral_Scores_16.tre #We don't want to keep that tree as it is confusing to the users and only contain some of the info in the tree Astral_Scores.tre
```

Detailed scripts for visualisations are available on request.

