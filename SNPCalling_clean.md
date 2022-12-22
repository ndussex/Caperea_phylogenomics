# Creating consensus sequences

### bcftools mpileup
The first step is to call variants using bcftools mpileup:

```
REF="GCF_002837175.2_ASM283717v2_genomic_headers_fixed.fa"
genome=$REF".fai"
module load BCFtools BEDTools SAMtools # specific to local server


First bcftools mpileup using the bam files:

#Call variants (making sure to output required field FORMAT/AD)
for b in 	*realn.bam ; do
  out=${b%.bam}.bcf.gz
  echo $out 
   bcftools mpileup -O u -Q 30 -q 30 -t AD -A -B -f $REF $b   | bcftools call -m -M -O b  --threads 8  -o $out  
  bcftools index $out
done
```
## Filtering indels and identify low quality sites

We then start filtering. We mark SNPs within 3 bp of indels with the SnpGap filter for exclusion. We get a list of all sites not in vcf, filter indels.

We then make sure we take  a list of all SNPs not in the ALLSITES vcf as low quality sites (excluded as indels or in the mpileup). They will be coded as Ns in the final consensus.

It is shown below as a small self-running script that can be applied to all species.


```
REF="GCF_002837175.2_ASM283717v2_genomic_headers_fixed.fa"
genome=$REF".fai"
module load BCFtools BEDTools SAMtools # specific to local server
##Usage
#for b in *.bcf.gz; do sbatch -A uoo00116 -t 1-00:00:00  --mem=8G addfiltersnpgapandcomplement.sh $b ; done # server specific
b=$1 # each bcf
echo $b
out=SnpGap.$b
bcftools filter --SnpGap 3    -O b -o $out $b
bcftools view  --exclude-types indels -O v -o $out".recode.vcf" $b
ln -s $out".recode.vcf" $out".vcf"
bedtools complement -i  $out".vcf" -g $genome > $out".bed"
newname=$(basename  $out ".bcf.gz")".vcf"
mv $out $newname
bgzip $newname
```


The script below first apply further filters. It considers all sites with less than 5 of coverage as low quality. Those sites are added to the complement created above which contained indels and mpileup ignored sites as low quality. They are all coded as N in consensus.


```
#####Written by Nic dussex for UPPMAX, adapted by Ludovic Dutoit for Mahuika.

###usage
#cd $vcf_dir
# for i in $(ls SnpGap*.vcf.gz); do sbatch -A uoo00116 --qos=short ./SNP_calling_part2_consensus.sh   $i; done

module load GATK VCFtools SAMtools # server specific


#copy files from home to tmp directory

vcf_dir='bams'
OUT='bams/'
REF='non_ambiguous/whale/whale/GCF_002837175.2_ASM283717v2_genomic_headers_fixed.fa'



########

vcfgz=$1
### 1. Generate VCF file of sites that will be masked by “N” in the consensus fasta: 1) sites of DP < 3 and 2) truly polymorphic sites
#picard CreateSequenceDictionary -R $REF
gatk IndexFeatureFile --input $vcfgz
gatk SelectVariants  -R $REF --variant $vcfgz -select "DP < 5" -O $vcfgz"_lowDP.vcf.gz"
vcftools --gzvcf $vcfgz --maf 0.00001 --max-maf 0.99999 --min-alleles 2 --remove-indels --recode --out $vcfgz"_trueSNPs"


# 2. Convert VCF files of sites to be masked into tab delimited files of the site positions, concatenate and sort them.
bgzip -cd $vcfgz"_lowDP.vcf.gz" | grep -v '^#' | awk '{print $1 "\t" $2}' > $vcfgz"_lowDP.list"
grep -v '^#' $vcfgz"_trueSNPs.recode.vcf" | awk '{print $1 "\t" $2}' > $vcfgz"_trueSNPs.list"
cat $vcfgz"_lowDP.list" $vcfgz"_trueSNPs.list" | sort -k1,1 -k2,2n -u > $vcfgz"_lowDP_trueSNPs.list"

# Zip and index the vcf file of true SNPs to save space
bgzip $vcfgz"_trueSNPs.recode.vcf"
tabix -p vcf $vcfgz"_trueSNPs.recode.vcf.gz"


output=$(basename $vcfgz ".vcf.gz")".unmaskedLOWCOV_consensus.fa"
bcftools consensus --mask $vcfgz"_lowDP_trueSNPs.list" --fasta-ref $REF $vcfgz > $output
#outputclean=$(basename $vcfgz ".vcf.gz")".maskedLOWCOV_consensus.fa"
#rm $output

echo "DONE: OUTPUT" $output

```

We now have unmasked consensus for each species.

### Masking


```
bedtools maskfasta -fi SnpGap.Blue.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.Blue.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/Blue_consensus.fa
bedtools maskfasta -fi SnpGap.Bow.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.Bow.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/Bow_consensus.fa
bedtools maskfasta -fi SnpGap.Fin.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.Fin.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/Fin_consensus.fa
bedtools maskfasta -fi SnpGap.Gray1.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.Gray1.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/Gray1_consensus.fa
bedtools maskfasta -fi SnpGap.Gray2.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.Gray2.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/Gray2_consensus.fa
bedtools maskfasta -fi SnpGap.Hump.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.Hump.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/Hump_consensus.fa
bedtools maskfasta -fi SnpGap.Minke.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.Minke.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/Minke_consensus.fa
bedtools maskfasta -fi SnpGap.NARight.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.NARight.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/NARight_consensus.fa
bedtools maskfasta -fi SnpGap.PygmyWhale.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.PygmyWhale.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/PygmyWhale_consensus.fa
bedtools maskfasta -fi SnpGap.Sei1.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.Sei1.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/Sei1_consensus.fa
bedtools maskfasta -fi SnpGap.Sei2.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.Sei2.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/Sei2_consensus.fa
bedtools maskfasta -fi SnpGap.Sperm.merged.rmdup.merged.realn.bcf.gz.unmaskedLOWCOV_consensus.fa -bed SnpGap.Sperm.merged.rmdup.merged.realn.bcf.gz.bed -fo ludo_calling/Sperm_consensus.fa
cd /home/ludovic.dutoit/nobackup/whale_clean/ludo_calling
# ln -s /home/ludovic.dutoit/nobackup/whale_clean/bams/whale_for_Kieren/ludo_calling*.fa .
```

All the positions that were not in the variant calling vcf (indels, low quality, etc) are now masked. We need to add the repeat masking the repeat regions of the Sperm on all the species.


First, I make fasta files for all the species but Sperm in upper case.

```python
# Python v3
import os
filenames= [ file for file in  os.listdir(".") if file.endswith(".fa") and not "novariants" in file and not "bed" in file and not "GCF" in file and not "Sperm" in file]
for filename in filenames:
	print(filename)
	output = open(filename.split(".fa")[0]+".correct_case.fa","w")
	from Bio import SeqIO
	import re
	for seq_record in SeqIO.parse(filename, "fasta"):
			print (seq_record.name)
			output.write(">"+seq_record.id+"\n")
			output.write(str(seq_record.seq).upper()+"\n")
	output.close()

#Just copy the sperm whale as is, it is the only one with correct case
filename= [ file for file in  os.listdir(".") if "fa" in file and not "novariants" in file and not "bed" in file and not "GCF" in file and  "Sperm" in file][0]
os.system(" cp "+ filename +" "+filename.split(".fa")[0]+".correct_case.fa")
```
```python
#python


from datetime import datetime
def timeit():![image](https://user-images.githubusercontent.com/4376065/181658443-b7aeab68-876c-4ea5-b7b7-b814f2446636.png)
![Uploading image.png…]()

	now = datetime.now()
	current_time = now.strftime("%H:%M:%S")
	print("Current Time =", current_time)

import os
from Bio import SeqIO
import re
filenames= [ file for file in  os.listdir(".") if "fa" in file and not "novariants" in file and not "bed" in file and not "GCF" in file and "correct" in file and "Sperm" in file]
for filename in filenames:
	print(filename)
	timeit()
	output = open("positions_to_mask"+filename+".bed","w")
	for seq_record in SeqIO.parse(filename, "fasta"):
		print (seq_record.name)
		#seq_record.seq = re.sub(r'[a-z]',"N",str(seq_record.seq)) #HARDM
		variants = [seq_record.name+"\t"+str(m.start())+"\t"+str(m.start()+1)+"\n" for m in re.finditer(r'[RYSWKMBDHVNatgcryswkmbdhv]', str(seq_record.seq))]
		output.write("".join(variants))
	output.close()
	timeit()
	print("merging")
	os.system("bedtools merge -i positions_to_mask"+filename+".bed > positions_to_mask_cleaned"+filename+".bed")
	timeit()

```




### Mask those positions


```python
#python
#module load BEDTools

import os,re
filenames =[file for file in os.listdir(".") if not "novariants" in file and not "GCF" in file and not "bed" in file and file.endswith(".fa") and "correct" in file]
for filename in filenames:
	print(filename)
	os.system("bedtools maskfasta -fullHeader -fi "+ filename+ " -fo "+ re.sub('\.fasta$', '',filename)+"novariantsHardMasked.fasta" +"  -bed positions_to_mask_cleanedSperm_consensus.correct_case.fa.bed")


```
# Splitting

The goal here is to obtain each window is one fasta file with all species involved.

### Identify windows

Obtain windows in pieces of 20kb. this is slightly more cumbersome that needs be. Windows could potentially be created in a single bash command but it is a remnant from trying to have windows with equal number of valid sites as per Arnason et al. We are not doing this.

```python
#python
import os
from Bio import SeqIO
length=20000
nwindows=0
windows = open("windows"+str(length)+".bed","w")
sequences = SeqIO.parse("Gray2_consensus.correct_case.fanovariantsHardMasked.fasta", "fasta")
for seq in sequences:
	print(seq)
	total_valid = 0
	start = 0
	for pos,char in enumerate(str(seq.seq)):
		#print(pos,char)
		total_valid+=1
		if total_valid ==length:
			window = str(seq.name)+"\t"+str(start)+"\t"+str(pos+1)+"\n"
			windows.write(window)
			nwindows+=1
			print(nwindows,window)
			##DO A CHECK HERE

			#print(window)
			start=pos+1
			total_valid =0
windows.close()
```
### Splitting each species into all the windows
one fasta sequence per window, one file per species

```python
#python
import os
#os.mkdir("20k_windows_nonmasked/")
filenames =[file for file in os.listdir(".") if  "HardMasked" in file and not file.endswith(".fai")][6:]
for file in filenames:
	print(file)
	os.system("bedtools getfasta -fi "+file+" -bed windows20000.bed -fo "+"20k_windows_nonmasked/"+file+"_splitted.fa" )
```

### Grouping sequence per window, not per species

one file per sequence:
```python
#python
import os
from Bio import SeqIO
#os.mkdir("20kwindow_by_window_nonmasked/")
#os.system("rm 20kwindow_by_window_nonmasked/* ")
for filename in os.listdir("20k_windows_nonmasked"):
	print(filename)
	species=filename.split(".")[0]
	sequences = SeqIO.parse("20k_windows_nonmasked/"+filename, "fasta")
	nseq=0
	for seq in sequences:
		nseq+=1
		if nseq%5000==0: print(nseq)#ca 65000 for 20k
		seqname = species #for now the seqname is only the species but that might have to change once I understand the phylogenetic software
		output=open("20kwindow_by_window_nonmasked/"+seq.name.replace(":","_")+".fasta","a")
		output.write(">"+seqname+"\n"+str(seq.seq)+"\n")
		#print(">"+seqname+"\n"+str(seq.seq)+"\n")
		output.close()
```

At this stage, each window is one fasta file with all species involved.
	

