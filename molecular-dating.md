# Prepartion of input files for molecular dating in mcmctree (PAML)

## Step 1: Extract and concatenate exons for each protein coding gene in FASTA format

This is done in the file. [Extracting_gene_sequences.md]( Extracting_gene_sequences.md). This creates one file per gene, with all the exons concatenated into one sequence per species with a single N separating exons.


## Step 2: Create a tree for each gene and filter based on minimum mean bootstrap support

This step removes genes with low phylogenetic information content.

```
for INPUT in *.fa

do

# create tree using RAxML

raxmlHPC-PTHREADS-AVX_v8.2.0 -f a -m GTRGAMMA -s $INPUT -o Sperm_consensus -T 4 -# 100 -x 34567 -p 34567 -n $INPUT

mv RAxML_bipartitions.$INPUT RAxML_bipartitions.$INPUT.tree

# calculate mean bootstrap support

phykit bipartition_support_stats RAxML_bipartitions.$INPUT.tree > $INPUT.support

# read mean bootstrap support from file

meansupport=$(awk 'NR == 1 {print $2}' $INPUT.support)
meansupportroundeddown=${meansupport%.*} 

echo $meansupport
echo $meansupportroundeddown

# move files for genes that pass filter into new folder

if [[ $meansupportroundeddown -gt 80 ]] # minimum mean bootstrap support
then
  cp $INPUT ../passed_bootstrap/
  cp RAxML_bipartitions.$INPUT.tree ../passed_bootstrap
fi

# tidy up

rm ./$INPUT.reduced
rm ./$INPUT.support
rm ./RAxML*$INPUT*

done
```

## Step 3: Filter genes that passed the previous step based on their Robinson-Foulds distance to the species tree

This step removes genes with histories that deviate substantially from the species tree.

### Species tree
```
(((NARight_consensus,Bow_consensus),(((((Gray2_consensus,Gray1_consensus),(Hump_consensus,Fin_consensus)),(Blue_consensus,(Sei1_consensus,Sei2_consensus))),Minke_consensus),PygmyWhale_consensus)),Sperm_consensus);
```

```
for INPUT in *.tree

do

phykit robinson_foulds_distance ./speciestree.txt $INPUT > $INPUT.rfdistance

rfdistance=$(awk 'NR == 1 {print $1}' $INPUT.rfdistance)

echo $rfdistance

tmp=${INPUT#*.}
b=${tmp%.*}

if [[ $rfdistance -lt 5 ]] # minimum RF distance
then
  cp $INPUT ../passed_RF/
  cp $b ../passed_RF/
fi

rm ./$INPUT.rfdistance

done
```

## Step 4: Load trees of those genes that passed Step 4 into R and calculate coefficient of variation of root-to-tip distance 

This step measures the extent to which the evolution of a gene is "clock-like".

```
# load in .tree files for each gene

files <- Sys.glob("*.tree")

# define ingroup to use in downstream calculations (to exclude sperm whale branch from calculation)

ingroup <- c("NARight_consensus", "Blue_consensus", "Fin_consensus", "Gray1_consensus", "PygmyWhale_consensus", "Minke_consensus", "Bow_consensus", "Hump_consensus", "Sei1_consensus", "Sei2_consensus", "Gray2_consensus")

# loop on .tree files to calculate coefficient of variation of root-to-tip distance (test.distroot.cov) for each gene and write to file

for( i in files )
{
    test <- ape::read.tree(file=i)
    test.distroot <- distRoot(test, tips=ingroup)
	  test.distroot.mean <- mean(test.distroot)
	  test.distroot.sd <- sd(test.distroot)
	  test.distroot.cov <- test.distroot.sd/test.distroot.mean*10
	  print(test.distroot.cov)
	  write(test.distroot.cov, file=paste0(i, ".cov"))
}
```

## Step 5: Filter genes that passed the Step 4 based on their clock-likeness

This step reomves those genes that do not display clock-like evolution.

```
for INPUT in *.cov

do

cov=$(awk 'NR == 1 {print $1}' $INPUT)
covroundeddown=${cov%.*} 

echo $cov
echo $covroundeddown

tmp=${INPUT%.*}
b=${tmp%.*}
bb=${b#*.}

if [[ $covroundeddown -lt 1 ]]
then
  cp $tmp ../passed_cov/
  cp $bb ../passed_cov/
fi

rm ./$INPUT

done
```

## Step 5: Estimate substitution rate for each gene that passed filters

In our case, 344 genes passed all filters. These were concatenated into a single PHYLIP file and used as the input for analysis in baseml (part of the PAML package).

### baseml.ctl
```
      seqfile = ./344g_110422.phy
     treefile = ./whales.tree

        ndata = 344
      outfile = mlb       * main result file
        noisy = 3   * 0,1,2,3: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output

        model = 7   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
                    * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

        clock = 1   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below
        kappa = 2.3   * initial or fixed kappa

    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0.5  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 4  * # of categories in the dG, AdG, or nparK models of rates

        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states
       method = 0  * 0: simultaneous; 1: one branch at a time
   Small_Diff = 1e-6
    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
```

The genes were then sorted by their estimate substitution rate and binned into eight equal size partitions.

## Step 6: Run molecular dating analysis in mcmctree

Run with "usedata = 3" first to perform the approximate likelihood calculation. Then run with "usedata = 2" for final analysis.

### mcmctree.ctl (autocorrelated rates model)

```
          seed = -1
       seqfile = ./344g_8p_110422.phy
      treefile = ./whales.tree
       outfile = out

         ndata = 8
       seqtype = 0  * 0: nucleotides; 1:codons; 2:AAs
       usedata = 3    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = <0.3748  * safe constraint on root age, used if no fossil for root.

         model = 7    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 1    * alpha for gamma rates at sites
         ncatG = 4    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1    * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 20   * gamma prior for overall rates for genes
  sigma2_gamma = 2 0.5    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = 1: 0.1  0.1  0.1  0.01 .5  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr

         print = 1
        burnin = 10000
      sampfreq = 100
       nsample = 10000
```

### mcmctree.ctl (independent rates model)

```
          seed = -1
       seqfile = ./344g_8p_110422.phy
      treefile = ./whales.tree
       outfile = out

         ndata = 8
       seqtype = 0  * 0: nucleotides; 1:codons; 2:AAs
       usedata = 3    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = 2    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = <0.3748  * safe constraint on root age, used if no fossil for root.

         model = 7    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 1    * alpha for gamma rates at sites
         ncatG = 4    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1    * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 20   * gamma prior for overall rates for genes
  sigma2_gamma = 2 0.5    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = 1: 0.1  0.1  0.1  0.01 .5  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr

         print = 1
        burnin = 10000
      sampfreq = 100
       nsample = 10000
```

### treefile (max. root age 37.48 Ma)
```
(((NARight_consensus,Bow_consensus),((((Gray2_consensus,(Hump_consensus,Fin_consensus))'>0.0837',(Blue_consensus,Sei2_consensus)),Minke_consensus),PygmyWhale_consensus))'>0.182',Sperm_consensus)'>0.364<0.3748';
```
### treefile (max. root age 41.2 Ma)
```
(((NARight_consensus,Bow_consensus),((((Gray2_consensus,(Hump_consensus,Fin_consensus))'>0.0837',(Blue_consensus,Sei2_consensus)),Minke_consensus),PygmyWhale_consensus))'>0.182',Sperm_consensus)'>0.364<0.412';
```
### treefile (max. root age 52.4 Ma)
```
(((NARight_consensus,Bow_consensus),((((Gray2_consensus,(Hump_consensus,Fin_consensus))'>0.0837',(Blue_consensus,Sei2_consensus)),Minke_consensus),PygmyWhale_consensus))'>0.182',Sperm_consensus)'>0.364<0.524';
```
