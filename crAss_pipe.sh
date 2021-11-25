ed "/^>/s/^>/>${fnaid};/g" $fna > ${fna%.fna}_edited.fna!/bin/bash

####################################################################################
DATE=`date '+%Y-%m-%d-%H:%M'`
exec > >(tee -i CRASS_PIPE_${DATE}.log)
exec 2>&1
####################################################################################

main=$PWD
bin=$PWD/bin
Metagenome=$PWD/Assembled_metagenomes
DonorB_reads="/dataone/shekas3/DonorB_metagenomics/Trimmed_interleave"
ucfmt1_reads="/dataone/shekas3/UCFMT1/Trim_interleaved"
comple_crAss=$PWD/crAss_genomes/complete_genomes.fa
blastOut=$PWD/crAss_in_metagenome
res=$PWD/results
coverage=/dataone/shekas3/crAssphage/crAss_from_reads
crAssFinder=$PWD/crAssFinder
phylo=$PWD/phylogeny
temp=$PWD/temp

mkdir -p $blastOut $res $coverage $temp

echo "#####################################"
echo "Finding crAss in metagenomes assembly"
echo "#####################################"

cd $Metagenome
ext="*_contigs.fasta"
for fna in $ext
do
        echo $fna
        fnaid=${fna%$ext}
        echo $fnaid
	
	echo "######### create blast dbs ################"	
	out=${fnaid}_db.ndb
        if [ -f $out ]; then
                echo "'$out' exists."
        else
                echo "'$out' is not found."
		makeblastdb -in $fna -dbtype nucl -out ${fnaid}_db
	fi

	echo "######### search crAss in metagenomes #########"
	out=$blastOut/${fnaid}.blastout
        if [ -f $out ]; then
                echo "'$out' exists."
        else
                echo "'$out' is not found."
                blastn -db ${fnaid}_db -query $comple_crAss -out $out -outfmt "6 std qcovhsp" -num_threads 25
        fi	

done
cd $main

echo "#####################################"
echo "Finding crAss in metagenomes reads"
echo "#####################################"

echo "in donor B samples:"
out=$coverage/DonorB_D_Oct17_inter.consensus.fa
if [ -f $out ]; then
	echo "'$out' exists."
else
	echo "'$out' is not found."
	$crAssFinder $DonorB_reads $coverage 20
fi

echo "in ucfmt1 samples:"
out=$coverage/PMCL657_S18_001_inter.consensus.fa
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        $crAssFinder $ucfmt1_reads $coverage 20
fi


echo "##########################################"
echo "Counting number of N's in consensus files"
echo "and select good assemblies for phylogeny"
echo "##########################################"

out=$res/consensus_quality.txt

if [ -f $out ]; then
        echo "'$out' exists."
else
	echo "'$out' is not found."
	cd $coverage
	ext="_inter.consensus.fa"
	for file in *${ext}
	do
		echo $file
		sample=${file%$ext}
		echo $sample
		python $bin/crAss_check.py $file ${sample}_inter.consensus.crAss.fa >> $res/consensus_quality.txt
	done
fi

echo "#####################################"
echo "building phylogenies"
echo "#####################################"

mkdir -p $phylo

echo "cat consensus crAss genomes:"
out=$phylo/complete_crAss.fa
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        cat $coverage/*.consensus.crAss.fa > $out
	# adding the crAss assembly to the alignments
	cat $comple_crAss >> $out
fi



out=$coverage/PMCL657_S18_001_inter.consensus.fa
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        $crAssFinder $ucfmt1_reads $coverage 20
fi


out=$phylo/pcrA_crAss.fa
if [ -f $out ]; then
	echo "'$out' exists."
else
	echo "'$out' is not found."

	pcr_dir=$phylo/pcr_amplicons
	mkdir -p $pcr_dir
	echo "First fix the fasta header name"
	cd $coverage
	for amplicon in *gretel/pcr*.fasta
	do
		echo $amplicon
		#removing the pattern + forward slash
		pcr=${amplicon/inter_gretel\//}
		echo $pcr
		cp $amplicon $pcr_dir/$pcr
		#change header name
		sed "/^>/s/^>/>${pcr%.fasta}_/g" $pcr_dir/$pcr > $pcr_dir/${pcr%.fasta}_fix.fa
	done
	echo "cat PCR amplicons:"
	cat $pcr_dir/*_pcrA_fix.fa > $phylo/pcrA_crAss.fa
	cat $pcr_dir/*_pcrB_fix.fa > $phylo/pcrB_crAss.fa
	cat $pcr_dir/*_pcrC_fix.fa > $phylo/pcrC_crAss.fa 
	cd $main
fi


out=$phylo/complete_crAss.aln
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
	mafft --auto $phylo/complete_crAss.fa > $out
        #muscle -in $phylo/complete_crAss.fa -out $out -maxiters 2 -diags
fi


out=$phylo/pcrA_crAss.aln
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        mafft --auto $phylo/pcrA_crAss.fa > $phylo/pcrA_crAss.aln
	mafft --auto $phylo/pcrB_crAss.fa > $phylo/pcrB_crAss.aln
	mafft --auto $phylo/pcrC_crAss.fa > $phylo/pcrC_crAss.aln
fi


out=$phylo/complete_crAss.tre
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        FastTree -nt $phylo/complete_crAss.aln > $out
fi



out=$phylo/pcrA_crAss.tre
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        FastTree -nt $phylo/pcrA_crAss.aln > $phylo/pcrA_crAss.tre
	FastTree -nt $phylo/pcrB_crAss.aln > $phylo/pcrB_crAss.tre
	FastTree -nt $phylo/pcrC_crAss.aln > $phylo/pcrC_crAss.tre
fi




echo "################### DONE ####################################"
date -u
echo "#############################################################"
