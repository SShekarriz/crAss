#!/bin/bash

####################################################################################
DATE=`date '+%Y-%m-%d-%H:%M'`
exec > >(tee -i CRASS_PIPE_${DATE}.log)
exec 2>&1
####################################################################################


Metagenome=$PWD/Assembled_metagenomes
DonorB_reads="/dataone/shekas3/DonorB_metagenomics/Trimmed_interleave"
ucfmt1_reads="/dataone/shekas3/UCFMT1/Trim_interleaved"
comple_crAss=$PWD/crAss_genomes/complete_genomes.fa
blastOut=$PWD/crAss_in_metagenome
res=$PWD/results
coverage=/dataone/shekas3/crAssphage/crAss_from_reads
crAssFinder=$PWD/crAssFinder


mkdir -p $blastOut $res $coverage

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


echo "################### DONE ####################################"
date -u
echo "#############################################################"
