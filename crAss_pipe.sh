#!/bin/bash

####################################################################################
DATE=`date '+%Y-%m-%d-%H:%M'`
exec > >(tee -i CRASS_PIPE_${DATE}.log)
exec 2>&1
####################################################################################


Metagenome=$PWD/Assembled_metagenomes
comple_crAss=$PWD/crAss_genomes/complete_genomes.fa
blastOut=$PWD/crAss_in_metagenome
res=$PWD/results

mkdir -p $blastOut $res

echo "############################"
echo "Finding crAss in metagenomes"
echo "############################"

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


echo "################### DONE ####################################"
date -u
echo "#############################################################"
