#!/bin/bash

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
coverage=$PWD/crAss_from_reads
cover=$PWD/crAss_coverage
crAssFinder=$PWD/crAssFinder
crAssMapper=$PWD/crAssMapper
phylo=$PWD/phylogeny
annot=$PWD/annotation
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
echo "####### parsing the blast results ####"
echo "#####################################"

out=$res/crass_contigs.txt
if [ -f $out ]; then
	echo "The file '$out' exists."
else
	echo "The file '$out' is not found."
	Rscript $bin/BlastnParser.R $blastOut $res
fi


echo "#####################################"
echo "## extracting crAss contigs ####"
echo "#####################################"

cd $Metagenome
ext="*_contigs.fasta"
for fna in $ext
do
        echo $fna
        fnaid=${fna%$ext}
        echo $fnaid

	out=$blastOut/${fnaid}_crass_contigs.fa
	if [ -f $out ]; then
        	echo "The file '$out' exists."
	else
        	echo "The file '$out' is not found."
		python $bin/extract_ids_fasta.py $fna $res/crass_contigs.txt \
		$blastOut/${fnaid}_crass_contigs.fa
		awk '/^>/{print ">'${fnaid}'__" ++i; next}{print}' < $blastOut/${fnaid}_crass_contigs.fa > $blastOut/${fnaid}_crass_edited_contigs.fa
	fi
done


echo "#####################################"
echo "## Manually fixing contigs into one ###"
echo "#####################################"


echo "DonorB_D_2012"
sample=DonorB_D_2012
out=$blastOut/${sample}_crass_edited_contigs_temp.fa.combine
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."

        input=$blastOut/${sample}_crass_edited_contigs.fa
        output=$blastOut/${sample}_crass_edited_contigs_temp.fa
        #perl $bin/fasta_reversecomplement $input > $output
        cp $input $output
	samtools faidx $output
        first=${output}.first
        second=${output}.second
	third=${output}.third
        combine=${output}.combine
        samtools faidx $output DonorB_D_2012__1 > $first
        samtools faidx $output DonorB_D_2012__2 > $second
	samtools faidx $output DonorB_D_2012__3 > $third
	#second contig need rev
        perl $bin/fasta_reversecomplement $second > ${second}.rev
	cat ${second}.rev $third $first > $combine
        rm $input
        union -filter $combine > $input
fi

echo "Re-assemble DonorB 2012 using mapped reads"
B2012=$main/fixing_donorB2012_assembly/DonorB_2012_smapped_assembly/Edited_contigs.fa
out=${B2012%.fa}_temp.fa.combine
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
	input=${B2012}
	output=${B2012%.fa}_temp.fa
	cp $input $output
	first=${output}.first
	second=${output}.second
	third=${output}.third
	combine=${output}.combin
	samtools faidx $output DonorB_D_2012__1 > $first
	samtools faidx $output DonorB_D_2012__2 > $second
	samtools faidx $output DonorB_D_2012__3 > $third	
	cat $first $second $third > $combine
	rm $input
	union -filter $combine > $input
fi


echo "DonorB_D_2013"
sample=DonorB_D_2013
out=$blastOut/${sample}_crass_edited_contigs_temp.fa.combine
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
        input=$blastOut/${sample}_crass_edited_contigs.fa
        output=$blastOut/${sample}_crass_edited_contigs_temp.fa
        #perl $bin/fasta_reversecomplement $input > $output
        cp $input $output
        samtools faidx $output
        first=${output}.first
        second=${output}.second
        third=${output}.third
        combine=${output}.combine
        samtools faidx $output DonorB_D_2013__1 > $first
        # second contig only has 1333 bp out of a 63489 as crAss match
	# look into this further later, is this a host dna?
	samtools faidx $output DonorB_D_2013__2 > $second
        cat $first > $combine
        rm $input
        cp $combine $input
fi


echo "SHCM2_D"
sample=SHCM2_D
out=$blastOut/${sample}_crass_edited_contigs_temp.fa.combine
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
        input=$blastOut/${sample}_crass_edited_contigs.fa
        output=$blastOut/${sample}_crass_edited_contigs_temp.fa
        #perl $bin/fasta_reversecomplement $input > $output
        cp $input $output
	samtools faidx $output
        first=${output}.first
        second=${output}.second
        combine=${output}.combine
        samtools faidx $output SHCM2_D__1 > $first
        samtools faidx $output SHCM2_D__2 > $second
        cat $first $second > $combine
        rm $input
        union -filter $combine > $input
fi


echo "SHCM4_D"
sample=SHCM4_D
out=$blastOut/${sample}_crass_edited_contigs_temp.fa.combine
if [ -f $out ]; then
	echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
	input=$blastOut/${sample}_crass_edited_contigs.fa
	output=$blastOut/${sample}_crass_edited_contigs_temp.fa
	perl $bin/fasta_reversecomplement $input > $output
	samtools faidx $output
	first=${output}.first
	second=${output}.second
	combine=${output}.combine
	samtools faidx $output SHCM4_D__2 > $first
	samtools faidx $output SHCM4_D__1 > $second
	cat $first $second > $combine
	rm $input
	union -filter $combine > $input
fi


echo "818_ACTGAT"
sample=818_ACTGAT
out=$blastOut/${sample}_crass_edited_contigs_temp.fa.combine
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
        input=$blastOut/${sample}_crass_edited_contigs.fa
        output=$blastOut/${sample}_crass_edited_contigs_temp.fa
        #perl $bin/fasta_reversecomplement $input > $output
        cp $input $output
	samtools faidx $output
        first=${output}.first
        second=${output}.second
        combine=${output}.combine
        samtools faidx $output 818_ACTGAT__1 > $first
        perl $bin/fasta_reversecomplement $first > ${first}.rev
	samtools faidx $output 818_ACTGAT__2 > $second
        cat $second ${first}.rev > $combine
        rm $input
        union -filter $combine > $input
fi


echo "823_TGACCA"
sample=823_TGACCA
out=$blastOut/${sample}_crass_edited_contigs_temp.fa.combine
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
        input=$blastOut/${sample}_crass_edited_contigs.fa
        output=$blastOut/${sample}_crass_edited_contigs_temp.fa
        perl $bin/fasta_reversecomplement $input > $output
        #cp $input $output
	samtools faidx $output
        first=${output}.first
        second=${output}.second
        combine=${output}.combine
        samtools faidx $output 823_TGACCA__1 > $first
        cat $first > $combine
        rm $input
        cp $combine $input
fi


echo "#############################################################"
echo " crAss genome annotations:"
echo " using assembled contigs from metagenomes"
echo "#############################################################"

# Multiple genome annotation
annot_M=$annot/Multiple
# Single genome annotation
annot_S=$annot/Single
mkdir -p $annot_M $annot_S

crass_gb=$main/crAss_genomes/NC_024711.gb
crass_org=$main/crAss_genomes/NC_024711.fasta
crass_ref=$annot_S/NC_024711_rev.fa


echo "############## 2. Multiple genome annotations ########################"

All=$annot_M/All_crass
DonB=$annot_M/DonorB_crass
FMT=$annot_M/FMT_crass
mkdir -p $All $DonB $FMT

suffix="_crass_edited_contigs.fa"
g0=$main/crAss_genomes/NC_024711_rev.fa
g1=$blastOut/DonorB_D_2012${suffix}
g2=$blastOut/DonorB_D_2013${suffix}
g3=$blastOut/DonorB_D_2016${suffix}
g4=$blastOut/DonorB_D_May17${suffix}
g5=$blastOut/818_ACTGAT${suffix}
g6=$blastOut/822_CGATGT${suffix}
g7=$blastOut/823_TGACCA${suffix}
g8=$blastOut/SHCM2_D${suffix}
g9=$blastOut/SHCM4_D${suffix}


echo "All crAss genomes:"
data=All
out=$All/${data}-cogs.tsv
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
	input=$All/${data}.fa
	echo "order of genome based on whole-genome alignment"
	cat $g2 $g3 $g1 $g0 $g4 $g5 $g9 $g8 $g7 $g6 > $input
	echo "########### circular alignment #################"
	/dataone/common/software/CSA/CSA R $input
	echo "########### annotation #################"
	prokka --cpus 15 --proteins $crass_gb --outdir $All/${data}_prokka \
	--prefix ${data} $All/${data}-Rotated.fasta
	echo "########### All-vs-all alignment #################"
	minimap2 -X -N 50 -p 0.1 -c $All/${data}_prokka/${data}.fna $All/${data}_prokka/${data}.fna > $All/${data}.paf
	echo "########### GC-content #################"
	$bin/seq-gc -Nbw 50 $All/${data}_prokka/${data}.fna > $All/${data}-gc.tsv
	echo "##### cluster protein sequences into orthogroups #######"
	mmseqs easy-cluster $All/${data}_prokka/${data}.faa $All/${data}-mmseqs /tmp -e 1e-5 -c 0.7
	$bin/cluster-ids -t "cog%03d" < $All/${data}-mmseqs_cluster.tsv > $All/${data}-cogs.tsv

fi


echo "########## using the aligned genomes for single genomes ###"
echo $All/${data}-Rotated.fasta
align_G=$All/${data}_aligned_genomes
mkdir -p $align_G

out=$align_G/DonorB_D_2012.fa
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
	input=$All/${data}-Rotated.fasta
	samtools faidx $input NC_024711.1 > $align_G/NC_024711.1_rev.fa
	samtools faidx $input DonorB_D_2013__1 > $align_G/DonorB_D_2013.fa
	samtools faidx $input DonorB_D_2016__1 > $align_G/DonorB_D_2016.fa
	samtools faidx $input DonorB_D_May17__1 > $align_G/DonorB_D_May17.fa
	samtools faidx $input 818_ACTGAT__2 > $align_G/818_ACTGAT.fa
	samtools faidx $input 822_CGATGT__1 > $align_G/822_CGATGT.fa
	samtools faidx $input 823_TGACCA__1 > $align_G/823_TGACCA.fa
	samtools faidx $input SHCM2_D__1 > $align_G/SHCM2_D.fa
	samtools faidx $input SHCM4_D__2 > $align_G/SHCM4_D.fa
	samtools faidx $input DonorB_D_2012__2 > $align_G/DonorB_D_2012.fa
fi


echo "DonorB crAss genomes:"
data=DonB
out=$DonB/${data}-cogs.tsv
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        All=$DonB
	input=$All/${data}.fa
        cat $g0 $g2 $g3 $g4 > $input
        echo "########### circular alignment #################"
        /dataone/common/software/CSA/CSA R $input
        echo "########### annotation #################"
        prokka --cpus 15 --proteins $crass_gb --outdir $All/${data}_prokka \
        --prefix ${data} $All/${data}-Rotated.fasta
        echo "########### All-vs-all alignment #################"
        minimap2 -X -N 50 -p 0.1 -c $All/${data}_prokka/${data}.fna $All/${data}_prokka/${data}.fna > $All/${data}.paf
        echo "########### GC-content #################"
        $bin/seq-gc -Nbw 50 $All/${data}_prokka/${data}.fna > $All/${data}-gc.tsv
        echo "##### cluster protein sequences into orthogroups #######"
        mmseqs easy-cluster $All/${data}_prokka/${data}.faa $All/${data}-mmseqs /tmp -e 1e-5 -c 0.7
        $bin/cluster-ids -t "cog%03d" < $All/${data}-mmseqs_cluster.tsv > $All/${data}-cogs.tsv

fi


echo "FMT crAss genomes:"
data=FMT
out=$FMT/${data}-cogs.tsv
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        All=$FMT
        input=$All/${data}.fa
        cat $g2 $g5 $g6 $g7 > $input
        echo "########### circular alignment #################"
        /dataone/common/software/CSA/CSA R $input
        echo "########### annotation #################"
        prokka --cpus 15 --proteins $crass_gb --outdir $All/${data}_prokka \
        --prefix ${data} $All/${data}-Rotated.fasta
        echo "########### All-vs-all alignment #################"
        minimap2 -X -N 50 -p 0.1 -c $All/${data}_prokka/${data}.fna $All/${data}_prokka/${data}.fna > $All/${data}.paf
        echo "########### GC-content #################"
        $bin/seq-gc -Nbw 50 $All/${data}_prokka/${data}.fna > $All/${data}-gc.tsv
        echo "##### cluster protein sequences into orthogroups #######"
        mmseqs easy-cluster $All/${data}_prokka/${data}.faa $All/${data}-mmseqs /tmp -e 1e-5 -c 0.7
        $bin/cluster-ids -t "cog%03d" < $All/${data}-mmseqs_cluster.tsv > $All/${data}-cogs.tsv

fi


echo "############## 1. Single genome annotations ########################"


out=$annot_S/DonorB_D_2013_prokka/DonorB_D_2013.faa
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        for sample in NC_024711.1_rev DonorB_D_2012 DonorB_D_2013 DonorB_D_2016 DonorB_D_May17 818_ACTGAT 822_CGATGT 823_TGACCA SHCM2_D SHCM4_D
        do
                echo $sample
                org=$align_G/${sample}.fa
                crass=$annot_S/${sample}.fa
                cp $org $crass
                prokka --cpus 15 --proteins $crass_gb --outdir $annot_S/${sample}_prokka \
                --prefix $sample $crass
        done
fi


echo "#####################################"
echo "## getting coverage information  ####"
echo "#####################################"


echo "1. for donorB samples:"
raw="/dataone/shekas3/DonorB_metagenomics/Trimmed_interleave"
crass_ref=$annot_S/NC_024711.1_rev.fa
echo $crass_ref
crAss_B2013=$annot_S/DonorB_D_2013.fa
echo $crAss_B2013
for sample in DonorB_D_2012 DonorB_D_2013 DonorB_D_2016 DonorB_D_May17 DonorB_D_Oct17
do
        echo $sample
        raw_sample=$raw/${sample}_inter.fastq
        echo "raw reads:$raw_sample"
        crass=$annot_S/${sample}.fa
        echo "assembled genome:$crass"
        $crAssMapper $raw_sample $cover 30 $crass $sample $crass_ref $crAss_B2013
done

echo "2. for FMT samples:"
raw="/dataone/shekas3/UCFMT1/Trim_interleaved"
for sample in 818_ACTGAT 822_CGATGT 823_TGACCA
do
        echo $sample
        raw_sample=$raw/${sample}_001_inter.fastq
        echo $raw_sample
        crass=$annot_S/${sample}.fa
        echo $crass
        $crAssMapper $raw_sample $cover 30 $crass $sample $crass_ref $crAss_B2013
done


echo "3. for SHCM samples:"
raw="/datastore/shekas3/SHCM/SHCMall/SHCM_All_trim_interleaved"
for sample in SHCM2_D SHCM4_D
do
        echo $sample
        raw_sample=$raw/${sample%_D}_Stool_001_inter.fastq
        echo $raw_sample
        crass=$annot_S/${sample}.fa
        echo $crass
        $crAssMapper $raw_sample $cover 30 $crass $sample $crass_ref $crAss_B2013
done


echo "#################################################################"
echo " SNP analysis- using assembled crAss for each sample via breseq:"
echo "#############################################################"


echo "########### breseq #################"
raw="/dataone/shekas3/DonorB_metagenomics/Trimmed_interleave"
out=$cover/DonorB_D_2012_breseq.tsv
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        for sample in DonorB_D_2012 DonorB_D_2016 DonorB_D_May17 DonorB_D_2013
        do
                echo $sample
                raw_sample=$raw/${sample}_inter.fastq
                echo $raw_sample
		output=$cover/${sample}_breseq
                echo $output
		gbf=$annot_S/${sample}_prokka/${sample}.gbf
                echo $gbf
		breseq -j 30 -p -o $output -r $gbf $raw_sample
		gd=$cover/${sample}_breseq.gd
		cp $output/output/output.gd $gd
		gdtools ANNOTATE -o $cover/${sample}_breseq.tsv -f TSV -r $gbf $gd 
        done

fi

out=$cover/DonorB_D_Oct17_B2013_breseq.tsv
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        for sample in DonorB_D_Oct17 DonorB_D_2012 DonorB_D_2016 DonorB_D_May17 DonorB_D_2013
        do
                echo $sample
                raw_sample=$raw/${sample}_inter.fastq
                echo $raw_sample
                output=$cover/${sample}_B2013_breseq
                echo $output
                gbf=$annot_S/DonorB_D_2013_prokka/DonorB_D_2013.gbf
                echo $gbf
                breseq -j 30 -p -o $output -r $gbf $raw_sample
                gd=$cover/${sample}_B2013_breseq.gd
                cp $output/output/output.gd $gd
                gdtools ANNOTATE -o $cover/${sample}_B2013_breseq.tsv -f TSV -r $gbf $gd
        done

fi


echo "########### breseq #################"
raw="/dataone/shekas3/UCFMT1/Trim_interleaved"
out=$cover/818_ACTGAT_breseq.tsv
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        for sample in 818_ACTGAT 822_CGATGT 823_TGACCA
        do
                echo $sample
                raw_sample=$raw/${sample}_001_inter.fastq
                echo $raw_sample
                output=$cover/${sample}_breseq
                echo $output
                gbf=$annot_S/${sample}_prokka/${sample}.gbf
                echo $gbf
                breseq -j 30 -p -o $output -r $gbf $raw_sample
                gd=$cover/${sample}_breseq.gd
                cp $output/output/output.gd $gd
                gdtools ANNOTATE -o $cover/${sample}_breseq.tsv -f TSV -r $gbf $gd
        done

fi


out=$cover/818_ACTGAT_B2013_breseq.tsv
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        for sample in 818_ACTGAT 822_CGATGT 823_TGACCA
        do
                echo $sample
                raw_sample=$raw/${sample}_001_inter.fastq
                echo $raw_sample
                output=$cover/${sample}_B2013_breseq
                echo $output
                gbf=$annot_S/DonorB_D_2013_prokka/DonorB_D_2013.gbf
                echo $gbf
                breseq -j 30 -p -o $output -r $gbf $raw_sample
                gd=$cover/${sample}_B2013_breseq.gd
                cp $output/output/output.gd $gd
                gdtools ANNOTATE -o $cover/${sample}_B2013_breseq.tsv -f TSV -r $gbf $gd
        done

fi

echo "########### breseq #################"
raw="/datastore/shekas3/SHCM/SHCMall/SHCM_All_trim_interleaved"
out=$cover/SHCM4_D_breseq.tsv
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        for sample in SHCM2_D SHCM4_D
        do
                echo $sample
                raw_sample=$raw/${sample%_D}_Stool_001_inter.fastq
                echo $raw_sample
                output=$cover/${sample}_breseq
                echo $output
                gbf=$annot_S/${sample}_prokka/${sample}.gbf
                echo $gbf
                breseq -j 30 -p -o $output -r $gbf $raw_sample
                gd=$cover/${sample}_breseq.gd
                cp $output/output/output.gd $gd
                gdtools ANNOTATE -o $cover/${sample}_breseq.tsv -f TSV -r $gbf $gd
        done

fi


out=$cover/SHCM4_D_B2013_breseq.tsv
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        for sample in SHCM2_D SHCM4_D
        do
                echo $sample
                raw_sample=$raw/${sample%_D}_Stool_001_inter.fastq
                echo $raw_sample
                output=$cover/${sample}_B2013_breseq
                echo $output
                gbf=$annot_S/DonorB_D_2013_prokka/DonorB_D_2013.gbf
                echo $gbf
                breseq -j 30 -p -o $output -r $gbf $raw_sample
                gd=$cover/${sample}_B2013_breseq.gd
                cp $output/output/output.gd $gd
                gdtools ANNOTATE -o $cover/${sample}_B2013_breseq.tsv -f TSV -r $gbf $gd
        done

fi


echo "#####################################"
echo "Mergin all _B2013_BRESEQ outputs"
echo "#####################################"

out=$cover/AllMerged_B2013_breseq.tsv
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
	gbf=$annot_S/DonorB_D_2013_prokka/DonorB_D_2013.gbf
	gdtools ANNOTATE -o $cover/AllMerged_B2013_breseq.tsv -f TSV -r $gbf $cover/*_B2013_breseq.gd
fi

echo "#####################################"
echo "building phylogenies"
echo "#####################################"

mkdir -p $phylo

echo "cat consensus crAss genomes:"
out=$phylo/complete_crAss.tre
if [ -f $out ]; then
        echo "'$out' exists."
else
        echo "'$out' is not found."
        echo $annot_S/*.fa
	cat $annot_S/*.fa > $phylo/complete_crAss.fa
	mafft --auto $phylo/complete_crAss.fa > $phylo/complete_crAss.aln
	FastTree -nt $phylo/complete_crAss.aln > $phylo/complete_crAss.tre
fi

:<<"END"
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


END


echo "################### DONE ####################################"
date -u
echo "#############################################################"
