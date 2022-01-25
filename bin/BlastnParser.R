#!/usr/bin/env Rscript
# SShekarriz Jul05,2021
# Making a Bin info file to be imported into anvio
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = "crAss_in_metagenome"
  args[2] = "results"
}

#################
library(tidyverse)
#################
# input directory
path=args[1]
# output directory
outpath=args[2] 



path="crAss_in_metagenome"
outpath="results"


patt=".blastout"
COLS <- c("Metagenome","qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
"qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp")
data.frame(sample.id = paste(dir(path, 
                                 pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(sample.id, ~ read_tsv(
  file.path(path, .), 
            col_names = F))) %>%
         unnest() %>%
  mutate(sample.id = gsub(patt, "", sample.id)) -> blastOuts
colnames(blastOuts) <- COLS

#clean up the result and add more var based on contig and crAss info
blastOuts %>%
  mutate(temp=sseqid) %>%
  separate(temp, c("NodeID", "temp2"), sep="_length_") %>%
  separate(temp2, c("ContigLength", "ContigCov"), sep = "_cov_") %>%
  mutate(ContigLength= as.numeric(ContigLength)) %>%
  # define crAss genomes
  mutate(crAssGenomes= case_when(
    qseqid == "MT074136.1" ~ paste("MT074136.1__Bacteroides_phage_DAC15"),
    qseqid == "MT074138.1" ~ paste("MT074138.1__Bacteroides_phage_DAC17"),
    qseqid == "NC_024711.1" ~ paste("NC_024711.1__Uncultured_crAssphage"),
    qseqid == "MH675552.1" ~ paste("MH675552.1__Bacteroides_phage_crAss001"),
    TRUE ~ paste("Others")))-> crAss_contigs


out=paste(outpath, "crAss_hits.txt", sep = "/")
crAss_contigs %>%
  group_by(Metagenome, crAssGenomes) %>% 
  summarise(total_length= sum(length)) %>%
  filter(total_length >= 40000) -> crAss_hits
write.table(crAss_hits, out,
          sep = "\t", row.names = F, quote = F)


crAss_contigs %>%
  left_join(crAss_hits %>% mutate(Hit_type= paste("crAss_hit")),
            by = c("Metagenome", "crAssGenomes")) %>%
  filter(Hit_type == "crAss_hit") %>%
  group_by(Metagenome, crAssGenomes, sseqid) %>% 
  summarise(total_length= sum(length)) %>%
  #remove any hit shorter than 1000bp in a contig
  filter(total_length >= 1000) %>%
  mutate(contig_type= case_when(total_length >= 80000 ~ paste("Complete"),
                                between(total_length, 10000, 80000) ~ paste("Semi_complete"),
                                total_length <= 10000 ~ paste("Fragmented"),
                                TRUE ~ paste("Other_cases")))-> contig_hits

out=paste(outpath, "contig_hits.txt", sep = "/")
contig_hits %>%
  rename(contig= sseqid) %>%
  group_by(Metagenome) %>%
  write.table(out,
          sep = "\t", row.names = F, quote = F)


out=paste(outpath, "crass_contigs.txt", sep = "/")
contig_hits %>%
  ungroup() %>%
  filter(crAssGenomes == "NC_024711.1__Uncultured_crAssphage") %>%
  filter(Metagenome != "SHCM1_D") %>%
  #filter(contig_type %in% c("Complete", "Semi_complete")) %>%
  select(sseqid) %>%
  write.table(out,
              sep = "\t", row.names = F, col.names = F, quote = F)



