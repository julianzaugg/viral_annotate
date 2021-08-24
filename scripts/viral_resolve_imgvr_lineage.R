#!/usr/bin/env Rscript


## Takes BLAST (format 6) results from BLASTing against IMG/VR and
## resolves, where possible, the consensus taxonomy string for
## each query

## $1: Format 6 BLAST results. Tab-delimited.
## $2: GFF file from Prodigal for BLASTed sequences
## $3: IMG/VR reference taxonomy file. Tab-delimited. Two columns: UViG ID, Taxonomic classification
## $4: Minimum fraction of genes from a viral sequence that are required to have a significant hit to resolve taxonomy for the sequence
## $5: Minimum fraction of genes WITH significant hit (passing $4) that need to agree at a taxonomy level to assign taxonomy at that level
## $6: Output directory

args = commandArgs(trailingOnly=TRUE)

MIN_FRACTION_GENES_HIT = as.numeric(args[4])
MIN_FRACTION_MAJORITY_RULE = as.numeric(args[5])
# MIN_FRACTION_GENES_HIT = 0.3
# MIN_FRACTION_MAJORITY_RULE = 0.5

library(tidyverse)

# -----------------------------------------------------------------------
# Load BLAST data
blast_results.df <- read.table(args[1], sep="\t", header=T)
# blast_results.df <- read.table("viral_proteins_imgvr_diamond_blast.tsv", sep = "\t", header = T)

# Load GFF file from prodigal
# Note - It is likely that not all genes will be in the BLAST results (no significant hit),
# indeed, may be the case no genes from a sequence will have a significant hit
gff.df <- read.delim(args[2], sep="\t", header=F, comment.char="#")
# gff.df <- read.delim("all_samples_viral_sequences.gff", header=F, comment.char="#")
names(gff.df) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

# Create protein ID column (seq_id + prot ID)
gff.df$protein_id <- with(gff.df, paste0(seqid,"_", gsub("ID=[0-9]{1,10}_(.*?);.*", "\\1", gff.df$attributes)))

# Load taxonomy reference file.
# Note - not all UViGs have a lineage
imgvr_taxonomy.df <- read.delim(args[3], sep="\t", header=T)
# imgvr_taxonomy.df <- read.delim("IMGVR_taxonomy.tsv", header = T, sep = "\t")
names(imgvr_taxonomy.df) <- c("UViG", "Taxonomic_classification")
rownames(imgvr_taxonomy.df) <- imgvr_taxonomy.df$UViG

# -----------------------------------------------------------------------
# Clean Subject_ID in the BLAST results to be able to search against taxonomy reference. This should be the UViG ID
# Everything after first '|' assumed to be specific reference protein ID
blast_results.df$Subject_ID_cleaned <- gsub("\\|.*","", blast_results.df$Subject_ID)

# Get list of subject IDs
subject_ids.v <- unique(blast_results.df$Subject_ID_cleaned)

# Subjects that were not found in the IMG/VR taxonomy reference
missing_subject_ids.v <- subject_ids.v[!subject_ids.v %in% rownames(imgvr_taxonomy.df)]

# Even after basic cleaning, some entries won't match exactly. We need to search for any
# partial match to the subject ID to resolve these
missing_subject_matches.l <- list()
for (sid in missing_subject_ids.v){
  missing_subject_matches.l[sid] <- list(grep(paste0(sid, "[\\|\\.]"), imgvr_taxonomy.df$UViG, value =T))
}
# FIXME - handle no matches? Simply remove from BLAST results?

# Assume the first entry/match is appropriate for lineage information (not ideal if more than one)
blast_results.df$Subject_ID_cleaned[blast_results.df$Subject_ID_cleaned %in% missing_subject_ids.v] <-
  as.character(lapply(blast_results.df$Subject_ID_cleaned[blast_results.df$Subject_ID_cleaned %in% missing_subject_ids.v],
       function(x) missing_subject_matches.l[[x]][[1]]))

# Get the lineage information for the final set of cleaned Subject IDs. These should match a UViG ID
subject_ids.v <- unique(blast_results.df$Subject_ID_cleaned)
# summary(subject_ids.v %in% imgvr_taxonomy.df$UViG)
subjects_imgvr_taxonomy.df <- imgvr_taxonomy.df[subject_ids.v,]
subjects_imgvr_taxonomy.df <- subjects_imgvr_taxonomy.df[!is.na(subjects_imgvr_taxonomy.df$UViG),]

# summary(blast_results.df$Subject_ID_cleaned %in% subjects_imgvr_taxonomy.df$UViG)
# ------------------------------------------
# Create lineage table for each gene. Sequence ID (virus) - Protein ID (query) - Subject ID (BLAST hit) - Subject ID  cleaned (BLAST hit/UViG ID) - lineage (IMG/VR)...
protein_lineages.df <- gff.df[,c("seqid","protein_id")]
# Combine with BLAST results (note - not all genes present)
protein_lineages.df <- left_join(protein_lineages.df, blast_results.df[,c("Query_ID", "Subject_ID", "Subject_ID_cleaned")], by= c("protein_id" = "Query_ID"))
# Combine with lineage information results
protein_lineages.df <- left_join(protein_lineages.df, subjects_imgvr_taxonomy.df, by = c("Subject_ID_cleaned" = "UViG"))
protein_lineages.df[protein_lineages.df == ""] <- NA
# protein_lineages.df[!complete.cases(protein_lineages.df),]

# Separate taxonomies
protein_lineages.df <- separate(protein_lineages.df, "Taxonomic_classification", into = c("Root", "Domain", "Phylum", "Class", "Order", "Family","Genus", "Species"), remove =F, sep = ";")

# Splitting taxa strings that are not specified at certain taxa levels will produce NA entries at those levels.
# NA entries should be changed to "Unassigned". Also account for empty "" levels/strings
protein_lineages.df[,c("Root", "Domain", "Phylum", "Class", "Order", "Family","Genus", "Species")][protein_lineages.df[,c("Root", "Domain", "Phylum", "Class", "Order", "Family","Genus", "Species")] == ""] <- "Unassigned"
protein_lineages.df[,c("Root", "Domain", "Phylum", "Class", "Order", "Family","Genus", "Species")][is.na(protein_lineages.df[,c("Root", "Domain", "Phylum", "Class", "Order", "Family","Genus", "Species")])] <- "Unassigned"

# Construct taxonomy strings
protein_lineages.df$taxonomy_species <- with(protein_lineages.df, paste(Root, Domain, Phylum, Class, Order, Family, Genus, Species, sep =";"))
protein_lineages.df$taxonomy_genus <- with(protein_lineages.df, paste(Root, Domain, Phylum, Class, Order, Family, Genus, sep =";"))
protein_lineages.df$taxonomy_family <- with(protein_lineages.df, paste(Root, Domain, Phylum, Class, Order, Family, sep =";"))
protein_lineages.df$taxonomy_order <- with(protein_lineages.df, paste(Root, Domain, Phylum, Class, Order, sep =";"))
protein_lineages.df$taxonomy_class <- with(protein_lineages.df, paste(Root, Domain, Phylum, Class, sep =";"))
protein_lineages.df$taxonomy_phylum <- with(protein_lineages.df, paste(Root, Domain, Phylum, sep =";"))
protein_lineages.df$taxonomy_domain <- with(protein_lineages.df, paste(Root, Domain, sep =";"))
protein_lineages.df$taxonomy_root <- with(protein_lineages.df, paste(Root, sep =";"))


# ------------------------------------------

# Count the number of genes per sequence
# protein_lineages.df will have the same entries as gff.df
Sequence_gene_count.df <-
  protein_lineages.df %>%
  dplyr::group_by(seqid) %>%
  dplyr::summarise(Number_of_Genes = n())

Sequence_gene_count_with_hit.df <-
  protein_lineages.df %>%
  dplyr::group_by(seqid) %>%
  dplyr::filter(!is.na(Subject_ID_cleaned)) %>%
  dplyr::summarise(Number_of_genes_with_hit = n())

Sequence_gene_count.df <- left_join(Sequence_gene_count.df, Sequence_gene_count_with_hit.df)
Sequence_gene_count.df <-
  Sequence_gene_count.df %>%
  mutate(Fraction_of_genes_with_hit = Number_of_genes_with_hit/Number_of_Genes)

# Get the sequences to resolve lineages for, i.e. those with the minimum fraction of genes with a hit
sequences_to_resolve.v <- Sequence_gene_count.df[Sequence_gene_count.df$Fraction_of_genes_with_hit >= MIN_FRACTION_GENES_HIT,]$seqid

# ------------------------------------------------------------------
# Now resolve the majority rules taxonomy

# For each sequencing, calculate the number of genes, and fraction of total,
# with the same taxonomy string at each taxonomy level
sequence_gene_lineage_fractions.df <-
  protein_lineages.df[,c("seqid","protein_id", "Subject_ID_cleaned",
                         "taxonomy_species","taxonomy_genus","taxonomy_family","taxonomy_order","taxonomy_class","taxonomy_phylum", "taxonomy_domain", "taxonomy_root")] %>%
  # dplyr::filter(seqid %in% sequences_to_resolve.v) %>% # Filter to those sequences we are resolving lineages for
  dplyr::group_by(seqid) %>%
  dplyr::mutate(Number_of_genes = n()) %>% # Count total number of genes
  dplyr::filter(!is.na(Subject_ID_cleaned)) %>% # Remove genes that have no corresponding UViG hit.
  # Count genes at each taxonomy level
  dplyr::group_by(seqid, taxonomy_species) %>%
  dplyr::mutate(N_genes_species = n(), Fraction_species = n()/Number_of_genes) %>%
  dplyr::group_by(seqid, taxonomy_genus) %>%
  dplyr::mutate(N_genes_genus = n(), Fraction_genus = n()/Number_of_genes) %>%
  dplyr::group_by(seqid, taxonomy_family) %>%
  dplyr::mutate(N_genes_family = n(), Fraction_family = n()/Number_of_genes) %>%
  dplyr::group_by(seqid, taxonomy_order) %>%
  dplyr::mutate(N_genes_order = n(), Fraction_order = n()/Number_of_genes) %>%
  dplyr::group_by(seqid, taxonomy_class) %>%
  dplyr::mutate(N_genes_class = n(), Fraction_class = n()/Number_of_genes) %>%
  dplyr::group_by(seqid, taxonomy_phylum) %>%
  dplyr::mutate(N_genes_phylum = n(), Fraction_phylum = n()/Number_of_genes) %>%
  dplyr::group_by(seqid, taxonomy_domain) %>%
  dplyr::mutate(N_genes_domain = n(), Fraction_domain = n()/Number_of_genes) %>%
  dplyr::group_by(seqid, taxonomy_root) %>%
  dplyr::mutate(N_genes_root = n(), Fraction_root = n()/Number_of_genes) %>%
  as.data.frame()

# Save this output
# TODO - also note those sequences with no BLAST/lineage results?
sequence_gene_lineage_fractions.df <- left_join(sequence_gene_lineage_fractions.df, Sequence_gene_count.df[,c("seqid", "Number_of_genes_with_hit", "Fraction_of_genes_with_hit")], by = "seqid")
# write.table(sequence_gene_lineage_fractions.df,
            # file="sequence_gene_lineage_fractions.tsv",
            # sep="\t", quote=F, row.names=F, col.names=T)
write.table(sequence_gene_lineage_fractions.df,
            file=paste0(args[6], "/sequence_gene_lineage_fractions.tsv"),
            sep="\t", quote=F, row.names=F, col.names=T)

# Sequences (seqid) that are missing from the results after the following
# did not have a majority of genes agreeing on a lineage, i.e. "No majority lineage" or simply "Virus"
sequence_resolved_lineage_gene_fractions.df <-
  sequence_gene_lineage_fractions.df %>%
  # dplyr::select(-protein_id, -Subject_ID_cleaned, -Number_of_genes, -Number_of_genes_with_hit, -Fraction_of_genes_with_hit) %>%
  dplyr::select(-protein_id, -Subject_ID_cleaned) %>%
  dplyr::select(-c(N_genes_species,N_genes_genus,N_genes_family,N_genes_order,N_genes_class,N_genes_phylum, N_genes_domain, N_genes_root)) %>%
  unique() %>%
  # gather(key = "N_gene_group", value = "N_gene_value", c(N_genes_species,N_genes_genus,N_genes_family,N_genes_order,N_genes_class,N_genes_phylum, N_genes_domain, N_genes_root)) %>%
  gather(key = "Fraction_group", value = "Fraction_value", c(Fraction_species,Fraction_genus,Fraction_family,Fraction_order,Fraction_class,Fraction_phylum, Fraction_domain, Fraction_root)) %>%
  # gather(key = "Taxonomy_group", value = "Taxonomy_value", c(taxonomy_species,taxonomy_genus,taxonomy_family,taxonomy_order,taxonomy_class,taxonomy_phylum, taxonomy_domain, taxonomy_root)) %>%
  dplyr::filter(Fraction_value > MIN_FRACTION_MAJORITY_RULE) %>%
  dplyr::group_by(seqid) %>%
  # dplyr::slice_min(n=1, order_by = Fraction_value) %>%
  dplyr::arrange(factor(Fraction_group, levels = c("Fraction_species","Fraction_genus","Fraction_family","Fraction_order","Fraction_class","Fraction_phylum", "Fraction_domain", "Fraction_root"))) %>%
  dplyr::slice(n=1) %>%
  as.data.frame()

# Clean up summary table and define the final majority taxonomy
sequence_resolved_lineage_gene_fractions.df$Majority_taxonomy_level <- gsub("Fraction","taxonomy", sequence_resolved_lineage_gene_fractions.df$Fraction_group)


sequence_resolved_lineage_gene_fractions.df$Majority_taxonomy <-
  unlist(lapply(sequence_resolved_lineage_gene_fractions.df$seqid,
       function(x) {
         ltl <- sequence_resolved_lineage_gene_fractions.df[sequence_resolved_lineage_gene_fractions.df$seqid == x,"Majority_taxonomy_level"]
         sequence_resolved_lineage_gene_fractions.df[sequence_resolved_lineage_gene_fractions.df$seqid == x,ltl]
       }
       )
       )
sequence_resolved_lineage_gene_fractions.df$Fraction_group <- NULL
names(sequence_resolved_lineage_gene_fractions.df)[names(sequence_resolved_lineage_gene_fractions.df) == "Fraction_value"] <- "Fraction_genes_at_majority_level"

# ------------------------------------------------------------
# Add sequences that are missing to the table

# Create table for missing sequences
temp <- data.frame("seqid" = unique(gff.df$seqid),"taxonomy_species" = NA,"taxonomy_genus" = NA,"taxonomy_family" = NA,"taxonomy_order" = NA,"taxonomy_class" = NA,"taxonomy_phylum" = NA,"taxonomy_domain" = NA,"taxonomy_root" = NA,
                   "Number_of_genes" = NA,"Number_of_genes_with_hit" = NA,"Fraction_of_genes_with_hit" = NA,"Fraction_genes_at_majority_level" = NA,"Majority_taxonomy_level" = NA,"Majority_taxonomy" = NA)
temp <- temp[!temp$seqid %in% sequence_resolved_lineage_gene_fractions.df$seqid,]
rownames(temp) <- temp$seqid

# Determine if a sequence was missing because not enough genes with hits
# or no majority reached
for (sid in temp$seqid){
  if (sid %in% Sequence_gene_count.df[Sequence_gene_count.df$Fraction_of_genes_with_hit < MIN_FRACTION_GENES_HIT,]$seqid){
    temp[sid,]$Majority_taxonomy_level <- "Not enough genes with hits to reference"
    # temp[sid,]$Majority_taxonomy <- "Not enough genes with hits to reference"
  } else if (! sid %in% sequence_resolved_lineage_gene_fractions.df$seqid){
    temp[sid,]$Majority_taxonomy_level <- "No majority lineage"
    # temp[sid,]$Majority_taxonomy <- "No majority lineage"
  }
  temp[sid,]$Number_of_genes <- Sequence_gene_count.df[Sequence_gene_count.df$seqid == sid,]$Number_of_Genes
  temp[sid,]$Number_of_genes_with_hit <- Sequence_gene_count.df[Sequence_gene_count.df$seqid == sid,]$Number_of_genes_with_hit
  temp[sid,]$Fraction_of_genes_with_hit <- Sequence_gene_count.df[Sequence_gene_count.df$seqid == sid,]$Fraction_of_genes_with_hit
}

sequence_resolved_lineage_gene_fractions.df <- rbind(sequence_resolved_lineage_gene_fractions.df, temp)
rownames(sequence_resolved_lineage_gene_fractions.df) <- NULL

sequence_resolved_lineage_gene_fractions.df <-
  sequence_resolved_lineage_gene_fractions.df[,c("seqid", "Number_of_genes", "Number_of_genes_with_hit",
                                                 "Fraction_of_genes_with_hit","Fraction_genes_at_majority_level",
                                                 "Majority_taxonomy_level", "Majority_taxonomy")]
names(sequence_resolved_lineage_gene_fractions.df)[1] <- "Sequence"



# Note regarding output/categories
# “Not enough genes with hits to reference”:
# Not enough genes hit reference database to try and classify sequence, no taxonomy string will be assigned
#
# “No majority lineage”:
# No agreement on majority taxonomy (> 50%), no taxonomy string will be assigned
#
# When “Majority_taxonomy_level” = “taxonomy_species” and where “Majority_taxonomy” is all “unassigned”:
# Majority of genes hit something in reference database, but the majority taxonomy is completely unassigned

# write.table(sequence_resolved_lineage_gene_fractions.df,
            # file="sequence_resolved_lineages.tsv",
            # sep="\t", quote=F, row.names=F, col.names=T)
write.table(sequence_resolved_lineage_gene_fractions.df,
            file=paste0(args[6], "/sequence_resolved_lineages.tsv"),
            sep="\t", quote=F, row.names=F, col.names=T)


