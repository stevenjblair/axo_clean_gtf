# Environment and Libraries -----------------------------------------------

setwd("/Users/sjblair/Projects/amex_CleanGTF")

library(pacman)
p_load(dplyr,rtracklayer,stringr)

# Import GTF --------------------------------------------------------------

gtf <- import.gff("./AmexT_v47-AmexG_v6.0-DD.gtf.gz", format ='gtf')

# Mutate GTF --------------------------------------------------------------

# First remove the unwanted columns
mcols(gtf) <- mcols(gtf)[, !colnames(mcols(gtf)) %in% c("gene_name", "homolog", "ORF_type", "CDS", "peptide")]

# Now, remove the annotation from the transcript_id to only be the transcript identifier
gtf$transcript_id <- sub(".*(AMEX60DD.*)$", "\\1", mcols(gtf)$transcript_id)

# Export GTF --------------------------------------------------------------

export.gff(gtf, "./AmexT_v47-AmexG_v6.0-DD_noNames.gtf") # Write the modified GTF file
# tr2g --------------------------------------------------------------------

# create a tr2g
df <- as.data.frame(gtf)
tr2g <- df %>% 
  filter(!is.na(gene_id) & !is.na(transcript_id)) %>% # removes the rows that would have an NA in one of the columns to keep
  distinct(gene_id, transcript_id, .keep_all = TRUE) %>% # remove redundant rows
  select(transcript_id, gene_id)

# Export tr2g -------------------------------------------------------------

write.table(tr2g, "AmexT_v47-AmexG_v6.0-DD_tr2g.tsv", row.names = F, quote = F,col.names = F)
