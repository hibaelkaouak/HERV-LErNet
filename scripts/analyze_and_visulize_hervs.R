library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(visNetwork)

# Load and prepare TAD (Topologically Associated Domain) data
# Reading domain list file from Arrowhead algorithm and formatting chromosome names to UCSC style
domains_df <- read.table("../data/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt", header = TRUE, stringsAsFactors = FALSE)[, c("chr1", "x1", "x2")]
names(domains_df) <- c("chr", "start", "end")
domains_df$chr <- paste0("chr", domains_df$chr)
domains_gr <- GRanges(
  seqnames = domains_df$chr,
  ranges = IRanges(start = domains_df$start, end = domains_df$end)
)

# Load chromatin loop data from HiCCUPS algorithm
# These represent long-range chromatin interactions between distant genomic regions
loops_df <- read.table("../data/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt", header = TRUE, stringsAsFactors = FALSE) [, c("chr1",	"x1",	"x2",	"chr2",	"y1",	"y2")]

# Create GRanges objects for both anchors of the chromatin loops
# Each loop has two anchors (regions that interact with each other)
loop_anchor1 <- GRanges(
  seqnames = loops_df$chr1,
  ranges = IRanges(start = loops_df$x1, end = loops_df$x2)
)

loop_anchor2 <- GRanges(
  seqnames = loops_df$chr2,
  ranges = IRanges(start = loops_df$y1, end = loops_df$y2)
)

# Set chromosome naming style to UCSC format (e.g., "chr1" instead of "1")
seqlevelsStyle(loop_anchor1) <- "UCSC"
seqlevelsStyle(loop_anchor2) <- "UCSC"

# Import HERV (Human Endogenous Retrovirus) data from BED file
# HERVs are remnants of ancient retroviral infections in the human genome
hervs_gr <- import("../data/package-entities-erv.bed", format = "BED")
hervs_gr$name <- sub("^[^:]+:(.*)", "\\1", hervs_gr$name)

# Find overlaps between HERVs and TADs to identify which HERVs are located within TADs
hits <- findOverlaps(hervs_gr, domains_gr)

# Extract indices of overlapping features for later reference
query_idx <- queryHits(hits)      
subject_idx <- subjectHits(hits) 

# Create a dataframe containing information about HERVs and their overlapping TADs
overlap_df <- data.frame(
  herv_chr = as.character(seqnames(hervs_gr)[query_idx]),
  herv_start = start(hervs_gr)[query_idx],
  herv_end = end(hervs_gr)[query_idx],
  herv_name = hervs_gr$name[query_idx],
  
  domain_chr = as.character(seqnames(domains_gr)[subject_idx]),
  domain_start = start(domains_gr)[subject_idx],
  domain_end = end(domains_gr)[subject_idx]
)

# Import gene annotation data from GENCODE
genes_gr <- import("../data/gencode.v47.annotation.gtf", format = "gtf")

# Filter to only include actual genes (not other feature types like exons, transcripts, etc.)
genes_gr <- genes_gr[genes_gr$type == "gene"]

# Create GRanges object for TADs that contain HERVs
# This will be used to find genes that are in the same TADs as HERVs
tads_from_hervs_gr <- GRanges(
  seqnames = overlap_df$domain_chr,
  ranges = IRanges(start = overlap_df$domain_start, end = overlap_df$domain_end)
)

# Find overlaps between genes and TADs containing HERVs
gene_hits <- findOverlaps(genes_gr, tads_from_hervs_gr)
gene_idx <- queryHits(gene_hits)
tad_idx <- subjectHits(gene_hits)

# Create a dataframe linking HERVs to genes through TAD co-localization
# This represents potential regulatory relationships where HERVs and genes share the same TAD
herv_gene_df <- data.frame(
  #herv_chr = overlap_df$herv_chr[tad_idx],
  #herv_start = overlap_df$herv_start[tad_idx],
  #herv_end = overlap_df$herv_end[tad_idx],
  herv_name = overlap_df$herv_name[tad_idx],
  
  #tad_chr = overlap_df$domain_chr[tad_idx],
  #tad_start = overlap_df$domain_start[tad_idx],
  #tad_end = overlap_df$domain_end[tad_idx],
  
  #gene_chr = as.character(seqnames(genes_gr)[gene_idx]),
  #gene_start = start(genes_gr)[gene_idx],
  #gene_end = end(genes_gr)[gene_idx],
  #gene_id = genes_gr$gene_id[gene_idx],
  gene_name = genes_gr$gene_name[gene_idx],
  #gene_type = genes_gr$gene_type[gene_idx]
  
  source = "TAD"
)

# Find HERVs that overlap with chromatin loop anchors
herv_in_anchor1 <- findOverlaps(hervs_gr, loop_anchor1)
herv_in_anchor2 <- findOverlaps(hervs_gr, loop_anchor2)

# Extract indices for HERVs in anchor1 and their corresponding loop indices
herv_idx1 <- queryHits(herv_in_anchor1)
loop_idx1 <- subjectHits(herv_in_anchor1)

# Get the partner anchor2 regions for loops where anchor1 contains HERVs
anchor2_partner1 <- loop_anchor2[loop_idx1]

# Extract indices for HERVs in anchor2 and their corresponding loop indices
herv_idx2 <- queryHits(herv_in_anchor2)
loop_idx2 <- subjectHits(herv_in_anchor2)

# Get the partner anchor1 regions for loops where anchor2 contains HERVs
anchor1_partner2 <- loop_anchor1[loop_idx2]

# Find genes that overlap with the partner anchors
gene_hits_partner1 <- findOverlaps(genes_gr, anchor2_partner1)
gene_hits_partner2 <- findOverlaps(genes_gr, anchor1_partner2)

# Create dataframe for HERVs in anchor1 connected to genes in anchor2
df1 <- data.frame(
  herv_name = hervs_gr$name[herv_idx1[subjectHits(gene_hits_partner1)]],
  gene_name = genes_gr$gene_name[queryHits(gene_hits_partner1)],
  source = "loop"
  #loop_source = "anchor1"
)

# Create dataframe for HERVs in anchor2 connected to genes in anchor1
df2 <- data.frame(
  herv_name = hervs_gr$name[herv_idx2[subjectHits(gene_hits_partner2)]],
  gene_name = genes_gr$gene_name[queryHits(gene_hits_partner2)],
  source = "loop"
  #loop_source = "anchor2"
)

# Combine both dataframes to get all HERV-gene connections through loops
herv_loop_gene_df <- bind_rows(df1, df2)

# Remove duplicates if needed
herv_loop_gene_df <- distinct(herv_loop_gene_df)

# Print summary statistics
#print(paste("Found", nrow(herv_loop_gene_df), "HERV-gene connections through loops"))
#print(paste("Number of unique HERVs involved:", length(unique(herv_loop_gene_df$herv_name))))
#print(paste("Number of unique genes involved:", length(unique(herv_loop_gene_df$gene_name))))

# Combine the HERV-gene associations found through both TAD co-localization and chromatin loops
combined_df <- rbind(herv_gene_df, herv_loop_gene_df)

# Extract a random sample of 10 unique HERVs from the combined dataset and display them
#unique_hervs_sample <- sample(unique(combined_df$herv_name), 10)
#print(unique_hervs_sample)

# Example: Focus on a specific HERV for network visualization
herv_id <- "ERV:MSTA" # You can change this to analyze different HERVs

# Filter data for that specific HERV
# This extracts only the connections (TAD and loop) involving the HERV of interest
subset_df <- combined_df %>% filter(herv_name == herv_id)

# Remove any duplicate entries that might have been introduced during data processing
subset_df <- unique(subset_df)

# Limit the number of genes to display in the network visualization
max_genes <- 20
subset_df <- subset_df %>% sample_n(min(nrow(.), max_genes))

# Check if we have any genes associated with the selected HERV
if (nrow(subset_df) > 0) {
  # Create nodes: HERV as green dot, genes as light green triangles
  # Each node has an ID, label, shape, color, and font size
  nodes <- data.frame(
    id = c(herv_id, paste(subset_df$gene_name, subset_df$source, sep = "_")),
    label = c(herv_id, subset_df$gene_name),
    shape = c("dot", rep("triangle", nrow(subset_df))),
    color = c("green", rep("lightgreen", nrow(subset_df))),
    font = list(size = 20),
    stringsAsFactors = FALSE
  )
  
  # Create edges (connections): TAD connections in green, loop connections in blue
  # Each edge connects from the HERV to a gene, with color indicating connection type
  edges <- data.frame(
    from = rep(herv_id, nrow(subset_df)),
    to = paste(subset_df$gene_name, subset_df$source, sep = "_"),
    color = ifelse(subset_df$source == "TAD", "green", "blue"),
    title = subset_df$source,
    stringsAsFactors = FALSE
  )
  
  # Visualize the network using the visNetwork package
  visNetwork(nodes, edges, height = "600px", width = "100%") %>%
    visNodes(borderWidth = 2) %>%
    visEdges(smooth = TRUE, arrows = "to") %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visLayout(randomSeed = 42)
  
} else {
  print(paste("No genes associated with", herv_id))
}
