library( data.table )

setwd("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/single_cell_integration_dna_rna/expr_chromothripsis")

# sc-rna bams
path_to_rna_bams <- "/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/GT_tumor/scRNA/aligment_star_2pass/bams/"
bamFiles <- list.files(path = path_to_rna_bams,pattern = ".Aligned.sortedByCoord.out.bam$",recursive = T,full.names = T)
baiFiles <- list.files(path = path_to_rna_bams,pattern = ".Aligned.sortedByCoord.out.bam.bai$",recursive = T,full.names = T)

# featureCounts_res 
mat <- read.delim("/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/GT_tumor/scRNA/aligment_star_2pass/featureCounts_res/mat.txt",check.names = F,stringsAsFactors = F)
rownames(mat) <- mat[,1]
mat <- mat[,-1]
colnames(mat) <- gsub(basename(colnames(mat)),pattern = ".Aligned.sortedByCoord.out.bam",replacement = "")

mat <- mat[ , colSums(mat) > 1000 ] # keep only cells with at least N detected genes
mat <- mat[ rowSums(mat) > 0,] # remove the genes that appear in none of the remaining cells

# Normalisation
nrm <- t( t(mat) / colSums(mat) )

# # keep only genes having UMI count not null in at leat x fraction of cells
# genes_umi <- names(which( (rowSums( nrm > 0 )/ncol(nrm)) >= genes_umi_threshold))
# nrm <- nrm[ genes_umi ,]
# counts <- counts[ genes_umi ,]

## get genes coordinates
gene_ordering <- read.delim(file = "/icgc/dkfzlsdf/analysis/B260/users/n790i/tools/binning_the_genome/humangenes_biomart_GRCh37p13.sort.bed",stringsAsFactors = F,skip = 1,header = F)
gene_ordering <- gene_ordering_file[which(gene_ordering$V1 != "MT"),]
colnames(gene_ordering) <- c("chr","start","end","GeneName","GeneID")





## plot nrm_counts per genomic segment
chrom <- "7"

keep_genes <- intersect(rownames(nrm),unique(gene_ordering$GeneID[grep(gene_ordering$chr,pattern = chrom)]))

expr <- nrm[keep_genes,]




