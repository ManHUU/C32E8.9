# Ensembl Release 105
# https://ftp.ensembl.org/pub/release-105/fasta/caenorhabditis_elegans/dna/
# rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-105/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz .
System(paste0("rsync -av ",
              "rsync://ftp.ensembl.org/ensembl/pub/release-105/",
              "fasta/caenorhabditis_elegans/dna/",
              "Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz ",
              "./resources/"), Log_name = "./logs/Create_Reference_index.log", Append = FALSE)

# Retrieve GTF file
# rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-105/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.105.gtf.gz
System(paste0("rsync -av ",
              "rsync://ftp.ensembl.org/ensembl/pub/release-105/gtf/caenorhabditis_elegans/",
              "Caenorhabditis_elegans.WBcel235.105.gtf.gz ", 
              "./resources/"), Log_name = "./logs/Create_Reference_index.log", Append = TRUE)

# NCBI Primary Assembly
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/

# Parameters
n_clust <- 6
fasta_ref <- "./resources/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
gtf_ref <- "./resources/Caenorhabditis_elegans.WBcel235.105.gtf"
star_genome <- "./data/STAR_genome"
STAR_dir <- "./tools/STAR-2.7.10b/bin/MacOSX_x86_64"
# ---------------------------------
# Build Reference Genome using STAR
# ---------------------------------
command <- paste0(STAR_dir, "/STAR --runThreadN ", n_clust, 
                  " --runMode genomeGenerate --genomeDir ", star_genome,
                  " --genomeFastaFiles ", fasta_ref,
                  " --genomeSAindexNbases 12",
                  " --sjdbGTFfile ", gtf_ref, 
                  " --sjdbOverhang 100 --limitGenomeGenerateRAM 35924399488")
System(command, Log_name = "./logs/Create_Reference_index.log", Append = TRUE)

# -----------------------
# Convert GTF to R-Object
# -----------------------
elegans_gtf <- as.data.frame(rtracklayer::import(gtf_ref))
# save R Object
saveRDS(elegans_gtf, file = "./data/R_data/elegans_gtf.RDS")
# select only type of gene
gtf <- readRDS("./data/R_data/elegans_gtf.RDS")
gtf_gene <- gtf[gtf$type == "gene", ]
rownames(gtf_gene) <- gtf_gene$gene_id
saveRDS(gtf_gene, file = "./data/R_data/c_elegans_gene_gtf.RDS")

