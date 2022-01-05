#!/usr/bin/env Rscript:
library(rtracklayer)
f = "/mnt/raid/454_data/cuscuta/RNA_seq/RNA_seq_mapping/Ceuropea_asm_v4_mapping/tmp/test2.gtf"
fout = "/mnt/raid/454_data/cuscuta/RNA_seq/RNA_seq_mapping/Ceuropea_asm_v4_mapping/tmp/test2.gff3"
f = commandArgs(T)[1]
fout = commandArgs(T)[2]
gtf = import(f, format = "gtf")

exons = gtf[gtf$type == "exon"]
exons$Parent = exons$transcript_id
exons$ID=paste0(exons$transcript_id, ".", exons$exon_number)
mRNA = gene = gtf[gtf$type == "transcript"]
mRNA$ID=mRNA$transcript_id
mRNA$Parent=mRNA$gene_id
mRNA$type = "mRNA"

gene$ID = mRNA$gene_id
gene$type = "gene"

gff = sort(c(gene, mRNA, exons))

gff$FPKM=NULL
gff$TPM=NULL
gff$cov=NULL
gff$gene_id=NULL
gff$transcript_id=NULL
gff$Name=gff$ID


export(gff, con = fout, format="gff3")
