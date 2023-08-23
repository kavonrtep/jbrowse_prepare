#!/usr/bin/env Rscript

# inputs BED1, BED2 , output is PAF, chromosome size1 and size2

library(optparse)
options_list <- list(
    make_option(c("-a", "--bed1"), action = "store", type = "character",
                help = "bed1", default = NA),
    make_option(c("-b", "--bed2"), action = "store", type = "character",
                help = "bed2", default = NA),
    make_option(c("-o", "--output"), action = "store", type = "character",
                help = "output", default = NA),
    make_option(c("-A", "--size1"), action = "store", type = "character", help=".fai for A genome"),
    make_option(c("-B", "--size2"), action = "store", type = "character", help=".fai for B genome")
    )

parser <- OptionParser(option_list = options_list)
opt <- parse_args(parser, args = commandArgs(TRUE))

a <- read.table(opt$bed1, header = FALSE, stringsAsFactors = FALSE)[,1:6] # only first 4 columns
colnames(a) <- c("qseqid", "qstart", "qend", "Name", "score", "strandA")
b <- read.table(opt$bed2, header = FALSE, stringsAsFactors = FALSE)[,1:6] # only first 4 columns
colnames(b) <- c("sseqid", "sstart", "send", "Name", "score", "strandA")

sa <- read.table(opt$size1, header = FALSE, stringsAsFactors = FALSE)[,1:2]
colnames(sa) <- c("qseqid", "qlen")
sb <- read.table(opt$size2, header = FALSE, stringsAsFactors = FALSE)[,1:2]
colnames(sb) <- c("sseqid", "slen")

a <- merge(a, sa, by = 1, all.x = TRUE)
b <- merge(b, sb, by = 1, all.x = TRUE)

# merge
ab <- merge(a, b, by = 4, all = FALSE)
ab$m <- ab$length <- ab$qend - ab$qstart
ab$mapq <- 255
ab$strand <- ifelse(ab$strandA.x == ab$strandA.y, "+", "-")

# convert  to PAF
paf_cols <- c("qseqid", "qlen", "qstart", "qend", "strand", "sseqid", "slen", "sstart", "send", "length", "m", "mapq")


paf <- ab[,paf_cols]

# export paf:
write.table(paf, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



