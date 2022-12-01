#!/usr/bin/env Rscript
library(optparse)
parser <- OptionParser()
option_list <- list(

  make_option(c("-i", "--input_table"), action = "store", type = "character",
              help = "input table", default = NULL),
  make_option(c("-w", "--window"), action = "store", type = "numeric",
              help = "window size for binning, default 100000 ", default = 100000)

)
description <- "Creates bigwig from gff listed in input table, density is calculated in provided window "

get_density <- function(x, tw=1000000){
  cvg <- coverage(x)
  bins <- tileGenome(seqlengths(x), tilewidth = tw)
  d <- binnedAverage(unlist(bins), cvg, "coverage")
  d
}



epilogue <- ""
parser <- OptionParser(option_list = option_list, epilogue = epilogue, description = description,
                       usage = "usage: %prog COMMAND [OPTIONS]")
opt <- parse_args(parser, args = commandArgs(TRUE))

if (is.null(opt$input_table)){
  print_help(parser)
  stop('specify input table!')
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(Biostrings)
})

# read input table:
t <- read.table(opt$input_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# referecne path is first line in input table:
ref_seq_path <-  paste(t[1,"dirname"],t[1,"filename"], sep = "/")
print(ref_seq_path)

# get list of gff files:
file_paths <-paste(t$dirname,t$filename, sep = "/")
is_gff <- t$format == "gff3"

gff_files <- file_paths[is_gff]


if (length(gff_files) == 0){
  stop("No gff files found in input table")
}else{
  print(paste("Found", length(gff_files), "gff files"))
  print(gff_files)

  # geting reference sequence:
  r <- readDNAStringSet(ref_seq_path)
  sl <- seqlengths(r)
  #rm(r)
  for (g in gff_files){
    print(paste("Processing", g))
    out_file <- paste(g, ".bw", sep = "")
    # check if output file already exists:
    if (file.exists(out_file)){
      print(paste("Output file", out_file, "already exists, skipping"))
    }else{
      # if this fails skip it and go to next file:
        tryCatch({
            gff <- import(g)
            # remove seqlevels that are not in reference sequence:
            gff <- gff[seqlevels(gff) %in% seqlevels(r)]
            L = sl[seqlevels(gff)]
            names(L) = seqlevels(gff)
            seqlengths(gff) <- L
            seqlengths(gff)[is.na(seqlengths(gff))] <- 100000
            d <- get_density(gff, tw=opt$window)
            d$score=d$coverage
            print(summary(d$score))
            export(d, format = "BigWig", con = out_file)
        }, error = function(e) {
            print(paste("Error processing", g, "skipping"))
            print(e)
            print("--------------------------------")
        })
    }
  }
}


