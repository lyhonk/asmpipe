#!/usr/bin/env Rscript

# Load libraries

suppressPackageStartupMessages(library("optparse", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("ggplot2", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("reshape2", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("data.table", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("tools", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))

# asm functions

#
# sliding windows
#

swsCalc<-function(x, win=list(L=75, R=75)){
  library(zoo)
  xLen <- length(x)
  winLen<-win$L + win$R + 1
  sws <- rollsum(x,winLen)
  sws_head<-tail(cumsum(x[1:(winLen-1)]),win$L)
  sws_tail<-rev(tail(cumsum(rev(x[(xLen - winLen + 2):xLen])),win$R))
  sws <- append(sws_head, sws)
  sws <- append(sws, sws_tail)
  return(sws)
}

#
# Generate tagmeth indexing for chr.tagmeth.csv
# Argument:
#   file.path - full path of chr.tagmeth.csv 
#   params - list(mscore, umscore, swsize)
# Return:
#   TagmethIndex Object

GenerateTagmethIndex <- function(genome, chr.tagmeth, chr.name, chr.pos, params = list(mscore = 0.8, umscore = 0.2, swsize = 50))
{
  perbase.tagcount <- data.frame(pos = chr.pos)
  perbase.tagcount$methylated <- rep(0, nrow(perbase.tagcount))
  perbase.tagcount$unmethylated <- rep(0, nrow(perbase.tagcount))
  
  chr.tagmeth.methylated <- chr.tagmeth[chr.tagmeth$score > params$mscore, ]
  chr.tagmeth.unmethylated <- chr.tagmeth[chr.tagmeth$score < params$umscore, ]
  
  stats <- as.data.frame(table(chr.tagmeth.methylated$mid))
  stats[, 1] <- as.integer(as.character(stats[, 1]))
  
  perbase.tagcount$methylated[stats[, 1]] <- stats[, 2]
  
  stats <- as.data.frame(table(chr.tagmeth.unmethylated$mid))
  stats[, 1] <- as.integer(as.character(stats[, 1]))
  
  perbase.tagcount$unmethylated[stats[, 1]] <- stats[, 2]
  
  chr.tagmeth.indexing.methylated <- swsCalc(perbase.tagcount$methylated, win=list(L = params$swsize, R = params$swsize))
  chr.tagmeth.indexing.unmethylated <- swsCalc(perbase.tagcount$unmethylated, win=list(L = params$swsize, R = params$swsize))
  return(data.frame(pos = chr.pos, midx = chr.tagmeth.indexing.methylated, umidx = chr.tagmeth.indexing.unmethylated))
}

# 
# Get ASMs of the tagmeth indexing
# RLE method is used
#

GetASMOfTagmethIndex <- function(tagmeth.indexing, range.size = 100, critera = list(meth, umeth))
{
  # generate tagmeth marker
  
  marker <- (tagmeth.indexing$midx > critera$meth) & (tagmeth.indexing$umidx > critera$umeth)
  
  # get the rle result
  
  idx.rle <- rle(marker)
  idx.ranges <- data.table(value = idx.rle$values, size = idx.rle$lengths, end = cumsum(idx.rle$lengths))
  idx.ranges$start <- idx.ranges$end - idx.ranges$size
  
  # filter out the required ranges
  
  tagmeth.asm <- idx.ranges[idx.ranges$value & idx.ranges$size > range.size,]
  
  return(tagmeth.asm)
}

# generate asm figures

PlotTagMethRange <- function(genome, chr.tagmeth, chr.name, range.start, range.end, plot.file.name) {
  
  # get cg positions
  
  refseq <- genome[[chr.name]][range.start:range.end]
  cg.pos <- range.start + start(matchPattern("CG", refseq)) - 2
  
  # get tags in the range
  
  tagmeth.range <- chr.tagmeth[chr.tagmeth$chr == chr.name & chr.tagmeth$mid > range.start & chr.tagmeth$mid < range.end]
  
  # convert tag cg positions
  
  tagmeth <- tagmeth.range
  tagmeth$mpos <- lapply(strsplit(tagmeth.range$mpos, ","), as.integer)
  tagmeth$upos <- lapply(strsplit(tagmeth.range$upos, ","), as.integer) 
  tagmeth$xpos <- lapply(strsplit(tagmeth.range$xpos, ","), as.integer)
  
  # sort by score
  
  tagmeth <- tagmeth[order(-tagmeth$score), ]
  tagmeth$tag <- NULL
  tagmeth$tag <- c(1:nrow(tagmeth))
  
  # prepare plot data
  
  tagmeth.plot.data <- tagmeth[, list(index = match(c(unlist(mpos), unlist(upos), unlist(xpos)) + pos, cg.pos), 
                                      type = c(rep("M", mcount), rep("U", ucount), rep("X", ucount))),
                               by = tag]
  tagmeth.plot.data <- tagmeth.plot.data[!is.na(tagmeth.plot.data$index)]
  
  # plot
  
  plot <- ggplot(tagmeth.plot.data, aes(index, tag, shape = type)) + geom_point() + scale_shape_manual(values=c(19, 1, 0))
  ggsave(plot.file.name, plot = plot)
}

##Specify desired options in a list

option_list <- list(
  make_option(c("-u","--unmethy-score"), help="methylation score threshhold for unmethlated tag", default = 0.2),
  make_option(c("-m","--methy-score"), help="methylation score threshhold for methylated tag", default = 0.8),
  make_option(c("-w","--sliding-window"), help="sliding window size", default = 50),
  make_option(c("-a","--asm-size"), help="asm range size", default = 500),
  make_option(c("-i","--methy-index"), help="methy index threshhold", default = 5),
  make_option(c("-x","--unmethy-index"), help="unmethy index threshold", default = 5),
  make_option(c("-l","--genome-library"), help="Bioconductor BSgenome library name", default = "BSgenome.Mmusculus.UCSC.mm9"),
  make_option(c("-n","--genome-name"), help="genome library object name. ex: \"Mmusculus\", \"Hsapiens\", \"Scerevisiae\"", default = "Mmusculus"),
  make_option(c("-t","--genome-type"), help="genome type , example mm9, mm10, hg19, hg18, default is NULL", default = "")
)

# Get command line options

arguments <- parse_args(OptionParser(usage = "%prog [options] tagmethPath", option_list = option_list), positional_arguments = 1)
opt <- arguments$options

kUnmethyScore <- opt$`unmethy-score`
kMethyScore <- opt$`methy-score`
kSlidingWindow <- opt$`sliding-window`
kAsmSize <- opt$`asm-size`
kMethyIndex <- opt$`methy-index`
kUnmethyIndex <- opt$`unmethy-index`
kGenomeLibrary <- opt$`genome-library`
kGenomeName <- opt$`genome-name`
kGenomeType <- opt$`genome-type`

kTagmethPath <- arguments$args

# load tagmeth files

if(!file.exists(kTagmethPath)){
  stop("tagmeth file path \"", kTagmethPath ,"\" does not exist.")
}

tagmeth.filenames <- list.files(kTagmethPath)
filename.base <- file_path_sans_ext(basename(kTagmethPath))
idx.meth.wig.filename <- paste(filename.base, ".idx.meth.wig", sep = "")
idx.umeth.wig.filename <- paste(filename.base, ".idx.umeth.wig", sep = "")
asm.bed.filename <- paste(filename.base, ".asm.bed", sep = "")
asm.csv.filename <- paste(filename.base, ".asm.csv", sep = "")

# load the genome library

kGenomeTypeList <- list(
  mm9  = list(genome.library="BSgenome.Mmusculus.UCSC.mm9",genome.name="Mmusculus"),
  mm10 = list(genome.library="BSgenome.Mmusculus.UCSC.mm10",genome.name="Mmusculus"),
  hg18 = list(genome.library="BSgenome.Hsapiens.UCSC.hg18",genome.name="Hsapiens"),
  hg19 = list(genome.library="BSgenome.Hsapiens.UCSC.hg19",genome.name="Hsapiens")
)
kGenome <- NULL

if ( kGenomeType %in% names(kGenomeTypeList) ){
  suppressPackageStartupMessages(library(kGenomeTypeList[[kGenomeType]][["genome.library"]], character.only = TRUE, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
  kGenome <- get(kGenomeTypeList[[kGenomeType]][["genome.name"]]) 
}else {
  suppressPackageStartupMessages(library(kGenomeLibrary, character.only = TRUE, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
  kGenome <- get(kGenomeName)
}

if ( is.null(kGenome)){
  stop( "Load Biocondutor Genome Library ERROR " )
}

# read and process tagmeth file by chrom

idx.meth.wig.lines <- ''
idx.umeth.wig.lines <- ''
asm.bed.lines <- ''
asm.csv.lines <- 'chr\tstart\tend'
for (file.name in tagmeth.filenames)
{
  file.path <- paste(kTagmethPath, "/", file.name, sep="")
  chr.name <- strsplit(file.name, "\\.")[[1]][1]
  
  if(file_test("-f", file.path))
  {
    message("    processing ", file.name, "\t", date())
      
    # Load tagmeth file
      
    chr.tagmeth <- fread(file.path, sep="\t",  
                         colClasses=c("factor","character","character", "character", 
                                      "integer","integer", "character", "character", "character",
                                      "integer", "integer", "integer"),
                         showProgress=FALSE, data.table=FALSE)
    # Assign colnames
      
    colnames(chr.tagmeth) <- c("tag", "type", "cvt", "chr", "pos", "size", "mpos", "upos", "xpos", "mcount", "ucount", "xcount")
    
    # Calculate the midpoint of tags
      
    chr.tagmeth$mid <- round(chr.tagmeth$pos + chr.tagmeth$size / 2)
      
    # Calculate methylation score
      
    chr.tagmeth$score <- chr.tagmeth$mcount / (chr.tagmeth$mcount + chr.tagmeth$ucount + chr.tagmeth$xcount)
    
    # build methylation virtual vectors
      
    chr.pos <- c(1: length(kGenome[[chr.name]]))
    
    # save tagmeth object
    
    saveRDS(chr.tagmeth, paste(chr.name, ".tagmeth.rds", sep = ""))
    
    # genrate tagmeth indexing
    
    chr.tagmeth.idx <- GenerateTagmethIndex(kGenome, chr.tagmeth, chr.name, chr.pos, params = list(mscore = kMethyScore, umscore = kUnmethyScore, swsize = kSlidingWindow))
    
    # generate asm files
    
    chr.tagmeth.asm <- GetASMOfTagmethIndex(chr.tagmeth.idx, kAsmSize, critera = list(meth = kMethyIndex, umeth = kUnmethyIndex))
    chr.tagmeth.asm$chr <- rep(chr.name, nrow(chr.tagmeth.asm))
    
    # write wigfiles & bedfiles
    
    chr.wig.title <- paste("fixedStep chrom=", chr.name, " start=", min(chr.tagmeth.idx$pos), " step=1 span=1", sep = "")
    write(chr.wig.title, idx.meth.wig.filename, append = TRUE)
    write(chr.wig.title, idx.umeth.wig.filename, append = TRUE)
    write(paste(chr.tagmeth.idx$midx, sep = '\n'), idx.meth.wig.filename, append = TRUE)
    write(paste(-chr.tagmeth.idx$umidx, sep = '\n'), idx.umeth.wig.filename, append = TRUE)
    
    write(paste(paste(chr.tagmeth.asm$chr, chr.tagmeth.asm$start, chr.tagmeth.asm$end, sep = '\t'), sep = '\n'), asm.csv.filename, append = TRUE)
    write(paste(paste(chr.tagmeth.asm$chr, chr.tagmeth.asm$start, chr.tagmeth.asm$end, rep(1, nrow(chr.tagmeth.asm)), sep = '\t'), sep = '\n'), asm.bed.filename, append = TRUE)
  }
}


