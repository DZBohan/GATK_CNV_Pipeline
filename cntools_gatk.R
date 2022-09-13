library(CNTools)

# set the arguments
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
gene.info.f <- args[2]
gatk_file <- args[3]
id <- args[4]

# read the geninfo file
geneinfo.hg38 = read.table(gene.info.f,header=T,sep="\t")

# read the results segment files of gatk
gatk_out <- read.table(gatk_file,header=T)

# standardiza the table of gatk output
names(gatk_out)[1] <- "ID"
names(gatk_out)[2] <- "chrom"
names(gatk_out)[3] <- "loc.start"
names(gatk_out)[4] <- "loc.end"
names(gatk_out)[5] <- "num.mark"
names(gatk_out)[6] <- "seg.mean"
gatk_out$ID <- c(id)

if(gatk_out[1,2]=="chr1"){gatk_out$chrom <- substring(gatk_out$chrom, 4)}
gatk_out <- na.omit(gatk_out)

# run cntools on the segments results of gatk
cnseq_gatk <- CNSeg(gatk_out)
rd_gatk <- getRS(cnseq_gatk, by="gene", imput=FALSE, XY=FALSE, geneMap=geneinfo.hg38, what = "mean")
rs_gatk <- rs(rd_gatk)

# remove 0 value in the table 
rs_gatk[rs_gatk==0] <- NA
rs_gatk <- na.omit(rs_gatk)

# output the txt file
write.table (rs_gatk, file = paste(id,"cntools","gatk","txt",sep = "."), quote =FALSE, sep ="\t", row.names=FALSE)
