###
### Read a table of counts of the three pipelines
### and an ICC file.
###

rm(list=ls())

### filename

fn.counts = "../Tables/16S/CMP_CCR116S-Normaux-norm-MALT-QIIME-DADA2_counts_rev-22-11-2019.csv"
fn.taxa   = "../Tables/16S/CMP_CCR116S-Normaux-norm-MALT-QIIME-DADA2_taxa_rev-22-11-2019.csv"
fn.icc    = "./CMP_CCR116S-Normaux-norm-MALT-QIIME-DADA2_ICC-export_rev-31-01-2020.csv"

###
### read files
###

dta  = read.csv( fn.counts, header=T, sep="\t", row.names=1 )
taxa = read.csv( fn.taxa  , header=T, sep="\t", row.names=1 )
icc  = read.csv( fn.icc   , header=T, sep="\t", stringsAsFactors=F ) 

counts = dta

###
### Five dataset where not computed using DADA2
### REmove them
###

remo = grep( "FR.116", colnames(counts))
counts = counts[, -remo ]
remo = grep( "FR.169", colnames(counts))
counts = counts[, -remo ]
remo = grep( "FR.276", colnames(counts))
counts = counts[, -remo ]
remo = grep( "FR.400", colnames(counts))
counts = counts[, -remo ]
remo = grep( "FR.716", colnames(counts))
counts = counts[, -remo ]

###
### No normalisation is mandatory (files are already norm)
###

counts.norm = counts

###
### Compute mean of abundance(normalisez) by taxid
###

counts.mean       = rowMeans( counts.norm )

###
### Dataframe for annotation of the heatmap
###

annot = data.frame( A=counts.mean, Phylum=taxa$Phylum, Genus=taxa$Genus )

###
### Add ICC computation
###

annot = merge( annot, icc, by="Genus" )

###
### Plot ICC versun mean abundance
###

RESX=4096
RESY=4096

RESX=2048
RESY=2048
png( "SuppFigure-resolution300dpi.png", width=RESX, height=RESY, res=300 )
#png( "icc-16S-versus-mean-abundance.png" )

par( cex.lab=1.6 )
par( cex.axis=1.6 )
plot( annot$A, annot$ICC3_single_fixed_raters, xlab="Mean normalized abundance", ylab="ICC3", pch=19, lwd=5, col="blue" ) 
abline( h=0.50, col="red", lty=2, lwd=3 )
abline( h=0.70, col="orange", lty=2, lwd=3 )
abline( h=0.90, col="green", lty=2, lwd=3 )
abline( v=2, col="grey", lty=3, lwd=2 )
abline( v=0.5, col="grey", lty=3, lwd=2 )
dev.off()
