###
### Read a table of counts for mOTU-v2 and 16S pipelines
### (MALT+MEGAN6/Qiime1/DADA2)
### + WGS: filtered mOTU == 0 for each control-samples
###
### author: Denis MESTIVIER (denis.mestivier@u-pec.fr)
###

rm(list=ls())

library(pheatmap)

###
### Filenames
###

fnWGS = "../Tables/WGS/EC2M3-WGS-motus_rev-04-11-2019_bc-taxo_sumGenus-controls-cleaningDM-GENUSNAME.csv"
fn16S.counts = "../Tables/16S/CMP_CCR116S-Normaux-norm-mergingALL-MALT-QIIME-DADA2_counts-genus_cleaningDM-GENUSNAME_rev-23-03-2020.csv"
outfnstats = "./checking-CMP_mean-abundance_Normaux_MOTU-MALT-QIIME-DADA2_counts-genus_cleaningDM-GENUSNAME.csv"

###
### read 16S (genus => col1)
### read WGS (genus => col51)
###

dta16S  = read.csv( fn16S.counts, header=T, sep="\t", row.names=1 )
dtaWGS  = read.csv( fnWGS, header=T, sep="\t", row.names=51 )

# WGS: filter mOTU for which no control had been identified
sr = apply( dtaWGS, 1, sum )
mOTUat0 = which( sr==0 )
dtaWGS = dtaWGS[ -mOTUat0, ]

###
### Merge WGS and 16S dataframes using genus 
###

dtaWGS$genus = rownames(dtaWGS)
dta16S$genus = rownames(dta16S)

z = merge( dtaWGS, dta16S, by="genus", all=T )
rownames(z)=z$genus
z = z[, -1]

###
### Five dataset where not computed using DADA2
### Remove them...
###

remo = grep( "FR.116", colnames(z))
z = z[, -remo ]
remo = grep( "FR.169", colnames(z))
z = z[, -remo ]
remo = grep( "FR.276", colnames(z))
z = z[, -remo ]
remo = grep( "FR.400", colnames(z))
z = z[, -remo ]
remo = grep( "FR.716", colnames(z))
z = z[, -remo ]

###
### Compute mean of abundance by genus 
### for each pipeline
###

iMOTU  = grep( ".mOTU$", colnames(z) )
iMALT  = grep( ".M$",  colnames(z) )  
iDADA2 = grep( ".D$",  colnames(z) )  
iQIIME = grep( ".Q$",  colnames(z) )  

# be carefull: some otu has 'NA' due to merging(all=T)
counts.mean.MOTU = rowMeans( z[ , iMOTU  ], na.rm=T ) 
counts.mean.MALT = rowMeans( z[ , iMALT  ], na.rm=T ) 
counts.mean.QIIME= rowMeans( z[ , iQIIME ], na.rm=T ) 
counts.mean.DADA2= rowMeans( z[ , iDADA2 ], na.rm=T ) 
counts.mean.16S  = rowMeans( z[ , c(iMALT, iQIIME, iDADA2) ] )

###
### Dataframe of mean abundances
###

counts.means = data.frame( Genus=rownames(z), M=counts.mean.MALT, Q=counts.mean.QIIME, D=counts.mean.DADA2, MOTU=counts.mean.MOTU )

# Add a column of the number of NA for 16S
# which is a basic evaluation of their concordance

counts.means$concor16S = apply( counts.means[, c("M","D","Q") ], 1, function(x){ sum(!is.na(x)) } )
counts.means$concorALL = apply( counts.means[, c("M","D","Q","MOTU") ], 1, function(x){ sum(!is.na(x)) } )

###
### Some Genus are "unusable" such as "g__uncultured..." or "s__" or "g__unknwon..."
### because they are not shared by the different databases used or represent some
### little variations between genus names (for example between QIIME1 and DADA2).
### They are filtered
###

L  = grep( "s__", counts.means$Genus )
counts.means = counts.means[ -L, ]
z = z[ -L, ]

L = grep( "g__uncultured", counts.means$Genus )
counts.means = counts.means[ -L, ]
z = z[ -L, ]

L = grep( "g__unknown", counts.means$Genus )
counts.means = counts.means[ -L, ]
z = z[ -L, ]

L = grep( "^[0-9]", counts.means$Genus )
counts.means = counts.means[ -L, ]
z = z[ -L, ]

###
### Should keep genus identified by at least WGS and one 16S 
### This will remove genus 'unknown' or 'unclassified'
### or 'exotic' genus specific to one pipeline.
###

L = which( counts.means$concorALL>=2 )
z = z[ L, ]
counts.means = counts.means[ L, ]

###
### redo a dataframe for 16S
###

z.16S = z[, c(iMALT, iDADA2, iQIIME ) ]

###
### Write table of abundances (not useful because it is done elsewhere)
### just for checking
###

#write.table( z, outfnstats, sep="\t", row.names=T )

###
### Ok let's draw
###

annotROWS = counts.means[ , c("MOTU", "M", "D", "Q" ) ] 

annotCOLS = data.frame( Pipeline=c( rep("M", length(iMALT)), rep( "Q  ", length(iQIIME)), rep( "D", length(iDADA2) ) ) )


# Filtering : keep most abundants
# and order by names (three pipeline/name)

L = order( counts.means$MOTU, decreasing=T )
K = order( colnames( z.16S ) )

### Zoom on the twenty-most abundant genus according
### to WGS abundance

L1 = L[1:51]
K1 = K[1:30]

zn = z.16S[ L1,K1 ]
zr = annotROWS[ L1, ]

za = sapply( strsplit( colnames(zn), "\\." ), "[[", 3 )
za = data.frame( Pipeline=za )
rownames(za) = colnames(zn)

collab = colnames( zn )
collab[ seq(from=1, to=length(collab), by=3) ] = ""
collab[ seq(from=3, to=length(collab), by=3) ] = ""
collab = gsub( ".M", "", collab )

breaksBetweenSample = grep( ".Q$", colnames(zn) )

# I had to guest these values
RESX=4096*0.75
RESY=4096

png( file="figure3-resolution300dpi.png", width = RESX, height = RESY, res = 300 )

#png( file="heatmap-WGS-16S-log2-highAbundanceWGS.png", width = 1024, height = 1024 )
pheatmap( log2(zn+0.001),  annotation_col=za, annotation_row=zr, cluster_cols=F, cluster_rows=F, gaps_col=breaksBetweenSample, labels_col=collab, angle_col=c("315"), fontsize_cols=10, fontsize_rows=8 )

dev.off()

### Zoom for publication for low-abundance genera
#
#L1 = L[20:40]
#K1 = K[1:30]
#
#zn = z.16S[ L1,K1 ]
#zr = annotROWS[ L1, ]
#breaksBetweenSample = grep( ".Q$", colnames(zn) )
#
#collab = colnames( zn )
#collab[ seq(from=1, to=length(collab), by=3) ] = ""
#collab[ seq(from=3, to=length(collab), by=3) ] = ""
#collab = gsub( ".M", "", collab )
#
#png( file="heatmap-16S-log2-zoom-low-abundance-sub-class_zoom.png", width = 480, height = 480 )
#pheatmap( log2( zn+0.001),  annotation_row=zr, cluster_cols=F, cluster_rows=F, gaps_col=breaksBetweenSample, labels_col=collab, angle_col=c("315"), na_col="grey" )
#dev.off()

### Zoom for publication for low-abundance genera

#L1 = L[20:60]
#K1 = K[1:30]
#
#zn = z.16S[ L1,K1 ]
#zr = annotROWS[ L1, ]
#breaksBetweenSample = grep( ".Q$", colnames(zn) )
#
#collab = colnames( zn )
#collab[ seq(from=1, to=length(collab), by=3) ] = ""
#collab[ seq(from=3, to=length(collab), by=3) ] = ""
#collab = gsub( ".M", "", collab )
#
#png( file="heatmap-16S-log2-zoom-low-abundance-sub-class-40_zoom.png", width = 480, height = 480 )
#pheatmap( log2( zn+0.001),  annotation_row=zr, cluster_cols=F, cluster_rows=F, gaps_col=breaksBetweenSample, labels_col=collab, angle_col=c("315") )
#dev.off()

### Complete heatmap for supp. figures

L1 = L
K1 = K

zn = z.16S[ L1,K1 ]
breaksBetweenSample = grep( ".Q$", colnames(zn) )

za = sapply( strsplit( colnames(zn), "\\." ), "[[", 3 )
za = data.frame( Pipeline=za )
rownames(za) = colnames(zn)

collab = colnames( zn )
collab[ seq(from=1, to=length(collab), by=3) ] = ""
collab[ seq(from=3, to=length(collab), by=3) ] = ""
collab = gsub( ".M", "", collab )

RESX=8192*0.7
RESY=8192*1.15
png( file="SuppFigureheatmap-resolution300dpi.png", width = RESX, height = RESY, res=300 )

#png( file="heatmap-WGS-16S-log2-global.png", width = 2000, height = 3000 )

pheatmap( log2( zn+0.001),  annotation_col=za, annotation_row=zr, cluster_cols=F, cluster_rows=F, gaps_col=breaksBetweenSample, labels_col=collab, angle_col=c("315"), cellwidth=5, cellheight=8 )
dev.off()
