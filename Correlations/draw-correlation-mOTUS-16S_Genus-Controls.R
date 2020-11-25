###
### Read a table of counts for mOTU-v2 and 16S pipelines
### (MALT+MEGAN6/Qiime1/DADA2)
### + WGS: filtered mOTU == 0 for each control-samples
### + filtered `unknown`/`unclassified`
###
### Draw correlation between WGS and 16S abudance by 16Spipeline
###

rm(list=ls())
library(VennDiagram)
library(RColorBrewer)

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

###
### Dataframe of mean abundances
###

counts.means = data.frame( Genus=rownames(z), M=counts.mean.MALT, Q=counts.mean.QIIME, D=counts.mean.DADA2, MOTU=counts.mean.MOTU )

# ajoute une colonne du nombre de NA sur les 16S
# donc une evaluation sommaire de leur concordance

counts.means$concor16S = apply( counts.means[, c("M","Q","D") ], 1, function(x){ sum(!is.na(x)) } )
counts.means$concorALL = apply( counts.means[, c("M","Q","D","MOTU") ], 1, function(x){ sum(!is.na(x)) } )

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
### Write table of abundances (no useful because it is done elsewhere)
### just for checking
###

write.table( counts.means, outfnstats, sep="\t", row.names=T )

###
### Filtering : keep most abundants
### and order by names (three pipeline/name)
###

L = order( counts.means$MOTU, decreasing=T )

counts.means = counts.means[ L, ]
z = z[L, ]

###
### Select only genera representing at 98% of abundance 
###

### Je fais la courbe d'abondance cumulee
### que je dois normaliser a 100%

w  = counts.means$MOTU
L  = which( !is.na( w ) )
w  = w[L]
w  = w * 100 / sum( w )
cs = cumsum( w )
L  = which( cs < 98.0 )

counts.means.filt = counts.means[ L, ]

write.table( counts.means.filt, "filt.csv", row.names=F )

#length( is.na( counts.means.filt$M ) )
#table( counts.means.filt$concor16S )

# On dessine les correlations entre pipelines

xl = c(0.1,15)
yl = c(0,15)

w = counts.means.filt

w$MOTU = 100.0 * w$MOTU # pour revenir aux memes unites

par(pch=16)
par(cex.lab=1.2)
plot( w$MOTU, w$M, col="black", xlim=xl, ylim=yl, xlab="MOTU" )
par(new=T)
plot( w$MOTU, w$Q, col="blue", xlim=xl, ylim=yl, xlab="MOTU" )
par(new=T)
plot( w$MOTU, w$D, col="red", xlim=xl, ylim=yl, xlab="MOTU" )

