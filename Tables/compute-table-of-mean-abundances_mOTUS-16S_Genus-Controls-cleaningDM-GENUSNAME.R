###
### Read table of counts for mOTU-v2 and 16S pipelines
### (MALT+MEGAN6/Qiime1/DADA2)
### Rq: gene-names have been cleaned
### + WGS: filtred mOTU == 0 for each control-samples
###
### Author: Denis Mestivier (denis.mestivier@u-pec.fr)
###

rm(list=ls())

###
### Filenames
###

fn16S = "./16S/CMP_CCR116S-Normaux-norm-mergingALL-MALT-QIIME-DADA2_counts-genus_cleaningDM-GENUSNAME_rev-23-03-2020.csv"
fnWGS = "./WGS/EC2M3-WGS-motus_rev-04-11-2019_bc-taxo_sumGenus-controls-cleaningDM-GENUSNAME.csv"
outfnstats = "./CMP_mean-abundance_Normaux_MOTU-MALT-QIIME-DADA2_counts-genus_cleaningDM-GENUSNAME_rev-01-06-2020.csv"

###
### read 16S
###

dta16S  = read.csv( fn16S, header=T, sep="\t" )
dtaWGS  = read.csv( fnWGS, header=T, sep="\t" )

# WGS: filter mOTU for which no control had been identified
#      last col (#51) => taxpath
sr = apply( dtaWGS[,1:50], 1, sum )
mOTUat0 = which( sr==0 )
dtaWGS = dtaWGS[ -mOTUat0, ]

###
### Merge WGS and 16S dataframes using taxpath
###

z = merge( dtaWGS, dta16S, by="genus", all=T )
rownames(z)=z$genus
z = z[, -1]

###
### Five dataset where not computed using DADA2
### Remove them
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
### Compute mean of abundance by taxid
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

counts.means = data.frame( Genus=rownames(z), MOTU=counts.mean.MOTU, M=counts.mean.MALT, Q=counts.mean.QIIME, D=counts.mean.DADA2 )

# ajoute une colonne du nombre de NA sur les 16S
# donc une evaluation sommaire de leur concordance

counts.means$concor16S = apply( counts.means[, c("M","Q","D") ], 1, function(x){ sum(!is.na(x)) } )
counts.means$concorALL = apply( counts.means[, c("M","Q","D","MOTU") ], 1, function(x){ sum(!is.na(x)) } )

###
### Write table of abundances
###

write.table( counts.means, outfnstats, sep="\t", row.names=F )
