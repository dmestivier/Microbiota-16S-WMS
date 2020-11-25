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

#write.table( counts.means, outfnstats, sep="\t", row.names=T )

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


### Maintenant on peut faire le diagramme de Venn

m = counts.means.filt

# Lignes ou des Genus ont ete detectes
L.motu = which( !is.na( m$MOTU ) )
L.dada2= which( !is.na( m$D) )
L.qiime= which( !is.na( m$Q) )
L.malt = which( !is.na( m$M) )

# Liste des genus
genus.motu  = m$Genus[ L.motu  ]
genus.dada2 = m$Genus[ L.dada2 ]
genus.qiime = m$Genus[ L.qiime ]
genus.malt  = m$Genus[ L.malt  ]

# Chart Produit un fichier
# si filename=NULL:
#    zz = venn.diagram(...)
#    grid.draw(zz)
# sinon il produit un fichier
myCol <- brewer.pal(4, "Pastel2")
#myCol <- brewer.pal(4, "Pastel2")

vd = venn.diagram( x = list( genus.motu, genus.dada2, genus.qiime, genus.malt),
	category.names = c( "mOTU-v2", "DADA2", "QIIME1", "MALT" ),
	#filename = NULL,
	filename = 'Figure-2B-VENN-resolution300dpi.png',
	output=FALSE,
	# Output features
	imagetype="png" ,
	height = 1024*1.2, 
	width = 1024*1.2, 
	resolution = 300,
	compression = "lzw",
	# Circles
	lwd = 1,
	#lty = 'blank',
	fill = myCol,
	# Numbers
	cex = 0.6,
	fontface = "italic",
	fontfamily = "sans",
	# Set names
	cat.cex = 0.8,
	cat.fontface = "bold",
	cat.default.pos = "outer",
	cat.fontfamily = "sans"
	#cat.pos = c(-27, 27, 135)
	#cat.dist = c(0.055, 0.055, 0.085), 
	#rotation = 1
)
