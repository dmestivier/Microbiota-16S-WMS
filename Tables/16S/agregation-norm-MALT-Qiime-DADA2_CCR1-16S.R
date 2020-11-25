### 
### Read three OTU table (+TAXA / CSV format)
### Do summarization at genus level
### Normalize to percent
### Merge dataframe.
###

rm( list=ls() )

###
### filename
### QIIME/MALT => Control/Adenoma/CCR
### DADA2      => Control
###

fnQIIME = "otu_taxa_EC2M3_16S_16S_clean_silva.csv"
fnMALT  = "Comparison-EC2M3-16S_MALT-SILVA128-fromSAM_rev-19-01-2018_Taxonomy-mod.csv"
fnDADA2 = "otu_taxa_EC2M3_16S_DADA2-Filt150_Silva128_counts-and-annotation.csv"
fnNormaux = "EC2M3_CCR1-16S_Normal.csv"

QIIME  = read.csv( fnQIIME, header=T, sep="\t" )
MALT   = read.csv( fnMALT, header=T, sep="\t" )
DADA2  = read.csv( fnDADA2, header=T, sep=" " )

normaux = read.csv( fnNormaux, header=T, sep="\t" )
normaux$sampleID = gsub( "-", ".", normaux$sampleID )

###
### Build taxpath (ie concatenation of rank name)
###

QIIME$taxa = paste( QIIME$Phylum, QIIME$Class, QIIME$Order, QIIME$Family, QIIME$Genus, sep="." )
MALT$taxa = paste( MALT$Phylum, MALT$Class, MALT$Order, MALT$Family, MALT$Genus, sep="." )
DADA2$taxa = paste( DADA2$Phylum, DADA2$Class, DADA2$Order, DADA2$Family, DADA2$Genus, sep="." )

###
### Deduplicate 
###

### MALT

uMALT = unique( MALT$taxa )
iMALT = 2:130
nMALT = as.data.frame( matrix( 0, nrow=length(uMALT), ncol=length(iMALT) ) )
colnames(nMALT) = colnames(MALT)[iMALT]
rownames(nMALT) = uMALT 

for( ti in uMALT )
{
	nMALT[ti, ] = colSums( MALT[ which( MALT$taxa == ti ), iMALT ] ) 
}

### QIIME 

uQIIME = unique( QIIME$taxa )
iQIIME = 2:130
nQIIME = as.data.frame( matrix( 0, nrow=length(uQIIME), ncol=length(iQIIME) ) )
colnames(nQIIME) = colnames(QIIME)[iQIIME]
rownames(nQIIME) = uQIIME 

for( ti in uQIIME )
{
	nQIIME[ti, ] = colSums( QIIME[ which( QIIME$taxa == ti ), iQIIME ] ) 
}

### DADA2

uDADA2 = unique( DADA2$taxa )
iDADA2 = 2:46
nDADA2 = as.data.frame( matrix( 0, nrow=length(uDADA2), ncol=length(iDADA2) ) )
colnames(nDADA2) = colnames(DADA2)[iDADA2]
rownames(nDADA2) = uDADA2 

for( ti in uDADA2 )
{
	nDADA2[ti, ] = colSums( DADA2[ which( DADA2$taxa == ti ), iDADA2 ] ) 
}

###
### Keep only Control sample
###

colsNormauxQIIME = match( normaux$sampleID, colnames(nQIIME) )
colsNormauxMALT  = match( normaux$sampleID, colnames(nMALT) )
colsNormauxDADA2 = match( normaux$sampleID, colnames(nDADA2) )
# DADA2 miss 5 Control / Will be fixe
colsNormauxDADA2 = na.omit(colsNormauxDADA2)

nQIIME = nQIIME[ , colsNormauxQIIME ]
nMALT  = nMALT [ , colsNormauxMALT  ]
nDADA2 = nDADA2[ , colsNormauxDADA2 ]

###
### Normalisation
###

nQIIME = apply( nQIIME, 2, function(x){ 100.0*x / sum(x) } )
nMALT  = apply( nMALT , 2, function(x){ 100.0*x / sum(x) } )
nDADA2 = apply( nDADA2, 2, function(x){ 100.0*x / sum(x) } )

nQIIME = as.data.frame( nQIIME )
nMALT  = as.data.frame( nMALT  )
nDADA2 = as.data.frame( nDADA2 )

###
### Suffix by pipeline name
### Add "taxa" column for merging
###

colnames(nMALT)  = paste0( colnames(nMALT),  ".M", sep="" )
colnames(nQIIME) = paste0( colnames(nQIIME), ".Q", sep="" )
colnames(nDADA2) = paste0( colnames(nDADA2), ".D", sep="" )

nMALT$taxa = rownames(nMALT)
nQIIME$taxa = rownames(nQIIME)
nDADA2$taxa = rownames(nDADA2)


###
### Merge  with only consensus genus
### - pour les counts
### - pour les annotations
###

z = merge( nMALT, nQIIME, by="taxa", all=F )
z = merge( z, nDADA2, by="taxa", all=F )

###
### Add the taxid/NCBI thank's to MALT pipeline
### and that we keep only commons taxa between the three pipeline
###

ik = match( z$taxa, MALT$taxa )
ta = MALT[ ik, c("taxa", "OTUID" ) ]
z  = merge( ta, z, by="taxa" )

# !!! can do this only because we merge *only* commons taxa
# !!! between the three dataframe
keepTaxa = match( z$taxa, MALT$taxa )
keepCols = c( "OTUID", "taxa", "Domain","Phylum","Class","Order","Family","Genus")
taxo = MALT[ keepTaxa, keepCols ]
colnames(taxo)[1] = "taxid"
write.table( taxo, "CMP_CCR116S-Normaux-norm-MALT-QIIME-DADA2_taxa_rev-22-11-2019.csv", row.names=F, sep="\t" )

###
### Export
### Table of taxa should be import in libreoffice
###

z = z[, -1] # taxa not necessary now
colnames(z)[1] = "taxid"
write.table( z, "CMP_CCR116S-Normaux-norm-MALT-QIIME-DADA2_counts_rev-22-11-2019.csv", row.names=F, sep="\t")
