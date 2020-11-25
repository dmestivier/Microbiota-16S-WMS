### 
### Read table OTU (+taxa / CSV format)
### and do genus summarization
###
### Summarization is done using genus-name and not
### the full taxpath
###
### Normalize data
###
### then we merge data.Frame
### and keep missing OTU
### => merge( ..., all=T )
###
### Rq: we are working on OTU table with
### genus-names cleaned by D. Mestivier
###
### Author : Denis Mestivier (denis.mestivier@u-pec.fr)
###

rm( list=ls() )

###
### filename.
### File = abundance + taxonomy
### !!! QIIME/MALT have the whole "French cohort" (Control/Adenoma/CCR)
### but DADA2 just have the Control sample
###

fnQIIME = "otu_taxa_EC2M3_16S_16S_clean_silva-cleaningDM.csv"
fnMALT  = "Comparison-EC2M3-16S_MALT-SILVA128-fromSAM_rev-19-01-2018_Taxonomy-mod.csv"
fnDADA2 = "otu_taxa_EC2M3_16S_DADA2-Filt150_Silva128_counts-and-annotation-cleaningDM.csv"
fnNormaux = "EC2M3_CCR1-16S_Normal.csv"

fnOUT = "CMP_CCR116S-Normaux-norm-mergingALL-MALT-QIIME-DADA2_counts-genus_cleaningDM-GENUSNAME_rev-23-03-2020.csv"

QIIME  = read.csv( fnQIIME, header=T, sep="\t", stringsAsFactors = F )
MALT   = read.csv( fnMALT, header=T, sep="\t", stringsAsFactors = F )
DADA2  = read.csv( fnDADA2, header=T, sep="\t", stringsAsFactors = F )

normaux = read.csv( fnNormaux, header=T, sep="\t" )
normaux$sampleID = gsub( "-", ".", normaux$sampleID )

###
### Build a taxa (= concatenation of rank names)
### NOT USED HERE
###

QIIME$taxa = paste( QIIME$Phylum, QIIME$Class, QIIME$Order, QIIME$Family, QIIME$Genus, sep="." )
MALT$taxa = paste( MALT$Phylum, MALT$Class, MALT$Order, MALT$Family, MALT$Genus, sep="." )
DADA2$taxa = paste( DADA2$Phylum, DADA2$Class, DADA2$Order, DADA2$Family, DADA2$Genus, sep="." )

###
### Some OTU have several row for the same genus.
### We had to summarized
###

### MALT

# trans-code NA to NONE (g__NA)
ii = which( is.na(MALT$Genus) )
MALT$Genus[ii] = "g__NA"

uMALT = unique( MALT$Genus )
iMALT = 2:130
nMALT = as.data.frame( matrix( 0, nrow=length(uMALT), ncol=length(iMALT) ) )
colnames(nMALT) = colnames(MALT)[iMALT]
rownames(nMALT) = uMALT 

for( ti in uMALT )
{
	nMALT[ti, ] = colSums( MALT[ which( MALT$Genus == ti ), iMALT ] ) 
}

### QIIME 
# trans-code NA to NONE (g__NA)
ii = which( is.na(QIIME$Genus) )
QIIME$Genus[ii] = "g__NA"

uQIIME = unique( QIIME$Genus )
iQIIME = 2:130
nQIIME = as.data.frame( matrix( 0, nrow=length(uQIIME), ncol=length(iQIIME) ) )
colnames(nQIIME) = colnames(QIIME)[iQIIME]
rownames(nQIIME) = uQIIME 

for( ti in uQIIME )
{
	nQIIME[ti, ] = colSums( QIIME[ which( QIIME$Genus == ti ), iQIIME ] ) 
}

### DADA2
# trans-code NA to NONE (g__NA)
ii = which( is.na(DADA2$Genus) )
DADA2$Genus[ii] = "g__NA"

uDADA2 = unique( DADA2$Genus )
iDADA2 = 2:46
nDADA2 = as.data.frame( matrix( 0, nrow=length(uDADA2), ncol=length(iDADA2) ) )
colnames(nDADA2) = colnames(DADA2)[iDADA2]
rownames(nDADA2) = uDADA2 

for( ti in uDADA2 )
{
	nDADA2[ti, ] = colSums( DADA2[ which( DADA2$Genus == ti ), iDADA2 ] ) 
}

###
### Keep only Control samples
###

colsNormauxQIIME = match( normaux$sampleID, colnames(nQIIME) )
colsNormauxMALT  = match( normaux$sampleID, colnames(nMALT) )
colsNormauxDADA2 = match( normaux$sampleID, colnames(nDADA2) )
# Be aware DADA2 miss 5 Control / will be update
colsNormauxDADA2 = na.omit(colsNormauxDADA2)

nQIIME = nQIIME[ , colsNormauxQIIME ]
nMALT  = nMALT [ , colsNormauxMALT  ]
nDADA2 = nDADA2[ , colsNormauxDADA2 ]

###
### Do normalisation => percentage
###

nQIIME = apply( nQIIME, 2, function(x){ 100.0*x / sum(x) } )
nMALT  = apply( nMALT , 2, function(x){ 100.0*x / sum(x) } )
nDADA2 = apply( nDADA2, 2, function(x){ 100.0*x / sum(x) } )

nQIIME = as.data.frame( nQIIME )
nMALT  = as.data.frame( nMALT  )
nDADA2 = as.data.frame( nDADA2 )

###
### Add suffix for 16S pipeline
### Add "taxa" column for merging
###

colnames(nMALT)  = paste0( colnames(nMALT),  ".M", sep="" )
colnames(nQIIME) = paste0( colnames(nQIIME), ".Q", sep="" )
colnames(nDADA2) = paste0( colnames(nDADA2), ".D", sep="" )

nMALT$genus = rownames(nMALT)
nQIIME$genus = rownames(nQIIME)
nDADA2$genus = rownames(nDADA2)

###
### Merge. Keep every genus
###

z = merge( nMALT, nQIIME, by="genus", all=T )
z = merge( z, nDADA2, by="genus", all=T )

###
### Export
###

write.table( z, fnOUT, row.names=F, sep="\t")
