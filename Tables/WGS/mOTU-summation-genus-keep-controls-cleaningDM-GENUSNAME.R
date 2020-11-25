###
### Read a table (normalised counts) from mOTU
### and do a summation up to Genus level.
### and keep only controls samples from CCR1
###

rm(list=ls())

###
### Files to use
###

fnMOTU = "./EC2M3-WGS-motus_rev-04-11-2019_bc-taxo-cleaningDM.csv"
fnMOTU.out = "./EC2M3-WGS-motus_rev-04-11-2019_bc-taxo_sumGenus-controls-cleaningDM-GENUSNAME.csv"
fnNormaux = "EC2M3_CCR1-16S_Normal.csv"

###
### WGS - Only mOTUS-v2 (already normalized)
###

MOTU = read.csv( fnMOTU, header=T, sep="\t", stringsAsFactors=F )

normaux = read.csv( fnNormaux, header=T, sep="\t" )
normaux$sampleID = gsub( "-", ".", normaux$sampleID )

###
### Now, we deduplicate (aka some mOTUS are associated to the same Genus name)
### Normalised abundance is summed (because all the data.Frame
### summed to 1 for each column)
###

# change NA to g__unknwon
ii = which( is.na(MOTU$Genus) )
MOTU$Genus[ii] = "g__unknown"

uMOTU = unique( MOTU$Genus )
iMOTU = 11:ncol(MOTU)
nMOTU = as.data.frame( matrix( 0, nrow=length(uMOTU), ncol=length(iMOTU) ) )
colnames(nMOTU) = colnames(MOTU)[iMOTU]
rownames(nMOTU) = uMOTU
for( ti in uMOTU )
{
	cat("Do: ", ti, "\n" )
	nMOTU[ti, ] = colSums( MOTU[ which( MOTU$Genus == ti ), iMOTU ] )
}

###
### Keep only Controls from 16S
###

colsNormauxMOTU = match( normaux$sampleID, colnames(nMOTU) )
nMOTU = nMOTU[ , colsNormauxMOTU ]
colnames(nMOTU) = paste0( colnames(nMOTU), ".mOTU", sep="" )
nMOTU$genus = rownames(nMOTU)

###
### Export
###

write.table( nMOTU, fnMOTU.out, sep="\t", row.names=F )
