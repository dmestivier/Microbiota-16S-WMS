###
### Lit une matrice de comptage avec 3 pipeline / sample
### pour toutes les bacteries communes aux 3 pipelines
### et calcule , pour chaque bacterie, un ICC
###

library(psych)

rm(list=ls())

###
### Filename of abundances for 16S
###

fn = "../Tables/16S/CMP_CCR116S-Normaux-norm-mergingALL-MALT-QIIME-DADA2_counts-genus_cleaningDM-GENUSNAME_rev-23-03-2020.csv"

###
### Load abundance table + genus (1st col)
###

cpt = read.csv( fn, h=T, sep="\t", row.names=1 )

###
### Five dataset where not computed using DADA2
### REmove them
###

remo = grep( "FR.116", colnames(cpt))
cpt  = cpt[, -remo ]
remo = grep( "FR.169", colnames(cpt))
cpt = cpt[, -remo ]
remo = grep( "FR.276", colnames(cpt))
cpt = cpt[, -remo ]
remo = grep( "FR.400", colnames(cpt))
cpt = cpt[, -remo ]
remo = grep( "FR.716", colnames(cpt))
cpt = cpt[, -remo ]

###
### Sort FR samplenames (put Pipelines together)
###

L = order( colnames(cpt) )
cpt = cpt[, L ]

###
### Il faudrait aussi ajouter l'abondance moyenne
### et les concordances... avant la "gestion" des is.na()
###

###
### Change NA to 0
###

cpt[ is.na(cpt) ] = 0 

###
### for each genus
###

for( bi in seq(from=1, to=nrow(cpt), by=1 ) )
{
	# extrait les comptages et Genus
	xi = cpt[ bi, ]
	#gi = tax[ bi, "Genus" ]
	gi = rownames(cpt)[bi]
	gi = as.character(gi)

	# les samples sont par ordre et les pipelines aussi
	# a l'interieur d'un individu. On peut donc remplir
	# directement une matrix (lig=indiv x col = pipeline)
	w = matrix( xi, ncol=3, byrow=T )

	# Calcule les indices ICC
	res = ICC(w)
	# pb !!! 	rrr = icc(w, model="oneway", type="consistency", unit="single" )
	#           il faut des nombres au lieu de reels ?
	
	# Output
	cat( "ID:\t", bi )
	cat( "\tGenus:\t", gi )
	cat( "\tICC1_single_raters_absolute\t", res$results$ICC[1] )
	cat( "\tICC2_single_random_raters\t", res$results$ICC[2] )
	cat( "\tICC3_single_fixed_raters\t", res$results$ICC[3] )
	cat( "\tICC1-p-value\t", res$results$p[1] )
	cat( "\tiCC2-p-value\t", res$results$p[2] )
	cat( "\tICC3-p-value\t", res$results$p[3] )
	cat( "\n" )
}
