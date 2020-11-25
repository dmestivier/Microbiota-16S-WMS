###
### Construit le diagramme de Venn de l'agregation WGS+16S
###


rm(list=ls())
library(VennDiagram)
library(RColorBrewer)

###
### Lecture des fichiers
###

fn = "../Tables/CMP_mean-abundance_Normaux_MOTU-MALT-QIIME-DADA2_counts-genus_cleaningDM-GENUSNAME_rev-01-06-2020.csv"

m = read.csv( fn, header=T, sep="\t", stringsAsFactor=F )

###
### Some Genus are "unusable" such as "g__uncultured..." or "s__" or "g__unknwon..."
### because they are not shared by the different databases used or represent some
### little variations between genus names (for example between QIIME1 and DADA2).
### They are filtered
###

L  = grep( "s__", m$Genus )
m = m[ -L, ]

L = grep( "g__uncultured", m$Genus )
m = m[ -L, ]

L = grep( "g__unknown", m$Genus )
m = m[ -L, ]

L = grep( "^[0-9]", m$Genus )
m = m[ -L, ]

###
### Generate sets of words (Genus) for each pipeline
###

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
	filename = 'Figure-2A-VENN-resolution300dpi.png',
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
