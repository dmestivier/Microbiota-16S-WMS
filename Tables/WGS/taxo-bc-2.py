###
### lit un fichier de taxo genere via un BIOM
### et s'assure que les niveaux taxo sont au 
### bon endroit, car on peut avoir des ;p__.._;o__;...
###
### Le fichier commence par K colonnes (tab-sep)
### qui ne sont pas touchees. La valeur de K est obtenue
### du header qui va recherche "Domain" ou "Kingdom"
###

import sys
import string

###
### new taxpath
###

def new_tp():
    n = {}
    n["d__"] = "NA"
    n["p__"] = "NA"
    n["c__"] = "NA"
    n["o__"] = "NA"
    n["f__"] = "NA"
    n["g__"] = "NA"
    n["s__"] = "NA"
    n["NA" ] = "NA" # junk case

    return n

###
### Rempli le tp avec les info de la ligne
### Les k premieres colonnes ne sont pas touchees
###

def fill_tp( cols, n, k ):

    # pour chaque level dans la ligne
    for lvl in cols[k:]:
        # Dans certains cas on peut avoir des espaces en debut !!!
        lvl = lvl.lstrip()
        k = lvl[0:3]   # pour "p__"
        #k = lvl[0:2]  # pour "p_"

        if k in n.keys():
            n[k] = lvl
        else:
            pass
            #print("Error with " + lvl + " and " + k )

    return n

###
### Affiche le taxpath
###

def print_tp2( x, n, k ):
    output = "\t".join( x[0:k] )
    output += "\t" + n["d__" ] 
    output += "\t" + n["p__" ] 
    output += "\t" + n["c__" ] 
    output += "\t" + n["o__" ] 
    output += "\t" + n["f__" ] 
    output += "\t" + n["g__" ] 
    output += "\t" + n["s__" ] 

    print(output)

def print_tp( x, n ):
    output = x
    output += "\t" + n["d_" ] 
    output += "\t" + n["p_" ] 
    output += "\t" + n["c_" ] 
    output += "\t" + n["o_" ] 
    output += "\t" + n["f_" ] 
    output += "\t" + n["g_" ] 
    output += "\t" + n["s_" ] 

    print(output)

####################################################
###
### ligne de commande
###
####################################################

if len(sys.argv)!=2:
    print("Syntaxe : %s export-BIOM.csv\n" % sys.argv[0] )
    sys.exit(1)

# mon fichier
fn = sys.argv[1]

# GO
fid = open( fn, "rt" )
header = fid.readline()
# Attention la ligne peut contenir des taxa entre '"'
header = header.replace('"','')
print(header[:-1])

# trouve la colonne de debut de la taxo
cols = header.split("\t")
if "Domain" in cols:
    k = cols.index( "Domain" )
elif "Kingdom" in cols:
    k = cols.index( "Kingdom" )
else:
    print("Error, aucun debut de taxonomie trouve !!!" )
    sys.exit(1)

sys.stderr.write("Fichier : " + fn + ". La taxonomie commence en colonne " + str(k) + "\n" )

###
### GO
###

for lig in fid:
    # Attention la ligne peut contenir des taxa entre '"'
    lig = lig.replace('"','')
    # decoupe a ligne
    cols = lig[:-1].split("\t")
    # reconstruit le taxpath
    tp = new_tp()
    tp = fill_tp( cols, tp, k )
    #print_tp( lig[:-1], tp )
    print_tp2( cols, tp, k )

fid.close()
