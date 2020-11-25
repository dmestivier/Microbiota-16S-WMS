# read a file (output from mOTUS (merged))
# some taxpath contain empty levels which cause problems
# when we try to import later.
# Fix this.

###
### new taxpath
###

def new_tp():
    n = {}
    n["k__"] = "NA"
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

def fill_tp( cols, n ):

    # pour chaque level dans la ligne
    for lvl in cols:
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
### (change k__ to d__ )

def print_tp( n ):
    output  = n["k__" ].replace( "k__", "d__" )
    output += "\t" + n["p__" ]
    output += "\t" + n["c__" ] 
    output += "\t" + n["o__" ] 
    output += "\t" + n["f__" ] 
    output += "\t" + n["g__" ] 
    output += "\t" + n["s__" ] 

    return output

################################################

import sys

# filename
fn = "EC2M3-WGS-motus_rev-04-11-2019.csv"

fid = open( fn, "rt" )
header = fid.readline()
header = header[:-1]
# Change colname n.2
header = header.split("\t")
header[1] = "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies"
print( "\t".join(header))

# for each line
for lig in fid:
    cols = lig[:-1].split("\t")
    taxpath = cols[1].split("|")
    
    tp = new_tp()
    tp = fill_tp( taxpath, tp )

    cols[1] = print_tp( tp ) 
    print( "\t".join(cols) )
# fin
fid.close()
