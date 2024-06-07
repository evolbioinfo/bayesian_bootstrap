params.msa = "data/alignment.fasta"
params.results="results"
params.collapse=0.1
params.collapseref=false
params.nboot = 200

nboot = params.nboot
results=params.results
msa = file(params.msa)
collapse=params.collapse
collapseref=params.collapseref

// Auto-detects alphabet (protein / nucleotide)
process alphabet {
    label 'goalign'

    tag "$id"
    publishDir "${results}${id}/", mode: 'link'

    input:
    tuple val(id), file(msa)

    output:
    tuple val(id), stdout
    tuple val(id), path("alphabet.txt")

    script:
    """
    goalign stats alphabet --auto-detect -i ${msa} > alphabet.txt
    cat alphabet.txt
    """
}

// Reformats input alignment into phylip and clean sequence names 
// (in case there are newick specific characters)
process rename {
    label 'goalign'

    tag "$id"
    publishDir "${results}${id}/", mode: 'link'

    input:
    tuple val(id), file(msa)

    output:
    tuple val(id), path("${msa.baseName}_renamed.phy")

    script:
    """
    goalign reformat phylip -i $msa --auto-detect | goalign rename --clean-names -o ${msa.baseName}_renamed.phy
    """
}

// Computes the length of the alignment
process alilen {
    label 'goalign'

    tag "$id"
    
    input:
    tuple val(id), file(msa)

    output:
    tuple val(id), stdout

    script:
    """
    printf \$(goalign stats length -p -i ${msa})
    """
}

// Infers reference tree
process inferRefTree {
    tag "$id"
    publishDir "${results}${id}/", mode: 'link'

    label 'phyml'

    input:
    tuple val(id), val(alphabet), file(msa)

    output:
    tuple val(id), file("${msa}_phyml_stats.txt")
    tuple val(id), file("${msa}_phyml_tree.txt")

    script:
    if( alphabet == 'nucleotide' )
    """
    phyml -i ${msa} -m GTR -c 4 -d nt -a e -o tlr -b 0 --r_seed 123456
    """
    else
    """
    phyml -i ${msa} -m LG -c 4 -d aa -a e -o tlr -b 0 --r_seed 123456
    """
}

// Generates Bayesian Bootstrap weights
process genWeightBoot {
    label 'goalign'

    tag "$id"
    
    publishDir "${results}${id}/", mode: 'link'

    input:
    tuple val(id), file(msa)
    val nboot

    output:
    tuple val(id), file("weights.txt")

    script:
    """
    goalign build weightboot -p -i ${msa} -n ${nboot} > weights.txt
    """
}

// Generate Bayesian Bootstrap trees
process inferWeightBootTrees {
     tag "$id"
     label 'phyml'

     input:
     tuple val(id), val(alphabet), file(w), file(msa)

     output:
     tuple val(id), file("${msa}_phyml_stats.txt")
     tuple val(id), file("${msa}_phyml_tree.txt")

     script:
     if( alphabet == 'nucleotide' )
     """
     phyml -i ${msa} --weights=${w} -m GTR -c 4 -d nt -a e -o tlr -b 0 --r_seed 123456
     """
     else
     """
     phyml -i ${msa} --weights=${w} -m LG -c 4 -d aa -a e -o tlr -b 0 --r_seed 123456
     """
}

// Computes Bayesian Bootstrap Supports
process computeWeightSupports {
     tag "$id"

     publishDir "${results}${id}/", mode: 'link'
     
     input:
     tuple val(id), val(length), file(ref), file('boot')
     val collapse
     val collapseref

     output:
     tuple val(id), file("reftree_weightsupport.nw")

     script:
     if( collapseref )
     """
     gotree collapse length -i <(cat $boot) -l ${collapse/length} |  gotree compute support fbp -i <(gotree collapse length -i ${ref} -l ${collapse/length}) -b - -o reftree_weightsupport.nw
     """
     else
     """
     gotree collapse length -i <(cat $boot) -l ${collapse/length} |  gotree compute support fbp -i ${ref} -b - -o reftree_weightsupport.nw
     """
}

process computeMetrics {
    label 'goalign'

    tag "$id"
    publishDir "${results}${id}/", mode: 'link'


    input:
    tuple val(id), val(alphabet), path(msa), path(supporttree)

    output:
    tuple val(id), path("*homoplasies.txt")

    script:
    if( alphabet == 'nucleotide' )
    """
    ALI=\$(goalign stats char -p --per-sites -i ${$msa} | pars.pl)
    LEN=\$(goalign stats length -p -i $msa)
    ML=\$(gotree stats -i $supporttree | cut -f 6| tail -n 1)
    MLHOMO=\$(awk -v ml=\$ML -v len=\$LEN -v ali=\$ALI 'BEGIN{print (ml*len-ali)*100/ali}')
    echo "ID\tAlphabet\tAliParsimony\tAliLength\tHomoplasy" > ${id}_homoplasies.txt
    echo -e "$id\t$alphabet\t\$ALI\t\$LEN\t\$MLHOMO" >> ${id}_homoplasies.txt
    """
    else
    """
    ALI=\$(goalign stats char -p --per-sites -i ${$msa} | pars_prot.pl)
    LEN=\$(goalign stats length -p -i $msa)
    ML=\$(gotree stats -i $supporttree | cut -f 6| tail -n 1)
    MLHOMO=\$(awk -v ml=\$ML -v len=\$LEN -v ali=\$ALI 'BEGIN{print (ml*len-ali)*100/ali}')
    echo "ID\tAlphabet\tAliParsimony\tAliLength\tHomoplasy" > ${id}_homoplasies.txt
    echo -e "$id\t$alphabet\t\$ALI\t\$LEN\t\$MLHOMO" >> ${id}_homoplasies.txt
    """
}

process treeStats {
    input:
    tuple val(id), path(supporttree)

    output:
    tuple val(id), path("*edges_ref.txt")

    script:
    """
    gotree stats edges -i $supporttree > ${id}_edges_ref.txt
    """
}

process drawFigures {
    tag "$id"

    publishDir "${results}${id}/", mode: 'link'

    input:
    tuple val(id), path(edges), val(alilen)

    output:
    tuple val(id), path("*.svg")

    script:
    """
#!/bin/env Rscript
library(ggplot2)
library(forcats)
library(zoo)
library(dplyr)
library(plyr)

reftree = read.table($edges,header=T,sep="\t",na.strings="N/A")

tmplen=reftree[reftree\$terminal=="false", ]
tmplen[order(tmplen\$length*$alilen),"sortidx"]=seq(1,length(tmplen\$length)
tmplen\$sortidx=tmplen\$sortidx/length(tmplen\$sortidx)
svg("$id_branch_lengths.svg",width=4.7,height=2.2)
ggplot(tmplenhomo,aes(y=sortidx,x=$alilen*length))+geom_point(size=0.10)+theme_bw()
dev.off()
    """
}


workflow {
    msa = Channel.fromPath(params.msa).map{it -> [it.baseName(), it]}
    renamed=rename(msa)
    len = alilen(renamed).map{it ->  [it[0],Integer.parseInt(it[1].trim())]}
    alpharaw = alphabet(renamed)
    alphastr=alpharaw[0].map{it -> [it[0],it[1].trim()]}

    reftree = inferRefTree(alphastr.combine(renamed, by:0))

    weights=genWeightBoot(renamed,nboot).splitText(by: 1, file:true)
    weightboot = inferWeightBootTrees(alphastr.combine(weights, by:0).combine(renamed, by:0))
    weightboottrees = weightboot[1].groupTuple()
    weightsupport = computeWeightSupports(len.combine(reftree[1], by:0).combine(weightboottrees, by:0), collapse, collapseref)

    metrics = computeMetrics(alphastr.combine(renamed, by:0).combine(weightsupport,by: 0)))
    stats = treeStats(weightsupport)
    figures = drawFigures(stats.combine(len, by: 0))
}
