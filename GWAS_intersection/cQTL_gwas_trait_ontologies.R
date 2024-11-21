library(data.table)
library(ontologyIndex)
library(webr)
library(ggplot2)

# GWAS catalog SNPs overlaping cQTLs accounting for LD (R2>0.8)
ol = fread("gwascat_consensus_LD_r2unphased_08_revision_phase3.txt")
length(unique(ol$`DISEASE/TRAIT`)) # 304
gwascat_ont = fread("gwas_catalog_v1.0.2-associations_e112_r2024-06-07.tsv")

ol1 = merge(ol, gwascat_ont[, c("SNPS", "DISEASE/TRAIT", "MAPPED_TRAIT_URI", "MAPPED_TRAIT")], 
            by.x=c("ID_B", "DISEASE/TRAIT"), by.y=c("SNPS", "DISEASE/TRAIT"))

gwas_traits = gsub("http://www.ebi.ac.uk/efo/", "", ol1$MAPPED_TRAIT_URI)
gwas_traits = unlist(strsplit(gwas_traits, ", "))

# Get ontologies for our cQTLs' GWAS traits
efo = get_ontology("http://www.ebi.ac.uk/efo/efo.obo")
ancestorsList = sapply(unique(gwas_traits) , function(x){
  if(length(efo$ancestors[[x]])>0) {  return(efo$ancestors[[x]]) }
  else{ 
    if(length(efo$ancestors[[gsub("_", ":", x)]])>0) { 
      return(efo$ancestors[[gsub("_", ":", x)]])  }
    else {return(NULL)}
  }
})
names(ancestorsList)<-sapply(unique(gwas_traits), 
    function(x)ifelse(length(efo$ancestors[[x]])>0, x, gsub("_", ":", x)))

ancestors=unlist(ancestorsList)

# Cut off at min 5 traits 
topEFO = sort(table(ancestors)[grep("EFO", names(table(ancestors)))],decreasing = T)[1:35] 

## For each "top" category,
## go through all GWAS trait categories and see if this top category is in the list of its ancestors
## if it is, get all cQTLs that overlap this GWAS trait category (note it can include >1 GWAS), accounting for LD (r2>0.8)
## then put cQTLs for this top category together across the GWAS traits and count their total

nSnps = sort(sapply(names(topEFO), function(broader_efo)
  length(unique(unlist(lapply(names(ancestorsList), function(gwas_trait)
    if(broader_efo%in%ancestorsList[[gwas_trait]])
      ol1[MAPPED_TRAIT_URI==paste0("http://www.ebi.ac.uk/efo/", gsub(":", "_", gwas_trait))]$ID_A)
  )))
), decreasing=T)


### Figure 7A (pie chart appearance finalised manually)
bq = fread("GUESS_BaseQTL_contactQTLs_hg19_rsids_LD08.txt")
sumTable = data.table(snp = unique(bq$ID_A))
sumTable[, isGWASld := snp %in% ol1$ID_A]
sumTable[, isGWASol := snp %in% ol1$ID_B]
PieDonut(sumTable, aes(isGWASld, isGWASol),col="cyan") # GWAS_stats_pie.pdf

### Figure 7B
namedNSnps = nSnps
names(namedNSnps) = sapply(names(nSnps), function(x)efo$name[x])
namedNSnps = namedNSnps[3:length(namedNSnps)]
namedNSnps = namedNSnps[namedNSnps>=10]
par(mar = c(5, 12, 4, 2) + 0.1)
barplot(sort(namedNSnps), horiz = T, las=1, cex.names=0.6, cex.lab=0.8,
        xlab="Number of cQTLs overlapping GWAS SNPs (R2>0.8)", col="cyan") 