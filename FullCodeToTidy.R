#load packages
library(stringr)
library(tidyr)
library(XML)
library("rentrez")
library(RColorBrewer)
library(ggplot2)
library(ggtree)
library(ape)
library(phangorn)
library(rgl)
library(caTools)
library(magick)


#sequences of B-tubulin gene (rep sequences) from all species in Epichloe Taxa
#Epichloe_taxa <- read.csv("~/Summer Scholarship 2019/21_species_tree/Epichloe_taxa.csv", stringsAsFactors = FALSE)
#all_recs <- entrez_fetch(db="nuccore", id = Epichloe_taxa$rep_seq, rettype="fasta", stringsAsFactors=FALSE)
#cat(strwrap(substr(all_recs, 1, 500)), sep="\n")
#write(all_recs, file="unaligned.fasta")


#make tree for ASTRAL comparison. Tidied up into a little function that allows 
#you to specify a path containing trees and the minimum number of tips to
#include.
read_trees <- function(tree_dir, min_tax = 16){
    #assumes directory _only_ includes trees
    all_path_trees <- list.files(tree_dir, full.names = TRUE) #filepath to best trees
    all_trees <- lapply(all_path_trees, read.tree)
    to_return <- all_trees[ sapply(all_trees, Ntip) == min_tax ]
    class(to_return) <- c("multiPhylo", class(to_return))
    to_return
}

all_sixteen_tips <- read_trees("~/Shared Folder/fifteen_trees/")
all_trees <- read_trees("~/Shared Folder/fifteen_trees/", min_tax=1)
#write.tree(all_sixteen_tips, file = "all_full_trees.phy") #newick format for ASTRAL (use write.nexus for nexus format)
write.tree(all_read_trees, file = "all_trees.phy") #all_trees file not used

#reinput ASTRAL tree to root + show node confidences
species_tree <- read.tree("~ASTRAL FILE")
rooted_tree <- root(species_tree, "Cpur") #root with tag name at Cpur
plot.phylo(rooted_tree, show.tip.label = TRUE, show.node.label = TRUE)

#make consensus net
cnet <- consensusNet(all_sixteen_tips, .2) #change for different relationships
plot(cnet, "2D") #can show.edge.label =TRUE for clarity
plot(cnet)
play3d(spin3d(axis=c(0,1,0), rpm=6), duration=10)
# create animated gif file - messy at current
movie3d(spin3d(axis=c(0,1,0), rpm=6), duration=10, dir = "~")

#produce dist_topo, see any obvious discrepancies / unpublishable
dist_topo <- dist.topo(all_sixteen_tips, method = "PH85")
tree <- upgma(dist_topo)
plot(tree)




#using B tub gene to produce tree and analysis

tree <- read.tree("~/Summer Scholarship 2019/21_species_tree/RAxML_bestTree.tree")
plot.phylo(tree)
rooted_tree <- root(tree, "JX083460.1")
old_names <- rooted_tree$tip.label
name_map <- c(L78276.1="glyceriae", L78290.1="bromicola", AF250743.1="elymi", GQ421707.1="sibirica", EF422757.1="ganuensis", AF457494.1="ganuensis_inebrian", JX083460.1="C_purpurea", L78271.1="brachyelyti", AF323371.1="aotearoae", AF250757.1="typhina_poae", AF176266.1="typhina_poae_canariensis", AY707694.1="typhina_poae_aonikenkana", AF457493.1="typhina_poae_huerfana", KC936108.1="sylvatica_sylvatica", L78281.1="typhina_clarkii", L78288.1="typhina_typhina", JF718439.1="sylvatica_pollinensis", L06961.1="baconii", EU526824.1="stromatolonga", KC936144.1="festucae_lolli", L06955.1="festucae", AF457469.1="mollis", L06958.1="amarillans")
new_names <- name_map[old_names]
plot.phylo(rooted_tree, show.tip.label = FALSE, x.lim = 0.6)
tiplabels(new_names, adj = c(-0.1), frame = "none")


tip_names <- read.csv("~/Summer Scholarship 2019/21_species_tree/Epichloe_taxa.csv", stringsAsFactors = FALSE)
tip_names$rep_seq <- paste0(tip_names$rep_seq, ".1") #change names to suit
name_map <- tip_names$tag_name
names(name_map) <- tip_names$rep_seq
# check all tips are in name map   #any(is.na(name_map[rooted_tree$tip.label]))
rooted_tree$tip.label <- name_map[rooted_tree$tip.label]
calib <- ape::makeChronosCalib(rooted_tree, age.min=58, age.max=59)
time_tree <- chronos(rooted_tree, lambda = 0.2, calibration = calib)
plot(time_tree)
axisPhylo()

epi_info <- read.csv("~/Summer Scholarship 2019/MetaDataTree1/Epichloe_taxon_data_updated.csv")
no_root_tree <- drop.tip(rooted_tree, " Claviceps")
sibiricaroot <- root(no_root_tree, "sibirica") #Murray Root Advice
gtree <- ggtree(sibiricaroot)
join <- gtree %<+% epi_info
x <- join + geom_tiplab(offset = 0, hjust = 0) + geom_tiplab(aes(color = sexual)) + theme(legend.position = "right", legend.title = element_text()) + labs(fill = "Sexual Reproduction") + scale_color_brewer(palette = "Set2") + geom_point(aes(color = genome))

geo_mat <- epi_info[,6:11]
row.names(geo_mat) <- epi_info$tag_name

cols=c("magenta2", "magenta4")
gheatmap(x, geo_mat, offset = 0.015, colnames_angle = 90, colnames_offset_y = -1, color = "black", width = 0.5) + scale_fill_manual(values=cols) + theme(legend.title = element_text()) + labs(fill = "Distribution") + labs(title = "Sexual reproduction and global distribtion of Epichloe")  
#heatmap of geography


#producing heatmap of host information
subset_meta <- subset(EpichloeTaxaComma, HybridStatus %in% "non-hybrid")
six <- EpichloeTaxaComma[6,]
eighteen <- EpichloeTaxaComma[18,] #adding 6 and 18 due to formatting of Epichloe data, doesn't automatically include as non-hybrids
addedsix <- rbind(subset_meta, six)
allnhy <- rbind(addedsix, eighteen)
allnhy <- allnhy[order(as.numeric(rownames(allnhy))),,]
HostNumber <- str_count(allnhy[,"KnownHostRange"], ",") + 1  #number of hosts
HostSubStr <- str_split(allnhy[,"KnownHostRange"], ",")
splited <- str_split(HostSubStr, boundary("word"))
HostGenus <- str_extract_all(splited, "\\b[A-Z]\\w+") #Host Genus as extracted with capital letters

ugen <- sapply(HostGenus, unique)
host_assoc <- data.frame(taxon = rep(allnhy$ï..Taxon, times=lengths(ugen)), plant = unlist(ugen)) 

"""
#produce tree of grass host species
taxon_names <- rotl::tnrs_match_names(names = unique(host_assoc$plant))
easy_ones <- taxon_names[taxon_names$number_matches == 1,]
grass_tree <- rotl::tol_induced_subtree(ott_ids = easy_ones$ott_id)
plot(grass_tree)
plot(grass_tree, type = 'u')
"""

host_assoc$observed <- "present"
host_assoc <- host_assoc[-c(46,47), ] #remving additional species shouldn't be included in analysis

wide_host <- pivot_wider(host_assoc, names_from=plant, values_from=observed, values_fill= list(observed="absent"))
wide_host[5,1] <- "EpichloÃ« bromicola                      Leuchtm. & Schardl"  #adding information not included

Tag <- Epichloe_taxon_data_updated[,c(1,2)] #dataframe with just tagname and scientific name
Host_Tags <- merge(Tag, HostDistribution1, by = "taxon") #HostDistribtion one is altered wide_host (removed spaces between taxon names to enable merging)
Host_Tags <- Host_Tags[,-3] #remove row number column from previous dataset
write.csv(Host_Tags, "MergeBromicola.csv")  #two bromicola rows, merging manually
Host_Tags <- read.csv("~/MergeBromicola.csv")

Hosts <- Host_Tags[,(4:30)] #just rows of host Genus
row.names(Hosts) <- Host_Tags$tag_name

#making heatmap with host Genus
host <- gheatmap(x, Hosts, offset = 0.01, colnames_angle = 90, colnames_offset_y = -2, color = "black", width = 0.5) + scale_fill_brewer(palette = "Paired") + theme(legend.title = element_text()) + labs(fill = "Host Presence") + labs(title = "Sexual reproduction and hosts of Epichloe") + geom_tiplab(align=TRUE)  
ggsave("Hosts.pdf", plot = last_plot(), height = 12, width = 28)

#making hosts into tribes (families and subfam identical)

Genera <- colnames(Hosts)
Generarep <- lapply(Genera, entrez_search, db="taxonomy")
Genera_rec <- entrez_fetch(db="taxonomy", id=Generarep, rettype="xml", parsed=TRUE)
family <- get_taxon(Genera_rec, "family")
ufam <- unique(family)

tag_names <- Host_Tags$tag_name

columns <- c("Genus", "Family", "tag_name")
GenusFamily <- data.frame(tag_names)
GenusFamily$Genus <- HostGenus

#making vectors of family, subfamily and tribe. Have to iterate over all Genera[] componenets individually due to memory restrictions
GenFam <- character()
GenFam <- c(GenFam, (get_family(Genera[27]))) #remove [6] as two families one Genus
GenFam

GenSubFam <- character()
GenSubFam <- c(GenSubFam, (get_subfamily(Genera[27])))
GenSubFam

GenTribe <- character()
GenTribe <- c(GenTribe, (get_tribe(Genera[27])))
GenTribe

Taxonomy <- data.frame(Genera, GenFam, GenSubFam, GenTribe)
HostT <- Hosts
colnames(HostT) <- Taxonomy$GenTribe

write.csv(HostT, file = "HostTribes.csv")

host <- gheatmap(z, HostT, offset = 0.01, colnames_angle = 90, colnames_offset_y = -2, color = "black", width = 0.5) + scale_fill_brewer(palette = "Paired") + theme(legend.title = element_text()) + labs(fill = "Host Presence") + labs(title = "Sexual reproduction and hosts of Epichloe") + geom_tiplab(align=TRUE)  
ggsave("HostTribe.pdf", plot = last_plot(), height = 12, width = 28)

get_tribe <- function(Genu){
  taxon_search <- entrez_search(db="taxonomy", term = Genu)
  taxon_rec <- entrez_fetch(db="taxonomy", id=taxon_search$ids, rettype="xml", parsed=TRUE)
  get_taxon(taxon_rec, "tribe")
}

get_subfamily <- function(Genu){
  taxon_search <- entrez_search(db="taxonomy", term = Genu)
  taxon_rec <- entrez_fetch(db="taxonomy", id=taxon_search$ids, rettype="xml", parsed=TRUE)
  get_taxon(taxon_rec, "subfamily")
}

get_family <- function(Genu){
  taxon_search <- entrez_search(db="taxonomy", term = Genu)
  taxon_rec <- entrez_fetch(db="taxonomy", id=taxon_search$ids, rettype="xml", parsed=TRUE)
  get_taxon(taxon_rec, "family")
}

get_taxon <- function(rec, rank){
  xp  <-paste0("//LineageEx/Taxon/Rank[.='", rank,"']/../ScientificName")
  res <- xpathSApply(rec, xp, xmlValue)
  if(is.null(res)){
    return(NA)
  }
  res
}
