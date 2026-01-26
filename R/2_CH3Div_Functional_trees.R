# ###################################################################################
# CH3Div: Generate functional trees
#
# $Date: 2024-08-14
#
# Author: Yohann Chauvier, yohann.chauvier@wsl.ch
# Aquatic Ecology Group
# Swiss Federal Research Institute EAWAG
# 
# Description: Based on our imputed traits, we are going to generate
# functional distance trees per taxa groups.
#
# ###################################################################################

### ==================================================================
### Initialise system
### ==================================================================


# Start with house cleaning
rm(list = ls()); graphics.off()

# max.print
options(max.print=500)

# R Library
library(ape)
library(cluster)
library(vegan)

# Setting working directory
setwd("./data")

# Functions
scr = list.files("../R/functions",full.names=TRUE)
invisible(lapply(scr, source))


### ==================================================================
### Pre-processing trait data
### ==================================================================


# EU species
EU.sp = c("EEA","FaunaEuropeae","EIVE","AraneaeEuropeae",
	"Loncarevic_et_al.2024","Middleton-Welling_et_al.2020")

# % missing trait values per data.frame (per traits and overall)
all.f = list.files("./x_TraitCH_1.0/raw_traits")
ex.col = c("Species","species","group","Genus","genus","Family","family","Order","order","Class","class",
	"Phylum","phylum","GBIF_accepted","GBIF_rank","gbif_accepted","gbif_rank","originCH",
	"originEUR","Swiss_status","IUCN_status","group")
df.perct = 
lapply(1:length(all.f),function(x){
	imp.f = read.table(paste0("./x_TraitCH_1.0/raw_traits/",all.f[x]),header=TRUE); table(imp.f$Source)
	imp.f = imp.f[,!colnames(imp.f)%in%ex.col]
	imp.f = imp.f[!imp.f$Source%in%EU.sp,]
	imp.f$Source = NULL					   
	o1 = length(imp.f[!is.na(imp.f)])*100/(nrow(imp.f)*ncol(imp.f))
	o2 = apply(imp.f,2,function(y) length(y[!is.na(y)])*100/length(y))
	return(list(o1,o2))
})
names(df.perct) = all.f ;df.perct

# Load imputed s2z trait files
s2z.ref = list.files("./x_TraitCH_1.0/missRanger_imputed_traits/")
s2z.t = s2z.ref[!grepl("_v2|S_EVAL_|all_compiled_|_IMP|focus",s2z.ref)]
trait.df = lapply(s2z.t,function(x) read.table(paste0("./x_TraitCH_1.0/missRanger_imputed_traits/",x),header=TRUE))

# Load evaluation files (same order)
s2z.e = list.files("./x_TraitCH_1.0/missRanger_evaluations/")
eval.df = lapply(s2z.e,function(x) {
	oo = read.table(paste0("./x_TraitCH_1.0/missRanger_evaluations/",x),header=TRUE)
	oo = oo[!oo$Traits%in%"Originch",]
	return(oo)
})


### ==================================================================
### Filter columns in our trait data.frame
### ==================================================================


# Loop over and keep only columns with (1) >= 90% data completeness OR
# (2) OOB errors < 0.2 + > 20% trait completeness (Soria & al. 2020)
toTree = list()
for (i in 1:length(trait.df))
{
	# col of traits to select from
	col.traits = trait.df[[i]][,c("Species","Source",eval.df[[i]]$Traits)]
	col.traits = col.traits[!col.traits$Source%in%EU.sp,] 
	col.id = col.traits$Species
	col.traits[,c("Species","Source")] = NULL

	# Which columns (eval) and (%) ?
	col.k0 = df.perct[[i]][[2]][eval.df[[i]]$Traits] >= 90
	col.k1 = eval.df[[i]]$OBB_pred_error < 0.2
	col.k2 = df.perct[[i]][[2]][eval.df[[i]]$Traits] > 20
	col.k1[is.na(col.k1)] = FALSE

	# Change into numeric OR factor
	trait.filter = col.k0 | (col.k1+col.k2)%in%2
	df.out = data.frame(Species=col.id,col.traits[,trait.filter])
	cat.cols = sapply(df.out, function(x) is.character(x) | all(na.omit(x)%in%c(0,1)))
	df.out[cat.cols] = lapply(df.out[cat.cols],as.factor)
	toTree[[i]] = df.out
}
names(toTree) = gsub("S_MEAN_s2z_raw_traits_|\\.txt","",s2z.t)

# Quick generate a summary table of traits used for calculating functional trees
df = as.data.frame(lapply(seq_along(toTree), function(x) {
	dd = names(toTree[[x]])
	length(dd) = max(lengths(lapply(toTree, names)))
	return(dd)
}))
names(df) = toupper(names(toTree)) # Make sure we have the same name as in CH3Div
names(df) = c("AMPHIBIANS","BIRDS","BEES","BRYOPHYTES","BEETLES","MAY-STONE–CADDISFLIES","FISHES",
	"FUNGI","BUTTERFLIES","LICHENS","MAMMALS","MOLLUSCS","DRAGONFLIES","GRASSHOPPERS","VASCULAR_PLANTS",
	"REPTILES","SPIDERS")
df[is.na(df)] = "_"
write.csv(df[-1,], "Trait_select_summary.csv", row.names = FALSE)


### ==================================================================
### Build the 16 functional trees (making sure all names are in phylo)
### ==================================================================


# List phylogenies
phylo.f = list.files("./x_OpenTreeOfLife/")

# Loop over taxa
lapply(1:length(toTree),function(x){

	cat(names(toTree)[x],"\n")

	# Next
	qd = "./x_functional_trees"
	#if (any(grepl(names(toTree)[x],list.files(qd)))) {return(NULL)}

	# Extract only traits
	otraits = toTree[[x]][,-1]
	row.names(otraits) = toTree[[x]][,1]

	# Transform data (here 0 to 1) & apply gower distance with the new phylopars filled trait matrix
	# ===> For gower & continous: Maire & al. 2015, Pavoine 2009
	dis.traits = daisy(otraits,metric="gower",stand=FALSE)

	# Test which dendrogram is the best (Mouchet & al. 2008)
	h.toTest = c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
	hclust.scores =
	sapply(h.toTest,function(y){
		f.tree = hclust(dis.traits,method=y)
		cor.ind = cor(cophenetic(f.tree),dis.traits)
		return(1-cor.ind^2)
	})

	# Apply a mantel test to compare/permut gower and winner's distance matrix
	# (Following Brun & al. 2019, Thuiller & al. 2015/2014, Petchey & Gaston 2006)
	hclust.win = hclust(dis.traits,method=h.toTest[which.min(hclust.scores)])
	m.test = mantel(dis.traits,cophenetic(hclust.win),permutations=9999)
	print(m.test$statistic)

	# Convert, check species names and save
	hclust.ape = as.phylo(hclust.win)
	phylo.t = readRDS(paste0("./x_OpenTreeOfLife/",phylo.f[grepl(paste0("_",names(toTree)[x]),phylo.f)]))
	print(all(phylo.t$tip.label%in%hclust.ape$tip.label))
	saveRDS(hclust.ape,file=paste0("./x_functional_trees/1sfunctionaltree_",names(toTree)[x]))
})
