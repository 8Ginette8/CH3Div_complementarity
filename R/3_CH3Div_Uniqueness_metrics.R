# ###################################################################################
# CH3Div: Uniqueness + diversity calcualtions
#
# $Date: 2024-08-23
#
# Author: Yohann Chauvier, yohann.chauvier@wsl.ch
# Aquatic Ecology Group
# Swiss Federal Research Institute EAWAG
# 
# Description: we want here to calculate per species their phylo-
# genetic and funtional uniqueness of Switzerland (additionally, we can infer the
# rarity range and distinctiveness).
#
# ###################################################################################

### ==================================================================
### Initialise system
### ==================================================================


# Cleaning
rm(list = ls()); graphics.off()

# max.print
options(max.print=500)

# R Library
library(terra)
library(funrar)
library(ape)
library(bigmemory)

# Setting working directory
setwd("./data/")

# Functions
scr = list.files("../R/functions",full.names=TRUE)
invisible(lapply(scr, source))


### ==================================================================
### Open and prepare files
### ==================================================================


# Infer the taxa table with GUILDS + JETZ habitats
imp.f = list.files("./x_TraitCH_1.0/missRanger_imputed_traits")
imp.f = imp.f[grepl("S_MEAN_s2z_",imp.f)]
imp.o = lapply(imp.f,function(x){
	oo = read.table(paste0("./x_TraitCH_1.0/missRanger_imputed_traits/",x),header=TRUE)
	oo = oo[,c("Species","Group","Ter","Aqu")]
	oo$Group = gsub("S_MEAN_s2z_raw_traits_|.txt","",x)
	return(oo)
})
imp.o = do.call("rbind",imp.o)

# Open all phylogenies
phylo.f = list.files("./x_OpenTreeOfLife/")
phylo.o = lapply(phylo.f,function(x) readRDS(paste0("./x_OpenTreeOfLife/",x)))
names(phylo.o) = phylo.f

# Open all functional trees
funct.f = list.files("./x_functional_trees")
funct.o = lapply(funct.f,function(x) readRDS(paste0("./x_functional_trees/",x)))
names(funct.o) = funct.f

# Create a 2nd name version of the raster files
files.m = list.files("./SDMapCHv1_10km")
files.m = gsub("10km_|_reg_.*","",files.m)
files.mm = files.m
sp.subsp = sapply(strsplit(files.m,"\\."),function(x) length(x)==2 | any(x%in%"x"))
files.mm[sp.subsp] = files.m[sp.subsp]
files.mm[!sp.subsp] = files.m[!sp.subsp]
files.mm = gsub("\\."," ",files.mm)
files.mm[!sp.subsp] = unlist(lapply(strsplit(files.mm[!sp.subsp]," "),function(x) paste(x[1],x[2],"subsp.",x[3])))
corsp.names = data.frame(names.maps=files.m,names.traits=files.mm) # check unique or not
corsp.names = merge(corsp.names,imp.o,by.x="names.traits",by.y="Species",all.x=TRUE)


### ==================================================================
### Create big.matrix community data (250m)
### ==================================================================


# Open a default file to know which dim should be our stored big.matrices
r.test = rast("./SDMapCHv1_10km/10km_Abax.baenningeri_reg_covariate_ensemble.tif")
nrow.d = dim(r.test)[1]*dim(r.test)[2]
ncol.d = length(list.files("./SDMapCHv1_10km"))

# Create empty big.matrix and fill in them with N-SDMs distribution data
com.mat = big.matrix(
	nrow.d,ncol.d,
	type = "integer",
	backingfile = paste0("NSDM_cm.bin"),
	descriptorfile = paste0("NSDM_cm.desc"),
	backingpath = "../outputs/community_matrix/",
	dimnames = list(NULL,corsp.names[,1]),
	init = NA
)
lfiles = list.files("./SDMapCHv1_10km/")
for (i in 1:length(lfiles)) {
	cat('\r',round(i*100/length(lfiles),2),"%")
	sp.ras = rast(paste0("./SDMapCHv1_10km/",lfiles[i]))
	str.tar = gsub("10km_|_reg_.*","",lfiles[i])
	com.mat[,colnames(com.mat)%in%gsub("\\."," ",str.tar)] = c(sp.ras[])
}
com.mat = attach.big.matrix("../outputs/community_matrix/NSDM_cm.desc")


### ==================================================================
### Calculate uniqueness indices (function and phylo) of each species
### ==================================================================


# Prepare file output
corsp.names$PhyloUniqueness = NA
corsp.names$FunctionUniqueness = NA
corsp.names$Restrictedness = NA

# For each taxa group we subset the current community matrix
gtaxa = names(table(corsp.names$Group))
for (i in 1:length(gtaxa))
{
	cat("\n","Processing","--",gtaxa[i],"--","\n")

	# Extract N-SDMs species from the target group
	tar.sp = corsp.names[corsp.names$Group%in%gtaxa[i],]

	# Target files
	tar.phylo = phylo.o[grepl(gtaxa[i],phylo.f)][[1]]
	tar.funct = funct.o[grepl(gtaxa[i],funct.f)][[1]]
	tar.comm = com.mat[,colnames(com.mat)%in%tar.sp$names.traits]

	# Quick coverage test
	print(all(tar.sp$names.traits%in%tar.phylo$tip)) # Normal as the phylo cannot include all names
	print(all(tar.sp$names.traits%in%tar.funct$tip))

	# Normalize distances from 0 to 1 + convert community matrix to relvative values (Grenier & al. 2018)
	phylo.distM = normalize(cophenetic.phylo(tar.phylo))
	funct.distM = normalize(cophenetic.phylo(tar.funct))
	tar.comm = make_relative(tar.comm)

	# Calculate overall uniqueness for each species (based on current "known" distribution)
	phyloU = uniqueness(tar.comm,phylo.distM)
	functU = uniqueness(tar.comm,funct.distM)
	restrict = restrictedness(tar.comm)

	# Add mean for missing uniqueness info
	sp.missing = tar.sp$names.traits[!tar.sp$names.traits%in%tar.phylo$tip]
	if (length(sp.missing)>0){
		phyloU = rbind(phyloU,data.frame(species=sp.missing,Ui=mean(phyloU$Ui)))
		phyloU = phyloU[order(phyloU$species),]
	}

	# Assign
	if (all(phyloU$species == tar.sp$names.traits)) {tar.sp$PhyloUniqueness = phyloU$Ui} else {print(FALSE)}
	if (all(functU$species == tar.sp$names.traits)) {tar.sp$FunctionUniqueness = functU$Ui} else {print(FALSE)}
	if (all(restrict$species == tar.sp$names.traits)) {tar.sp$Restrictedness = restrict$Ri} else {print(FALSE)}

	# Save
	write.table(tar.sp,paste0("../outputs/sp_uniqueness/Ui_indices_",gtaxa[i],".txt"),row.names=FALSE)
}