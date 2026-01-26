# ###################################################################################
# CH3Div: Conservation prioritizations (run on LINUX ONLY)
#
# $Date: 2024-10-31
#
# Author: Yohann Chauvier, yohann.chauvier@wsl.ch
# Aquatic Ecology Group
# Swiss Federal Research Institute EAWAG
# 
# Description: We have the N-SDMs and the uniqueness indices for each species! We
# can now set up our conservation prioritization with Zonation5
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
library(doMC)
library(foreach)
library(stringr)

# Setting working directory
setwd("./data/")

# Functions
scr = list.files("../R/functions",full.names=TRUE)
invisible(lapply(scr, source))


### ==================================================================
### Preparing folders, files and names
### ==================================================================


# Extract PAs and mask path
wpa = list.files("./templates",full.names=TRUE)
Z5.wpa = normalizePath(wpa[grepl("Z5_PAs",wpa)])
Z5.mask = normalizePath(wpa[grepl("Z5_mask",wpa)])

# Load Z5 templates
features.t = read.table("../outputs/zonation/feature_list.txt",header=TRUE)
setting.f = readLines("../outputs/zonation/example_dhW.z5")
setting.f = setting.f[1:2] # Need 3 arguments

# Antoine's file (temporara, to remove)
o.a = read.csv("1_SDMapCHv1_CH3Div_TerAqu.csv")[,-1]

# Load uniqueness-weight files + normalize from 1 to 2 + save new version
Uf = list.files("../outputs/sp_uniqueness",full.names=TRUE)
Uf = lapply(Uf[!grepl("non_zero_",Uf)],function(x) {
	#
	u.o = read.table(x,header=TRUE)
	u.o$PhyloUniqueness = normalize(u.o$PhyloUniqueness)
	u.o$FunctionUniqueness = normalize(u.o$FunctionUniqueness)
	u.o$SummedUniqueness = u.o$PhyloUniqueness + u.o$FunctionUniqueness
	u.o$SummedUniqueness[u.o$SummedUniqueness==0] = min(u.o$SummedUniqueness[u.o$SummedUniqueness!=0])
	# Add Ter/Aqu with Antoine's information (temporary, to remove these lines when new traits are ready)
	merge.temp = merge(u.o,o.a,by.x="names.traits",by.y="species",all.x=TRUE)
	merge.temp[merge.temp$CH3Div.scheme%in%"terrestrial","Ter"] = 1
	merge.temp[merge.temp$CH3Div.scheme%in%"aquatic","Aqu"] = 1
	#
	out.path = paste0(sub("[^/]*$","",x),"non_zero_version/",gsub(".*/","",x))
	write.table(merge.temp[,names(u.o)],out.path,row.names=FALSE)
	return(u.o)
})
names(Uf) = sapply(Uf,function(x) x$Group[1])

# Unique scenario declinations
f.taxa.config = names(Uf)


### ==================================================================
### Setting up all Z5 feature files
### ==================================================================


# All unique scenario paths
all.f = levels(interaction("../outputs/zonation",f.taxa.config,sep="/"))
all.f = normalizePath(all.f)

# N-SDMs path
nsdm.f = list.files("./SDMapCHv1_10km",full.names=TRUE)
nsdm.tif = normalizePath(nsdm.f)

# Loop over taxa
for (j in 1:length(f.taxa.config))
{
	# Subset of nsdm per taxa + weights combine
	nsdm.ptaxa = sapply(paste0(Uf[[j]]$names.maps,"_"),function(x) nsdm.tif[grepl(x,nsdm.tif)])
	Ui.ptaxa = Uf[[j]][,c("names.maps","SummedUniqueness")]
	print(all(Ui.ptaxa$names.maps == gsub("_$", "", names(nsdm.ptaxa))))

	# Setting up the z5 features file
	features.f = as.data.frame(matrix(NA,length(Ui.ptaxa$SummedUniqueness),2))
	names(features.f) = names(features.t)[3:4]
	features.f$weight = Ui.ptaxa$SummedUniqueness
	features.f$filename = nsdm.ptaxa
	row.names(features.f) = Uf[[j]]$names.maps
		
	# Subset aquatic and terrestrial
	if (any(Uf[[j]]$Aqu%in%1)) {feat.aqu = features.f[Uf[[j]]$Aqu%in%1,]} else {feat.aqu = data.frame()}
	feat.ter = features.f[!Uf[[j]]$names.maps%in%row.names(feat.aqu),]

	# Save in target folders
	folders.tar = all.f[grepl(paste0("\\\\",f.taxa.config[j]),all.f)]
	if (nrow(feat.ter)!=0) {
		write.table(feat.ter,paste0(folders.tar,"/feature_list_ter.txt"),row.names=FALSE)
	}
	if (nrow(feat.aqu)!=0) {
		write.table(feat.aqu,paste0(folders.tar,"/feature_list_aqu.txt"),row.names=FALSE)
	}
}


### ==================================================================
### Setting up all Z5 setting files
### ==================================================================


# Setting up the setting files per scenario folder
for (i in 1:length(all.f))
{
	# Default
	setting.out = setting.f

	# Feature file
 	setting.out[1] = "feature list file = feature_list.txt"

	# Mask layer
	setting.out[2] = paste("analysis area mask layer =",Z5.mask)

	# Save file (ter & aqua by default)
	hab = c("_ter","_aqu")
	lapply(1:length(hab),function(x) {
		setting.out[1] = paste0("feature list file = feature_list",hab[x],".txt")
		if (any(grepl(hab[x],list.files(all.f[i])))){
			writeLines(setting.out,paste0(all.f[i],"/settings",hab[x],".z5"))
		}
	})
}


### ==================================================================
### Run zonation in parallel for z5 is found
### ==================================================================


# Zonation 5 (parallel x3)
registerDoMC(cores=3)
foreach (i=1:length(all.f),.packages=c("doMC","foreach")) %dopar%
{
	# Cat
	cat(all.f[i],"\n")

	# Next
	#if (all(c("out_Aqu","out_Ter")%in%list.files(all.f[i]))) {next}

	# Scan files
	zfiles = list.files(all.f[i],full.name=TRUE)
	setting.f = zfiles[grepl("settings_",zfiles)]

	# Label -wha
	z.param = "-wa"

	# Run
	lab = str_to_title(gsub(".*_|.z5","",setting.f))
	for (j in 1:length(setting.f)) {

		# Run zonation
		out.f = paste0(all.f[i],paste0("/out_",lab[j]))
		system(command = paste("../zonation5 --mode=CAZ2",paste(z.param,setting.f[j]),out.f))
	}
}