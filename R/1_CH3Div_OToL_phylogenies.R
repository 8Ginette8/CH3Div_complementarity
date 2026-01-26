# ###################################################################################
# CH3Div: Open Tree of Life to generate group phylogenies
#
# $Date: 2024-07-23
#
# Author: Yohann Chauvier, yohann.chauvier@wsl.ch
# Aquatic Ecology Group
# Swiss Federal Research Institute EAWAG
# 
# Description: Here we want to use the Open Tree of Life package to generate for
# each taxo group a phylogeny. We could have done it with more elaborated phylogenies,
# but given the number of taxonomic group we have, we cannot ensure such process for
# every groups. So simpler to apply OToL for everything...
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
library(rotl)
library(ape)
library(stringr)
library(httr)
library(phytools)

# Setting working directory
setwd("./data/")

# Functions
scr = list.files("../R/functions",full.names=TRUE)
invisible(lapply(scr, source))


### ==================================================================
### Process priliminary files
### ==================================================================


# EU species
EU.sp = c("EEA","FaunaEuropeae","EIVE","AraneaeEuropeae",
	"Loncarevic_et_al.2024","Middleton-Welling_et_al.2020")

# Open all s2z (swiss) trait files (contain all species names)
s2z.l = list.files("./x_TraitCH_1.0/raw_traits")
s2z.txt = lapply(s2z.l,function(x){
	oo = read.table(paste0("./x_TraitCH_1.0/raw_traits/",x),header=TRUE)
	oo = oo[!oo$Source%in%EU.sp,]
	return(oo)
})


### ==================================================================
### Get the phylogenies
### ==================================================================


# Deactivate secure connection (ROTL bugg) (FORCE when does not work)
#set_config(config(ssl_verifypeer = 0L))

# Prepare group IDs and OToL groups
tnrs_contexts()
g.IDs = gsub("s2z_raw_traits_|\\.txt","",s2z.l)
g.otol = c("Amphibians","Birds","Insects","Land plants","Insects","Insects","Vertebrates","Fungi","Insects",
	"Fungi","Mammals","Molluscs","Insects","Insects","Land plants","Tetrapods","Arachnids")

# Loop over files to extract species names to generate the phylogenies
for (i in 1:length(g.IDs))
{
	cat ("Generating '",g.IDs[i],"' ( n =",i,")","phylogeny...","\n")

	# Next
	#if (any(grepl(g.IDs[i],list.files("./x_OpenTreeOfLife/")))) {next}

	# Species + rm potential "subsp." to enlarge the match
	taxa = s2z.txt[[i]][,c("Species","GBIF_rank")]
	taxa.sp = taxa[!taxa$GBIF_rank%in%c("SUBSPECIES","VARIETY"),]
	taxa.subsp = taxa[taxa$GBIF_rank%in%c("SUBSPECIES","VARIETY"),]
	#
	taxa.subsp1 = taxa.subsp
	taxa.subsp1[,"Species"] = gsub("subsp\\. |var.\\ ","",taxa.subsp1[,"Species"])
	#
	taxa.subsp2 = taxa.subsp1
	c.subsp = taxa.subsp2$GBIF_rank%in%"SUBSPECIES" & !grepl(" x ",taxa.subsp2$Species)
	c.var = taxa.subsp2$GBIF_rank%in%"VARIETY" & !grepl(" x ",taxa.subsp2$GBIF_rank)
	taxa.subsp2$Species[c.subsp] = unlist(sapply(strsplit(taxa.subsp2$Species[c.subsp]," "),
		function(x) paste(x[1],x[2],"subsp.",x[3])))
	taxa.subsp2$Species[c.var] = unlist(sapply(strsplit(taxa.subsp2$Species[c.var]," "),
		function(x) paste(x[1],x[2],"var.",x[3])))
	#
	taxa = rbind(taxa.sp,taxa.subsp1,taxa.subsp2)
	taxa = taxa[!duplicated(taxa$Species),]

	# Check the names in the database (some won't be found --> e.g., aggregates or hybrids)
	resolved.n = tnrs_match_names(taxa$Species,context_name=g.otol[i])
	resolved.n$S2Z.rank = taxa$GBIF_rank
	resolved.n = resolved.n[!is.na(resolved.n$unique_name),]
	resolved.n$search_string = str_to_sentence(resolved.n$search_string)

	# Enlarge results for better matching with S2Z
	r.sp.k = resolved.n[!resolved.n$S2Z.rank%in%c("SUBSPECIES","VARIETY"),]
	r.subsp.k = resolved.n[resolved.n$S2Z.rank%in%c("SUBSPECIES","VARIETY"),]
	#
	r.subsp1 = r.subsp.k
	r.subsp1[,"search_string"] = gsub("subsp\\. |var.\\ ","",r.subsp1[,"search_string"])
	r.subsp1.k = r.subsp1[!r.subsp1$search_string%in%r.subsp.k$search_string,]
	#
	r.subsp2 = r.subsp1
	c.subsp = r.subsp2$S2Z.rank%in%"SUBSPECIES"
	c.var = r.subsp2$S2Z.rank%in%"VARIETY"
	r.subsp2$search_string[c.subsp] = unlist(sapply(strsplit(r.subsp2$search_string[c.subsp]," "),
		function(x) paste(x[1],x[2],"subsp.",x[3])))
	r.subsp2$search_string[c.var] = unlist(sapply(strsplit(r.subsp2$search_string[c.var]," "),
		function(x) paste(x[1],x[2],"var.",x[3])))
	r.subsp2.k = r.subsp2[!r.subsp2$search_string%in%r.subsp.k$search_string,]
	#
	resolved.n = rbind(r.sp.k,r.subsp.k,r.subsp1.k,r.subsp2.k)

	# Remove duplicated that are not S2Z or if no S2Z keep relevant names to keep the correspondance
	no.dup = resolved.n[!all_duplicated(resolved.n$unique_name),]
	all.dup = resolved.n[all_duplicated(resolved.n$unique_name),]
	u.dup = unique(all.dup$unique_name)
	if (!nrow(all.dup)%in%0)
	{
		s2z.keep = data.frame(DUP=u.dup,S2Z=NA,S2Z.rank=NA)
		for (j in 1:nrow(s2z.keep))
		{
			# Extract s2z relevant names
			tar.dup = all.dup[all.dup$unique_name%in%u.dup[j],]
			o.s2z = s2z.txt[[i]][s2z.txt[[i]]$Species%in%tar.dup$search_string,]
			# Case where the related names have already been captured in "no_dup" (synonym anomaly)
			if (nrow(o.s2z)==0) {next} 
			ms2z = merge(tar.dup,o.s2z[,c("Species","Source","GBIF_rank")],by.x="search_string",by.y="Species")
			if (any(ms2z$Source%in%"SPEED2ZERO")) {
				s2z.k = ms2z[ms2z$Source%in%"SPEED2ZERO",]
			} else {
				if (any(ms2z$is_synonym%in%FALSE))
				{
					s2z.k = ms2z[ms2z$is_synonym%in%FALSE,]
				}
				else if (identical(s2z.k,character(0))) {
					s2z.k = ms2z[ms2z$GBIF_rank%in%"SPECIES",]
				} else {
					s2z.k = ms2z
				}
			}

			# Just keep the highest score in case we obtain several
			s2z.f = s2z.k[which.max(s2z.k$score),c("search_string","S2Z.rank")]
			s2z.keep[j,c("S2Z","S2Z.rank")] = s2z.f
		}

		# Remove duplicated unique_name + anomalies + assign correct search_string
		all.dup2 = all.dup[!duplicated(all.dup$unique_name),]
		all.dup2$S2Z.rank = NULL
		s2z.keep = s2z.keep[!is.na(s2z.keep$S2Z),]
		mdup = merge(all.dup2,s2z.keep,by.x="unique_name",by.y="DUP")
		mdup$search_string = mdup$S2Z
		mdup = mdup[,names(resolved.n)]
		no.dup = rbind(no.dup,mdup)
	} 

	# Check potential errors in the OTT id
	ott.id = no.dup$ott_id
	isT = is_in_tree(ott.id)
	qfound = lapply(ott.id[!isT],function(x) try(tol_node_info(x)))
	it.found = sapply(qfound,function(x) !grepl("Error",x))
	toTree = no.dup[isT,]

	# Correct some anomalies
	#if (any(it.found%in%TRUE)){
	#	stop("Check manually!")
	#}

	# Generate tree and find missing names
	my.tree = tol_induced_subtree(ott_ids=toTree$ott_id)
	tree.label = gsub("_ott.*","",my.tree$tip.label)
	tree.label = gsub("_"," ",tree.label)
	q.missing = !no.dup$unique_name%in%tree.label

	# If there are, we need to add missing tips to the "Most Recent Common Ancestor" of the closest genus
	if (any(q.missing))
	{
		sp.missing = paste0(gsub(" ","_",no.dup$unique_name[q.missing]),"_ott",no.dup$ott_id[q.missing])
		tar.genus = sapply(sp.missing, function(x) strsplit(x,"_")[[1]][1])
		for (j in 1:length(tar.genus))
		{
			genus.tips = my.tree$tip.label[grep(paste0(tar.genus[j],"_"),my.tree$tip.label)]
			if (length(genus.tips)>0)
			{
				mrca.node = getMRCA(my.tree, genus.tips)
				my.tree = bind.tip(my.tree,sp.missing[j],NULL,where=mrca.node)
			} else {
				print("Genus could not be found...")
			}
		}
	}

	# Change tip labels to S2Z names
		#
	my.tree$tip.label = gsub(".*_","",my.tree$tip.label)
	my.tree$tip.label = gsub("ott","",my.tree$tip.label)
	tip.merge = data.frame(tip=my.tree$tip.label)
	f.merge = merge(tip.merge,no.dup,by.x="tip",by.y="ott_id",sort=FALSE,all.x=TRUE)
		
		# PROBABLY not useful anymore but it doesn't do any harm to keep it in case of debugging...
	tar.subsp = f.merge$S2Z.rank%in%"SUBSPECIES"
	c.subsp = sapply(f.merge$search_string,function(x) !any(strsplit(x," ")[[1]]%in%"subsp."))
	splsubsp = lapply(f.merge$search_string,function(x) strsplit(x," ")[[1]])
	if (any(tar.subsp)){
		f.merge[tar.subsp & c.subsp,"search_string"] = sapply(splsubsp[tar.subsp & c.subsp],
			function(x) paste(x[1],x[2],"subsp.",x[3]))
	}
		#
	tar.var = f.merge$S2Z.rank%in%"VARIETY"
	c.var = sapply(f.merge$search_string,function(x) !any(strsplit(x," ")[[1]]%in%"var."))
	splvar = lapply(f.merge$search_string,function(x) strsplit(x," ")[[1]])
	if (any(tar.var)){
		f.merge[tar.var & c.var,"search_string"] = sapply(splvar[tar.var & c.var],
			function(x) paste(x[1],x[2],"var.",x[3]))
	}
		#
	tar.aggr = f.merge$S2Z.rank%in%"AGGREGATE"
	c.aggr = sapply(f.merge$search_string,function(x) !any(strsplit(x," ")[[1]]%in%"aggr."))
	splaggr = lapply(f.merge$search_string,function(x) strsplit(x," ")[[1]])
	if (any(tar.aggr)){
		f.merge[tar.aggr & c.aggr,"search_string"] = sapply(splaggr[tar.aggr & c.aggr],
			function(x) paste(x[1],x[2],"aggr.",x[3]))
	}

	# Save name changes
	print(length(my.tree$tip.label) == length(f.merge$search_string))
	my.tree$tip.label = f.merge$search_string

	# Reformtaing text for bryophytes...(some odd name writing in all our databases)
	if (g.IDs[i]%in%"bryophytes"){
		my.tree$tip.label[grepl("Sciuro-hypnum",my.tree$tip.label)] =
			gsub("h","H",my.tree$tip.label[grepl("Sciuro-hypnum",my.tree$tip.label)])
	}

	# Generate branch length (Grafen's (1989) computation of branch lengths)
	# Balances power of 1 by default
	my.tree.brlen = compute.brlen(my.tree,power=1)

	# Save tree
	print(all(my.tree.brlen$tip%in%s2z.txt[[i]]$Species)) # ToCheck if all S2Z names are in the tree
	saveRDS(my.tree.brlen,file=paste0("./x_OpenTreeOfLife/sphylotreeCH_",g.IDs[i]))
}

# Reactivate secure connection
#set_config(config(ssl_verifypeer = 1L))