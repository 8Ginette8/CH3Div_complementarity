# Small function to convert a data.frame using three columns to inform on sp. traits (sp,traitID,value) to
# a data.frem per unique species names and unique trait columns
unique_to_df = function(x,na.rm=TRUE){

	# Create a new matrix comprising unique sp. name + traits
	out.df = as.data.frame(matrix(NA,nrow=length(unique(x[,1])),ncol=length(unique(x[,2]))+1))
	out.df[,1] = unique(x[,1])
	names(out.df) = c("species",unique(x[,2]))

	# Fill out the columns with loops
	for (i in 1:nrow(out.df))
	{
		# Extract and summarize traits
		dat.tar = x[x[,1]%in%out.df$species[i],]
		traits = unique(dat.tar$trait_name)

		# Loop over to calculate for each traits the mean / modal value
		for (j in 1:length(traits))
		{
			# Extract and calculate mean or mode
			dat.tar2 = dat.tar[dat.tar$trait_name%in%traits[j],"value"]
			dat.test = gsub("[0-9]|\\.","",dat.tar2)
			if (all(dat.test=="")){
				out.df[i,traits[j]] = mean(as.numeric(dat.tar2),na.rm=na.rm)
			} else {
				out.df[i,traits[j]] = Mode(dat.tar2)
			}
		}
	}
	return(out.df)
}