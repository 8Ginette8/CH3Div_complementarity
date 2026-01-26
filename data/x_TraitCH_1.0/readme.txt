@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ TraitsCH
@
@ Date: 2025/03/18
@
@ Author: Yohann Chauvier, yohann.chauvier@wsl.ch
@ Aquatic Ecology Group
@ Swiss Federal Research Institute EAWAG
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


1) Folder 'raw_traits' --> compiled trait tables
	Format: Space-delimited text file: 17 taxonomic groups, non-imputed)

2) Folder 'missRanger_imputed_traits' --> imputed trait tables
	Format: Space-delimited text file: 17 groups * (25 imputation reps + average 'MEAN')

3) Foler 'missRanger_evaluations' --> detailed evaluation of the imputation models
	Format: Traits * (25 imputation reps + average ('OBB_pred_error') + Data coverage ('data_completeness_prct')

4) 'Delarze_bryo_class_desc.csv' --> detailed class descritpion for Delarze habitats (bryophytes)

5) 'Flora-indicativa_plant_class_desc.csv' --> detailed class description for FloraIndicativa habitats (vascular plants)

6) 'metadata.xlsx' --> detailed description of all compiled traits and references for 17 taxonomic groups

7) 'InfoGuildes.csv' --> detailed description of Swiss Habitat Guilds ('GUILDE.x')

8) 'to_GenusSpecies.R' --> Function to remove/correct scientific names duplicates, incomplete entries, and unknown names