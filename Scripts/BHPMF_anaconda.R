#### Global path ####
#setwd("/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/Stages/Stage_M2_Mongolie/R_stats") 
#setwd("/media/lucas.dugerdil/Samsung_T5/Documents/Recherche/Stages/Stage_M2_Mongolie/R_stats") 
setwd("/home/lucas.dugerdil/Documents/Recherche/R_stats") 
#setwd("/media/lucas.dugerdil/Samsung_T5/Documents/Recherche/R_stats") 
# DB.path = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/Data_Bases/"
DB.path = "/home/lucas.dugerdil/Documents/Recherche/Data_Bases/"

# 
library("BHPMF") # en cas de probleme réinstaller le package library(devtools) install_github("fisw10/BHPMF")

Gap.filling.Uz = T
if(Gap.filling.Uz == T){
  MT.Uz <- readRDS("Resultats/Uzbekistan/Traits/Traits_values/MT_TUSD_scaled.Rds")
  TBT <- data.frame(read.csv(file="Import/Uzbekistan/Traits/Trait_ID_names.csv",sep=",",dec=".", header=T))
  MT.Uz <- MT.Uz[which(rowSums(MT.Uz[names(MT.Uz)[grep("TRY", names(MT.Uz))]], na.rm = T) != 0),]
  
  Only.trait <- as.matrix(MT.Uz[names(MT.Uz)[grep("TRY", names(MT.Uz))]])
  # Only.trait <- as.matrix(Only.trait[names(Only.trait) %in% paste("TRY_", TBT$TRY_lab[TBT$Trait_continu == T], sep = "")])
  Only.hiera <- MT.Uz[names(MT.Uz)[!grepl("TRY", names(MT.Uz))]]
  Only.hiera <- subset(Only.hiera, select = -c(GrowthForm, Other.clade))
  Only.hiera <- Only.hiera[rev(names(Only.hiera))]
  Only.hiera <- data.frame(lapply(Only.hiera, as.factor))
  Only.hiera <- cbind(plant_id = row.names(MT.Uz), Only.hiera)

  new.folder <- "/home/lucas.dugerdil/Documents/Recherche/R_stats/Resultats/Uzbekistan/Traits/Traits_values/"
  saveRDS(Only.hiera, paste(new.folder, "Hierarchie_gap-filling.Rds", sep = ""))
  saveRDS(Only.trait, paste(new.folder, "TM_before_gap-filling.Rds", sep = ""))
  
  #### Gapfilling model ####
  GapFilling(Only.trait, Only.hiera, verbose = T,
             mean.gap.filled.output.path = "/tmp/MT_TUSD_scale_gf.txt",
             std.gap.filled.output.path = "/tmp/MT_TUSD_scale_gf_sd.txt")
  
  list.of.files <- list.files("/tmp", "*txt", full.names = T)
  file.copy(list.of.files, new.folder, overwrite = T)
  
  MT.ACA.sd <- read.table(file = "Resultats/Uzbekistan/Traits/Traits_values/MT_TUSD_scale_gf.txt", sep="\t",dec=".", header=T, stringsAsFactors = F)
  MT.ACA.gf <- read.table(file = "Resultats/Traits/Traits_values/MT_TUSD_gf_scale_sd.txt", sep = "\t", dec=".", header=T)
  }
  


