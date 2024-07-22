#### Global path ####
#setwd("/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/Stages/Stage_M2_Mongolie/R_stats") 
#setwd("/media/lucas.dugerdil/Samsung_T5/Documents/Recherche/Stages/Stage_M2_Mongolie/R_stats") 
setwd("/home/lucas.dugerdil/Documents/Recherche/R_stats") 
#setwd("/media/lucas.dugerdil/Samsung_T5/Documents/Recherche/R_stats") 
# DB.path = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/Data_Bases/"
DB.path = "/home/lucas.dugerdil/Documents/Recherche/Data_Bases/"
source("Scripts/TNRS.R")

#### Libraries ####
# suppressPackageStartupMessages()
library(ggplot2)
library(ggpp)
library(ggpubr) # stat_cor
library(ggpmisc) # stat_poly_line
library(RColorBrewer)
library(patchwork)
library(dplyr) # mutate
library(maps)
library(mapproj)
library(gstat)
library(reshape)
library(ggnewscale) # function new_scale_color
library("readr")
library(ggthemes)
library(tibble) # as_tibble
library(scales)       # quoi du type d'echelle sur les graphs ggplot
library(lubridate)
library(rgdal)

# usethis::use_git_config(user.name="Lucas Dugerdil", user.email="lucas.dugerdil@lebib.org")
# library(usethis)

#### Functions ####
Clean.trait.discon <- function(MD, Trait.ID, Keep.string, No.aggregation = F){
  if(missing(Keep.string)){Keep.string <- NULL}
  MD <- MD[MD$TraitID %in% Trait.ID, c("AccSpeciesName", "TraitID", "View()")]
  names(MD)[1] <- c("species")
  MD$ID <- seq(1:nrow(MD))
  #MD$TraitID <- as.numeric(MC$TraitID)
  MD <- as_tibble(reshape2::dcast(MD, species + ID ~ TraitID, value.var = "OrigValueStr"))
  
  #### Photosynthesis pathway ####
  if(is.null(MD[["22"]]) == F){
    MD[["22"]] <- gsub("c3", "C3", MD[["22"]])
    MD[["22"]] <- gsub("C3\\?", "C3", MD[["22"]])
    MD[["22"]] <- gsub("C4\\?", "C4", MD[["22"]])
    MD[["22"]] <- gsub("c4", "C4", MD[["22"]])
    
    MD[["22"]] <- gsub("C3/C4", 0.5, MD[["22"]])
    MD[["22"]] <- gsub("C3/CAM", -0.5, MD[["22"]])
    MD[["22"]] <- gsub("C3", 0, MD[["22"]])
    MD[["22"]] <- gsub("C4", 1, MD[["22"]])
    MD[["22"]] <- gsub("CAM", -1, MD[["22"]])
    MD[["22"]] <- as.numeric(MD[["22"]])
    }
  
  #### Evergreenness ####
  if(is.null(MD[["37"]]) == F){
    MD[["37"]][which(MD[["37"]] == "Evergreen")] <- "evergreen"
    MD[["37"]][which(MD[["37"]] == "exchanger")] <- "evergreen"
    MD[["37"]][which(MD[["37"]] == "yes")] <- "evergreen"
    MD[["37"]][which(MD[["37"]] == "3")] <- "evergreen"
    MD[["37"]][which(MD[["37"]] == "?")] <- NA
    MD[["37"]][which(MD[["37"]] == "n.d.")] <- NA
    MD[["37"]][which(MD[["37"]] == "aphyllous")] <- NA
    MD[["37"]][which(MD[["37"]] == "E")] <- "evergreen"
    MD[["37"]][which(MD[["37"]] == "evergreen type 2")] <- "evergreen"
    MD[["37"]][which(MD[["37"]] == "evergreen  type 1")] <- "evergreen"
    MD[["37"]][which(MD[["37"]] == "always persistent green")] <- "evergreen"
    
    MD[["37"]][which(MD[["37"]] == "Deciduous")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "aestival")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "vernal")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "hibernal")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "aestival")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "no")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "D")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "Db")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "W")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "winter deciduous")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "deciduous type 1")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "deciduous type 2")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "deciduous type 3")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "5")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "1")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "megaphanerophyte")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "Nonevergreen")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "always summer green")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "always overwintering green")] <- "deciduous"
    MD[["37"]][which(MD[["37"]] == "always spring green")] <- "deciduous"
    
    MD[["37"]][which(MD[["37"]] == "2")] <- "semi-evergreen"
    MD[["37"]][which(MD[["37"]] == "SEMI")] <- "semi-evergreen"
    MD[["37"]][which(MD[["37"]] == "ED")] <- "semi-evergreen"
    MD[["37"]][which(MD[["37"]] == "deciduous/evergreen")] <- "semi-evergreen"
    MD[["37"]][which(MD[["37"]] == "semi-deciduous")] <- "semi-evergreen"
    MD[["37"]][which(MD[["37"]] == "semideciduous")] <- "semi-evergreen"
    MD[["37"]][which(MD[["37"]] == "winter semi-deciduous")] <- "semi-evergreen"
    MD[["37"]][which(MD[["37"]] == "drought semi-deciduous")] <- "semi-evergreen"
    MD[["37"]][which(MD[["37"]] == "drought-deciduous")] <- "semi-evergreen"
    MD[["37"]][which(MD[["37"]] == "variable")] <- 0.5
    MD[["37"]][which(MD[["37"]] == "semi-evergreen")] <- 0.5
    MD[["37"]][which(MD[["37"]] == "deciduous")] <- 0
    MD[["37"]][which(MD[["37"]] == "evergreen")] <- 1
    MD[["37"]] <- as.numeric(MD[["37"]])}
  
  #### Woodiness ####
  if(is.null(MD[["38"]]) == F){
    MD[["38"]][which(MD[["38"]] == "Semi-woody")] <- 0.5
    MD[["38"]][which(MD[["38"]] == "semi-woody")] <- 0.5
    MD[["38"]][which(MD[["38"]] == "wood at base")] <- 0.5
    MD[["38"]][which(MD[["38"]] == "woody at base")] <- 0.5
    MD[["38"]][which(MD[["38"]] == "Suffrutex")] <- 0.5
    MD[["38"]][which(MD[["38"]] == "non-woody/woody")] <- 0.5
    MD[["38"]][which(MD[["38"]] == "Variable")] <- 0.5
    MD[["38"]][which(MD[["38"]] == "variable")] <- 0.5
    MD[["38"]][which(MD[["38"]] == "Herb")] <- 0
    MD[["38"]][which(MD[["38"]] == "h")] <- 0
    MD[["38"]][which(MD[["38"]] == "H")] <- 0
    MD[["38"]][which(MD[["38"]] == "Herbaceous")] <- 0
    MD[["38"]][which(MD[["38"]] == "Grass&Sedges")] <- 0
    MD[["38"]][which(MD[["38"]] == "non-woody")] <- 0
    MD[["38"]][which(MD[["38"]] == "non woody")] <- 0
    MD[["38"]][which(MD[["38"]] == "Woody")] <- 1
    MD[["38"]][which(MD[["38"]] == "W")] <- 1
    MD[["38"]][which(MD[["38"]] == "Y")] <- 1
    MD[["38"]][which(MD[["38"]] == "w")] <- 1
    MD[["38"]][which(MD[["38"]] == "3")] <- 1
    MD[["38"]][which(MD[["38"]] == "2")] <- 1
    MD[["38"]][which(MD[["38"]] == "woody")] <- 1
    MD[["38"]] <- as.numeric(MD[["38"]])
    }
  
  #### Growth forme ####
  if(is.null(MD[["48"]]) == F){
    MD[["48"]][which(MD[["48"]] == "Shrub")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == " shrub")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "shrug")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "shrub*")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "shrub*****")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "twining shrub")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "twiner")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "twiner*")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "woody*")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "woody")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "Woody")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "succulent shrub")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "climbing shrub")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "climbing shrub*")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "semi-shrub")] <- "subshrub"
    MD[["48"]][which(MD[["48"]] == "sub-shrub")] <- "subshrub"
    MD[["48"]][which(MD[["48"]] == "sub shrub")] <- "subshrub"
    MD[["48"]][which(MD[["48"]] == "sub shrub*")] <- "subshrub"
    MD[["48"]][which(MD[["48"]] == "woody herb")] <- "subshrub"
    MD[["48"]][which(MD[["48"]] == "cushion grass")] <- "subshrub"
    MD[["48"]][which(MD[["48"]] == "cushion plant")] <- "subshrub"
    MD[["48"]][which(MD[["48"]] == "creeper")] <- "subshrub"
    MD[["48"]][which(MD[["48"]] == "creeper*")] <- "subshrub"
    MD[["48"]][which(MD[["48"]] == "woody creeper")] <- "subshrub"
    MD[["48"]][which(MD[["48"]] == "forb*")] <- "forb"
    MD[["48"]][which(MD[["48"]] == "crucifer")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "herb*")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "Herb")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "Herb, Pere")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "per grass")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "climbing herb*")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "climbing herb")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "twining herb*")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "twining herb")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "herbaceous twiner")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "herbaceous climber*")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "tussock grass")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "Grass")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "bunchgrass")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "moss")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "rosette")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "succulent*")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "leaf succulent")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "succulent")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "succulent herb")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "neophyte")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "neophyte*")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "halophyte")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "epiphyte")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "Non-woody epiphyte")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "non-woody epiphyte")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "geophyte")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "geophyte")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "geophyte*")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "hemicryptophyte")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "therophyte")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "cyperoid")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "cyperoid*")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "grass")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "grass*")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "Graminoid")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "bamboo")] <- "tree"
    MD[["48"]][which(MD[["48"]] == "shrubby bamboo")] <- "tree"
    MD[["48"]][which(MD[["48"]] == "Tree")] <- "tree"
    MD[["48"]][which(MD[["48"]] == "small tree")] <- "tree"
    MD[["48"]][which(MD[["48"]] == "Small_Tree")] <- "tree"
    MD[["48"]][which(MD[["48"]] == "tree ")] <- "tree"
    MD[["48"]][which(MD[["48"]] == "treelet*")] <- "tree"
    MD[["48"]][which(MD[["48"]] == "tree*")] <- "tree"
    MD[["48"]][which(MD[["48"]] == "tree**** (means find other one in original, add *)")] <- "tree"
    MD[["48"]][which(MD[["48"]] == "woody vine")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "woody vine*")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "Liana")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "liana")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "liana*")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "climber*")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "Vine")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "vine*")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "vine")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "climbing legume")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "herbaceous vine")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "herbaceous climber")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "hook climber")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "leaning climber")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "tendril climber")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "shrubby climber")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "woody climber")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "woody climber*")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "hemiepiphyte")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "epiphyte")] <- "climber"
    MD[["48"]][which(MD[["48"]] == "epiphyte*")] <- "climber"
    
    MD[["48"]][which(MD[["48"]] == "hemiparasite")] <- "parasite"
    MD[["48"]][which(MD[["48"]] == "hemiparasite*")] <- "parasite"
    MD[["48"]][which(MD[["48"]] == "Parasite")] <- "parasite"
    MD[["48"]][which(MD[["48"]] == "parasite*")] <- "parasite"
    MD[["48"]][which(MD[["48"]] == "Sedge")] <- "sedge"
    MD[["48"]][which(MD[["48"]] == "Aquatic")] <- "aquatic"
    MD[["48"]][which(MD[["48"]] == "aquatic sedge")] <- "aquatic"
    MD[["48"]][which(MD[["48"]] == "aquatic sedge*")] <- "aquatic"
    MD[["48"]][which(MD[["48"]] == "aquatic herb")] <- "aquatic"
    MD[["48"]][which(MD[["48"]] == "aquatic herb*")] <- "aquatic"
    MD[["48"]][which(MD[["48"]] == "aquatic*")] <- "aquatic"
    MD[["48"]][which(MD[["48"]] == "hydrophyte")] <- "aquatic"
    MD[["48"]][which(MD[["48"]] == "hydrophyte")] <- "aquatic"
    MD[["48"]][which(MD[["48"]] == "marsh")] <- "aquatic"
    
    
    MD[["48"]][which(MD[["48"]] == "1")] <- NA
    MD[["48"]][which(MD[["48"]] == "10")] <- NA
    MD[["48"]][which(MD[["48"]] == "2")] <- NA
    MD[["48"]][which(MD[["48"]] == "3")] <- NA
    MD[["48"]][which(MD[["48"]] == "4")] <- NA
    MD[["48"]][which(MD[["48"]] == "5")] <- NA
    MD[["48"]][which(MD[["48"]] == "Reed")] <- NA
    MD[["48"]][which(MD[["48"]] == "Pam")] <- NA
    MD[["48"]][which(MD[["48"]] == "scrambler")] <- NA
    MD[["48"]][which(MD[["48"]] == "scrambling herb")] <- NA
    
    MD[["48"]][which(MD[["48"]] == "subshrub")] <- "shrub"
    MD[["48"]][which(MD[["48"]] == "forb")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "sedge")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "graminoid")] <- "herb"
    MD[["48"]][which(MD[["48"]] == "parasite")] <- "Other"
    MD[["48"]][which(MD[["48"]] == "climber")] <- "Other"
    MD[["48"]][which(MD[["48"]] == "aquatic")] <- "Other"
    MD[["48"]][which(MD[["48"]] == "other")] <- "Other"
    MD[["48"]][which(MD[["48"]] == "herb")] <- "Herb"
    MD[["48"]][which(MD[["48"]] == "tree")] <- "Tree"
    MD[["48"]][which(MD[["48"]] == "shrub")] <- "Shrub"
    }

  #### Select unique value of GF ####  
  if(is.null(Keep.string) == F & any(names(MD) %in% as.character(Keep.string)) == T){
    MD.str <- MD[names(MD) %in% c("species", as.character(Keep.string))]
    names(MD.str)[2] <- "To.agg"
    library(plyr)
    MD.str <- plyr::count(MD.str, c("species", "To.agg"))
    MD.str <- MD.str[is.na(MD.str$To.agg)==F,]
    Duplicated.species <- MD.str$species[duplicated(MD.str$species)]
    # Selected.type <- MD.str[0]
    for(i in Duplicated.species){
      Sub.df <- MD.str[MD.str$species == i,]
      
      # print(Sub.df$To.agg[Sub.df$freq == max(Sub.df$freq)][1])
      MD.str$To.agg[MD.str$species == i] <- Sub.df$To.agg[Sub.df$freq == max(Sub.df$freq)][1]
      # print(MD.str$To.agg[MD.str$species == i])
    }
    MD.str <- MD.str[-c(3)]
    MD.str <- MD.str[!duplicated(MD.str),]
    names(MD.str)[2] <- as.character(Keep.string)
    Go = T
    }
  else{Go = F}
    
  #### Merge MD/MC ####
  MD <- MD[-2]
  MD <- MD[!names(MD) %in% as.character(Keep.string)]
  if(No.aggregation == F){
    MD <- aggregate(MD, list(MD[["species"]]), FUN = mean, na.action = na.omit)
    MD <- MD[-2]
    names(MD)[1] <- c("species")
    }
  
  if(is.null(Keep.string) == F & Go == T){
    MD <- full_join(MD, MD.str, by = "species")
  }
  
  return(MD)
}

Clean.trait.con <- function(MC, Trait.ID, Average){
  if(missing(Average)){Average = F}
  MC <- MC[MC$TraitID %in% Trait.ID, c("AccSpeciesName", "TraitID", "StdValue")]#, "UnitName", "Comment")]
  names(MC)[1] <- c("species")
  MC <- MC[which(!is.na(MC$species)),]
  
  # PB <- unique(MC[c(2,4,5)])
  # PB <- PB[sort(PB$TraitID),]
  # View(PB)
  
  MC$StdValue <- as.numeric(MC$StdValue)
  MC$TraitID <- as.numeric(MC$TraitID)
  MC <- aggregate(MC, list(MC$species, MC$TraitID), FUN = mean)
  MC <- MC[-c(2,3)]
  names(MC)[c(1,2)] <- c("species", "TraitID")
  MC <- reshape2::dcast(MC[-4], species ~ TraitID, na.rm = F)
  
  #### Valeur moyenne des traits équivalents ####
  if(Average == T){
    MC[["3108"]] <- rowMeans(MC[as.character(3108:3114)], na.rm = T)
    MC <- MC[-match(as.character(3109:3114), names(MC))]  
    
    MC[["3115"]] <- rowMeans(MC[as.character(3115:3117)], na.rm = T)
    MC <- MC[-match(as.character(3116:3117), names(MC))]  
    # print(MC)
    }
  
  return(MC)
}

PCA.bioclim <- function(MP, transp_OK, Site.name, Type.samples, Ellipse, Shape, Show.centroid, Show.arrow,
                        Csv.sep, Scale.PCA, Groupes, Cluster.core, Cluster.core.lab, return.pick, Contrib,
                        Save.path, Manu.lim.x, Manu.lim.y, Dot.size, Dot.opac, Ellipse.opa, Density.contour, 
                        Opa.range, Reverse.dim, Show.annot, Show.Plotly, PCA.site, Marg.density.plot,
                        Symbol.path = NULL, Symbol.pos = NULL, Legend.position, Density.type, Save.plot, H, W, Num.facet, Legend.size){
  #### Settings ####
  library(vegan)
  library(ggrepel)
  library("FactoMineR") # FactoMineR pour la PCA (Le et al. 2008)
  library("factoextra")
  if(missing(Csv.sep)){Csv.sep = "\t"}
  if(missing(Site.name)){Site.name = "Site1"}
  if(missing(Cluster.core)){Cluster.core = NULL}
  if(missing(Num.facet)){Num.facet = NULL}
  if(missing(Ellipse.opa)){Ellipse.opa = 0.4}
  if(missing(Shape)){Shape = 21}
  if(missing(Cluster.core.lab)){Cluster.core.lab = "Biome (Dinerstein et al., 2017)"}
  if(missing(PCA.site)){PCA.site = F}
  if(missing(Marg.density.plot)){Marg.density.plot = F}
  if(missing(Reverse.dim)){Reverse.dim = F}
  if(missing(Show.arrow)){Show.arrow = T}
  if(missing(Show.annot)){Show.annot = T}
  if(missing(Ellipse)){Ellipse = F}
  if(missing(Contrib)){Contrib = T}
  if(missing(return.pick)){return.pick = F}
  if(missing(Show.centroid)){Show.centroid = F}
  if(missing(transp_OK)){transp_OK = T}
  if(missing(Legend.position)){Legend.position = "right"}
  if(Show.centroid == F){Centroide = "quali"}
  if(Show.centroid == T){Centroide = NULL}
  if(missing(Legend.size)){Legend.size = 1}
  if(missing(Scale.PCA)){Scale.PCA = 1}
  if(missing(Save.path)){Save.path = NULL}
  if(missing(Groupes)){Groupes = NULL}
  if(missing(Manu.lim.x)){Manu.lim.x = NULL}
  if(missing(Manu.lim.y)){Manu.lim.y = NULL}
  if(missing(Type.samples)){Type.samples = NULL}
  if(missing(Show.Plotly)){Show.Plotly = F}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(Density.contour)){Density.contour = F}
  if(missing(Density.type)){Density.type = "polygon"}
  if(missing(Dot.size)){Dot.size = NULL}
  if(missing(Dot.opac)){Dot.opac = NULL}
  if(missing(Opa.range)){Opa.range = c(0.01,0.1)}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  # if(missing(Leg.row)){Leg.row = 1}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Data pulishing ####
  print("Lets PCA bioclimatic !")
  Remove.name <- c("Biome", "Biome.no", "Ecosystem", "Latitude", "Longitude", "Bioclim",
                  "Aridity", "Type", "AP_NAP", "PFT", "GrowthForm", "Country")
  Keep.xdata <- MP[intersect(names(MP), Remove.name)]
  MP <- MP[setdiff(names(MP), Remove.name)]
  names(MP) <- gsub("TRY_", "", names(MP))
  
  #### Groupes ####
  if(is.null(Groupes) == F){
    Groupes <- melt(Groupes)
    Groupes <- Groupes$L1[match(names(MP), Groupes$value)]
    Groupes[is.na(Groupes)] <- "Unknown"
    Groupes <- factor(Groupes)
    # print(Groupes)
  }
  # else{Groupes <- as.factor("royalblue")}
  # else{
   # }
  
  #### Transforming the data + PCA calcul ####
  if(nlevels(as.factor(is.na(MP))) >= 2){
    library(missMDA)
    nb <- estim_ncpPCA(MP, ncp.max=5) ## Time consuming, nb = 2
    MP.comp <- imputePCA(MP, ncp = nb[[1]])
    MP <- MP.comp$completeObs
    }

  MP.pca <- PCA(MP, graph = FALSE, scale.unit = transp_OK)
  PCA <- data.frame(MP.pca$ind$coord)
  PCA <- PCA[,1:2]
  names(PCA) <- c("PC1", "PC2")
  
  if(is.null(Keep.xdata$Type) == F){
    PCA <- cbind(PCA, Type = Keep.xdata$Type, Biome = Keep.xdata[[Cluster.core]])
    PCA <- PCA[PCA$Type == "MPS",]
    Point.surf <- geom_point(inherit.aes = F, PCA, mapping = aes(x = PC1, y = PC2, fill = Biome), shape = Shape, colour = "grey10", alpha = 0.5, na.rm = T)
    Opac <- 0.3
    Ptaille <- 2
  }
  else{
    Point.surf <- NULL
    Opac <- 0.7
    Ptaille <- 1.5 
    }
  if(is.null(Dot.size) == F){Ptaille <- Dot.size}
  if(is.null(Dot.opac) == F){Opac <- Dot.opac}
    
  #### Color settings ####
  if(is.null(Groupes) == F){
    my_orange = c("Water" = "darkblue",
                  "Altitude" = "brown",
                  "Temperature" = "darkred")

    Scale.color.vec <- scale_color_manual(values = my_orange, name = "Proxies")
    
    }
  else{
    Groupes <- "Unique"
    my_orange <- data.frame(Unique = "royalblue")
    Scale.color.vec <- scale_color_manual(values = my_orange, guide = "none")
    }
  
  if(is.null(Cluster.core) == F){
    if(is.numeric(Keep.xdata[[Cluster.core]]) == T){
      my_orange2 = brewer.pal(n = 11, "RdYlBu")[Keep.col2[-c(3,4,5,7,8,9)]] 
      orange_palette2 = colorRampPalette(my_orange2)
      my_orange2 = rev(orange_palette2(length(seq(min(Keep.xdata[[Cluster.core]]), max(Keep.xdata[[Cluster.core]]), by = 200))))
      Scale.fill <- scale_fill_gradientn(colours = my_orange2, guide = "colourbar", 
                                         name = Cluster.core.lab,
                                         breaks = seq(round(min(Keep.xdata[[Cluster.core]]),digits = -3), max(Keep.xdata[[Cluster.core]]), by = 1000),
                                         na.value = "white")}
    else{
      values.bi = c("Deserts & Xeric Shrublands" = "#C88282",
                    "Temperate Grasslands, Savannas & Shrublands" = "#ECED8A",
                    "Montane Grasslands & Shrublands" = "#D0C3A7",
                    "Temperate Conifer Forests" = "#6B9A88",
                    "Temperate Broadleaf & Mixed Forests" = "#3E8A70",
                    "N/A" = "#FFEAAF",
                    "Tundra" = "#A9D1C2",
                    "Boreal Forests/Taiga" = "#8FB8E6",
                    "Tropical & Subtropical Coniferous Forests" = "#99CA81",
                    "Mangroves" = "#FE01C4",
                    "Flooded Grasslands & Savannas" = "#BEE7FF",
                    "Tropical & Subtropical Moist Broadleaf Forests" = "#38A700",
                    "Plant_height" = "royalblue",
                    "Leaf_thickness" = "darkorange",
                    "Photosynthesis_pathway" = "purple",
                    "Woodyness" = "#323232",
                    "Leaf_size" = "darkred",
                    "Variable" = "#aa373aff",
                    "Algal" = "#4666E9",
                    "NAP" = "#b5ab32ff",
                    "Herb" = "#b5ab32ff",
                    "Shrub" = "#aa373aff",
                    "Other" = "grey90",
                    "Unknown" = "grey90",
                    "AP" = "#0f6b31ff",
                    "Tree" = "#0f6b31ff",
                    "1_Hyper-arid" = "#8c510a",
                    "2_Arid" = "#bf812e", 
                    "3_Semi-Arid" = "#dfc27e",
                    "4_Dry sub-humid" = "#f5e9bf", 
                    "5_Humid" = "#80cec1",
                    Mongolia = "#3e96bdff", Chine = "#f02a26", Uzbekistan = "#6fb440", 
                    Tajikistan = "#e4af08", Russia = "#0035a9", Azerbaijan = "#094227",
                   "Chol cold desert-steppes" = "#7916C4", 
                   "Tugai riparian forest" = "#BB0268", 
                   "Chol warm deserts" = "#bb0202", 
                   "Adyr desert-steppes" = "#ff5400", 
                   "Adyr steppes" = "#e6c607", 
                   "Tau riparian forest" = "#2C9740", 
                   "Tau thermophilous woodlands" = "#85682D", 
                   "Tau juniper steppe-forest" = "#176E5B",
                   "Tau steppes" = "#bab133",
                   "Alau cryophilous steppe-forest" = "#54a697",
                   "Alau meadows" = "#197CDA" 
                    )
      
      
      
      values.bi <- values.bi[which(names(values.bi) %in% unique(Keep.xdata[[Cluster.core]]))]
      Scale.fill <- scale_fill_manual(values = values.bi, name = Cluster.core.lab, drop = T)
      Scale.color <- scale_color_manual(values = values.bi, name = Cluster.core.lab, drop = T)
    }
    Col.select <- Keep.xdata[[Cluster.core]]
    Fill.select <- Keep.xdata[[Cluster.core]]
    }
  else{
    Scale.fill <- scale_fill_manual(values = "brown", name = "Species density")
    Scale.color <- scale_color_manual(values = "brown", name = "Species density")
    Opac = 0.2
    Col.select = "brown"
    Fill.select = "brown"
    }
  
  if(Marg.density.plot == F){
    My_title <- paste(Num.facet, Site.name, Type.samples, sep = " ")
    }
  else{My_title <- NULL}
  
  #### Density contours ####
  if(Density.contour == T){
    Data.contour <- data.frame(MP.pca$ind$coord)
    Data.contour$Col.select <- Col.select
    Data.contour <- unique(Data.contour)
    if(Density.type == "contour"){
      Density.color.line <- "grey20"
      Scale.opa <- NULL
      Point.surf <- NULL
      Dot.up <- "point"
      }
    if(Density.type == "polygon"){
      Density.color.line <- NA
      Scale.opa <- scale_alpha_continuous(range = Opa.range)
      Point.surf <- geom_point(inherit.aes = F, Data.contour, mapping = aes(x = Dim.1, y = Dim.2, colour = Col.select), shape = 16, alpha = Opac, na.rm = T, size = Ptaille)
      Dot.up <- "none"
    }
    Density.contour <- stat_density_2d(data = Data.contour, mapping = aes(x = Dim.1, y = Dim.2, alpha = ..level..),
                                       geom = Density.type, colour = Density.color.line, show.legend = F,
                                       bins = 6)

    }
  else{
    Dot.up <- "point"
    Scale.opa <- NULL
    Density.contour <- NULL}
  
  #### Reverse dim ####
  A = round(MP.pca$eig[1,2], digits = 0)
  B = round(MP.pca$eig[2,2], digits = 0)
  
  if(Reverse.dim == T){
    if(Show.annot == T){
      Note.n <- annotate("text", x = min(MP.pca$ind$coord[,2]), y = min(MP.pca$ind$coord[,1]), label = paste("n = ", nrow(MP), sep = ""), size = 5, hjust = 0)}
    else{Note.n <- NULL}
    Axes <- c(2,1)
    X.arrow <- c(2)
    Y.arrow <- c(1)
    Xlab <- labs(x = substitute(paste("PCA"[2], ~ "(", B, " %)", sep = " " )),
                 y = substitute(paste("PCA"[1], ~ "(", A, " %)", sep = " " )))
      }
  else{
    if(Show.annot == T){
      Note.n <- annotate("text", x = min(MP.pca$ind$coord[,1]), y = min(MP.pca$ind$coord[,2]), label = paste("n = ", nrow(MP), sep = ""), size = 5, hjust = 0, vjust = 0.5)}
    else{Note.n <- NULL}
    
    Axes = c(1,2)
    X.arrow <- c(1)
    Y.arrow <- c(2)
    Xlab <- labs(x = substitute(paste("PCA"[1], ~ "(", A, " %)", sep = " " )),
                 y = substitute(paste("PCA"[2], ~ "(", B, " %)", sep = " " )))
    }
  if(is.null(Manu.lim.x)==F){
    Lim.x <- xlim(Manu.lim.x)
    if(Show.annot == T){Note.n <- annotate("text", x = Manu.lim.x[1], y = min(MP.pca$ind$coord[,2]), label = paste("n = ", nrow(MP), sep = ""), size = 5, hjust = 0, vjust = 0.5)}
    else{Note.n <- NULL}  
    }
  else{Lim.x <- NULL}
  if(is.null(Manu.lim.y)==F){Lim.y <- ylim(Manu.lim.y)}
  else{Lim.y <- NULL}
  
  #### Ajout vectors ####
  Pouet <- data.frame(X0 = 0,
                      Y0 = 0,
                      X1 = MP.pca$var$coord[,X.arrow]*Scale.PCA,
                      Y1 = MP.pca$var$coord[,Y.arrow]*Scale.PCA,
                      Lab = row.names(MP.pca$var$coord),
                      Contrib = MP.pca$var$contrib[,c(1)] + MP.pca$var$contrib[,c(2)] , 
                      Groupes = Groupes)
  Pouet$Contrib <- Pouet$Contrib/sum(Pouet$Contrib)
  
  Arrow.lab <- geom_text_repel(data = Pouet, aes(x = X1, y = Y1, label = Lab, color = Groupes), min.segment.length = 10, force = 20)
  if(Contrib == T & Show.arrow == T){
    Arrow <- geom_segment(data = Pouet, aes(x = X0, y = Y0, xend = X1, yend = Y1, color = Groupes, size = Contrib), arrow = arrow(length=unit(0.2,"cm")))}
  if(Contrib == F & Show.arrow == T){
    Arrow <- geom_segment(data = Pouet, aes(x = X0, y = Y0, xend = X1, yend = Y1, color = Groupes, size = Contrib), arrow = arrow(length=unit(0.2,"cm")))}
  if(Show.arrow == F){
    Arrow <- NULL
    Arrow.lab <- NULL
    }
  Scale.size <- scale_size_continuous(range = c(0.2,1), guide = "none")
  
  #### Ajout symbole ####
  if(is.null(Symbol.path) == F){
    if(is.null(Symbol.pos)== T){Symbol.pos <- c(.9, .9, .16)}
    library(png)
    library(grid)
    img <- readPNG(Symbol.path)
    g <- rasterGrob(x = Symbol.pos[1], y = Symbol.pos[2], width = Symbol.pos[3], height = Symbol.pos[3], img, interpolate = T)
    Logo <- annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
  }
  else{Logo <- NULL}
  
  #### PLOT  ####
  p <- fviz_pca_biplot(MP.pca, na.rm = T,
                       geom.ind = Dot.up, axes = Axes,
                       pointshape = 16,
                       pointsize = Ptaille,
                       fill.ind = Fill.select,
                       alpha.ind = Opac, #alpha.var ="contrib",
                       col.ind = Col.select,
                       repel = T, arrowsize = 0.5,
                       invisible = Centroide, # enlève ou ajoute le centroïde
                       addEllipses = Ellipse, ellipse.level = 0.75, ellipse.type = "norm", # convex or norm
                       ellipse.alpha = Ellipse.opa, 
                       title = My_title,
                       # col.var = Groupes #a enleve parfois
                       col.var = NA
                       ) +
    guides(color = guide_legend(override.aes = list(size = Ptaille+2)),
           # fill = guide_legend(override.aes = list(size = Ptaille+2))
           # fill = "none"
           )+
    #  guides(color = guide_legend(override.aes = list(size = Ptaille+2, nrow = Leg.row)),
    # fill = guide_legend(override.aes = list(size = Ptaille+2, nrow = Leg.row)
    # ))+
    # guides(fill = guide_legend(nrow = Leg.row))+
    Xlab+ #Ylab+ 
    Lim.x+ Lim.y+
    # Axes+
    Scale.opa + 
    Density.contour +
    Scale.fill +
    Scale.color +
    Point.surf +
    new_scale_color() + 
    Arrow +
    Arrow.lab + #, nudge_y = 0.5, nudge_x = 0.5,
    Scale.size +
    Scale.color.vec +
    Note.n + Logo +
    theme(plot.background = element_blank(),
          legend.position = Legend.position,
          legend.title = element_text(size = (Legend.size+2)),
          legend.text = element_text(size = Legend.size),
          legend.key = element_blank(),
          panel.border = element_rect(NA, "black", linewidth = 1),
          panel.grid = element_line(linetype = "dashed"),
          axis.line = element_blank())
  
  
  #### Add margin density ####
  if(Marg.density.plot == T){
    MMarg <- data.frame(MP.pca$ind$coord)
    MMarg <- cbind(MMarg, Biome = Col.select)
    
    List.of.NA <- which(MMarg$Dim.1 < 1e-12 & MMarg$Dim.1 > -1e-12 & MMarg$Dim.2 < 1e-12 & MMarg$Dim.2 > -1e-12)
    if(length(List.of.NA) > 0){
      print("Remove NA from density.")
      MMarg <- MMarg[-List.of.NA,]}
    
    #### Density plots up ####
    plot_top <- ggplot(MMarg, aes(x = Dim.1, fill = Biome)) + 
      geom_density(alpha = 0.6, size = 0.1) + Scale.fill + 
      ggtitle(paste(Num.facet, Site.name, Type.samples, sep = " "))+
      #### Theme ####
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(), axis.text.y = element_text(hjust = 1, size = 6),
      axis.ticks.x.bottom = element_blank(),
      axis.title = element_blank(),
      axis.line.y = element_line(colour = "grey"),
      axis.ticks.y = element_line(colour = "grey"),
      #panel.border = element_rect(fill = NA, colour = "grey"),
      legend.title = element_text(),
      legend.key = element_blank(),
      legend.justification = c("center"),               # left, top, right, bottom
      legend.text = element_text(size = 8),
      panel.background = element_blank(),
      panel.spacing = unit(0.7, "lines"),
      legend.position = "none",
      strip.text.x = element_text(size = 12, angle = 0, face = "bold"),
      strip.placement = "outside",
      strip.background = element_rect(color = "white", fill = "white"),
      plot.margin=unit(c(0,0,0,0),"cm")
    )
    
    #### Density plots right ####
    plot_right <- ggplot(MMarg, aes(x = Dim.2, fill = Biome)) + 
      geom_density(alpha = 0.6, size = 0.1) + Scale.fill +
      coord_flip() +
      #### Theme ####
    theme(
      axis.line.y = element_blank(),
      axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      axis.line.x = element_line(colour = "grey"),
      axis.ticks.x = element_line(colour = "grey"),
      #panel.border = element_rect(fill = NA),
      legend.title = element_text(),
      legend.key = element_blank(),
      legend.justification = c("center"),               # left, top, right, bottom
      legend.position = "none",
      legend.text = element_text(size = 8),
      panel.background = element_blank(),
      panel.spacing = unit(0.7, "lines"),
      strip.text.x = element_text(size = 12, angle = 0, face = "bold"),
      strip.placement = "outside",
      strip.background = element_rect(color = "white", fill = "white"),
      plot.margin=unit(c(0,0,0,0),"cm")
    )
    
    
    #### Export ####
    layout <- "AAAAAAA#
             CCCCCCCB
             CCCCCCCB
             CCCCCCCB
             CCCCCCCB
             CCCCCCCB
             "
    pfull <- plot_top + plot_right + p + plot_layout(design = layout)
    print(pfull)
    }
  else{print(p)}
  #### Save html ####
  if(Show.Plotly == T){
    library(plotly)
    library(htmlwidgets)
  
    #### 3D PCA chart ####
    Pouet <- data.frame(MP.pca$ind$coord)
    Pouet$Col <- Col.select
    fig <- plot_ly(Pouet, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3, color = ~Col, 
                   colors = values.bi, text = ~paste('Site: ', row.names(Pouet), "\n Biome:", Pouet$Col, sep = "")) #%>%
    options(warn = - 1)
      # add_markers(size = 60)
    # 
    # fig <- fig %>%
    #   layout(
    #     # title = tit,
    #     xaxis = list(title = "PCA1"),
    #     scene = list(bgcolor = "#e5ecf6")
    #   )
        
    #### Export ####
    Save.plot.html <- gsub("pdf", "html", Save.plot)
    Keep.name <- gsub(".*\\/", "", Save.plot.html)
    Path.root <- paste(gsub(Keep.name, "", Save.plot.html), "HTML_files/", sep = "")
    if(file.exists(Path.root) == F){dir.create(Path.root)}
    Save.plot.html <- paste(Path.root, Keep.name, sep = "")
    # # p1_ly <- ggplotly(p)
    # # p1_ly <- p1_ly %>% layout()
    options(warn = - 1)
    saveWidget(fig, file = Save.plot.html)
    options(warn = - 1)
    # saveWidget(p1_ly, file = Save.plot.html)
  }
  
  #### Export data ####
  if(is.null(Save.path) == F){
    Site.name <- gsub(" ","_",Site.name)
    Save.path.Site <- gsub("\\.csv", "_PCA_site.csv", Save.path)
    Save.path.Taxon <- gsub("\\.csv", "_PCA_clim_param.csv", Save.path)
    write.table(data.frame(MP.pca$ind$coord), file = Save.path.Site, col.names = T, sep = ",")
    write.table(data.frame(MP.pca$var$coord), file = Save.path.Taxon, col.names = T, sep = ",")}
  if(is.null(Save.plot) == F){
    dev.off()}
  
  if(return.pick == T){return(p)}
  if(return.pick == F){
    if(PCA.site == T){
      Mat.ret <- data.frame(MP.pca$ind$coord)
      # print(Mat.ret)
      return(Mat.ret)
      }
    if(Contrib == T){return(MP.pca$var$contrib)}
    }
}

PCoIA.plot <- function(DF1, DF2, transp_OK, H, W, Save.plot, Show.outliers, Ylim, return.pick){ 
  #### Settings ####
  library(ade4)
  if(missing(Ylim)){Ylim = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(transp_OK)){transp_OK = T}
  if(missing(return.pick)){return.pick = F}
  if(missing(Show.outliers)){Show.outliers = T}
  
  #### Calcul PCA and PCoIA ####
  PCA.Pss <- PCA.bioclim(DF1, transp_OK = transp_OK, return.pick = F, Contrib = F, PCA.site = T)
  PCA.Psl <- PCA.bioclim(DF2, transp_OK = transp_OK, return.pick = F, Contrib = F, PCA.site = T)
  pro <- procrustes(X = PCA.Pss, Y = PCA.Psl, symmetric = T)
  PROTEST <- protest(PCA.Pss, PCA.Psl)
  ctest <- data.frame(rda1=pro$Yrot[,1],
                      rda2=pro$Yrot[,2], xrda1=pro$X[,1],
                      xrda2=pro$X[,2],biome=as.factor(MCWT.clim.PT_ss$Biome))
  
  ctest$Procrustes_error <- sqrt((ctest$rda1-ctest$xrda1)^2+(ctest$rda2-ctest$xrda2)^2)
  ctest <- ctest[order(ctest$Procrustes_error,decreasing = T),]
  row.names(ctest) <- seq(1:nrow(ctest))
  ctest$ID <- as.numeric(row.names(ctest))

  #### Colors settings ####
  values.bi = c("Deserts & Xeric Shrublands" = "#C88282",
                "Temperate Grasslands, Savannas & Shrublands" = "#ECED8A",
                "Montane Grasslands & Shrublands" = "#D0C3A7",
                "Temperate Conifer Forests" = "#6B9A88",
                "Temperate Broadleaf & Mixed Forests" = "#3E8A70",
                "N/A" = "#FFEAAF",
                "Tundra" = "#A9D1C2",
                "Boreal Forests/Taiga" = "#8FB8E6",
                "Tropical & Subtropical Coniferous Forests" = "#99CA81",
                "Mangroves" = "#FE01C4",
                "Flooded Grasslands & Savannas" = "#BEE7FF",
                "Tropical & Subtropical Moist Broadleaf Forests" = "#38A700",
                "Water" = "darkblue",
                "Temperature" = "darkred",
                "Plant_height" = "royalblue",
                "Leaf_thickness" = "darkorange",
                "Photosynthesis_pathway" = "purple",
                "Woodyness" = "#323232",
                "Leaf_size" = "darkred",
                "Variable" = "#aa373aff",
                #"Unknown" = "#aa373aff",
                "Algal" = "#4666E9",
                "NAP" = "#b5ab32ff",
                "Herb" = "#b5ab32ff",
                "Shrub" = "#aa373aff",
                "Other" = "grey90",
                "Unknown" = "grey90",
                "AP" = "#0f6b31ff",
                "Tree" = "#0f6b31ff"
  )
  Scale.fill <- scale_fill_manual(values = values.bi, guide = "none")
  Scale.color <- scale_color_manual(values = values.bi, guide = "none")
  
  if(Show.outliers == T){Col.out = "grey60"}
  if(Show.outliers == F){Col.out = NA}
  #### Other graphic param ####
  # A = round(MP.pca$eig[1,2], digits = 2)
  # B = round(MP.pca$eig[2,2], digits = 2)
  
  # Xlab <- labs(x = substitute(paste("RDA"[1], ~ "(", A, " %)", sep = " " )),
  #              y = substitute(paste("RDA"[2], ~ "(", B, " %)", sep = " " )))
  Xlab <- labs(x = substitute(paste("RDA"[1])),
               y = substitute(paste("RDA"[2])))
  
  My_theme <- theme(axis.line = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_rect(fill = NA),
                    plot.background = element_blank(), 
                    panel.grid = element_blank(), 
                    # strip.text = Strip.lab.disp,
                    strip.background = element_blank(),
                    legend.position = "none"
                    )
  
  if(is.null(Ylim) == F){Ylim = xlim(Ylim)}
  
  Note.n <- annotate("text", x = 0, y = max(ctest$Procrustes_error), 
                     label = paste("m12^2 = ", round(pro$ss, digits = 3), 
                                   ",\nPROTEST correlation = ", round(PROTEST$t0, digits = 3), 
                                   ",\nPROTEST p-value > ", round(PROTEST$signif, digits = 3), 
                                   sep = ""), size = 4, color = "grey30", hjust = 0, vjust = 1)
  #### Biomes labs ####
  ctest$Biom.lab <- as.factor(ctest$biome)
  levels(ctest$Biom.lab) <- c("Boreal Forests/Taiga"                        = "Taiga",# "TAIG",                       
                              "Deserts & Xeric Shrublands"                  = "Desert",#  "DESE",                 
                              "Montane Grasslands & Shrublands"             = "Mont. grass.",#  "MONT",            
                              "Temperate Broadleaf & Mixed Forests"         = "Mix. Forest",#  "TEBF",        
                              "Temperate Conifer Forests"                   = "Con. Forest",#  "TECF",                  
                              "Temperate Grasslands, Savannas & Shrublands" = "Steppe",# "STEP",
                              "Tundra"                                      = "TUND") 
  
  ctest <- ctest[order(ctest$Biom.lab),]
  ctest$ID_by_biome <- seq(1, nrow(ctest))
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Plots ####
  p <- ggplot(ctest) +
    geom_point(aes(x = rda1, y = rda2, color = biome), size = 1.3, na.rm = T) +
    # geom_point(aes(x=xrda1, y=xrda2, color=biome), size = 0) +
    geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2), color = "grey30", alpha = 0.4, arrow = arrow(length=unit(0.2,"cm")))+
    facet_wrap(vars(biome), scales = "free")+
    Xlab+
    # guides(color = NA)+
    Scale.color+
    My_theme
  
  p2 <- ggplot(ctest, aes(x = ID_by_biome, y = Procrustes_error)) +
    # stat_summary(aes(Procrustes_error), geom="line", lwd=2)+
    # geom_segment(aes(x = ID, y = 0, xend = ID, yend = Procrustes_error, color = biome))+
    geom_bar(aes(fill = biome), position = 'dodge', stat='identity')+
    geom_hline(yintercept = mean(ctest$Procrustes_error), color="grey30", lwd = 0.2)+
    geom_hline(yintercept = quantile(ctest$Procrustes_error)[c(2,4)], color="grey30", lwd = 0.2, linetype = "longdash")+
    xlab("Samples ID")+
    ylab("Procrustes errors")+
    Note.n+
    Scale.color+
    Scale.fill+
    theme_classic()+
    My_theme
    
  
  p3 <- ggplot(ctest) +
    geom_vline(xintercept = quantile(ctest$Procrustes_error)[c(2,4)], color="grey30", lwd = 0.2, linetype = "longdash")+
    geom_vline(xintercept = mean(ctest$Procrustes_error), color="grey10", lwd = 0.4)+
    geom_boxplot(aes(Procrustes_error, fill = biome, y = Biom.lab), outlier.colour = Col.out, notch = T, outlier.size = 0.8)+
    xlab("Procrustes errors")+
    ylab("Biomes")+
    coord_flip()+
    Ylim +
    Scale.fill+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    My_theme
  
  p2 <- p/(p2 + p3 + plot_layout(widths = c(3/4, 1/4))) + plot_layout(heights = c(2/3, 1/3))
  
  if(is.null(Save.plot) == F){
    print(p2)
    dev.off()}
  
  if(return.pick == T){
    return(p2)}
  else{return(ctest)}
  }

Biplot.bioclim <- function(MC.area, MC.samp, Same.mat, PClim1, PClim2, Show.density, Emprise.alpha, Add.reg, R2.pos,
                           Emprise.bin, Emprise.size, return.pick, Save.plot, H, W, Leg.pos, Mean.cal, Limites){
  #### Settings ####
  if(missing(Add.reg)){Add.reg = F}
  if(missing(Emprise.bin)){Emprise.bin = F}
  if(missing(return.pick)){return.pick = F}
  if(missing(Emprise.alpha)){Emprise.alpha = 1}
  if(missing(Emprise.size)){Emprise.size = 3}
  if(missing(Leg.pos)){Leg.pos = c(0.15, 0.85)}
  if(missing(Limites)){Limites = NULL}
  if(missing(Show.density)){Show.density = T}
  if(missing(Same.mat)){Same.mat = F}
  if(missing(Mean.cal)){Mean.cal = F}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(R2.pos)){R2.pos = "bottomleft"}
  
  
  MC.area <- MC.area[!is.na(MC.area$Biome),] 
  MC.samp <- MC.samp[!is.na(MC.samp$Biome),] 
  # MC.samp <- MC.samp[!MC.samp$Biome == "N/A",] 
  # MC.samp <- MC.samp[!MC.samp$Biome == "Tundra",] 
  # MC.samp <- MC.samp[!MC.samp$Biome == "Temperate Conifer Forests",]
  # MC.samp <- MC.samp[!MC.samp$Biome == "Temperate Broadleaf & Mixed Forests",]
  # MC.samp <- MC.samp[!MC.samp$Biome == "Tropical & Subtropical Coniferous Forests",]
  # MC.samp <- MC.samp[!MC.samp$Biome == "Boreal Forests/Taiga",]
  # MC.samp <- MC.samp[!MC.samp$Biome == "Deserts & Xeric Shrublands",]
  # 
  # MC.samp[MC.samp$Biome == "Montane Grasslands & Shrublands", "Biome"] <- "A"
  # MC.samp[MC.samp$Biome == "Temperate Grasslands, Savannas & Shrublands", "Biome"] <- "B"
  # 
  #### Stats ####
  if(Mean.cal == T){
    mu2 <- aggregate(MC.samp[[PClim2]], list(MC.samp$Biome), FUN=mean) 
    names(mu2)[1]<- "Biome"
    mu1 <- aggregate(MC.samp[[PClim1]], list(MC.samp$Biome), FUN=mean) 
    names(mu1)[1]<- "Biome"
    Mean.line1 <- geom_vline(data = mu1, aes(xintercept = x, color = Biome), linetype="dashed")
    Mean.line2 <- geom_vline(data = mu2, aes(xintercept = x, color = Biome), linetype="dashed")
    }
  else{
    Mean.line1 <- NULL
    Mean.line2 <- NULL
  }
  
  #### Colours settings ####
  Col.vec <- c("Deserts & Xeric Shrublands" = "#C88282",
               "Temperate Grasslands, Savannas & Shrublands" = "#ECED8A",
               "Montane Grasslands & Shrublands" = "#D0C3A7",
               "Temperate Conifer Forests" = "#6B9A88",
               "Temperate Broadleaf & Mixed Forests" = "#3E8A70",
               "N/A" = "#FFEAAF",
               "Tundra" = "#A9D1C2",
               "Boreal Forests/Taiga" = "#8FB8E6",
               "Tropical & Subtropical Coniferous Forests" = "#99CA81",
               "Mangroves" = "#FE01C4",
               "Flooded Grasslands & Savannas" = "#BEE7FF",
               "Tropical & Subtropical Moist Broadleaf Forests" = "#38A700"
               )
  
  #### Limits ####
  if(is.null(Limites) == F){
    Lim.x <- xlim(c(Limites[1], Limites[2]))
    Lim.y <- ylim(c(Limites[3], Limites[4]))
    Lim.up <- Lim.x 
    Lim.right <- xlim(c(Limites[3], Limites[4])) 
    }
  else{
    Lim.x <- NULL
    Lim.y <- NULL
    Lim.up <- NULL
    Lim.right <- NULL
    }
  
  #### Rasterization test ####
  # library(MASS)
  # library(tidyr)
  # xlim <- range(MC.samp[[PClim1]])
  # ylim <-range(MC.samp[[PClim2]])
  # newplot_data <- MC.samp %>% group_by(Biome) %>% do(Dens=kde2d(.[[PClim1]], .[[PClim2]], n=100, lims=c(xlim,ylim)))
  # newplot_data  %<>%  do(Biome=.$Biome, V=expand.grid(.$Dens$x,.$Dens$y), Value=c(.$Dens$z)) %>% do(data.frame(Biome=.$Biome,x=.$V$Var1, y=.$V$Var2, Value=.$Value))
  # newplot_data  %<>%   spread(Biome, value=Value) %>%
  #   mutate(Level = if_else(A > B, A, B), Biome = if_else(A > B,"Montane Grasslands & Shrublands", "Temperate Grasslands, Savannas & Shrublands"))
  
  #### Add regression lines ####
  if(Add.reg == T){print(R2.pos)
    if(R2.pos == "bottomleft"){
      R2.y = "bottom"
      R2.x = "left"}
    if(R2.pos == "bottomright"){
      R2.y = "bottom"
      R2.x = "right"} 
    if(R2.pos == "topright"){
      R2.y = "top"
      R2.x = "right"}
    if(R2.pos == "none"){
      R2.y = "none"
      R2.x = "none"}
    print(R2.x)
    Add.reg.line <- geom_smooth(data = MC.samp, mapping = aes(x = eval(parse(text = PClim1)), y = eval(parse(text = PClim2))),
                                method = "lm", se = F, span = 1000, size = 0.7, linetype = "dashed", colour = "grey20",
                                formula = y ~ x)  
    
    Add.r2 <- stat_poly_eq(data = MC.samp, 
                           label.y = R2.y, label.x = R2.x, 
                           size = 3, small.r = F, vstep = 0.07, p.digits = 1, na.rm = T,  
                           aes(x = eval(parse(text = PClim1)), y = eval(parse(text = PClim2)),
                             label =  sprintf("%s*\", \"*%s" ,
                                                after_stat(rr.label),
                                                # after_stat(r.squared),
                                                after_stat(p.value.label)
                           )))
  }
  else{
    Add.r2 <- NULL
    Add.reg.line <- NULL
  
  }
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  #### Same matrice 1 et 2 ####
  if(Same.mat == F){
    if(Emprise.bin == F){
      Point.grey <- geom_point(data = MC.area, mapping = aes(x = eval(parse(text = PClim1)), y = eval(parse(text = PClim2))), 
                               size = Emprise.size, shape = 15, alpha = Emprise.alpha, colour = "grey70")}
    else{Point.grey <- geom_hex(data = MC.area, mapping = aes(x = eval(parse(text = PClim1)), y = eval(parse(text = PClim2)), 
                                # alpha = ..density.., 
                                fill = Biome),
             bins = Emprise.size, 
                                  # size = Emprise.size, shape = 15, 
             alpha = Emprise.alpha, fill = "grey70"
    )}
    
    } 
  else{Point.grey <- NULL}
  
  #### PLOT ####
  p <- ggplot()+
    Point.grey + 
    # geom_hex(data = MC.samp, mapping = aes(x = eval(parse(text = PClim1)), y = eval(parse(text = PClim2)), alpha = ..density.., fill = Biome),
    #          bins = 30#, colour = "grey60"
    #          ) +
    stat_density_2d(data = MC.samp, mapping = aes(x = eval(parse(text = PClim1)), y = eval(parse(text = PClim2)), alpha = ..level.., fill = Biome),
                    geom = "polygon", #colour = "grey90",
                    bins = 8) +
    geom_point(data = MC.samp, mapping = aes(x = eval(parse(text = PClim1)), y = eval(parse(text = PClim2)), fill = Biome), size = 1.5, shape = 21, alpha = 1, colour = "grey40")+
    #geom_point(data = MC.samp, mapping = aes(x = eval(parse(text = PClim1)), y = eval(parse(text = PClim2))), fill = NA, size = .7, shape = 21, alpha = 0.6, colour = "black")+
    scale_fill_manual(values = Col.vec, name = "Biomes (Dinerstein et al., 2017)")+
    Add.r2 + Add.reg.line + 
    # scale_alpha_continuous(limits=c(0,0.02)) +
    scale_alpha_continuous(limits=c(0,0.6)) +
    scale_color_manual(values = Col.vec, name = "Biomes (Dinerstein et al., 2017)")+
    guides(colour = FALSE, alpha = F)+
    xlab(PClim1)+
    ylab(PClim2)+
    Lim.x +
    Lim.y +
    #### Theme ####
  theme(
    axis.line= element_blank(),
    axis.ticks.x.bottom = element_line(colour = "grey"),
    panel.border = element_rect(fill = NA, colour = "grey"),
    legend.title = element_text(),
    legend.key = element_blank(),
    legend.justification = c("center"),               # left, top, right, bottom
    legend.text = element_text(size = 8),
    panel.background = element_blank(),
    panel.spacing = unit(0.7, "lines"),
    strip.text.x = element_text(size = 12, angle = 0, face = "bold"),
    strip.placement = "outside",
    # legend.position = "none",
    legend.position = Leg.pos,
    strip.background = element_rect(color = "white", fill = "white"),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")
  )
  
  
  #### Density plots up ####
  plot_top <- ggplot(MC.samp, aes(x = eval(parse(text = PClim1)), fill = Biome)) + 
    geom_density(alpha = 0.6, size = 0.2) +
    Mean.line1 +
    Lim.up +
    scale_color_manual(values = Col.vec) + 
    scale_fill_manual(values = Col.vec)+
    #### Theme ####
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(), axis.text.y = element_text(hjust = 1, size = 6),
      axis.ticks.x.bottom = element_blank(),
      axis.title = element_blank(),
      axis.line.y = element_line(colour = "grey"),
      axis.ticks.y = element_line(colour = "grey"),
      #panel.border = element_rect(fill = NA, colour = "grey"),
      legend.title = element_text(),
      legend.key = element_blank(),
      legend.justification = c("center"),               # left, top, right, bottom
      legend.text = element_text(size = 8),
      panel.background = element_blank(),
      panel.spacing = unit(0.7, "lines"),
      legend.position = "none",
      strip.text.x = element_text(size = 12, angle = 0, face = "bold"),
      strip.placement = "outside",
      strip.background = element_rect(color = "white", fill = "white"),
      plot.margin=unit(c(0,0,0,0),"cm")
    )
    
  #### Density plots right ####
  plot_right <- ggplot(MC.samp, aes(x = eval(parse(text = PClim2)), fill = Biome)) + 
    geom_density(alpha = 0.6, size = 0.2) +
    Mean.line2 +
    Lim.right +
    coord_flip() + 
    scale_color_manual(values = Col.vec) + 
    scale_fill_manual(values = Col.vec)+
    #### Theme ####
    theme(
      axis.line.y = element_blank(),
      axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      axis.line.x = element_line(colour = "grey"),
      axis.ticks.x = element_line(colour = "grey"),
      #panel.border = element_rect(fill = NA),
      legend.title = element_text(),
      legend.key = element_blank(),
      legend.justification = c("center"),               # left, top, right, bottom
      legend.position = "none",
      legend.text = element_text(size = 8),
      panel.background = element_blank(),
      panel.spacing = unit(0.7, "lines"),
      strip.text.x = element_text(size = 12, angle = 0, face = "bold"),
      strip.placement = "outside",
      strip.background = element_rect(color = "white", fill = "white"),
      plot.margin=unit(c(0,0,0,0),"cm")
    )
    
 
  #### Export ####
  if(Show.density == T){
    print("here")
    layout <- "AAAAAAA#
               CCCCCCCB
               CCCCCCCB
               CCCCCCCB
               CCCCCCCB
               CCCCCCCB
               "
    pfull <- plot_top + plot_right + p + plot_layout(design = layout)
    print(pfull)
    }
  if(Show.density == F){
    pfull <- p
    print(pfull)}
    
  
  if(is.null(Save.plot) == F){dev.off()}
  return(pfull)
}

Trait.distribution <- function(MT, MP1, MP2, Selec.tax.to.shown, Select.trait, Plot.gaussian,
                               Show.var, Scale.var = F, Boundaries, Diaz.compa, Log.scale, W, H, Save.plot){
  #### Settings ####
  if(missing(Diaz.compa)){Diaz.compa = F}
  if(missing(Show.var)){Show.var = F}
  if(missing(Log.scale)){Log.scale = F}
  if(missing(Plot.gaussian)){Plot.gaussian = F}
  if(missing(Selec.tax.to.shown)){Selec.tax.to.shown = c("Poaceae", "Amaranthaceae")}
  if(missing(Select.trait)){Select.trait = c("Height", "SLA")}
  if(missing(Boundaries)){Boundaries = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(MP2)){MP2 = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  Trait.labels <- readRDS("Resultats/ACA/Traits/Corresp_tax/Trait_ID_names.Rds")
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Calculation (new) ####
  # names(MT)[grep("X", names(MT))] <- as.vector(Trait.labels[match(gsub("X", "", names(MT)[grep("X", names(MT))]), names(Trait.labels))])
  # 
  # Common.sp = intersect(MP1$species, MT$species)
  # MT <- MT[which(Common.sp %in% MT$species),]
  # MP1 <- MP1[which(Common.sp %in% MP1$species),]
  # MP1 <- MP1[match(MP1$species, MT$species),]
  # MT1 <- cbind(MT, MP1)
  # 
  # # Ptype <- gsub(".*\\.", "-", deparse(substitute(MP1)))
  # Ptype <- "-ss"
  # names(MP1) <- paste(names(MP1), Ptype, sep = "")
  # print(names((MP1)))
  # MT1 <- cbind(MT, MP1)
  # # View(MT1)
  # # MT1 <- left_join(MT, MP1, by = "species")
  # Keep.col <- match(c(paste(Selec.tax.to.shown, Ptype, sep = ""), Select.trait), names(MT1))
  # Keep.col <- Keep.col[!is.na(Keep.col)]
  # MT1 <- MT1[,Keep.col]
  # 
  # A <- melt(MT1, id = Select.trait)
  # A <- A[A$value == T,]
  # A <- A[rowSums(is.na(A[,Select.trait])) != ncol(A[,Select.trait]), ]
  # A$Type <- paste("ACASP", Ptype, sep = "")
  
  #### Calculation (old) ####
  names(MT)[grep("X", names(MT))] <- as.vector(Trait.labels[match(gsub("X", "", names(MT)[grep("X", names(MT))]), names(Trait.labels))])
  
  Ptype <- gsub(".*\\.", "-", deparse(substitute(MP1)))
  names(MP1) <- paste(names(MP1), Ptype, sep = "")
  MT1 <- cbind(MT, MP1)
  Keep.col <- match(c(paste(Selec.tax.to.shown, Ptype, sep = ""), Select.trait), names(MT1))
  Keep.col <- Keep.col[!is.na(Keep.col)]
  MT1 <- MT1[,Keep.col]
  
  A <- melt(MT1, id = Select.trait)
  A <- A[A$value == T,]
  A <- A[rowSums(is.na(A[,Select.trait])) != ncol(A[,Select.trait]), ]
  A$Type <- paste("ACASP", Ptype, sep = "")
  
  if(is.null(MP2) == F){
    MT <- MT[which(MT$species %in% intersect(MT$species, MP2$species)),] 
    MP2 <- MP2[which(MT$species %in% intersect(MT$species, MP2$species)),]
    
    Ptype2 <- gsub(".*\\.", "-", deparse(substitute(MP2)))
    names(MP2) <- paste(names(MP2), Ptype2, sep = "")
    MT2 <- cbind(MT, MP2)
    Keep.col <- match(c(paste(Selec.tax.to.shown, Ptype2, sep = ""), Select.trait), names(MT2))
    
    Keep.col <- Keep.col[!is.na(Keep.col)]
    MT2 <- MT2[,Keep.col]
    
    MT2 <- melt(MT2, id = Select.trait)
    MT2 <- MT2[MT2$value == T,]
    MT2 <- MT2[rowSums(is.na(MT2[,Select.trait])) != ncol(MT2[,Select.trait]), ]
    MT2$Type <- paste("ACASP", Ptype2, sep = "")
    A <- rbind(A, MT2)
    A$variable <- gsub(Ptype2, "", A$variable)
    }
  
  A$variable <- gsub(Ptype, "", A$variable)
  A <- A[,!(names(A) %in% c("value"))]
  A <- melt(A, id = c("variable", "Type"))
  A <- A[!is.na(A$value),]
  names(A)[1] <- "Taxa"
  A$Taxa <- factor(A$Taxa, levels = Selec.tax.to.shown, ordered = T)
  
  
  #### Add Diaz et al., 2016 ####
  if(Diaz.compa == T){
    DB.diaz <- data.frame(read.csv("Import/World_DB/Traits/TRY/Try202337115611480_Dataset/Dataset_Global_Spectrum/Species_mean_traits.csv", sep = "\t", dec = "."))
    DB.diaz <- DB.diaz[,!grepl("n.o", names(DB.diaz))]
    DB.diaz <- DB.diaz[,c(2,16:20,22,24)]
    names(DB.diaz) <- c("Species", "LeafArea", "LeafN", "SLA", "Height", "SeedMass", "LDMC", "SSD")
    # DB.diaz$SLA <- 1000/(DB.diaz$SLA)
    DB.diaz$SLA <- 1/(DB.diaz$SLA)
    DB.diaz <- DB.diaz[,match(c("Species", Select.trait), names(DB.diaz))]
    DB.diaz <- melt(DB.diaz, id = "Species")
    DB.diaz$Taxa <- "Global Spectrum"
    DB.diaz$Type <- "Global Spectrum"
    
    # print(names(DB.diaz))
    
    A <- full_join(A, DB.diaz, by = c("Taxa", "Type", "variable", "value"))
    Level_diaz <-  c("Global Spectrum", Selec.tax.to.shown)
    A$Taxa <- factor(A$Taxa, levels = Level_diaz, ordered = T)
    Barre.global <- geom_vline(xintercept = 1.5, lty="dashed")
    B0 <- geom_rect(data = data.frame(xmin = 0, xmax = 1.5), inherit.aes = F, mapping = aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill = "grey50", alpha = .5, color = "grey50", size = 0, linetype = 0)
    
    }
    else{
      Level_diaz <-  Selec.tax.to.shown
      Barre.global <- NULL
      B0 <- NULL
      }
  
  A$variable <- factor(A$variable, levels = unique(A$variable), ordered = T)
  
  #### Parameters label ####
  A$Lab <- as.character(A$variable)
  A$Lab[A$Lab == "SLA"] = "SLA~log(mm2.mg^-1)"
  A$Lab[A$Lab == "SSD"] = "SSD~(mg.mm^-3)"
  A$Lab[A$Lab == "SeedMass"] = "SeedMass~log(mg)"
  A$Lab[A$Lab == "Height"] = "Height~log(m)"
  A$Lab[A$Lab == "LeafArea"] = "LeafArea~log(mm^2)"
  A$Lab[A$Lab == "LeafN"] = "N[leaf]~(mg.g^-1)"
  
  #### Log et non log ####
  if(Log.scale == T){
    Trait.log.dist2 <- c("SeedMass", "Height", "LeafArea", "SLA")
    A[A$variable %in% Trait.log.dist2,4] <- log(A[A$variable %in% Trait.log.dist2,4])
    # MT.hist.log <- A[A$variable %in% Trait.log.dist2,]
    # MT.hist <- A[!A$variable %in% Trait.log.dist2,]
    # X.axe <- element_blank()
    
    # MT.hist.log <- NULL
    MT.hist <- A
    X.axe <- element_text(size = 10, angle = 45, hjust = 1, vjust = 1)
    
    }
  else{
    # Trait.log.dist2 <- c("SeedMass", "Height", "LeafArea")
    # A[A$variable %in% Trait.log.dist2,4] <- log(A[A$variable %in% Trait.log.dist2,4])
    
    # MT.hist.log <- NULL
    MT.hist <- A
    X.axe <- element_text(size = 10, angle = 45, hjust = 1, vjust = 1)
    }
  
  #### Calculate variance ####
  if(Show.var == T){
    Var.matrice <- function(A, Log, Scale.var = F){
      A <- A[which((!is.na(A$value))),]
      MVar <- setNames(data.frame(matrix(nrow = length(unique(A$Taxa)), ncol = length(unique(A$variable)))), unique(A$variable))
      row.names(MVar) <- unique(A$Taxa)
      for(i in unique(A$variable)){
        B <- A[A$variable == i,]
        for(j in unique(A$Taxa)){
          C <- B[B$Taxa == j,]
          if(Log == F){My_var <- var(C$value, na.rm = T)}
          # if(Log == T){My_var <- var(log(C$value))}
          # if(Log == F){My_var <- sd(C$value)}
          # if(Log == T){My_var <- sd(log(C$value))}
          MVar[j,i] <- My_var
        }
      }
      if(Scale.var == T){MVar <- scale(MVar, center = F)}
      else{
        MVar <- data.frame(MVar)
        MVar$Taxa <- row.names(MVar)
        }
      
      MVar <- melt(MVar)
      names(MVar) <- c("Taxa", "Trait", "value")
      MVar$Taxa <- factor(MVar$Taxa, levels = Level_diaz, ordered = T)
      
      MVar$Lab <- as.character(MVar$Trait)
      MVar$Lab[MVar$Lab == "SLA"] = "SLA~(mm2.mg-1)"
      MVar$Lab[MVar$Lab == "SSD"] = "SSD~(mg.mm^-3)"
      MVar$Lab[MVar$Lab == "SeedMass"] = "SeedMass~log(mg)"
      MVar$Lab[MVar$Lab == "Height"] = "Height~log(m)"
      MVar$Lab[MVar$Lab == "LeafArea"] = "LeafArea~log(mm^2)"
      MVar$Lab[MVar$Lab == "LeafN"] = "N[leaf]~(mg.g^-1)"
      return(MVar)}
    
    MVar <- Var.matrice(MT.hist, Log = F, Scale.var = Scale.var)
    # print(MVar)
    Taxa.lab <- element_blank()
    Ticks <- element_blank()
    Tit.lab <- element_blank()
    Title <- NULL
    Legend.pos <- "none"
    
    #### Anova ? ####
    print(unique(MT.hist$Type))
    Two.way <- aov(value ~ Type, data = MT.hist)
    print(summary(Two.way))
    }
  else{
    Title <- xlab("Pollen-types")
    Ticks <- NULL
    Tit.lab <- NULL
    Legend.pos <- "bottom"
    Taxa.lab <- element_text(size = 10, angle = 45, hjust = 1, vjust = 1)
    }
  
  
  #### Add vertical lines ####
  if(is.null(Boundaries) == T){
    Barre1 = max(grep("eae", Selec.tax.to.shown))+0.5
    Barre2 = min(grep("\\s", Selec.tax.to.shown))-0.5}
  else{
    Barre1 = Boundaries[1]
    Barre2 = Boundaries[2]
    Barre3 = Boundaries[3]
    
  }
  # print(typeof(Boundaries[1]))
  
  B1 <- geom_rect(data = data.frame(xmin = 0, xmax = Barre1), inherit.aes = F, mapping = aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill = "grey50", alpha = .5, color = "grey50", size = 0, linetype = 0)
  B2 <- geom_rect(data = data.frame(xmin = Barre1, xmax = Barre2), inherit.aes = F, mapping = aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill = "grey65", alpha = .5, color = "grey65", size = 0, linetype = 0)
  B3 <- geom_rect(data = data.frame(xmin = Barre2, xmax = Barre3), inherit.aes = F, mapping = aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill = "grey85", alpha = .5, color = "grey85", size = 0, linetype = 0)
  B4 <- geom_rect(data = data.frame(xmin = Barre3, xmax = length(Selec.tax.to.shown)+0.5), inherit.aes = F, mapping = aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill = "white", alpha = .5, color = "white", size = 0, linetype = 0)
  
  #### Plot type boite à moustache ####
  if(Plot.gaussian == F){
    #### Color ####
    Col.vec <- c("ACASP-sl" = "#ff7f2aff",
                 "ACASP-ss" = "#ffcc00ff",
                 "Global Spectrum" = "#ff6e91"
                 )
    
    Col.vec <- Col.vec[which(names(Col.vec) %in% unique(MT.hist$Type))]
    
    #### Plot 1 (normal) ####
    p1 <- ggplot(data = MT.hist, aes(x = Taxa, y = value))+
      facet_wrap(~Lab, scales = "free_y", labeller = label_parsed)+
      geom_violin(alpha = 1, fill = "grey40", na.rm = F, show.legend = F, trim = F, scale = "width", color = NA) +
      B0 + B1 + B2 + B3 + B4 +
      Barre.global + geom_vline(xintercept = Barre1, lty="dashed")+
      geom_vline(xintercept = Barre2, lty="dashed")+
      geom_vline(xintercept = Barre3, lty="dashed")+
      geom_boxplot(aes(fill = Type), outlier.color = "grey20", outlier.shape = 16,  outlier.size = .8, outlier.alpha = 0.3, alpha = 0.7, varwidth = F, na.rm = T, show.legend = T) +
      # geom_jitter(shape = 16, size = 1, color = "grey30", alpha = 0.35, position = position_jitter(width = .15), na.rm = F) +
      # stat_compare_means(method = "anova", label.y = 100)+
      scale_fill_manual(values = Col.vec)+
      ylab("Trait values")+
      theme(
        panel.background = element_blank(), panel.grid = element_blank(), 
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = Taxa.lab, axis.ticks.x = Ticks,
        panel.border = element_rect(fill = NA), legend.background = element_blank(), strip.background = element_blank(),
        strip.text = element_text(size = 12), plot.margin = margin(0,0,0,0)
      )

    }
  
  if(Show.var == T){
  p1.var <- ggplot(data = MVar, aes(x = Taxa, y = value))+
    facet_wrap(~Lab, scales = "free_y", labeller = label_parsed)+
    ylab("Variance")+
    geom_segment(aes(xend = Taxa, yend = 0))+
    geom_point()+
    geom_vline(xintercept = Barre1, lty="dashed")+
    geom_vline(xintercept = Barre2, lty="dashed")+
    geom_vline(xintercept = Barre3, lty="dashed")+ Barre.global +
    B0 + B1 + B2 + B3 + B4 +
    theme(panel.background = element_blank(), panel.grid = element_blank(),
          axis.title.x = element_blank(),
          panel.border = element_rect(fill = NA), plot.margin = margin(0,0,0,0),
          strip.text = element_blank(),
          axis.text.y = element_text(size = 7),
          legend.background = element_blank(), strip.background = element_blank(),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1)
          # axis.text.x = element_blank()
          )
  }
    
  if(Plot.gaussian == T){
    #### Color ####
    Col.vec <- c("Artemisia" = "#C88282",
                 "Brassicaceae" = "#ECED8A",
                 "Amaranthaceae" = "#D0C3A7",
                 "Pinus haploxylon" = "#6B9A88",
                 "Pinus diploxylon" = "#3E8A70",
                 # "N/A" = "#FFEAAF",
                 # "Tundra" = "#A9D1C2",
                 "Fabaceae" = "#8FB8E6",
                 # "Tropical & Subtropical Coniferous Forests" = "#99CA81",
                 # "Mangroves" = "#FE01C4",
                 "Ericaceae" = "#BEE7FF"#,
                 # "Tropical & Subtropical Moist Broadleaf Forests" = "#38A700"
    )
    
    Theme.dist <- theme(
      panel.background = element_blank(),
      # legend.position = "none",
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
      panel.border = element_rect(fill = NA),
      strip.text = element_text(size = 12),
      strip.background = element_rect(color = "white", fill = "white"),
    )
    
    #### Plot 2 (normal) ####
    p1 <- ggplot(MT.hist.log, aes(x = value, fill = Taxa)) + 
      facet_wrap(vars(variable), scales = "free")+
      scale_x_continuous(trans = 'log10',
                         breaks = trans_breaks('log10', function(x) 10^x),
                         labels = trans_format('log10', math_format(10^.x)))+
      # scale_fill_manual(values = Col.vec)+
      geom_density(alpha = 0.7, size = 0.1, na.rm = T, outline.type = "full")+
      Theme.dist
    
    #### Plot 2 (log) ####
    p1.log <- ggplot(MT.hist, aes(x = value, fill = Taxa)) + 
      facet_wrap(vars(variable), scales = "free")+
      geom_density(alpha = 0.7, size = 0.1, na.rm = F)+
      # geom_density(alpha = 0.7, size = 0.1, na.rm = T)+
      # scale_fill_manual(values = Col.vec) +
      Theme.dist
    }
  #### Save and export ####
    if(Show.var == F){p <- p1}
    if(Show.var == T){
      p <- p1/p1.var + plot_layout(heights = c(8/10,2/10))
      }
  print(p)
  if(is.null(Save.plot) == F){dev.off()}
  return(p)
  }

LRelation.CWT.clim <- function(CWT, Select.trait, Select.Pclim, Select.eco, Trait.lim, Add.n,
                               Strip.lab, Bit.map, Leg.pos, Facet.scale, Add.linear, Alpha, Save.plot, H, W){
  #### Settings ####
  if(missing(Alpha)){Alpha = 1}
  if(missing(Trait.lim)){Trait.lim = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Select.trait)){Select.trait = NULL}
  if(missing(Select.Pclim)){Select.Pclim = NULL}
  if(missing(Select.eco)){Select.eco = NULL}
  if(missing(Add.linear)){Add.linear = NULL}
  if(missing(Add.linear)){Add.linear = NULL}
  if(missing(Strip.lab)){Strip.lab = T}
  if(missing(Bit.map)){Bit.map = F}
  if(missing(Add.n)){Add.n = F}
  if(missing(Leg.pos)){Leg.pos = "right"}
  if(missing(Facet.scale)){Facet.scale = "free_x"}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
 
  #### Select facets ####
  if(is.null(Select.trait) == T){Select.trait <- grep("TRY", names(CWT))}
  else{Select.trait <- which(names(CWT) %in% Select.trait)}
  
  if(is.null(Select.Pclim) == T){Select.Pclim <- setdiff(seq(1, length(names(CWT))), Trait)}
  else{Select.Pclim <- which(names(CWT) %in% Select.Pclim)}
  
  CWT.m1 <- melt(CWT[c(1, Select.trait)], id = "Site")
  names(CWT.m1) <- c("Site", "Trait", "Trait.Val")
  CWT.m2 <- melt(CWT[c(1, Select.Pclim)], id = "Site")
  names(CWT.m2) <- c("Site", "Clim", "Clim.Val")
  CWT.m <- merge(CWT.m1, CWT.m2, by = "Site")
  
  if(is.null(Select.eco) == T){CWT.m$Eco.Val <- "OK"}
  else{Select.eco <- which(names(CWT) %in% Select.eco)
    CWT.m3 <- melt(CWT[c(1, Select.eco)], id = "Site")
    names(CWT.m3) <- c("Site", "Eco", "Eco.Val")
    CWT.m <- merge(CWT.m, CWT.m3, by = "Site")
    }
  
  n.trait <- nlevels(CWT.m$Trait)
  CWT.m$Trait <- as.factor(gsub("TRY_", "", CWT.m$Trait))
  CWT.m$Clim <- as.factor(gsub("_wc", "", CWT.m$Clim))
  CWT.m$Clim <- as.factor(gsub("_chel", "", CWT.m$Clim))
  
  #### Regression ####
  if(is.null(Add.linear) == F){
    # Add.linear <- geom_smooth(method='lm', se = F, color='turquoise4', size = .7, linetype = "dashed")
    # Add.r2 <- stat_regline_equation(label.y = -2.7, aes(label = ..rr.label..), color = "turquoise4", size = 3.5)
    Add.linear <- stat_poly_line(method='lm', se = T, color='turquoise4', fill = "turquoise4", size = .7, linetype = "dashed")
    # Add.r2 <- stat_poly_eq(label.y = "bottom", label.x = "right", aes(label = ..rr.label..), color = "turquoise4", size = 3.5)
    
    if(Add.n == F){
      Add.r2 <- stat_poly_eq(label.y = "bottom", label.x = "right", color = "turquoise4", size = 3.5, small.r = F,
                             aes(label =  sprintf("%s*\", \"*%s" ,
                                                  after_stat(rr.label),
                                                  # after_stat(r.squared),
                                                  after_stat(p.value.label)
                             )))}
    else{
      Add.r2 <- stat_poly_eq(label.y = "bottom", label.x = "right", color = "turquoise4", size = 3, small.r = F,
                                aes(label =  sprintf("%s*\", \"*%s*\", \"*%s" ,
                                                     after_stat(rr.label),
                                                     # after_stat(r.squared),
                                                     after_stat(p.value.label),
                                                     after_stat(n.label)
                                )))}
    }
  
  #### Color def + lab ####
  Colors.biomes <- c("Deserts & Xeric Shrublands" = "#C88282",
                     "Temperate Grasslands, Savannas & Shrublands" = "#ECED8A",
                     "Montane Grasslands & Shrublands" = "#D0C3A7",
                     "Temperate Conifer Forests" = "#6B9A88",
                     "Temperate Broadleaf & Mixed Forests" = "#3E8A70",
                     "N/A" = "#FFEAAF",
                     "Tundra" = "#A9D1C2",
                     "Boreal Forests/Taiga" = "#8FB8E6",
                     "Tropical & Subtropical Coniferous Forests" = "#99CA81",
                     "Mangroves" = "#FE01C4",
                     "Flooded Grasslands & Savannas" = "#BEE7FF",
                     "Tropical & Subtropical Moist Broadleaf Forests" = "#38A700",
                     "Pollen" = "#bea33a",
                     "Vegetation" = "#6d956f")
  if(Bit.map == F){
    Plot.dot <- geom_point(aes(fill = Eco.Val, color = Eco.Val), size = 1, shape = 16, alpha = Alpha)
    
    Fill.scale <- scale_fill_manual(values = Colors.biomes, name = "Biomes (Dinerstein et al., 2017)")
    Col.scale <- scale_color_manual(values = Colors.biomes, name = "Biomes (Dinerstein et al., 2017)")}
  else{
    Plot.dot <- geom_hex(data = CWT.m, mapping = aes(x = Clim.Val, y = Trait.Val, fill = ..count..), bins = 60, color = NA)
    if(unique(CWT.m$Eco.Val) == "Pollen"){col.scale.hex <- c("#bea03733", "#a07e09f6")}
    if(unique(CWT.m$Eco.Val) == "Vegetation"){col.scale.hex <- c("#6289620d", "#3a723af6")}
    Fill.scale <- scale_fill_gradientn(colors = col.scale.hex,
                                       guide = "legend",
                                       values = scales::rescale(c(1,5,10,30,100,1000), from = c(1,1000)),
                                       # values = scales::rescale(c(1,1000), from = c(1,1000)),
                                       name = "Count")
    Col.scale <- NULL
    }
  
  if(is.null(Trait.lim) == F){Trait.lim <- ylim(Trait.lim)}
  else{Trait.lim <- NULL}
  
  #### Annotations names Strig.lab = F ####
  if(Strip.lab == F){
    Strip.lab.disp <- element_blank()
    S.trait <- setNames(data.frame(as.factor(unique(CWT.m$Trait)), rep(1,nlevels(CWT.m$Trait)), rep(1,nlevels(CWT.m$Trait))), c("Lab", "x","y"))
    S.clim <- setNames(data.frame(as.factor(unique(CWT.m$Clim)), rep(1,nlevels(CWT.m$Clim)), rep(1,nlevels(CWT.m$Clim))), c("Lab", "x","y"))
    
    Theme.null <- theme(axis.line = element_blank(), axis.title = element_blank(),
                        strip.text = element_blank(), axis.text = element_blank(),
                        axis.ticks = element_blank(), plot.background = element_blank(),
                        panel.grid = element_blank(), panel.background = element_blank())
    
    p.up <- ggplot(S.trait, mapping = aes(x = x, y = y))+
      facet_wrap(vars(Lab), scales = "free_x", ncol = n.trait)+
      geom_text(aes(label = Lab))+ Theme.null
    
    p.right <- ggplot(S.clim, mapping = aes(x = x, y = y))+
      facet_wrap(vars(Lab), scales = "free_x", nrow = n.trait)+
      geom_text(aes(label = Lab), angle = 270,  hjust=0.5, vjust=1)+ Theme.null
    }
  else{Strip.lab.disp <- element_text(hjust = 0)}
  
  #### PLOT ####
  p <- ggplot(data = CWT.m, mapping = aes(x = Clim.Val, y = Trait.Val))+
    Plot.dot +
    Add.linear +
    Add.r2 +
    xlab("Climate parameters")+
    ylab("CWM-Traits (z-score)")+
    facet_wrap(Clim ~ Trait, scales = Facet.scale, ncol = n.trait)+ #, ncol=2, scales="free_y", switch = "y") +
    Trait.lim +
    Fill.scale +
    Col.scale +
    #### Theme ####
  theme(
    axis.line = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.background = element_blank(), 
    legend.position = Leg.pos,
    panel.grid = element_blank(), 
    strip.text = Strip.lab.disp,
    strip.background = element_blank(),
    )
  
  if(Strip.lab == F){p <- p.up + plot_spacer() + p + p.right + plot_layout(nrow = 2, heights = c(1/40,39/40), widths = c(39/40,1/40))}
  #### Export ####
  print(p)
  if(is.null(Save.plot) == F){dev.off()}
  return(p)
}

BPlot.CWT.clim <- function(CWT, Select.trait, Select.Pclim, Select.eco, Subgroup, Nuage.point, Only.one.col = F,
                           Alpha.point, Annotate, ANOVA.letters, Save.anova, Save.plot, Legend.size, H, W){
  #### Settings ####
  library(multcompView) # multcompLetters
  if(missing(Save.anova)){Save.anova = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Select.trait)){Select.trait = NULL}
  if(missing(Select.Pclim)){Select.Pclim = NULL}
  if(missing(Select.eco)){Select.eco = NULL}
  if(missing(Subgroup)){Subgroup = NULL}
  if(missing(Annotate)){Annotate = F}
  if(missing(ANOVA.letters)){ANOVA.letters = F}
  if(missing(Nuage.point)){Nuage.point = F}
  if(missing(Legend.size)){Legend.size = 8}
  if(missing(Alpha.point)){Alpha.point = 0.1}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  #### Select facets ####
  if(is.null(Select.trait) == T){Select.trait <- grep("TRY", names(CWT))}
  else{
    A <- list()
    for(i in 1:length(Select.trait)){
      A[[i]] <- grep(Select.trait[i], names(CWT))  
    }
    
    Select.trait <- unlist(A)
    # Select.trait <- c(Select.trait, which(names(CWT) == "Site"))  
    # print(which(names(CWT) == "Site"))
    # Select.trait <- which(names(CWT) %in% Select.trait)
    }
  
  # if(is.null(Select.Pclim) == T){Select.Pclim <- setdiff(seq(1, length(names(CWT))), Trait)}
  if(is.null(Select.Pclim) == T){Select.Pclim <- 1}
  else{Select.Pclim <- which(names(CWT) %in% Select.Pclim)}
  if(is.null(Select.eco) == T){Select.eco <- NULL}
  else{Select.eco <- which(names(CWT) %in% Select.eco)}
  # print(CWT[c(1, Select.trait)])
  CWT.m1 <- melt(CWT[c(1, Select.trait)], id = "Site")
  names(CWT.m1) <- c("Site", "Trait", "Trait.Val")
  CWT.m2 <- melt(CWT[c(1, Select.Pclim)], id = "Site")
  names(CWT.m2) <- c("Site", "Clim", "Clim.Val")
  CWT.m3 <- melt(CWT[c(1, Select.eco)], id = "Site")
  names(CWT.m3) <- c("Site", "Eco", "Eco.Val")
  
  CWT.m <- left_join(CWT.m1, CWT.m2, by = "Site")
  CWT.m <- left_join(CWT.m, CWT.m3, by = "Site")
  CWT.m <- CWT.m[!is.na(CWT.m$Eco.Val),] 
  CWT.m$Subgroup <- NA
  if(is.null(Subgroup) == F){
    CWT.m$Trait <- as.character(CWT.m$Trait)
    for(i in 1:length(Subgroup)){
      CWT.m$Subgroup[grep(Subgroup[i], CWT.m$Trait)] <- Subgroup[i]
      CWT.m$Trait[grep(Subgroup[i], CWT.m$Trait)] <- gsub(Subgroup[i], "", CWT.m$Trait[grep(Subgroup[i], CWT.m$Trait)])
      }
    CWT.m <- CWT.m[!is.na(CWT.m$Subgroup),]
  }
  
  #### Labels ####
  CWT.m$Eco.val.lab <- as.factor(CWT.m$Eco.Val)
  levels(CWT.m$Eco.val.lab) <- c("Boreal Forests/Taiga"                        = "Taiga",# "TAIG",                       
                                 "Deserts & Xeric Shrublands"                  = "Desert",#  "DESE",                 
                                 "Montane Grasslands & Shrublands"             = "Mont. grass.",#  "MONT",            
                                 "Temperate Broadleaf & Mixed Forests"         = "Mix. Forest",#  "TEBF",        
                                 "Temperate Conifer Forests"                   = "Con. Forest",#  "TECF",                  
                                 "Temperate Grasslands, Savannas & Shrublands" = "Steppe",# "STEP",
                                 "Tundra"                                      = "TUND") 
    
  CWT.m$Trait <- gsub("TRY_", "", CWT.m$Trait)
  
  #### Parameters label ####
  My_lab <- c("Height[CWM]", "LeafArea[CWM]", "N[leaf-CWM]", "SeedMass[CWM]", "SLA[CWM]", "SSD[CWM]")
  CWT.m$Trait[CWT.m$Trait == "Height"] = My_lab[1]
  CWT.m$Trait[CWT.m$Trait == "LeafArea"] = My_lab[2]
  CWT.m$Trait[CWT.m$Trait == "LeafN"] = My_lab[3]
  CWT.m$Trait[CWT.m$Trait == "SeedMass"] = My_lab[4]
  CWT.m$Trait[CWT.m$Trait == "SLA"] = My_lab[5]
  CWT.m$Trait[CWT.m$Trait == "SSD"] = My_lab[6]
  # CWT.m$Trait <- factor(CWT.m$Trait, levels = My_lab, ordered = T)
  
  #### Color and fill settings ####
  Value.bi <- c(
                 "Boreal Forests/Taiga" = "#8FB8E6",
                 "Deserts & Xeric Shrublands" = "#C88282",
                 "Montane Grasslands & Shrublands" = "#D0C3A7",
                 "Temperate Broadleaf & Mixed Forests" = "#3E8A70",
                 "Temperate Conifer Forests" = "#6B9A88",
                 "Temperate Grasslands, Savannas & Shrublands" = "#ECED8A",
                 "N/A" = "#FFEAAF",
                 "Tundra" = "#A9D1C2",
                 "Tropical & Subtropical Coniferous Forests" = "#99CA81",
                 "Mangroves" = "#FE01C4",
                 "Flooded Grasslands & Savannas" = "#BEE7FF",
                 "Tropical & Subtropical Moist Broadleaf Forests" = "#38A700")
  
  if(Only.one.col == F){
    Value.bi2 <- c("_ACAV" = "#e2e2bd81", "_ACAV_gf" = "#baba7250",
                    "_sl" = "grey50", "_sl_gf" = "black", "_coarse_gf" = "black",
                    "_ss" = "#af8282ff", "_fine_gf" = "#640f0fff", "_ss_gf" = "#640f0fff")
    show.lab = NULL}
  else{
    One.col1 = "grey10"
    One.col2 = "grey10"
    Value.bi2 <- c("_ACAV" = One.col1, "_ACAV_gf" = One.col1, "_1_ACAV_gf" = One.col1,
                    "_sl" = One.col2, "_sl_gf" = One.col2, "_coarse_gf" = One.col2,  "_3_coarse_gf" = One.col2,
                    "_ss" = One.col2, "_fine_gf" = One.col2,  "_2_fine_gf" = One.col2, "_ss_gf" = One.col2)
    show.lab = guides(color = "none")}
  
  
  
  Value.bi <- Value.bi[which(names(Value.bi) %in% unique(CWT.m$Eco.Val))]
  Value.bi2 <- Value.bi2[which(names(Value.bi2) %in% unique(CWT.m$Subgroup))]
  
  
  My_fill <-  scale_fill_manual(values = Value.bi, name = "") # "Biomes (Dinerstein et al., 2017)"#, labels = levels(ACA.biom.proj$BIOME_NAME)
  My_col <- scale_colour_manual(values = Value.bi2, name = "", label = c("_sl" = "ACASP-sl",
                                                                       "_sl_gf" = "ACASP-sl (gap-filled)",
                                                                       "_coarse_gf" = "ACASP-coarse (gap-filled)",
                                                                       "_ss" = "ACASP-ss",
                                                                       "_fine_gf" = "ACASP-fine (gap-filled)",
                                                                       "_ss_gf" = "ACASP-ss (gap-filled)",
                                                                       "_ACAV" = "ACAV", "_1_ACAV" = "ACAV",
                                                                       "_ACAV_gf" = "ACAV (gap-filled)"))
  #### Annotate ####
    library(plyr)
    CWT.m.stat <- CWT.m[!is.na(CWT.m$Trait.Val),]
    # CWT.m.stat$Divis <- length(Select.Pclim)+length(Subgroup)
    mms.cor <- ddply(.data=CWT.m.stat, 
                     .(Trait), 
                     summarize, 
                     # n = paste("n =", round(length(Trait.Val)/Divis, digit = 1)),
                     n = paste("n =", length(Trait.Val)),
                     max = max(Trait.Val))
    mms.cor$num <- paste("(", toupper(letters[1:nrow(mms.cor)]) , ")", sep = "")
    
  if(Annotate == T){Note.n <- geom_text(inherit.aes=FALSE, data = mms.cor, aes(x = length(unique(CWT.m$Eco.Val)), y = max, label=n), colour="grey30", parse=FALSE, hjust = 0.5)}
  else{Note.n <- NULL}
    
  if(Nuage.point == T){
    Nuage.point <- geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = Alpha.point, na.rm = T)
    # Nuage.point <- stat_binhex(position = position_jitterdodge(jitter.width = 0.1), bins = 200, fill = NA, alpha = Alpha.point)
    # Nuage.point <- geom_rug(position = position_jitterdodge(jitter.width = 0.1), bins = 200, fill = NA, alpha = Alpha.point)
    }
  else{Nuage.point <- NULL}
  
    
  #### PLOT ####
  # print(CWT.m$Trait)
  p <- ggplot(data = CWT.m, aes(x = Eco.val.lab, y = Trait.Val, fill = Eco.Val, colour = Subgroup))+
    Nuage.point +
    geom_boxplot(outlier.colour = NA, outlier.alpha = 0.1, outlier.shape = 16, outlier.size = 1,
                 alpha = 0.9, notch = T, notchwidth = 0.7, varwidth = F, na.rm = T, show.legend = T) +
    Note.n +
    geom_text(data=mms.cor, aes(x = 1, y = max, label=num), colour = "black", size = 4, inherit.aes=FALSE, parse=FALSE, hjust = 1, vjust = 1)+
    xlab("Biomes")+
    ylab("Trait values (z-scores)")+
    facet_wrap(.~Trait, scales = "free_y", labeller = label_parsed)+#, ncol=2, scales="free_y", switch = "y") +
    My_fill + My_col +
    guides(fill = guide_legend(nrow = 2), colour = guide_legend(nrow = 2))+
    #### Theme ####
  theme(
    panel.background = element_blank(),
    legend.position = "bottom", panel.grid = element_blank(),
    axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1),
    legend.text = element_text(size = Legend.size), legend.background = element_blank(),
    legend.key = element_blank(),
    panel.border = element_rect(fill = NA),
    strip.text = element_text(size = 12), 
    strip.background = element_blank()
  )
  
    #### Anova ####
    if(is.null(Save.anova) == F){
      #### Graphs Settings plot ANOVA #### 
      # W = 800
      # H = 1000
      # Save.plot.anova <- gsub("\\.pdf", "_anova.pdf", Save.plot) 
      # pdf(file = Save.plot.anova, width = W*0.01041666666667, height = H*0.01041666666667)
      
      #### ANOVA calculation #### 
      Save.two.ways <- setNames(data.frame(matrix(NA, nrow = length(Subgroup), ncol = length(unique(CWT.m$Trait)))), unique(CWT.m$Trait))
      # ATTENTION PROBLEM: si le nombre de biomes change, il faudra changer le 15
      Save.two.ways.biom <- setNames(data.frame(matrix(NA, nrow = 15, ncol = length(unique(CWT.m$Trait)))), unique(CWT.m$Trait))
      Save.anova.res <- setNames(data.frame(matrix(ncol = length(unique(CWT.m$Trait)), nrow = 2)), unique(CWT.m$Trait))
      row.names(Save.anova.res) <- c("Sample type", "Biomes")
      Boxes <- ggplot_build(p)$data[[2]]
      
      for(i in unique(CWT.m$Trait)){
        #### Tukey test #### 
        A <- CWT.m[CWT.m$Trait == i,]
        A <- A[c(3,7,8,9)]
        # Two.way <- aov(Trait.Val ~ Subgroup + Eco.val.lab, data = A)
        Two.way <- aov(Trait.Val ~ Subgroup*Eco.val.lab, data = A)
        tukey.two.way <- TukeyHSD(Two.way, conf.level=.95)
        # tukey.two.way2 <- TukeyHSD(Two.way2, conf.level=.95)
        print(paste("Trait étudié:", i))
        # print(names(tukey.two.way))
        Anov.res <- data.frame(summary(Two.way)[[1]])[1:2,1:5]
        # print(Anov.res)
        Anov.res$ANOVA <- paste("DF = ", Anov.res$Df, 
                                ", F = ", round(Anov.res$F.value, digits = 2),
                                ", p = ", round(Anov.res$Pr..F., digits = 4), sep = "")
        Anov.res$ANOVA[grepl("p = 0$*", Anov.res$ANOVA)] <- gsub("p = 0", "p < 0.001", Anov.res$ANOVA[grepl("p = 0$*", Anov.res$ANOVA)])
        
        #### Tukey plot ####
        # par(mfrow=c(1,2))
        # # plot(Two.way$)
        # plot(tukey.two.way, las = 2, cex.axis = 0.6)
        # title(i)
        # par(mfrow=c(1,1))
        
        
        #### Generation des lettres ####
        if(ANOVA.letters == T){
          #### Generate letters ####
          generate_label_df <- function(TUKEY, variable){
            # Extract labels and factor levels from Tukey post-hoc 
            Tukey.levels <- TUKEY[[variable]][,4]
            Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
            #I need to put the labels in the same order as in the boxplot :
            Tukey.labels$treatment=rownames(Tukey.labels)
            Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
            return(Tukey.labels)
          }
          My_lab <- generate_label_df(tukey.two.way, "Subgroup:Eco.val.lab")
          names(My_lab) <- c('Letters','Eco.val.lab')
          yvalue <- aggregate(Trait.Val ~ Subgroup*Eco.val.lab, data = A, quantile, probs=1)
          row.names(yvalue) <- paste(yvalue$Subgroup, yvalue$Eco.val.lab, sep = ":")
          final <- merge(My_lab, yvalue, by = "row.names")
          final$Trait <-  i
          Corresp.eco.val <- unique(A[c("Eco.Val", "Eco.val.lab")])
          final$Eco.Val <- Corresp.eco.val$Eco.Val[match(final$Eco.val.lab.y, Corresp.eco.val$Eco.val.lab)]
          Boxes.i <- Boxes[Boxes$PANEL == which(unique(CWT.m$Trait) %in% i), c(11,3,7,15,6)]
          final <- cbind(final[order(final[[8]]),], Boxes.i)
    
          #### plot letters ####
          p <- p + geom_text(data = final, aes(x = x,y = ymax,label = Letters),
                             size = 3.5,vjust = -.5,show.legend = F) + show.lab
      }
      #### Save ####
      Save.anova.res[names(Save.two.ways.biom) == i] <- Anov.res$ANOVA
      Save.two.ways.biom[names(Save.two.ways.biom) == i] <- round(as.vector(tukey.two.way$Eco.val.lab[,4]), digits = 5)
      Save.two.ways[names(Save.two.ways) == i] <- round(as.vector(tukey.two.way$Subgroup[,4]), digits = 5)
      }
    Save.anova.res <- data.frame(t(Save.anova.res))
    
    #### Export ANOVA results & plots ####
    row.names(Save.two.ways.biom) <- row.names(tukey.two.way$Eco.val.lab)
    row.names(Save.two.ways) <- row.names(tukey.two.way$Subgroup)
    Save.two.ways.full <- rbind(Save.two.ways, Save.two.ways.biom)
    
    Path.to.create <- gsub("(.*/).*\\.csv.*","\\1",Save.anova)
    dir.create(file.path(Path.to.create), showWarnings = F)
    Save.path.Mfull <- gsub("\\.csv", "_tukey_subgroups.csv", Save.anova)
    Save.path.Mfull2 <- gsub("\\.csv", "_tukey_biomes.csv", Save.anova)
    Save.path.Mfull3 <- gsub("\\.csv", "_tukey_full.csv", Save.anova)
    write.table(Save.anova.res, file = Save.anova, row.names=T, col.names=NA, sep=",", dec = ".")
    write.table(Save.two.ways, file = Save.path.Mfull, row.names=T, col.names=NA, sep=",", dec = ".")
    write.table(Save.two.ways.biom, file = Save.path.Mfull2, row.names=T, col.names=NA, sep=",", dec = ".")
    write.table(Save.two.ways.full, file = Save.path.Mfull3, row.names=T, col.names=NA, sep=",", dec = ".")
    # dev.off()
  }

  #### Export ####
  print(p)
  if(is.null(Save.plot) == F){dev.off()}
  }

Mat.corel.CWT.clim <- function(Mclim, Mtrait, I.confiance, Display.pval, Disp.R, Display,
                      Title, Label, Save.path, return.pick, Bar.pos, my_color, Use.cor, Xlab.pos,
                      Permutation.test = F, Nb.permutations = 1000, Save.pval = NULL,
                      Method.cor, Save.plot, W, H, Order.by.alphabet, Average){
  #### Settings ####
  Save.tab = T
  library(corrplot)
  if(missing(Bar.pos)){Bar.pos = "b"}
  if(missing(Display)){Display = "full"}
  if(missing(I.confiance)){I.confiance = 0}
  if(missing(return.pick)){return.pick = F}
  if(missing(Save.path)){Save.tab = F}
  if(missing(Order.by.alphabet)){Order.by.alphabet = F}
  if(missing(Label)){Label = T}
  if(missing(Average)){Average = T}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Xlab.pos)){Xlab.pos = NULL}
  if(missing(Method.cor)){Method.cor = "pearson"}
  if(missing(Use.cor)){Use.cor = "pairwise.complete.obs"}
  if(missing(Display.pval)){Display.pval = "pch"}
  if(missing(Disp.R)){Disp.R = "pie"}
  if(missing(Mtrait)){Mtrait = NULL}
  if(missing(Title)){Title = NULL}
  if(missing(my_color)){my_color = colorRampPalette(c("royalblue", "white", "darkorange"))(50)}
  if(Display.pval == "blank"){Subt = paste("*The cells are blanked when ever the p-value is over the", I.confiance, "significant level.")}
  if(Display.pval == "pch"){Subt = paste("*The cells are crossed when ever the p-value is over the", I.confiance, "significant level.")}
  if(Label == F){Subt <- NULL}
  #### Calculations ####
  names(Mtrait) <- gsub("TRY_", "", names(Mtrait))
  
  if(Order.by.alphabet == T){Mtrait <- Mtrait[order(names(Mtrait))]}
  
  Mcor <- cor(Mclim, Mtrait, use = Use.cor, method = Method.cor)
  Mvar <- var(Mclim, Mtrait)
  Mtot <- cbind(Mclim, Mtrait)
  Resid <- cor.mtest(Mtot, conf.level = I.confiance)
  PV <- Resid$p
  
  if(is.null(Mtrait) == F){PV <- PV[(ncol(Mtrait)+1):(ncol(Mtrait)+ncol(Mclim)), 1:ncol(Mtrait)]} 
  row.names(PV) <- row.names(Mcor)
  colnames(PV) <- colnames(Mcor)
  
  #### Permutations test ####
  if(Permutation.test == T){
    library("jmuOutlier")
    print(paste("p-value calculated with permutation test made on ", Nb.permutations, " permutations.", sep = ""))
    PV.perm <- setNames(data.frame(matrix(NA, nrow = length(Mclim), length(Mtrait))), names(Mtrait))
    row.names(PV.perm) <- names(Mclim) 
    for(i in 1:length(names(Mclim))){
      for(j in 1:length(names(Mtrait))){
        set.seed(0)
        M.perm <- na.omit(cbind(Mclim[i], Mtrait[j]))
        PV.perm[i, j] <- perm.cor.test(M.perm[[1]], M.perm[[2]], num.sim = Nb.permutations)$p.value
        }
      }
    PV <- as.matrix(PV.perm)
    }

  #### Save p-values ####
  if(is.null(Save.pval) == F){
    Path.to.create <- gsub("(.*/).*\\.csv.*","\\1", Save.pval)
    dir.create(file.path(Path.to.create), showWarnings = F)
    write.table(PV, file = Save.pval, row.names=T, col.names = NA, sep=",", dec = ".")
    Save.path.RDS <- gsub("\\.csv", "\\.Rds", Save.pval)
    saveRDS(PV, Save.path.RDS)
    }
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Graphic param ####
  # my_color <- colorRampPalette(c("royalblue","white", "#ee4700ff"))(50)
  #### Plots ####
  CP <- corrplot(Mcor, type = Display, col = my_color, tl.pos = Xlab.pos,
                 tl.col="black", tl.srt=45, tl.cex = .7, method = Disp.R, cl.align.text = "l", cl.pos = Bar.pos, 
                 p.mat = PV, sig.level = (1-I.confiance), insig = Display.pval, pch.cex = 2)
  
  title(main = Title, sub = Subt)
  
  #### Average values ####
  if(Average == T){
    type.correl <- gsub(".*ACA", "", Title)
    type.correl <- gsub("\\s.*", "", type.correl)
    
    MeanR2 <- round((Mcor)^2, digits = 2)
    MeanR2 <- cbind(MeanR2, MeansCol = rowMeans(MeanR2, na.rm = T))
    MeanR2 <- rbind(MeanR2, MeansRow = colMeans(MeanR2, na.rm = T))
    
    print(MeanR2)
    
    R.mean.ACAV.t.chel <- round(mean(MeanR2["MeansRow", grepl("_chel", colnames(MeanR2))], na.rm = T), digits = 2)
    R.mean.ACAV.t.wc <- round(mean(MeanR2["MeansRow", grepl("_wc", colnames(MeanR2))]), digits = 2)
    R.mean.ACAV.gl.clim <- round(mean(MeanR2[grepl("_gf", rownames(MeanR2)), "MeansCol"]), digits = 2)
    R.mean.ACAV.clim <- round(mean(MeanR2[!grepl("_gf", rownames(MeanR2)), "MeansCol"]), digits = 2)

    print(paste("R² mean ACA-", type.correl, " vs. Climat-CHELSA:", R.mean.ACAV.t.chel, sep = ""))
    print(paste("R² mean ACA-", type.correl, " vs. Climat-WORLDCLIM:", R.mean.ACAV.t.wc, sep = ""))
    print(paste("R² mean ACA-", type.correl, "-NORMAL vs. Climat :", R.mean.ACAV.clim, sep = ""))
    print(paste("R² mean ACA-", type.correl, "-GAPFILLING vs. Climat :", R.mean.ACAV.gl.clim, sep = ""))
    }
  
  #### Export datas ####
  Mfull <- list(R2 = Mcor, Var = Mvar, p.val = PV)
  if(Save.tab == T){
    Path.to.create <- gsub("(.*/).*\\.csv.*","\\1", Save.path)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    write.table(Mfull, file = Save.path, row.names=T, col.names=NA, sep=",", dec = ".")}
  
  #### Return ####
  if(is.null(Save.plot) == F){dev.off()}
  if(return.pick == T){return(CP)}
  else{return(Mfull)}
}

Mat.corel.network <- function(M, Method, Min.R, Save.plot, Type.graph, Color.scale.bar, Title, W, H){
  #### Library ####
  if(missing(Method)){Method <- "pearson"}
  if(missing(Type.graph)){Type.graph <- "kk"}
  if(missing(Min.R)){Min.R <- 0.2}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(Color.scale.bar)){Color.scale.bar = T}
  if(missing(Title)){Title = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  library(RColorBrewer)
  library(igraph)
  library(corrr)
  library(ggraph)
  #library(hrbrthemes)
  #library(dplyr)
  #library(tidyverse) # mutate
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  
  #### Vertices tentation ####
  print(M)
  Edge <- list(CWM = c("SLA", "SeedMass", "Height", "LeafArea", "LeafN", "SSD"),
                     Temp = c("MAAT_wc", "MTCOQ_wc", "MTWAQ_wc"),
                     Prec = c("AI", "MAP_wc", "MPCOQ_wc", "MPWAQ_wc"))
  Edge <- melt(unlist(Edge))
  Edge <- data.frame(cbind(from = row.names(Edge), to = Edge$value))
  Edge$from <- gsub("[0-9]", "", Edge$from)
  print(Edge)
  
  #### Correlation calculation ####
  # res.cor <- correlate(M, use = "complete.obs") 
  res.cor <- correlate(M, use = "pairwise.complete.obs") 
  res.cor <- reshape2::melt(res.cor, id = "term")
  names(res.cor)[3] <- "r"
  
  # print(res.cor)
  
  graph_cors <- res.cor %>%
    filter(abs(r) > Min.R) %>%
    graph_from_data_frame(directed = FALSE)

  
  # graph_cors <- res.cor %>%
  #   filter(abs(r) > Min.R) #%>%
  #   #graph_from_data_frame(directed = T)
  # 
  # print(graph_cors)
  # graph_cors <- graph_from_data_frame(Edge, vertice = graph_cors)
  # print(graph_cors)
  
  #### Color settings ####
  Param.type <- levels(res.cor$variable)
  Param.type[grep("TRY", Param.type)] <- "Trait"
  print(Param.type)
  Param.type[setdiff(seq(1, length(Param.type)),grep("Trait", Param.type))] <- "clim"
  Param.type <- as.factor(Param.type)
  coul <- c("grey20", "grey60")
  coul2 <- c("grey80", "grey10")
  my_color <- coul[as.numeric(Param.type)]
  my_color_lab <- coul2[as.numeric(Param.type)]
  if(is.null(Title) == F){Title <- ggtitle(Title)}
  
  if(Color.scale.bar == T){Color.scale.bar.pos <- "bottom"}
  else{Color.scale.bar.pos <- "none"}
  
  #### Plot ####
  my_layout <- create_layout(graph_cors, layout = Type.graph, circular = T) # (kk) Kamada-Kawai force-directed algorithm
  my_layout$name <- gsub("TRY_", "", my_layout$name)
  my_layout$name <- gsub("_wc", "", my_layout$name)
  my_layout$name <- gsub("_chel", "", my_layout$name)
  my_layout$col <- my_layout$name
  my_layout$col[which(my_layout$col %in% c("Altitude", "Latitude", "Longitude"))] <- "#5b5b5b"
  my_layout$col[which(my_layout$col %in% c("SLA", "SeedMass", "Height", "LeafArea", "LeafN", "SSD"))] <- "#0e5a25cc"
  my_layout$col[which(my_layout$col %in% c("AI", "MAP", "MPCOQ", "MPWAQ", "MPCOQ_chel", "MPWAQ_chel", "MAP_chel"))] <- "#181650cc"
  my_layout$col[which(my_layout$col %in% c("MAAT", "MTCOQ", "MTWAQ", "MAAT_chel", "MTWAQ_chel", "MTCOQ_chel"))] <- "#761023cc"
  # print(my_layout)
  
  p <- ggraph(my_layout) +           
  # p <- ggraph(graph_cors, layout = 'eigen') +               # layout = "linear" / "kk", circular = T/F 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl', 'lgl'.
    # geom_edge_density(aes(fill = r)) +
    geom_node_voronoi(aes(fill = col), max.radius = 0.5, colour = 'white', alpha = 0.2) + 
    # geom_edge_arc(aes(colour = r, edge_width = abs(r), alpha = ..index..))+
    geom_edge_hive(aes(edge_alpha = abs(r), edge_width = abs(r), color = r, alpha = ..index..), strength = 0.5) +
    scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("royalblue", "white", "darkorange"), name = "Pearson R") +
    scale_edge_width_continuous(range = c(0,12))+
    scale_color_manual(limits = as.factor(my_layout$col), values = my_layout$col, name = NULL, labels = NULL, breaks = NULL) +
    scale_fill_manual(limits = as.factor(my_layout$col), values = my_layout$col, name = NULL, labels = NULL, breaks = NULL) +
    geom_node_point(aes(color = col), alpha = 0.8, size = 30, fill = "grey") +
    geom_node_text(aes(label = name), repel = F, color = "grey90", size = 5, fontface = "bold") +
    guides(edge_alpha = "none", edge_width = "none", node_point = NULL) +
    # guides(edge_alpha = "none", edge_width = "none", color = T) +
    # coord_fixed() +
    Title +
    theme_void()+
    theme(title = element_text(size = 20),
          legend.position = Color.scale.bar.pos,
          # plot.margin=unit(c(10,10,10,10),"pt"),
          # panel.background = element_blank(),
          legend.direction = "horizontal")
  
  #### Export ####
  print(p)
  if(is.null(Save.plot) == F){dev.off()}
  return(p)
}

SSM.CWT.clim <- function(MG, MC, Full.stats, Nb.max, Adj.R, Save.path){
  #### Settles ####
  if(missing(Save.path)){Save.path = NULL}
  if(missing(Full.stats)){Full.stats = F}
  if(missing(Nb.max)){Nb.max = 3}
  if(missing(Adj.R)){Adj.R = F}
  x <- names(MG)
  y <- names(MC)
  
  Export.param <- list()
  Tab.best <- data.frame(matrix(0, ncol = 6, nrow = 0))
  names(Tab.best) = c("Param.Clim", "NB.param", "R2", "RMSE", "AIC", "Formula")
  if(Adj.R == T){
    Tab.best <- data.frame(matrix(0, ncol = 7, nrow = 0))
    names(Tab.best) = c("Param.Clim", "NB.param", "R2", "RMSE", "AIC", "Formula", "Adj_R2")}
  
  #### Counting number of parameters for each models ####
  Nb.param.SSM <- function(M){
    chars <- unlist(strsplit(M,""))
    N <- length(chars[chars == "+"])
    N <- N+1
    return(N)
  }
  
  #### Loops on climate parameters ####
  for(i in 1:length(y)){
    #### Extract stats param ####
    print(paste("SSM calcul", round(i/length(y), digits = 2)*100, "%"))
    MG <- cbind(MG, MC[i])
    comb <- lapply(1:length(x), function(nvar) do.call(rbind, combn(x, nvar, simplify = F)))
    comb <- lapply(comb, function(..) apply(.., 1, function(.) paste(paste(y[i], "~"), paste(., collapse = " + "))))
    comb <- do.call(c, comb)
    
    ols  <- lapply(comb, function(.) lm(., data = MG))
    r2   <- sapply(ols, function(.) summary(.)$r.squared)
    if(Adj.R == T){Adj_R   <- sapply(ols, function(.) summary(.)$adj.r.squared)}
    All.occ <- sapply(comb, function(.) Nb.param.SSM(.))
    AIC.model <- sapply(ols, function(.) AIC(., k = 2))
    
    #### Calcul RMSE depuis le Linear Model ####
    Res   <- sapply(ols, function(.) summary(.)$residuals)
    RSS   <- sapply(ols, function(.) crossprod(summary(.)$residuals))
    MSE   <- RSS / nrow(Res)
    RMSE <- sqrt(MSE)        # The Root
    
    #### Export all param ####
    if(Full.stats == T){
      Full.s <- list(Models = comb, R2 = r2, RMSE = RMSE, Nparam = All.occ, AIC = AIC.model)
      if(Adj.R == T){Full.s <- list(Models = comb, R2 = r2, Adj_R = Adj_R, RMSE = RMSE, Nparam = All.occ, AIC = AIC.model)}
      Export.param <- append(Export.param, list(Full.s))
      names(Export.param)[i] = names(MC)[i]
    }
    
    #### Export best models ####
    if(Full.stats == F){
      #### Selection of the best models for all the number of parameters ####
      for(j in 1:length(x)){
        TT <- nrow(Tab.best)
        Mat.param.equal <- which(All.occ == j)
        R2.PE <- r2[Mat.param.equal]
        if(Adj.R == T){Adj_R.PE <- Adj_R[Mat.param.equal]}
        AIC.PE <- AIC.model[Mat.param.equal]
        RMSE.PE <- RMSE[Mat.param.equal]
        NAMES.PE <- All.occ[Mat.param.equal]
        
        #### Selection on the Nb.max number of best ####
        Best.R2 <- tail(sort(R2.PE), Nb.max)
        Index.best.R2 <- which(R2.PE %in% Best.R2)
        Best.RMSE <- RMSE.PE[Index.best.R2]
        Best.AIC <- AIC.PE[Index.best.R2]
        Best.formula <- NAMES.PE[Index.best.R2]
        if(Adj.R == T){Best.Adj_R2 <- Adj_R.PE[Index.best.R2]}
        
        #### Tableau synthese ####
        for(l in 1:length(Best.R2)){
          k <- l + TT
          Tab.best[k,1] <- y[i]
          Tab.best[k,2] <- j
          Tab.best[k,3] <- Best.R2[l]
          Tab.best[k,4] <- Best.RMSE[l]
          Tab.best[k,5] <- Best.AIC[l]
          Tab.best[k,6] <- names(Best.formula)[l]
          if(Adj.R == T){Tab.best[k,7] <- Best.Adj_R2[l]}
        }
      }
      
      Export.param <- Tab.best  
    }}
  
  #### Export datas ####
  if(is.null(Save.path) == F){
    Path.to.create <- gsub("(.*/).*\\.csv.*","\\1",Save.path)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    Save.path.Mfull <- gsub("\\.csv", "_SSM.csv", Save.path)
    write.table(Export.param, file = Save.path.Mfull, row.names=T, col.names=NA, sep=",", dec = ".")}
  
  return(Export.param)
}

SSM.CWT.biom <- function(MG, MB, Method, Full.stats, Nb.max, Adj.R, Save.path){
  #### Settles ####
  library(randomForest)
  if(missing(Method)){Method = "RF"}
  if(missing(Save.path)){Save.path = NULL}
  if(missing(Full.stats)){Full.stats = F}
  if(missing(Nb.max)){Nb.max = 3}
  if(missing(Adj.R)){Adj.R = F}
  names(MB)[ncol(MB)] <- "Species"
  y <- names(MB)
  
  
  #### Save list settle ####
  Export.param <- list()
  Tab.best <- data.frame(matrix(0, ncol = 6, nrow = 0))
  names(Tab.best) = c("Param.Clim", "NB.param", "R2", "RMSE", "AIC", "Formula")
  if(Adj.R == T){
    Tab.best <- data.frame(matrix(0, ncol = 7, nrow = 0))
    names(Tab.best) = c("Param.Clim", "NB.param", "R2", "RMSE", "AIC", "Formula", "Adj_R2")}
  
  #### Counting number of parameters for each models ####
  Nb.param.SSM <- function(M){
    chars <- unlist(strsplit(M,""))
    N <- length(chars[chars == "+"])
    N <- N+1
    return(N)
  }
  
  #### Loops on climate parameters ####
  for(i in 1:length(y)){
    #### Extract stats param ####
    print(paste("SSM calcul", round(i/length(y), digits = 2)*100, "%"))
    
    x <- names(MG)
    MG <- cbind(MG, MB[i])
    Keep.species <- names(which(table(MG$Species)>5))
    MG <- MG[which(MG$Species %in% Keep.species),]
    MG$Species = as.factor(MG$Species)
    MG$Species <- droplevels(MG$Species)
    
    #### Splitting train and test data ####
    set.seed(1234)
    data1 <- sample(2, nrow(MG), 
                    replace = T, 
                    prob = c(0.8, 0.2))
    train_cl <- MG[data1 == 1,]
    test_cl <- MG[data1 == 2,]
    
    #### Feature Scaling / standardization #### 
    train_cl[, 1:(ncol(train_cl)-1)] <- scale(train_cl[, 1:(ncol(train_cl)-1)])
    test_cl[, 1:(ncol(test_cl)-1)] <- scale(test_cl[, 1:(ncol(test_cl)-1)])
    
    #### Comb ####
    comb <- lapply(1:length(x), function(nvar) do.call(rbind, combn(x, nvar, simplify = F)))
    comb <- lapply(comb, function(..) apply(.., 1, function(.) paste(paste(y[i], "~"), paste(., collapse = " + "))))
    comb <- do.call(c, comb)
    
    # ols  <- 
    # ols  <- lapply(comb, function(.) lm(., data = MG))
    if(Method == "RF"){
      ols <- lapply(comb, function(.) randomForest(Species ~ ., data = train_cl, ntree = 1000, na.action = na.omit))
      y_pred   <- lapply(ols, function(.) predict(., newdata = test_cl))
      } # Random Forest
    
    #### Confusion Matrix ####
    cm <- lapply(y_pred, function(.) table(test_cl$Species, .))
    df.cm <- lapply(cm, function(.) confusionMatrix(.)$table)
    
    #convert confusion matrices to tables, and binding them together
    df.cm.col <- lapply(df.cm, function(.) . / rowSums(.))
    
    df.table <- lapply(df.cm, function(.) reshape2::melt(.))
    df.table.col <- lapply(df.cm.col, function(.) reshape2::melt(.))
    
    # /!\ bloqué 
    df.table <- lapply(df.table, function(.) left_join(., df.table.col, by =c("Var1", "y_pred")))
    return(df.table)
    
    #calculate accuracy and class accuracy
    acc.vector <- c(diag(df.cm)) / c(rowSums(df.cm))
    class.acc <- data.frame(y_pred = "Biome Acc.", Var1 = names(acc.vector), value = acc.vector)
    acc <- sum(diag(df.cm)) / sum(df.cm)
    
    r2   <- sapply(ols, function(.) summary(.)$r.squared)
    if(Adj.R == T){Adj_R   <- sapply(ols, function(.) summary(.)$adj.r.squared)}
    All.occ <- sapply(comb, function(.) Nb.param.SSM(.))
    AIC.model <- sapply(ols, function(.) AIC(., k = 2))
    
    #### Calcul RMSE depuis le Linear Model ####
    Res   <- sapply(ols, function(.) summary(.)$residuals)
    RSS   <- sapply(ols, function(.) crossprod(summary(.)$residuals))
    MSE   <- RSS / nrow(Res)
    RMSE <- sqrt(MSE)        # The Root
    
    #### Export all param ####
    if(Full.stats == T){
      Full.s <- list(Models = comb, R2 = r2, RMSE = RMSE, Nparam = All.occ, AIC = AIC.model)
      if(Adj.R == T){Full.s <- list(Models = comb, R2 = r2, Adj_R = Adj_R, RMSE = RMSE, Nparam = All.occ, AIC = AIC.model)}
      Export.param <- append(Export.param, list(Full.s))
      names(Export.param)[i] = names(MB)[i]
    }
    
    #### Export best models ####
    if(Full.stats == F){
      #### Selection of the best models for all the number of parameters ####
      for(j in 1:length(x)){
        TT <- nrow(Tab.best)
        Mat.param.equal <- which(All.occ == j)
        R2.PE <- r2[Mat.param.equal]
        if(Adj.R == T){Adj_R.PE <- Adj_R[Mat.param.equal]}
        AIC.PE <- AIC.model[Mat.param.equal]
        RMSE.PE <- RMSE[Mat.param.equal]
        NAMES.PE <- All.occ[Mat.param.equal]
        
        #### Selection on the Nb.max number of best ####
        Best.R2 <- tail(sort(R2.PE), Nb.max)
        Index.best.R2 <- which(R2.PE %in% Best.R2)
        Best.RMSE <- RMSE.PE[Index.best.R2]
        Best.AIC <- AIC.PE[Index.best.R2]
        Best.formula <- NAMES.PE[Index.best.R2]
        if(Adj.R == T){Best.Adj_R2 <- Adj_R.PE[Index.best.R2]}
        
        #### Tableau synthese ####
        for(l in 1:length(Best.R2)){
          k <- l + TT
          Tab.best[k,1] <- y[i]
          Tab.best[k,2] <- j
          Tab.best[k,3] <- Best.R2[l]
          Tab.best[k,4] <- Best.RMSE[l]
          Tab.best[k,5] <- Best.AIC[l]
          Tab.best[k,6] <- names(Best.formula)[l]
          if(Adj.R == T){Tab.best[k,7] <- Best.Adj_R2[l]}
        }
      }
      
      Export.param <- Tab.best  
    }}
  
  #### Export datas ####
  if(is.null(Save.path) == F){
    Path.to.create <- gsub("(.*/).*\\.csv.*","\\1",Save.path)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    Save.path.Mfull <- gsub("\\.csv", "_SSM.csv", Save.path)
    write.table(Export.param, file = Save.path.Mfull, row.names=T, col.names=NA, sep=",", dec = ".")}
  
  return(Export.param)
}

Map.biogeo.CWM <- function(MCWT, MCWT2, Select.trait, Type1, Type2, Crop.zone, Leg.pos,
                           Vertical, Show.diff, Strip.lab, Hex.size, H, W, Save.plot){
  #### Settings ####
  library(patchwork)
  # library(hexbin)
  # library(ggraph)
  if(missing(Hex.size)){Hex.size = 40}
  if(missing(MCWT2)){MCWT2 = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Select.trait)){Select.trait = NULL}
  if(missing(Crop.zone)){Crop.zone = NULL}
  if(missing(Leg.pos)){Leg.pos = "bottom"}
  if(missing(Type1)){Type1 = "Type1"}
  if(missing(Type2)){Type2 = "Type2"}
  # if(missing(Select.eco)){Select.eco = NULL}
  # if(missing(Add.linear)){Add.linear = NULL}
  # if(missing(Add.linear)){Add.linear = NULL}
  if(missing(Strip.lab)){Strip.lab = T}
  if(missing(Show.diff)){Show.diff = F}
  if(missing(Vertical)){Vertical = F}
  # if(missing(Bit.map)){Bit.map = F}
  # if(missing(Add.n)){Add.n = F}
  # if(missing(Leg.pos)){Leg.pos = "right"}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Crop area Map ####
  if(is.null(Crop.zone) == F){
    crs.select <- 3857
    Select.map <- st_as_sf(map('world', fill = TRUE, plot = FALSE, region = Crop.zone))
    Select.map <- st_transform(Select.map, crs = crs.select) #32647, 3857
    M.surf <- st_as_sf(MCWT, coords = 3:2, crs = 4326) #4326
    M.surf <- st_transform(M.surf, crs = crs.select) # CRS pseudo-mercador -> coord plane
    M.surf <- st_intersection(M.surf, Select.map)
    M.surf <- st_transform(M.surf, crs = crs.select) # CRS Mongolia 
    M.surf$Longitude <- st_coordinates(M.surf)[,1]
    M.surf$Latitude <- st_coordinates(M.surf)[,2]
    MCWT <- as.data.frame(M.surf)
    
    if(is.null(MCWT2) == F){
      M.surf2 <- st_as_sf(MCWT2, coords = 3:2, crs = 4326) #4326
      M.surf2 <- st_transform(M.surf2, crs = crs.select) # CRS pseudo-mercador -> coord plane
      M.surf2 <- st_intersection(M.surf2, Select.map)
      print(M.surf2)
      M.surf2 <- st_transform(M.surf2, crs = crs.select) # CRS Mongolia 
      M.surf2$Longitude <- st_coordinates(M.surf2)[,1]
      M.surf2$Latitude <- st_coordinates(M.surf2)[,2]
      MCWT2 <- as.data.frame(M.surf2)
    }
    Select.map <- st_transform(Select.map, crs = crs.select) 
    Fond.carte <- geom_sf(data = Select.map, alpha = 0, color = "black", size = .5)
    Fond.carte.1 <- NULL
  }
  else{
    Fond.carte.1 <- geom_polygon(data = ACA.bo.proj, aes(x=long, y=lat),colour="black", fill = NA, size = .5)
    Fond.carte <- geom_polygon(data = ACA.bo.co.proj, aes(x=long, y=lat, group = group), colour = "grey40", fill = NA, size = .25, linetype = 2)}
  
  #### Plot maps CWT trait biogeography ####
  MCWT <- MCWT[c(grep("Lat", names(MCWT)), grep("Lon", names(MCWT)), which(names(MCWT) %in% Select.trait))]
  # print(Select.trait)
  MCWT <- melt(MCWT, id = c("Latitude","Longitude"))
  names(MCWT)<- c("Latitude", "Longitude", "CWM_Traits", "Value")
  MCWT$Type <- Type1
  
  if(is.null(MCWT2) == F){
    MCWT2 <- MCWT2[c(grep("Lat", names(MCWT2)), grep("Lon", names(MCWT2)), which(names(MCWT2) %in% Select.trait))]
    MCWT2 <- melt(MCWT2, id = c("Latitude","Longitude"))
    names(MCWT2)<- c("Latitude", "Longitude", "CWM_Traits", "Value")
    MCWT2$Type <- Type2
    MCWT <- rbind(MCWT, MCWT2)}
  
  #### ICIIII
  if(Show.diff == T){
    library(hexbin)
    library(ggplot2)
    
    set.seed(2)
    xA <- rnorm(1000)
    set.seed(3)
    yA <- rnorm(1000)
    set.seed(4)
    zA <- sample(c(1, 0), 20, replace = TRUE, prob = c(0.2, 0.8))
    hbinA <- hexbin(xA, yA, xbins = 40, IDs = TRUE)
    
    A <- data.frame(x = xA, y = yA, z = zA)
    
    set.seed(5)
    xB <- rnorm(1000)
    set.seed(6)
    yB <- rnorm(1000)
    set.seed(7)
    zB <- sample(c(1, 0), 20, replace = TRUE, prob = c(0.4, 0.6))
    hbinB <- hexbin(xB, yB, xbins = 40, IDs = TRUE)
    
    B <- data.frame(x = xB, y = yB, z = zB)
    
    
    ggplot(A, aes(x, y, z = z)) + stat_summary_hex(fun = function(z) sum(z)/length(z), alpha = 0.8) +
      scale_fill_gradientn(colours = c("blue","red")) +
      guides(alpha = FALSE, size = FALSE)
    
    ggplot(B, aes(x, y, z = z)) + stat_summary_hex(fun = function(z) sum(z)/length(z), alpha = 0.8) +
      scale_fill_gradientn (colours = c("blue","red")) +
      guides(alpha = FALSE, size = FALSE)
    
    ## find the bounds for the complete data 
    xbnds <- range(c(A$x, B$x))
    ybnds <- range(c(A$y, B$y))
    nbins <- 30
    
    #  function to make a data.frame for geom_hex that can be used with stat_identity
    makeHexData <- function(df) {
      h <- hexbin(df$x, df$y, nbins, xbnds = xbnds, ybnds = ybnds, IDs = TRUE)
      data.frame(hcell2xy(h),
                 z = tapply(df$z, h@cID, FUN = function(z) sum(z)/length(z)),
                 cid = h@cell)
    }
    
    Ahex <- makeHexData(A)
    Bhex <- makeHexData(B)
    
    ##  not all cells are present in each binning, we need to merge by cellID
    byCell <- merge(Ahex, Bhex, by = "cid", all = T)
    
    ##  when calculating the difference empty cells should count as 0
    byCell$z.x[is.na(byCell$z.x)] <- 0
    byCell$z.y[is.na(byCell$z.y)] <- 0
    
    ##  make a "difference" data.frame
    Diff <- data.frame(x = ifelse(is.na(byCell$x.x), byCell$x.y, byCell$x.x),
                       y = ifelse(is.na(byCell$y.x), byCell$y.y, byCell$y.x),
                       z = byCell$z.x - byCell$z.y)
    
    ##  plot the results
    
    ggplot(Ahex) +
      geom_hex(aes(x = x, y = y, fill = z),
               stat = "identity", alpha = 0.8) +
      scale_fill_gradientn (colours = c("blue","red")) +
      guides(alpha = FALSE, size = FALSE)
    
    ggplot(Bhex) +
      geom_hex(aes(x = x, y = y, fill = z),
               stat = "identity", alpha = 0.8) +
      scale_fill_gradientn (colours = c("blue","red")) +
      guides(alpha = FALSE, size = FALSE)
    
    ggplot(Diff) +
      geom_hex(aes(x = x, y = y, fill = z),
               stat = "identity", alpha = 0.8) +
      scale_fill_gradientn (colours = c("blue","red")) +
      guides(alpha = FALSE, size = FALSE)
  }
  # print(MCWT)
  
  n.trait <- nlevels(MCWT$CWM_Traits)
  MCWT$CWM_Traits <- as.factor(gsub("TRY_", "", MCWT$CWM_Traits))
  MCWT$Type <- as.factor(MCWT$Type)
  n.type <- nlevels(MCWT$Type)
  
  MCWT$Longitude <- round(MCWT$Longitude, digits = 0)
  MCWT$Latitude <- round(MCWT$Latitude, digits = 0)
  
  #### Annotations names Strig.lab = F ####
  if(Strip.lab == F){
    Strip.lab.disp <- element_blank()
    S.trait <- setNames(data.frame(as.factor(unique(MCWT$CWM_Traits)), rep(1,nlevels(MCWT$CWM_Traits)), rep(1,nlevels(MCWT$CWM_Traits))), c("Lab", "x","y"))
    S.clim <- setNames(data.frame(as.factor(unique(MCWT$Type)), rep(1,nlevels(MCWT$Type)), rep(1,nlevels(MCWT$Type))), c("Lab", "x","y"))
    
    Theme.null <- theme(axis.line = element_blank(), axis.title = element_blank(),
                        strip.text = element_blank(), axis.text = element_blank(),
                        axis.ticks = element_blank(), plot.background = element_blank(),
                        panel.grid = element_blank(), panel.background = element_blank())
    
    if(Vertical == T){print("Verticalisation !")
      p.up <- ggplot(S.clim, mapping = aes(x = x, y = y))+ 
        facet_wrap(vars(Lab), scales = "free_x", ncol = n.trait) + 
        geom_text(aes(label = Lab))+ Theme.null
      
      p.right <- ggplot(S.trait, mapping = aes(x = x, y = y))+ 
        facet_wrap(vars(Lab), scales = "free_x", nrow = n.trait) +
        geom_text(aes(label = Lab), angle = 270,  hjust=0.5, vjust=1)+ Theme.null}
    else{
      p.up <- ggplot(S.trait, mapping = aes(x = x, y = y))+ 
        facet_wrap(vars(Lab), scales = "free_x", ncol = n.trait) + 
        geom_text(aes(label = Lab))+ Theme.null
      
      p.right <- ggplot(S.clim, mapping = aes(x = x, y = y))+ 
        facet_wrap(vars(Lab), scales = "free_x", nrow = n.trait) +
        geom_text(aes(label = Lab), angle = 270,  hjust=0.5, vjust=1)+ Theme.null
    }
    
    
    
  }
  else{Strip.lab.disp <- element_text(hjust = 0)}
  
  
  #### Graphical settings ####
  Map_theme <- ggplot2::theme(
    plot.background = element_blank(), panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.title = element_text(),
    legend.key = element_blank(),
    legend.justification = c("center"),               # left, top, right, bottom
    legend.text = element_text(size = 8),
    legend.position = Leg.pos,
    # legend.position = "none",
    panel.background = element_blank(),
    strip.text = Strip.lab.disp,
    # panel.spacing = unit(0.7, "lines"),
    # plot.margin=unit(c(0,0,0,0),"cm")
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")
  )
  
  if(Vertical == T){
    My_facet <- facet_wrap(CWM_Traits ~ Type, ncol = n.type)
    }
  else{My_facet <- facet_wrap(Type ~ CWM_Traits, nrow = n.type)}
  
  #### Plots ####
  p <- ggplot() +
    #### Map settings ####
  # geom_tile(data = DEM.low.ACA.df, aes(x = x, y = y, fill = DEM.low, color = DEM.low), alpha = 0.1) +
  # scale_color_gradientn(colours = c("#1a9641", "#ccea8f", "#ffffc0", "#fed18a", "#f89957", "#ed6e43", "#ececec"), guide = "legend", name = "Elevation (m a.s.l.)")+
  # scale_fill_gradientn(colours = c("#1a9641", "#ccea8f", "#ffffc0", "#fed18a", "#f89957", "#ed6e43", "#ececec"), guide = "legend", name = "Elevation (m a.s.l.)")+
  # new_scale_fill()+
    scale_fill_gradientn(
                       # colors = c("#40004b", "#c2a5cf", "#a6dba0","#00441b"),# label = c("a", "b", "c", "d", "e"),
                       # colors = c("royalblue", "grey80", "darkorange"),# label = c("a", "b", "c", "d", "e"),
                       colors = c("#40004b", "grey80", "#00441b"),# label = c("a", "b", "c", "d", "e"),
                       # values = scales::rescale(c(summary(MCWT$Value)[[1]],summary(MCWT$Value)[[2]],summary(MCWT$Value)[[4]],summary(MCWT$Value)[[5]],summary(MCWT$Value)[[6]]), from = c(summary(MCWT$Value)[[1]],summary(MCWT$Value)[[6]])),
                       # values = scales::rescale(c(summary(MCWT$Value)[[1]],-0.2,0.4,summary(MCWT$Value)[[6]]), from = c(summary(MCWT$Value)[[1]],summary(MCWT$Value)[[6]])),
                       values = scales::rescale(c(-5,-0.7,0.5,1.3,5), from = c(-5,5)),
                       name = "CWM-traits (z-scores)")+
    # scale_fill_distiller(palette = "PRGn", direction = 1, na.value = NA, values = c(0,0.4999,0.5,0.5001,1))+#, guide = T, name = NULL, labels = NULL, breaks = NULL)+    # geom_point(data = MCWT, mapping = aes(x = Longitude, y = Latitude, color = Value), shape = 16, size = 3, na.rm = T)+
    Fond.carte +
    Fond.carte.1 +
    stat_summary_hex(data = MCWT, mapping = aes(x = Longitude, y = Latitude, z = Value), 
                     fun = function(x) mean(x), position = "jitter",
                     bins = Hex.size, color = NA) +
    # geom_point(data = MCWT, mapping = aes(x = Longitude, y = Latitude, color = Value), shape = 15) +
    My_facet +
    Map_theme
    # ggtitle("(B) occurence density used \n for the ACA plant checklist")+
    # scalebar(ACA.bo.proj, dist = 1000, dist_unit = "km", st.size = 1.7, border.size = 0.2,
    # transform = TRUE, model = "WGS84", st.bottom = FALSE, st.dist = 0.03,
    # facet.lev = c("TRY_SSD")) +
    # north(ACA.bo.proj, symbol = 6, location = "topright", scale = 0.1) +
  
 
  
  #### Export plot map ####
  if(Show.diff == T){
    #### Methode 1 ####  
    # print(summary(MCWT$Value)[[1]])
    # print(summary(MCWT$Value)[[6]])
    d1 <- ggplot_build(p)$data[[3]][, 2:7]
    d1$x <- round(d1$x, digits = 0)
    d1$y <- round(d1$y, digits = 0)
    d1 <- reshape2::dcast(d1, x + y ~ PANEL, fun.aggregate = mean)
    names(d1) <- c("Longitude", "Latitude", "Veg", "Pol")
    d1 <- na.omit(d1)                                 # Apply na.omit function
    d1$Veg <- (d1$Veg-min(d1$Veg))/(max(d1$Veg)-min(d1$Veg))
    d1$Pol <- (d1$Pol-min(d1$Pol))/(max(d1$Pol)-min(d1$Pol))
    d1$diff <-abs(d1$Veg - d1$Pol)
    print(d1)
    
    p2 <- ggplot() +
      geom_polygon(data = ACA.bo.co.proj, aes(x=long, y=lat, group = group), colour = "grey40", fill = NA, size = .25, linetype = 2) +
      geom_polygon(data = ACA.bo.proj, aes(x=long, y=lat),colour="black", fill = NA, size = .5) +
      stat_summary_hex(data = d1, mapping = aes(x = Longitude, y = Latitude, z = diff), fun = function(x) mean(x), bins = Hex.size, color = NA) +
      scalebar(ACA.bo.proj, dist = 1000, dist_unit = "km", st.size = 1.7, border.size = 0.2,
               transform = TRUE, model = "WGS84", st.bottom = FALSE, st.dist = 0.03,
               facet.lev = c("TRY_SSD")) +
      north(ACA.bo.proj, symbol = 6, location = "topright", scale = 0.1)+ 
      Map_theme
    
    #### Methode 2 ####
    Save.fit <- ggplot_build(p)$data[[3]][c(2:7)]
    my_digits = 0
    Trait.1 <- Save.fit[Save.fit$PANEL == 1,] 
    Trait.1p <- Save.fit[Save.fit$PANEL == 4,]
    
    Trait.1$x <- round(Trait.1$x, digits = my_digits)
    Trait.1$y <- round(Trait.1$y, digits = my_digits)
    Trait.1p$x <- round(Trait.1p$x, digits = my_digits)
    Trait.1p$y <- round(Trait.1p$y, digits = my_digits)
    
    names(Trait.1p)[c(3:6)] <- paste(names(Trait.1p)[c(3:6)], "p", sep = "_")
    Trait.1 <-full_join(Trait.1, Trait.1p, by = c("x", "y"))
    # Trait.1$Erreur <- abs(Trait.1$value_p - Trait.1$value)/Trait.1$value*100
    Trait.1$Erreur <- abs(Trait.1$value_p - Trait.1$value)
    
    # p2 <- ggplot(Trait.1, aes(x = x, y = y, color = Erreur)) +
    #             geom_point(na.rm = T, size = 3)+
    #             scale_color_gradientn(
    #               # colors = c("#40004b", "#c2a5cf", "#a6dba0","#00441b"),# label = c("a", "b", "c", "d", "e"),
    #               # colors = c("royalblue", "grey80", "darkorange"),# label = c("a", "b", "c", "d", "e"),
    #               colors = c("#40004b", "grey80", "#00441b"),# label = c("a", "b", "c", "d", "e"),
    #               # values = scales::rescale(c(summary(MCWT$Value)[[1]],summary(MCWT$Value)[[2]],summary(MCWT$Value)[[4]],summary(MCWT$Value)[[5]],summary(MCWT$Value)[[6]]), from = c(summary(MCWT$Value)[[1]],summary(MCWT$Value)[[6]])),
    #               # values = scales::rescale(c(summary(MCWT$Value)[[1]],-0.2,0.4,summary(MCWT$Value)[[6]]), from = c(summary(MCWT$Value)[[1]],summary(MCWT$Value)[[6]])),
    #               values = scales::rescale(c(-5,-0.7,0.5,1.3,5), from = c(-5,5)), na.value = NA,
    #               name = "CWM-traits (z-scores)")+ Map_theme
    # 
  
    
    }
  # print(names(Trait.1p))
  # print(Trait.1)
  
  if(Strip.lab == F){
    p <- p.up + plot_spacer() + p + p.right + 
      plot_layout(nrow = 2, heights = c(1/40,39/40), widths = c(39/40,1/40))}
  
  if(Show.diff == T & Strip.lab == T){
    p <- p + p2
    print("pouet")
    
  }
  
  # p <- p+p2+plot_layout(nrow = 2, heights = c(2/3, 1/3))
  #### Export ####
  print(p)
  if(is.null(Save.plot) == F){dev.off()}
  # return(Trait.1)
  return(p)
}

R2.compar <- function(MP, MV, Xlab, Ylab, Bisectrice.area, Leg.pos, shade.areas, Pourc.coerrent = T,
                      p.val.show = T, Panel.lab, Show.Plotly, Repel.outliers, Nb.labs, H, W, Save.plot){
  #### Settings ####
  library(ggrepel)
  if(missing(Leg.pos)){Leg.pos = "top"}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(Panel.lab)){Panel.lab = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Xlab)){Xlab = NULL}
  if(missing(Ylab)){Ylab = NULL}
  if(missing(Nb.labs)){Nb.labs = NULL}
  if(missing(Bisectrice.area)){Bisectrice.area = NULL}
  if(missing(shade.areas)){shade.areas = F}
  if(missing(Repel.outliers)){Repel.outliers = T}
  if(missing(Show.Plotly)){Show.Plotly = F}
  library(RColorBrewer)
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Calc ####
  P1 <- suppressWarnings(melt(MP$R2))
  names(P1) <- c("X1", "X2", "R2_pol") 
  P2 <- suppressWarnings(melt(MV$R2))
  names(P2) <- c("X1", "X2", "R2_veg") 
  P1 <- left_join(P1, P2, by = c("X1", "X2"))
  P1$Func <- paste(P1$X2, " vs. ", P1$X1, "", sep = "")
  P1$Mean <- rowMeans(P1[c(3,4)])
  
  #### Keep only one combination possible ####
  Comb <- data.frame(gtools::combinations(nrow(MP$R2), 2, row.names(MP$R2)))
  Comb$Func <- paste(Comb$X2, " vs. ", Comb$X1, "", sep = "")
  # Dup <- duplicated(P1[c(3,4)]) #& duplicated(P1[c(1,2)], fromLast = T)
  # P1 <- P1[!Dup,]
  P1 <- P1[P1$Func %in% Comb$Func,]
  P1 <- P1[P1$R2_pol !=1,]
  
  #### p-value ####
  if(p.val.show == T){
    Mpval1 <- suppressWarnings(melt(MP$p.val))
    names(Mpval1) <- c("X1", "X2", "pval_pol") 
    Mpval2 <- suppressWarnings(melt(MV$p.val))
    names(Mpval2) <- c("X1", "X2", "pval_veg") 
    Mpval1 <- left_join(Mpval1, Mpval2, by = c("X1", "X2"))
    Mpval1$Func <- paste(Mpval1$X2, " vs. ", Mpval1$X1, "", sep = "")
    Mpval1 <- Mpval1[Mpval1$Func %in% Comb$Func,]
    Mpval1 <- Mpval1[Mpval1$pval_veg !=1,]
    Mpval1$Significant <- Mpval1$pval_pol < 0.001 & Mpval1$pval_veg < 0.001
    Mpval1$Significant[Mpval1$Significant == T] <- "p < 0.0001 (significant)"
    Mpval1$Significant[Mpval1$Significant == F] <- "p > 0.0001 (insignificant)"
    P1 = cbind(P1, Signific = Mpval1$Significant)
    }
  else{P1$Signific = "All dots"}
  
  #### Signe coerence ####
  if(Pourc.coerrent == T){
    print("**** r direction consistence between vegetation and pollen : ****")
    Tab.coer <- P1
    Tab.coer$Pos.pol <- P1$R2_pol >= 0
    Tab.coer$Pos.veg <- P1$R2_veg >= 0
    Tab.coer$Coerence <- Tab.coer$Pos.pol == Tab.coer$Pos.veg
    print(Tab.coer[c(5,9)])
    Pourc.coerence <- 100 - round(length(which(Tab.coer$Coerence == F))/nrow(Tab.coer)*100, digits = 0)
    print(paste(Pourc.coerence, "% of the relations are consistent in direction."))
    }
  
  #### Graph settings ####
  if(is.null(Xlab) == F){Xlab <- xlab(Xlab)}
  if(is.null(Ylab) == F){Ylab <- ylab(Ylab)}
  # if(is.null(Panel.lab) == F){Panel.name <- geom_text(inherit.aes=FALSE, aes(x = -1, y = 1, label = Panel.lab), colour="grey30", parse=T, hjust = 0.5)}
  if(is.null(Panel.lab) == F){Panel.name <- geom_text(inherit.aes=F, aes(x = -1, y = 1, label = Panel.lab), 
                                                      colour="grey30", hjust = 0.5, na.rm = T)}
  else{Panel.name <- NULL}
  Add.r2 <- stat_poly_eq(label.y = "bottom", label.x = "right", color = "turquoise4", size = 3.5, small.r = F, small.p = T,
                         aes(label =  sprintf("%s*\", \"*%s" ,
                                              after_stat(rr.label),
                                              # after_stat(r.squared),
                                              after_stat(p.value.label)
                         )))
  if(is.null(Bisectrice.area) == F){Bisectrice.area <- geom_abline(aes(slope = 1, intercept = 0), linewidth = Bisectrice.area, color = "grey70", alpha = 0.2, linetype = "solid")}
  # else{Bisectrice.area <-}
  
  if(shade.areas == T){
    Shade.box <- geom_rect(inherit.aes = F, mapping = aes(xmin = 0, xmax = 1 , ymin = -1, ymax = 0), fill = "grey90", colour = NA, na.rm = T)
    Shade.box2 <- geom_rect(inherit.aes = F, mapping = aes(xmin = 0, xmax = -1 , ymin = 1, ymax = 0), fill = "grey90", colour = NA, na.rm = T)
  }
  else{Shade.box = Shade.box2 <- NULL}
  #### Repel selection ####
  if(is.null(Nb.labs) == F){
    dist2d <- function(a) {
      b = c(0,0)
      c = c(1,1)
      
      v1 <- b - c
      v2 <- a - b
      m <- cbind(v1,v2)
      d <- abs(det(m))/sqrt(sum(v1*v1))
      } 
    
    P1$Dist <- NA
    for(i in 1:nrow(P1)){
      a <- c(P1$R2_pol[i], P1$R2_veg[i])
      P1$Dist[i] <- dist2d(a)
      }
    
    
    if(Repel.outliers == T){
      P1 <- P1[!is.na(P1$Dist),]
      Val4 <- sort(P1$Dist)[length(P1$Dist)-Nb.labs]
      Sub.DB <- subset(P1, Dist > Val4)
      }
    else{
      Val4 <- sort(P1$Dist)[Nb.labs]
      Sub.DB <- subset(P1, Dist <= Val4)
      # print(Sub.DB)
      }
    # print(Val4)
    Annot <- geom_text_repel(data = Sub.DB, aes(label = Func), size = 2.5, force = 100, max.overlaps = 100,
                             force_pull = 10,
                             nudge_x = 0.3,
                             # box.padding = 2.5,
                             nudge_y = 0.3,
                             # min.segment.length = 0, # draw all lines no matter how short
                             segment.size = 0.3,
                             # segment.curvature = 0.35,
                             segment.curvature = 0.4,
                             # segment.ncp = 0.2,
                             segment.angle = 35,
                             # label.size=NA, #no border/box
                             # fill = NA, #
                             # position = position_dodge(1),
                             # min.segment.length = 1,
                             show.legend = F)
    }
  else{Annot <- NULL}
  
  #### Best table of r ####
  P1$Pond <- 1-abs(P1$Dist)
  P1$Pond <- abs(P1$R2_pol)*P1$Pond
  P1 <- P1[order(P1$Pond, decreasing = T),]
  Table.best <- P1[c(1:10), c(5,3,4)]
  # print(Table.best)
  
  #### Plot ####
  p <- ggplot(P1, aes(x = R2_pol, y = R2_veg, color = Mean))+
    Xlab + Ylab + Shade.box + Shade.box2 +
    geom_abline(aes(slope = 1, intercept = 0), linewidth = .8, color = "grey60", linetype = "dashed")+
    Bisectrice.area +
    geom_point(aes(shape = Signific))+
    lims(x = c(-1,1), y = c(-1,1))+
    geom_smooth(method = "lm", se = F, span = 1000, linetype = "longdash", linewidth = 0.5, color = "turquoise4",
                formula = y ~ x)+    
    Add.r2 +
    scale_color_gradientn(colors = c("royalblue", "grey50", "darkorange", "red"), 
                          # values = scales::rescale(c(-1,0,0.1,1), from = c(-1,1)),
                          values = c(0,0.26,0.38,0.5,1),
                          name = "r-value")+
    Annot +
    Panel.name +
    theme(axis.line = element_line(colour = "grey30"), 
          plot.background = element_blank(), panel.background = element_blank(),
          panel.grid = element_blank(), legend.direction = "horizontal", legend.position = Leg.pos,
          legend.key = element_blank(), legend.text = element_text(size = 4, angle = 0),
          legend.key.size = unit(4, "mm"), legend.title = element_text(size = 8, vjust = 0.75),
          panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 1)
          
    )
  
  #### Save html ####
  if(Show.Plotly == T){
    library(plotly)
    library(htmlwidgets)
    # print(names(Mtern))
    Save.plot.html <- gsub("pdf", "html", Save.plot)
    Keep.name <- gsub(".*\\/", "", Save.plot.html)
    Path.root <- paste(gsub(Keep.name, "", Save.plot.html), "HTML_files/", sep = "")
    if(file.exists(Path.root) == F){dir.create(Path.root)}
    Save.plot.html <- paste(Path.root, Keep.name, sep = "")
    p1_ly <- ggplotly(p)
    p1_ly <- p1_ly %>% layout(boxmode = "group", boxpoints = F)
    options(warn = - 1) 
    saveWidget(p1_ly, file = Save.plot.html)
  }
  #### Export ####
  print(p)
  if(is.null(Save.plot) == F){dev.off()}
  return(p)
}

CWT.calculation <- function(MT, MP, Mclim, MPS.ACA.Biom, Accep.seuil, Remove.biom){
  #### Settings ####
  if(missing(Remove.biom)){Remove.biom = NULL}
  #### Clean matrix ####
  Keep.real.names <- row.names(MP)
  MP <- data.frame(t(MP))
  names(MP) <- Keep.real.names
  MP$id <- row.names(MP)
  names(MT) <- gsub("X", "", names(MT))
  names(MT) <- paste("TRY", names(MT), sep = "_")
  names(MT)[1] <- "id"
  MT <- MT[which(MT$id %in% MP$id),]
  # /!\ missing taxa ! 
  Missing.taxa <- setdiff(MP$id, MT$id)
  print("The following taxa are missing from the trait matrix:")
  print(Missing.taxa)
  # return(Missing.taxa)
  Missing.raw <- MT[setdiff(MP$id, MT$id),]
  Missing.raw$id <- setdiff(MP$id, MT$id)
  MT <- rbind(MT, Missing.raw)
  MT <- as_tibble(lapply(MT, function(x){x[is.nan(x)] <- NA ; x}))
  
  #### Merge MP + MT ####
  MPT <- merge(MT, MP, by = "id", all = T)
  MPT <- melt(MPT, id = names(MT))
  names(MPT)[names(MPT) == "variable"] <- "Site"
  names(MPT)[names(MPT) == "value"] <- "FA"
  
  #### Calculation of the CWT ####
  MCWT <- data.frame(Site = (unique(MPT["Site"])))
  MCWT.stat <- setNames(data.frame(NA,NA,NA), c("Trait", "Pour.sites.kept", "N.site"))
  MCWT.stat <- MCWT.stat[-1,]
  Tot.site.nb <- length(levels(MPT$Site))
  
  #### Verbose ####
  pb = txtProgressBar(min = 1, 
                      max = length(grep("TRY", names(MPT))),
                      width = 40,
                      initial = 0,  style = 3) 
  
  init <- numeric(length(grep("TRY", names(MPT))))
  end <- numeric(length(grep("TRY", names(MPT))))
  
  #### Main loop ####
  for(i in grep("TRY", names(MPT))){
    init[i] <- Sys.time()
    Trait.treat.i <- names(MPT)[i]
    A <- MPT[c("id", Trait.treat.i, "Site", "FA")]
    A$FA[is.na(A[[Trait.treat.i]])] <- NA
    all_abund = aggregate(FA ~ Site, A, sum)
    colnames(all_abund)[2] = "tot_abund"
    Keep.sites <- all_abund[all_abund[2] > Accep.seuil,]
    # print(Trait.treat.i)
    # print(nrow(Keep.sites))
    if(nrow(Keep.sites) > 0){
      N.site <- length(Keep.sites$tot_abund)
      Pourc.site.up.seuil <- round(length(Keep.sites$tot_abund)/Tot.site.nb, digits = 2)
      MCWT.stat[i,] <- c(Trait.treat.i, Pourc.site.up.seuil, N.site)
      A <- A[which(A$Site %in% Keep.sites$Site),]
      A = merge(A, all_abund, by = "Site")
      A$FA = A$FA/A$tot_abund
      XX <- aggregate(FA * eval(parse(text = Trait.treat.i)) ~ Site, A, sum, na.rm = T)
      # print(typeof(XX$Site))
    }
    else{print("zob")
      XX <- data.frame(Site = NA, X = "")
    }
    names(XX)[2] <- Trait.treat.i
    # XX$Site <- as.character(XX$Site)
    # print(typeof(XX$Site))
    MCWT <- left_join(MCWT, XX, by = "Site")
    
    end[i] <- Sys.time()
    setTxtProgressBar(pb, i)
    time <- round(seconds_to_period(sum(end - init)), 0)
    est <- length(grep("TRY", names(MPT))) * (mean(end[end != 0] - init[init != 0])) - time
    remainining <- round(seconds_to_period(est), 0)
    cat(paste(" // Execution time:", time,
              " // Estimated time remaining:", remainining), "")
  }
  close(pb)
  
  MCWT.stat <- MCWT.stat[-1,]
  MCWT.stat[nrow(MCWT.stat)+1,] <- c("Average", round(mean(as.numeric(MCWT.stat$Pour.sites.kept)), digits = 2), round(mean(as.numeric(MCWT.stat$N.site)), digits = 0))
  
  #### Merge CWT + climat ####
  Mclim[["Site"]] <- row.names(Mclim)
  MPS.ACA.Biom[["Site"]] <- row.names(MPS.ACA.Biom)
  MCWT.clim = merge(MCWT, Mclim, by = "Site")    # On fusionne les matrices CWM et CLIMAT
  MCWT.clim = merge(MCWT.clim, MPS.ACA.Biom, by = c("Site", "Latitude", "Longitude"))    # On fusionne les matrices CWM et CLIMAT
  if(is.null(Remove.biom) == F){MCWT.clim <- MCWT.clim[setdiff(seq(1,nrow(MCWT.clim)), which(MCWT.clim$Biome %in% Remove.biom)),]}
  
  #### Export ####
  Lexport <- list(MCWT = MCWT.clim, MCWT.stat = MCWT.stat, Missing.taxon = Missing.taxa)
  return(Lexport)
}


#### Import data ####
Import.data = F
if(Import.data == T){
  #### SIG ####
  SIG = F
  if(SIG == T){
    if(exists("Path.ACA.border") == F){
      library(raster) # processing spatial raster data. !!!overwrites dplyr::select!!!
      library(rgdal)
      library(stars)
      library(maptools)
      library(ggforce) # geom_arc_bar
      library(ggsn) # north and scalebar
      library(sp)
      library(sf)
      
      Path.ACA.border = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/Projets/ACA/Borders_ACA/Extern_border/Out_border_ACA.shp"
      # Path.ACA.border = "/media/lucas.dugerdil/Climacoptera/Documents/Recherche/SIG/Projets/ACA/Borders_ACA/Extern_border/Out_border_ACA.shp"
      ACA.bo = readOGR(Path.ACA.border)
      proj4string(ACA.bo) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
      # proj4string(ACA.bo) <- CRS("+proj=lcc +lat_0=30 +lon_0=95 +lat_1=15 +lat_2=65 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
      ACA.bo.proj = fortify(ACA.bo)
      
      Path.ACA.border.co = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/Projets/ACA/Borders_ACA/Full/ACA_border.shp"
      ACA.bo.co = readOGR(Path.ACA.border.co)
      proj4string(ACA.bo.co) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
      ACA.bo.co.proj = fortify(ACA.bo.co)
      
      Crop.biome = F
      if(Crop.biome == T){
        Path.ACA.biom = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/Data_vegetation/Ecoregions2017/Ecoregions2017.shp"
        ACA.biom = readOGR(Path.ACA.biom)
        proj4string(ACA.biom) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
        ACA.biom <- intersect(ACA.bo, ACA.biom)
        writeOGR(ACA.biom, dsn = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/Projets/ACA/Ecosystem2017_ACA_biome_R", layer = "ACA.biom", driver="ESRI Shapefile") #also you were missing the driver argument
      }
      else{
        Path.ACA.biom = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/Projets/ACA/Ecosystem2017_ACA_biome_R/ACA.biom.shp"
        ACA.biom = readOGR(Path.ACA.biom)
        ACA.biom@data$id <- rownames(ACA.biom@data)
        ACA.biom.proj = fortify(ACA.biom)
        ACA.biom.proj <- left_join(ACA.biom.proj, ACA.biom@data)
      }
      
      Crop.DEM = F
      if(Crop.DEM == T){
        Path.DEM = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/MNT/WorldClim2_MNT/wc2.1_2.5m_elev.tif"
        DEM = raster(Path.DEM)
        DEM.ACA <- crop(DEM, extent(ACA.bo))
        DEM.ACA <- mask(DEM.ACA, ACA.bo)
        writeRaster(DEM.ACA, "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/Projets/ACA/DEM/DEM_ACA_wc2.5m.tif", format = "GTiff")
  
        Path.DEM.low = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/MNT/WorldClim2_MNT/wc2.1_10m_elev.tif"
        DEM.low = raster(Path.DEM.low)
        DEM.low.ACA <- crop(DEM.low, extent(ACA.bo))
        DEM.low.ACA <- mask(DEM.low.ACA, ACA.bo)
        writeRaster(DEM.low.ACA, "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/Projets/ACA/DEM/DEM_ACA_wc10m.tif", format = "GTiff")
        
        }
      else{
        Path.DEM = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/Projets/ACA/DEM/DEM_ACA_wc2.5m.tif"
        DEM.ACA = raster(Path.DEM)
        DEM.ACA.df <- as(DEM.ACA, "SpatialPixelsDataFrame")
        DEM.ACA.df <- as.data.frame(DEM.ACA.df)
        colnames(DEM.ACA.df) <- c("DEM", "x", "y")
        Path.DEM.low = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/Projets/ACA/DEM/DEM_ACA_wc10m.tif"
        DEM.low.ACA = raster(Path.DEM.low)
        DEM.low.ACA.df <- as(DEM.low.ACA, "SpatialPixelsDataFrame")
        DEM.low.ACA.df <- as.data.frame(DEM.low.ACA.df)
        colnames(DEM.low.ACA.df) <- c("DEM.low", "x", "y")
      }}}
  
  #### GBIF ####
  GBIF = F
  if(GBIF == T){
    library("rgbif")
    
    if(file.exists(DB.path) == F){warning("The disk dure n'a pas ete branche.")}
    GBIF.oc.iran <- data.frame(read.csv(file = paste(DB.path, "Vegetation/Occurences/GBIF_v2.1.0/GBIF_Iran_plantlist.csv", sep = "") ,sep="\t",dec=".",header=T, row.names = 1, stringsAsFactors = FALSE))
    GBIF.pl.ACA <- data.frame(read.csv(file = paste(DB.path, "Vegetation/Occurences/GBIF_v2.1.0/GBIF_ACA_plantlist.csv", sep = "") ,sep="\t",dec=".",header=T, row.names = 1, stringsAsFactors = FALSE))
    #GBIF.oc.ACA <- data.frame(read.csv(file = paste(DB.path, "Vegetation/Occurences/GBIF_v2.1.0/GBIF_ACA_occurrence.csv", sep = "") ,sep="\t",dec=".",header=T, row.names = 1, stringsAsFactors = FALSE))
    GBIF.oc.ACA <- vroom::vroom(file = paste(DB.path, "Vegetation/Occurences/GBIF_v2.1.0/GBIF_ACA_occurrence.csv", sep = ""))
    
    GBIF.oc.ACA.spdf <- SpatialPointsDataFrame(coords = GBIF.oc.ACA[,c(23,22)], data = GBIF.oc.ACA, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    GBIF.oc.ACA <- GBIF.oc.ACA.spdf[ACA.bo,]
    GBIF.oc.ACA <- data.frame(GBIF.oc.ACA)
    saveRDS(GBIF.oc.ACA,  paste(DB.path, "Vegetation/Occurences/GBIF_v2.1.0/GBIF_ACA_occ_V2.Rds", sep = ""))
    
    #A <- occ_data(country = c("IR"), kingdomKey = 6, limit = 100000)
    #B <- occ_data(country = c("IR"), kingdomKey = 6, start = 50000, limit = 40000)
    #A <- occ_download(pred("country", "IR"), pred("kingdomKey", 6))
    
    #GBIF.oc.Iran <- as.data.frame(A$data)    
    }
  
  #### TRY ####
  TRY.imp = F
  if(TRY.imp == T){
    TRY <- vroom::vroom(file = paste(DB.path, "Traits/TRY_dec_2019/8066.csv", sep = ""))
    List.tax.TRY <- unique(TRY[,c("AccSpeciesName", "AccSpeciesID")])
    saveRDS(List.tax.TRY,  paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_plantlist.Rds", sep = ""))
    
    #### TNRS ####
    TNRS.check == F
    if(TNRS.check == T){
      #### function ####
      TNRS.check <- function(M) {
        print("Checklist cleaned by the TNRS API request")
        library(TNRS)
        results <- data.frame()
        Step = 1000
        
        #### Loops pour bp de taxa ####
        for(i in 1:ceiling(length(M)/Step)){
          Min <- 1 + (i - 1)*Step
          Max <- i*Step
          if(Max > length(M)){Max = length(M)}
          print(paste("From species", Min, "to", Max))
          Add.list <- TNRS(M[Min:Max])
          results <- rbind(results, Add.list)
        }
        print(paste("Number of species cleaned by the TNRS:", nrow(results), sep = " "))
        results$Accepted_family[results$Taxonomic_status == "No opinion"] <- results$Name_matched_accepted_family[results$Taxonomic_status == "No opinion"] 
        results$Accepted_name[results$Taxonomic_status == "No opinion"] <- results$Name_matched[results$Taxonomic_status == "No opinion"] 
        results <- results[c(2,3,15,14,21,41,33,34,7)]
        names(results) <- c("AccSpeciesName", "Score", "Score_G", "TNRS_Genus", "Score_F", "TNRS_Family", "Status_taxo", "Accepted_name", "Rank")
        return(results)
      }
      
      #### Clean data ####
      Import.checklist <- List.tax.TRY[["AccSpeciesName"]]
      Import.checklist <- Import.checklist[!duplicated(Import.checklist)]
      result <- TNRS.check(Import.checklist)
      TL.clean <- left_join(List.tax.TRY, results, by = "AccSpeciesName")
      Second.chance <- TL.clean$AccSpeciesName[is.na(TL.clean$Score)]
      TL.clean <- TL.clean[!is.na(TL.clean$Score),]
      result2 <- TNRS.check(Second.chance)
      TL.clean <- left_join(List.tax.TRY, result2, by = "AccSpeciesName")
      troisiem.chance <- TL.clean$AccSpeciesName[is.na(TL.clean$Score)]
      TL.clean <- TL.clean[!is.na(TL.clean$Score),]
      result3 <- TNRS.check(troisiem.chance)
      TL.clean <- left_join(TL.clean, result3, by = "AccSpeciesName")
      TL.clean <- TL.clean[!duplicated(TL.clean),]
      TL.clean <- TL.clean[TL.clean$Score > 0.9,]
      TL.clean$Genera <- gsub("\\s.*", "", TL.clean$Accepted_name)
      
      result.f <- rbind(result, result2, result3)
      result.f <- result.f[!duplicated(result.f),]
      A <- result.f[duplicated(result.f$AccSpeciesName)|duplicated(result.f$AccSpeciesName, fromLast = T),]
      A <- A[A$Score < 1,]
      result.f <- result.f[!(result.f$AccSpeciesName  %in% A$AccSpeciesName & result.f$Score < 1),]
      TL.clean <- left_join(List.tax.TRY, result.f, by = "AccSpeciesName")
      TL.clean$Accepted_name[which(TL.clean$Score < 0.91)] <- TL.clean$AccSpeciesName[which(TL.clean$Score < 0.91)]
      TL.clean$Genera <- gsub("\\s.*", "", TL.clean$Accepted_name)
      TL.clean <- TL.clean[c(1,2,7,11,9,10)]
      saveRDS(TL.clean,  paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_plantlist_TNRS.Rds", sep = ""))
      }
    else{saveRDS(TL.clean,  paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_plantlist.Rds", sep = ""))}    
    
    List.trait.TRY <- unique(TRY[,c("TraitName", "TraitID")])
    saveRDS(List.trait.TRY,  paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_traitlist.Rds", sep = ""))
    }
  
  #### MPS.ACA ####
  MPS.crop.ACA = F
  if(MPS.crop.ACA == T){
    #### Coord EPDB ####
    SIG.pollen = F
    if(SIG.pollen == T){
      MPS.13549.Cord <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DB13549_Coord.Rds")
      MPS.13549.Cord <- MPS.13549.Cord[which(is.na(MPS.13549.Cord$Long) == F),]
      # MPS.13549.Cord.spdf <- SpatialPointsDataFrame(coords = MPS.13549.Cord[,1:2], data = MPS.13549.Cord, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      MPS.13549.Cord.spdf <- SpatialPointsDataFrame(coords = MPS.13549.Cord[,1:2], data = MPS.13549.Cord, proj4string = CRS("+proj=lcc +lat_0=30 +lon_0=95 +lat_1=15 +lat_2=65 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
      MPS.ACA.Cord <- MPS.13549.Cord.spdf[ACA.bo,]
      MPS.ACA.Cord <- data.frame(MPS.ACA.Cord)
      saveRDS(MPS.ACA.Cord, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Coord.Rds")
      write.csv(MPS.ACA.Cord, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Coord.csv")
    }
    else{
      MPS.ACA.Cord <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Coord.Rds")
    }
    
    MPS.13549.PoP_sl <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DB13549_Pol_PT_sl.Rds")
    MPS.13549.PoP_ss <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DB13549_Pol_PT_ss.Rds")
    MPS.ACA.PoP_sl <- MPS.13549.PoP_sl[row.names(MPS.ACA.Cord),]
    MPS.ACA.PoP_sl <- MPS.ACA.PoP_sl[,!colSums(MPS.ACA.PoP_sl) == 0]
    MPS.ACA.PoP_ss <- MPS.13549.PoP_ss[row.names(MPS.ACA.Cord),]
    MPS.ACA.PoP_ss <- MPS.ACA.PoP_ss[,!colSums(MPS.ACA.PoP_ss) == 0]
    
    #### Functions clean TNRS pollen taxa ####
    Check.TNRS.APG.pol <- function(MP, Seuil, Table.conv.taxa){
      if(missing(Table.conv.taxa)){Table.conv.taxa = F}
      APG_diff <- data.frame(read.csv("Import/World_DB/Taxonomie/APG_II_APG_IV_diff.csv" ,sep=",",dec=".",header=T, stringsAsFactors = F))
      if(missing(Seuil)){Seuil = 0.90}
      
      #### Transpose ####
      t_df <- as.data.frame(t(MP),row.names = NA)
      t_df <- cbind(names(MP),t_df)
      names(t_df)[1] <- "ID"
      
      #### TNRS ####
      Tax.TNRS.pol <- TNRS.taxa.accepted(Taxa.vector = names(MP), Taxonstand.test = T)
      print(Tax.TNRS.pol[order(as.numeric(Tax.TNRS.pol$match.score)),])
      # return(Tax.TNRS.pol)
      Tax.TNRS.pol <- Tax.TNRS.pol[Tax.TNRS.pol$match.score >= Seuil,]
      # t_df$ID[which(t_df$ID %in% Tax.TNRS.pol$Name_submitted)] <- Tax.TNRS.pol$Accepted_name[which(Tax.TNRS.pol$Name_submitted %in% t_df$ID)]
      t_df$ID[which(t_df$ID %in% Tax.TNRS.pol$Name_submitted)] <- Tax.TNRS.pol$Accepted_name_TXSD[which(Tax.TNRS.pol$Name_submitted %in% t_df$ID)]
      t_df$ID <- gsub(" NA", "", t_df$ID)
      t_df$ID[which(t_df$ID %in% APG_diff$Old_name)] <- APG_diff$New_name[match(t_df$ID[which(t_df$ID %in% APG_diff$Old_name)], APG_diff$Old_name)]
      t_df$ID <- gsub(" subsp.\\s.*", "", t_df$ID)
      t_df$ID[t_df$ID == "Bignonia"] <- "Bignoniaceae"
      t_df$ID[t_df$ID == "Keteleeria"] <- "Pinaceae"
      t_df$ID[t_df$ID == "Buddleia"] <- "Buddleja"
      t_df$ID[t_df$ID == "Phyteuma"] <- "Ranunculaceae"
      t_df$ID[t_df$ID == "Anemoclema"] <- "Ranunculaceaee"
      t_df$ID[t_df$ID == "Chosenia"] <- "Salix"
      t_df$ID[t_df$ID == "Flacourtia"] <- "Salicaceae"
      t_df$ID[t_df$ID == "Lannea"] <- "Anacardiaceae"
      t_df$ID[t_df$ID == "Pteroceltis"] <- "Ulmus"
        
      Table.taxa <- data.frame(Old.name = names(MP), New.name = t_df$ID)
      #### Aggregate ####
      new_df=aggregate(t_df[-1], by=list(t_df[["ID"]]), sum)
      t_df2 <- as.data.frame(t(new_df[-1]),row.names = NA)
      names(t_df2) <- new_df$Group.1
      row.names(t_df2) <- names(new_df)[-1]
      MP <- t_df2
      
      #### Export ####
      if(Table.conv.taxa == F){return(MP)}
      if(Table.conv.taxa == T){return(Table.taxa)}
      }
     
    #### Application #### 
    MPS.ACA.PoP_sl.clean <- Check.TNRS.APG.pol(MPS.ACA.PoP_sl, Seuil = 0.90)
    MPS.ACA.PoP_ss.clean <- Check.TNRS.APG.pol(MPS.ACA.PoP_ss, Seuil = 0.78)
    Table.conv.taxa_sl <- Check.TNRS.APG.pol(MPS.ACA.PoP_sl, Seuil = 0.90, Table.conv.taxa = T)
    Table.conv.taxa_ss <- Check.TNRS.APG.pol(MPS.ACA.PoP_ss, Seuil = 0.78, Table.conv.taxa = T)
    names(MPS.ACA.PoP_ss.clean)[names(MPS.ACA.PoP_ss.clean) == "NA"] <- "Anemoclema"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$New.name == "NA"] <- "Anemoclema"
    names(MPS.ACA.PoP_ss.clean)[names(MPS.ACA.PoP_ss.clean) == "Asteridea"] <- "Asteroideae"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$New.name == "Asteridea"] <- "Asteroideae"
    
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$Old.name == "Aellenia"] <- "Halothamnus"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$Old.name == "Astralagus"] <- "Astragalus"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$Old.name == "Capparidaceae"] <- "Capparaceae"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_sl$Old.name == "Capparidaceae"] <- "Capparaceae"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$Old.name == "Dianthus petrorhagia"] <- "Dianthus"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$Old.name == "Gramineae"] <- "Poaceae"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$Old.name == "Hololachna"] <- "Reaumuria"
    Table.conv.taxa_sl$New.name[Table.conv.taxa_sl$Old.name == "Hololachna"] <- "Reaumuria"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$Old.name == "Onopordon"] <- "Onopordum"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$Old.name == "Knorringia"] <- "Polygonum"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$Old.name == "Ophelia"] <- "Swertia"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$Old.name == "Helianthemum"] <- "Helianthemum"
    Table.conv.taxa_ss$New.name[Table.conv.taxa_ss$Old.name == "Rhamnus"] <- "Rhamnus"
    
    #### Export ####
    saveRDS(MPS.ACA.PoP_sl, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_sl.Rds")
    saveRDS(MPS.ACA.PoP_ss, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_ss.Rds")
    saveRDS(MPS.ACA.PoP_sl.clean, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_sl_clean.Rds")
    saveRDS(MPS.ACA.PoP_ss.clean, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_ss_clean.Rds")
    saveRDS(Table.conv.taxa_sl, "Resultats/World_DB/Pollen/Merging_DB/DB13549/Table.conv.taxa_sl.Rds")
    saveRDS(Table.conv.taxa_ss, "Resultats/World_DB/Pollen/Merging_DB/DB13549/Table.conv.taxa_ss.Rds")
    
    }
  
  #### Vegetation plot ####
  Veget.import = F
  if(Veget.import == T){
    #### Jamsranjav_et_al_2018 ####
    Veg.Mong = F
    if(Veg.Mong == T){
      VPmong <- tibble::as_tibble(read.csv(file="Import/Mongolia/Vegetation/Jamsranjav_et_al_2018/Jamsranjav_et_al_2018_data.csv",sep=",", dec=".", header=T, stringsAsFactors = F))
      VPmong <- as_tibble(cbind(ID = paste("ACAV-", sprintf("%05d", as.numeric(seq(1, nrow(VPmong)))), sep = ""), Dataset = "Mongolia_Jamsranjav_2018", Country = "Mongolia",  VPmong))
      Taxon.VPmong <- tibble::as_tibble(read.csv(file="Import/Mongolia/Vegetation/Jamsranjav_et_al_2018/Taxon_correspo_Jam_2018.csv",sep=",", dec=".", header=T, stringsAsFactors = F))
      Clean.names <- strsplit(Taxon.VPmong$Species.Name, split = "\\s")
      Clean.names <- sapply(Clean.names, function(x) paste(x[][1], tolower(x[][2])))
      
      Taxon.VPmong.tnrs <- TNRS.taxa.accepted(Taxa.vector = Clean.names)
      Taxon.VPmong.tnrs <- Taxon.VPmong.tnrs[Taxon.VPmong.tnrs$match.score >= 0.95,]
      Clean.names <- left_join(data.frame(Name_submitted = Clean.names), Taxon.VPmong.tnrs, by = "Name_submitted")
      Taxon.VPmong <- as_tibble(cbind(Taxon.VPmong, Clean.names[c(2,5)]))
      Taxon.VPmong[is.na(Taxon.VPmong$Accepted_name),"Accepted_name"] <- Taxon.VPmong[is.na(Taxon.VPmong$Accepted_name),"Species.Name"] 
      Duplicated.taxa <- Taxon.VPmong[duplicated(Taxon.VPmong$Accepted_name)|duplicated(Taxon.VPmong$Accepted_name, fromLast = T),]
      
      VPmong.plot <- VPmong[c(1,42:386)]
      VPmong.plot$ARDR <- VPmong.plot$ARDR + VPmong.plot$ARGL
      VPmong.plot$GEDE <- VPmong.plot$GEDE + VPmong.plot$GEMA
      VPmong.plot$KOMA <- VPmong.plot$KOMA + VPmong.plot$KOMU
      VPmong.plot$SIAD <- VPmong.plot$SIAD + VPmong.plot$SIJE
      VPmong.plot <- VPmong.plot[setdiff(names(VPmong.plot), c("ARGL", "GEMA", "KOMU", "SIJE")) ]
      Taxon.VPmong <- Taxon.VPmong[-c(which(Taxon.VPmong$Species.Code %in% c("ARGL", "GEMA", "KOMU", "SIJE"))),]
      Taxon.VPmong <- Taxon.VPmong[!grepl("Unknown", Taxon.VPmong$Accepted_name),]
      VPmong.plot <- VPmong.plot[,!grepl("Unknown", names(VPmong.plot))]
      VPmong.plot <- VPmong.plot[,!grepl("UK", names(VPmong.plot))]
      
      names(VPmong.plot)[which(names(VPmong.plot) %in% Taxon.VPmong$Species.Code)] <- Taxon.VPmong$Accepted_name[which(names(VPmong.plot)[-1] %in% Taxon.VPmong$Species.Code)]
      VPmong.plot[-1] <- VPmong.plot[-1]/rowSums(VPmong.plot[-1])

      #### Aggregate duplicated taxa ###      
      t_df <- as.data.frame(t(VPmong.plot[-1]),row.names = NA)
      names(t_df) <- VPmong.plot$ID
      t_df <- cbind(names(VPmong.plot)[-1],t_df)
      names(t_df)[1] <- "ID"
      new_df=aggregate(t_df[-1], by=list(t_df[["ID"]]), sum)
      t_df2 <- as.data.frame(t(new_df[-1]),row.names = NA)
      names(t_df2) <- new_df$Group.1
      t_df2 <- cbind(names(new_df)[-1],t_df2)
      names(t_df2)[1] <- "ID"
      VPmong.plot <- t_df2
      saveRDS(VPmong, "Import/ACA/Vegetation/VPmong_Jamsranjav_metadata.Rds")
      saveRDS(VPmong.plot, "Import/ACA/Vegetation/VPmong_Jamsranjav.Rds")
      saveRDS(Taxon.VPmong, "Import/ACA/Vegetation/VPmong_Jamsranjav_taxalist.Rds")
      }
    else{
      VPmong <- readRDS("Import/ACA/Vegetation/VPmong_Jamsranjav_metadata.Rds")
      Taxon.VPmong <- readRDS("Import/ACA/Vegetation/VPmong_Jamsranjav_taxalist.Rds")
      VPmong.plot <- readRDS("Import/ACA/Vegetation/VPmong_Jamsranjav.Rds")}    
    
    #### sPlot in the ACA area ####
    Coordonee.splot = F
    if(Coordonee.splot == T){
      load("Import/ACA/Vegetation/sPlot_ACA.RData")
      header.ACA.spdf <- SpatialPointsDataFrame(coords = header.ACA[,7:8], data = header.ACA, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      sPlot.ACA.Cord <- header.ACA.spdf[ACA.bo,]
      sPlot.ACA.Cord <- as_tibble(sPlot.ACA.Cord)
      sPlot.ACA.Cord <- as_tibble(cbind(ID = paste("ACAV-", sprintf("%05d", as.numeric(seq((nrow(VPmong)+1), (nrow(sPlot.ACA.Cord)+nrow(VPmong))))), sep = ""), sPlot.ACA.Cord))
      
      header.ACA <- header.ACA[match(sPlot.ACA.Cord$PlotObservationID, header.ACA$PlotObservationID),]
      header.ACA <- as_tibble(cbind(ID = paste("ACAV-", sprintf("%05d", as.numeric(seq((nrow(VPmong)+1), (nrow(header.ACA)+nrow(VPmong))))), sep = ""), header.ACA))
      
      # List potential co-authors
      List.Pot.co.autaut <- db.out[match(unique(sPlot.ACA.Cord$GIVD.ID), db.out$`GIVD ID`),]
      
      write.csv(List.Pot.co.autaut, "Resultats/ACA/Vegetation/List_ACA_contributors.csv")
      
      saveRDS(sPlot.ACA.Cord, "Import/ACA/Vegetation/sPlot.ACA.Cord.Rds")
      saveRDS(header.ACA, "Import/ACA/Vegetation/sPlot.ACA.metadata.Rds")
      }
    else{sPlot.ACA.Cord <- readRDS("Import/ACA/Vegetation/sPlot.ACA.Cord.Rds")}
    
    #### Pulishing the relevés ####
    Pulish.taxa.TNRS = F
    if(Pulish.taxa.TNRS == T){
      splot.CWM.ACA <- CWM.ACA[which(CWM.ACA$PlotObservationID %in% sPlot.ACA.Cord$PlotObservationID),]
      splot.CWM.ACA$ID <- sPlot.ACA.Cord$ID[match(splot.CWM.ACA$PlotObservationID, sPlot.ACA.Cord$PlotObservationID)]
      splot.soil.clim <- soilclim.ACA[which(soilclim.ACA$PlotObservationID %in% sPlot.ACA.Cord$PlotObservationID),]
      
      splot.veget <- DT2.ACA[which(DT2.ACA$PlotObservationID %in% sPlot.ACA.Cord$PlotObservationID),]
      splot.veget$ID <- sPlot.ACA.Cord$ID[match(splot.veget$PlotObservationID, sPlot.ACA.Cord$PlotObservationID)]
      Taxon.splot <- unique(splot.veget$Species) 
      
      Taxon.splot.tnrs1 <- TNRS.taxa.accepted(Taxa.vector = Taxon.splot[1:3000])
      Taxon.splot.tnrs2 <- TNRS.taxa.accepted(Taxa.vector = Taxon.splot[3001:length(Taxon.splot)])
      Taxon.splot.tnrs <- rbind(Taxon.splot.tnrs1, Taxon.splot.tnrs2)
      splot.trait.ACA <- traits.eurasia[which(traits.eurasia$Species %in% splot.veget$Species),]
      
      saveRDS(Taxon.splot.tnrs, "Import/ACA/Vegetation/sPlot_taxa_TNRS_raw.Rds")
      saveRDS(splot.CWM.ACA, "Import/ACA/Vegetation/sPlot_ACA_CWM.Rds")
      saveRDS(splot.veget, "Import/ACA/Vegetation/sPlot_ACA_veget.Rds")
      saveRDS(splot.soil.clim, "Import/ACA/Vegetation/sPlot_ACA_soilclim.Rds")
      saveRDS(splot.trait.ACA, "Import/ACA/Vegetation/sPlot_ACA_trait.Rds")
      }
    else{
      Taxon.splot.tnrs <- readRDS("Import/ACA/Vegetation/sPlot_taxa_TNRS_raw.Rds")
      splot.CWM.ACA <- readRDS("Import/ACA/Vegetation/sPlot_ACA_CWM.Rds")
      splot.veget <- readRDS("Import/ACA/Vegetation/sPlot_ACA_veget.Rds")
      splot.soil.clim <- readRDS("Import/ACA/Vegetation/sPlot_ACA_soilclim.Rds")
      splot.trait.ACA <- readRDS("Import/ACA/Vegetation/sPlot_ACA_trait.Rds")
      }
    
    #### Taxa list of MV ####
    Full.taxa.list.MV = F
    if(Full.taxa.list.MV == T){
      MV_TL <- full_join(Taxon.VPmong, Taxon.splot.tnrs, by= c("Accepted_name", "Accepted_family"))
      MV_TL <-MV_TL[c(4,5)]
      MV_TL <-MV_TL[-1,]
      MV_TL$Accepted_name[grep(".*subsp.", MV_TL$Accepted_name)] <- gsub(" subsp.\\s.*", "", MV_TL$Accepted_name[grep(".*subsp.", MV_TL$Accepted_name)])
      MV_TL$Accepted_name[grep(".*var.", MV_TL$Accepted_name)] <- gsub(" var.\\s.*", "", MV_TL$Accepted_name[grep(".*var.", MV_TL$Accepted_name)])
      MV_TL$genus <- gsub("\\s.*", "", MV_TL$Accepted_name)
      MV_TL$Taxa_level <- unlist(purrr::map(strsplit(MV_TL$Accepted_name, split = "\\s"), length))
      MV_TL <- MV_TL[!MV_TL$Taxa_level == 0,]   
      
      MV_TL$Taxa_level[MV_TL$Taxa_level == 3] <- 2
      MV_TL$Accepted_name[MV_TL$Accepted_name == "Leguminosae"] <- "Fabaceae"
      MV_TL$Taxa_level[MV_TL$Taxa_level == 1][grep("eae", MV_TL$Accepted_name[MV_TL$Taxa_level == 1])] <- 3
      
      Only.genus <- MV_TL$Accepted_name[MV_TL$Taxa_level == 1]
      Species.genus <- MV_TL$genus[MV_TL$Taxa_level == 2]
      MV_TL <- MV_TL[!MV_TL$Accepted_name %in% intersect(Only.genus, Species.genus),]   
      MV_TL[which(MV_TL$Accepted_family == ""),2] <- c("Lecideaceae", "Thelotremataceae", "Lichinaceae", "Peltulaceae", "Pertusariaceae", "Trapeliaceae", "Pertusariaceae")
      
      saveRDS(MV_TL, "Resultats/ACA/Vegetation/MV_ACA_taxonlist.Rds")
      }
    else{
      MV_TL <- readRDS("Resultats/ACA/Vegetation/MV_ACA_taxonlist.Rds")
    }
    
    #### Test CWM de splot ####
    test.extract.CWM.splot = T
    if(test.extract.CWM.splot == T){
      List.trait.splot <- unique(splot.CWM.ACA$variable)
      Keep.trait <- List.trait.splot[c(4, 7, 8, 14, 16, 17, 20, 24, 26, 30)]
      splot.CWM.ACA <- splot.CWM.ACA[which(splot.CWM.ACA$variable %in% Keep.trait),c(8,2,6)]
      
      Trait.table <- c("SLA_mean","SeedMass_mean", "PlantHeight_mean", "LeafN_mean","LeafArea.undef.undef_mean","StemDens_mean", "LeafCN.ratio_mean", "LeafThickness_mean", "LDMC_mean", "LeafP_mean")
      Trait.table.2 = c("TRY_SLA","TRY_SeedMass","TRY_Height","TRY_LeafN","TRY_LeafArea","TRY_SSD", "TRY_CNRatio", "TRY_LeafThick", "TRY_LDMC", "TRY_LeafP")
      Trait.table <- tibble(Trait.table, Trait.table.2)
      splot.CWM.ACA$variable <- Trait.table$Trait.table.2[match(splot.CWM.ACA$variable, Trait.table$Trait.table)]
      splot.CWM.ACA <- tidyr::spread(splot.CWM.ACA, variable, CWM, drop = F, fill = NA)
      
      #MT.ACA.no.scale <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_trait_ACA_noscale.Rds", sep = ""))
      
      # test <- as_tibble(scale(MT.ACA.no.scale[5:8]))
      # test.2 <- as_tibble(scale(MCWT.clim.splot[4:ncol(MCWT.clim.splot)], scale = T, center = F))
      
      MCWT.clim.splot <- sPlot.ACA.Cord[c(1,9,8)]
      MCWT.clim.splot <- dplyr::full_join(MCWT.clim.splot, splot.CWM.ACA, by = "ID")
      MCWT.clim.splot[4:ncol(MCWT.clim.splot)] <- scale(MCWT.clim.splot[4:ncol(MCWT.clim.splot)])
      
      MV_climtot <- readRDS("Resultats/ACA/Vegetation/DBACA_Clim_tot.Rds")
      MV_climtot$ID <- row.names(MV_climtot)
      MCWT.clim.splot <- dplyr::full_join(MCWT.clim.splot, MV_climtot, by = "ID")
      
      MV.ACA.Biom <- readRDS("Resultats/ACA/Vegetation/DBACA_Biom.Rds")
      MV.ACA.Biom$ID <- row.names(MV.ACA.Biom)
      MCWT.clim.splot <- dplyr::full_join(MCWT.clim.splot, MV.ACA.Biom, by = "ID")
      
      MV.ACA.Land <- readRDS("Resultats/ACA/Vegetation/DBACA_LandCover.Rds")
      MV.ACA.Land$ID <- row.names(MV.ACA.Land)
      MCWT.clim.splot <- dplyr::full_join(MCWT.clim.splot, MV.ACA.Land, by = "ID")
      
      names(MCWT.clim.splot) <- gsub("\\.x", "", names(MCWT.clim.splot))
      names(MCWT.clim.splot) <- gsub("\\.y", "", names(MCWT.clim.splot))
      names(MCWT.clim.splot)[1] <- "Site"
      MCWT.clim.splot <- MCWT.clim.splot[!duplicated(names(MCWT.clim.splot))]
      MCWT.clim.splot <- MCWT.clim.splot[c(1:3,13,8,9,11,10,6,4,5,7,12,14:27,31,32,28,29,30)]
      
      saveRDS(MCWT.clim.splot, "Import/ACA/Vegetation/sPlot_ACA_CWM_clim.Rds")
    }
    Taiga.test = F
    if(Taiga.test == T){
      Remove.taig <- readRDS("Import/ACA/Vegetation/sPlot_remove_fake_taiga.Rds")
      
      sPlot.ACA.Cord.taig <- sPlot.ACA.Cord[sPlot.ACA.Cord$ID %in% Remove.taig,]
      
    }
    #### Merge veget data ####
    Merge.data.veg = F
    if(Merge.data.veg == T){
      Select.header <- c("ID","Dataset", "Country", "Latitude", "Longitude")
      MV.Cord <- as_tibble(rbind(VPmong[Select.header], sPlot.ACA.Cord[Select.header]))
      
      Taxon.to.keep <- unique(c(MV_TL$Accepted_name, MV_TL$genus))
      splot.veget$Accepted_species <- Taxon.splot.tnrs$Accepted_name[match(splot.veget$Species, Taxon.splot.tnrs$Name_submitted)]
      splot.veget <- splot.veget[splot.veget$Accepted_species %in% Taxon.to.keep,]
      splot.veget.unnest <- splot.veget[c("ID", "Relative_cover", "Accepted_species")]
      names(splot.veget.unnest) <- c("ID", "FA", "Species")
      splot.veget.unnest <- aggregate(FA ~ ID + Species, splot.veget.unnest, sum)
      splot.veget.unnest <- tidyr::spread(splot.veget.unnest, Species, FA, drop = T, fill = 0.00)
      splot.veget.unnest <- splot.veget.unnest[-2]
      
      MV.releve <- dplyr::full_join(VPmong.plot, splot.veget.unnest, by = "ID")
      names(MV.releve) <- gsub("\\.x", "", names(MV.releve))
      MV.releve[is.na(MV.releve)]<- 0
      MV.releve <- reshape2::melt(MV.releve, by = "ID")
      names(MV.releve) <- c("ID", "Species", "FA")
      all_abund = aggregate(FA ~ ID, MV.releve, sum)
      colnames(all_abund)[2] = "FA_tot"
      MV.releve = merge(MV.releve, all_abund, by = "ID")
      MV.releve$FA = MV.releve$FA / MV.releve$FA_tot
      MV.releve <- MV.releve[-c(4)]
      MV.releve$Species <- gsub("\\.y", "", MV.releve$Species)
      MV.releve$Species <- gsub("\\.y", "", MV.releve$Species)
      
      MV.releve$Species[grep(".*subsp.", MV.releve$Species)] <- gsub(" subsp.\\s.*", "", MV.releve$Species[grep(".*subsp.", MV.releve$Species)])
      MV.releve$Species[grep(".*var.", MV.releve$Species)] <- gsub(" var.\\s.*", "", MV.releve$Species[grep(".*var.", MV.releve$Species)])
      
      MV.releve[MV.releve$Species == "Leguminosae",2] <- "Fabaceae"
      
      saveRDS(MV.releve, "Resultats/ACA/Vegetation/MV_ACA_plot.Rds")
      saveRDS(MV.Cord, "Resultats/ACA/Vegetation/MV_ACA_cord.Rds")
      write.csv(MV.Cord, "Resultats/ACA/Vegetation/MV_ACA_cord.csv")
      }
    }
  
  
  }

#### API request ####
BIEN.request = F
if(BIEN.request == T){
  #### Occurrence ####
  #BIEN.oc.ACA <- BIEN_occurrence_spatialpolygons(spatialpolygons = ACA.bo) # Ne fonctionne pas.
  BIEN.oc.ACA <- BIEN_occurrence_country(c("Armenia", "Georgia", "Azerbaijan", "Mongolia", "Iran", "Afghanistan", "Kazakhstan", "Kyrgyzstan", "Pakistan", "Tajikistan", "Turkmenistan", "Uzbekistan"))
  
  BIEN.oc.Sib <- BIEN_occurrence_country("Russia")
  BIEN.oc.Sib <- BIEN.oc.Sib[which(is.na(BIEN.oc.Sib$longitude) == F),]
  BIEN.oc.Sib.spdf <- SpatialPointsDataFrame(coords = BIEN.oc.Sib[,c(3,2)], data = BIEN.oc.Sib, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  BIEN.oc.Sib <- BIEN.oc.Sib.spdf[ACA.bo,]
  BIEN.oc.Sib <- data.frame(BIEN.oc.Sib)
  
  BIEN.oc.Chi <- BIEN_occurrence_country(country = "China")
  BIEN.oc.Chi <- BIEN.oc.Chi[which(is.na(BIEN.oc.Chi$longitude) == F),]
  BIEN.oc.Chi.spdf <- SpatialPointsDataFrame(coords = BIEN.oc.Chi[,c(3,2)], data = BIEN.oc.Chi, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  BIEN.oc.Chi <- BIEN.oc.Chi.spdf[ACA.bo,]
  BIEN.oc.Chi <- data.frame(BIEN.oc.Chi)
  
  BIEN.oc.Ind <- BIEN_occurrence_country(country = "India")
  BIEN.oc.Ind <- BIEN.oc.Ind[which(is.na(BIEN.oc.Ind$longitude) == F),]
  BIEN.oc.Ind.spdf <- SpatialPointsDataFrame(coords = BIEN.oc.Ind[,c(3,2)], data = BIEN.oc.Ind, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  BIEN.oc.Ind <- BIEN.oc.Ind.spdf[ACA.bo,]
  BIEN.oc.Ind <- data.frame(BIEN.oc.Ind)
  
  BIEN.oc.ACA <- rbind(BIEN.oc.ACA[,1:6], BIEN.oc.Chi[,1:6], BIEN.oc.Sib[,1:6], BIEN.oc.Ind[,1:6])
  saveRDS(BIEN.oc.ACA,  paste(DB.path, "Vegetation/Occurences/BIEN/BIEN_ACA_plantlist.Rds", sep = ""))
  }

#### Plot maps ####
Plot.occur = F
if(Plot.occur == T){
  #### Check if map imported ####
  GBIF.oc.ACA <- readRDS(paste(DB.path, "Vegetation/Occurences/GBIF_v2.1.0/GBIF_ACA_occ_V2.Rds", sep = ""))
  BIEN.oc.ACA <-   readRDS(paste(DB.path, "Vegetation/Occurences/BIEN/BIEN_ACA_plantlist.Rds", sep = ""))
  MPS.ACA.Cord <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Coord.Rds")
  MV.Cord <- readRDS("Resultats/ACA/Vegetation/MV_ACA_cord.Rds")
  
  #### Map 1 : veget ####
  pmap <- ggplot() +
    #### Map settings ####
    geom_polygon(data = ACA.bo.co.proj, aes(x=long, y=lat, group = group), colour = "grey45", fill = "grey90", size = .25, linetype = 2) +
    geom_polygon(data = ACA.biom.proj, aes(x=long, y=lat, group = group, fill = factor(BIOME_NAME)), colour = NA, size = .05, alpha = 1) +
    ggtitle("(A) ACAV and ACAP")+
    scalebar(ACA.bo.proj, dist = 500, dist_unit = "km", st.size = 3, border.size = 0.2,
             transform = TRUE, model = "WGS84", st.bottom = FALSE, st.dist = 0.03,
             facet.lev = c("4.2ky")) +
    north(ACA.bo.proj, symbol = 6, location = "topright", scale = 0.1) +
    scale_x_continuous(name = "Longitude (°)") +
    scale_y_continuous(name = "Latitude (°)") +
    scale_fill_manual(values = c("Deserts & Xeric Shrublands" = "#C88282",
                                 "Temperate Grasslands, Savannas & Shrublands" = "#ECED8A",
                                 "Montane Grasslands & Shrublands" = "#D0C3A7",
                                 "Temperate Conifer Forests" = "#6B9A88",
                                 "Temperate Broadleaf & Mixed Forests" = "#3E8A70",
                                 "N/A" = "#FFEAAF",
                                 "Tundra" = "#A9D1C2",
                                 "Boreal Forests/Taiga" = "#8FB8E6",
                                 "Tropical & Subtropical Coniferous Forests" = "#99CA81",
                                 "Mangroves" = "#FE01C4",
                                 "Flooded Grasslands & Savannas" = "#BEE7FF",
                                 "Tropical & Subtropical Moist Broadleaf Forests" = "#38A700"
                                  ), name = "Biomes (Dinerstein et al., 2017)"#, labels = levels(ACA.biom.proj$BIOME_NAME)
                       ) +
    geom_polygon(data = ACA.bo.co.proj, aes(x=long, y=lat, group = group), colour = "grey20", fill = NA, size = .25, linetype = 2) +
    geom_polygon(data = ACA.bo.proj, aes(x=long, y=lat),colour="black", fill = NA, size = .5) +
    #### Points /lines ####
    geom_point(data = MV.Cord, mapping = aes(x = Longitude, y = Latitude), shape = 19, alpha = 0.7, fill = "white", colour = "darkred", size = 1)+ #
    geom_point(data = MPS.ACA.Cord, mapping = aes(x = Long, y = Lat), shape = 17, alpha = 0.7, fill = "white", colour = "darkblue", size = 1)+
    # geom_path(data=EASML, aes(x=long, y=lat), colour="royalblue", linetype = 2, size = .6) +
    coord_quickmap(xlim = c(43, 130), ylim = c(25, 64)) +
    #### Theme ####
  theme(
    # plot.background = element_rect(colour = "grey"), 
    panel.border = element_rect(colour = "grey", fill = NA),
    axis.line= element_line(colour = "grey", lineend = "butt"),
    axis.ticks.x.bottom = element_line(colour = "grey"),
    legend.title = element_text(),
    legend.justification = c("center"),               # left, top, right, bottom
    legend.text = element_text(size = 8),
    panel.background = element_blank(),
    # panel.spacing = unit(0.7, "lines"),
    plot.margin=unit(c(0,0,0,0),"cm")
  )
  
  
  
  ##### Map DEM ####
  pmap.dem <- ggplot() +
    #### Map settings ####
    ggtitle("(C) DEM")+
    scalebar(ACA.bo.proj, dist = 500, dist_unit = "km", st.size = 3, border.size = 0.2,
             transform = TRUE, model = "WGS84", st.bottom = FALSE, st.dist = 0.03,
             facet.lev = c("4.2ky")) +
    north(ACA.bo.proj, symbol = 6, location = "topright", scale = 0.1) +
    geom_tile(data = DEM.low.ACA.df, aes(x = x, y = y, fill = DEM.low, color = DEM.low), alpha = 1) +
    scale_fill_gradientn(colours = c("#1a9641", "#ccea8f", "#ffffc0", "#fed18a", "#f89957", "#ed6e43", "#ececec"), guide = "legend", name = "Elevation (m a.s.l.)")+
    scale_color_gradientn(colours = c("#1a9641", "#ccea8f", "#ffffc0", "#fed18a", "#f89957", "#ed6e43", "#ececec"), guide = "legend", name = "Elevation (m a.s.l.)")+
    geom_polygon(data = ACA.bo.co.proj, aes(x=long, y=lat, group = group), colour = "grey20", fill = NA, size = .25, linetype = 2) +
    geom_polygon(data = ACA.bo.proj, aes(x=long, y=lat),colour="black", fill = NA, size = .5) +
    #### Points ####
  geom_point(data = MV.Cord, mapping = aes(x = Longitude, y = Latitude), shape = 1, alpha = 0.6, fill = NA, colour = "grey30", size = 0.2)+ #
  geom_point(data = MPS.ACA.Cord, mapping = aes(x = Long, y = Lat), shape = 1, alpha = 0.6, fill = NA, colour = "grey30", size = 0.2)+
  coord_quickmap(xlim = c(43, 130), ylim = c(25, 64)) +
    #### Theme ####
  theme(
    # plot.background = element_rect(colour = "grey"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.title = element_text(),
    legend.key = element_blank(),
    legend.justification = c("center"),               # left, top, right, bottom
    legend.position = "right",
    # legend.position = "none",
    legend.text = element_text(size = 8),
    panel.background = element_blank(),
    # panel.spacing = unit(0.7, "lines"),
    # plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")
    plot.margin=unit(c(0,0,0,0),"cm")
  )
  
  ##### Map occurences ####
  pmap.occur <- ggplot() +
    #### Map settings ####
    scale_fill_gradientn(colors = c("royalblue", "darkorange", "darkred"),
                         labels = c("1-10", "10-100", "100-1000", "1000-10000", ">10000"),
                         breaks = c(1, 10, 100, 1000, 10000),
                       guide = "legend",
                       values = scales::rescale(c(1,10,100,1000,10000), from = c(1,10000)),
                       name = "Occurance (#.cell-1)")+
    geom_hex(data = GBIF.oc.ACA, mapping = aes(x = decimalLongitude, y = decimalLatitude),  alpha = 0.8, bins = 120, color = NA)+
    geom_hex(data = BIEN.oc.ACA, mapping = aes(x = longitude, y = latitude),  alpha = 0.8, bins = 120, color = NA)+
    #ggtitle("(B) occurence density used \n for the ACA plant checklist")+
    # scalebar(ACA.bo.proj, dist = 500, dist_unit = "km", st.size = 3, border.size = 0.2,
    #          transform = TRUE, model = "WGS84", st.bottom = FALSE, st.dist = 0.03,
    #          facet.lev = c("4.2ky")) +
    # north(ACA.bo.proj, symbol = 6, location = "topright", scale = 0.1) +
    # geom_polygon(data = ACA.bo.co.proj, aes(x=long, y=lat, group = group), colour = "grey20", fill = NA, size = .25, linetype = 2) +
    geom_polygon(data = ACA.bo.proj, aes(x=long, y=lat),colour="black", fill = NA, size = .5) +
    # coord_quickmap(xlim = c(43, 130), ylim = c(25, 64)) +
    # coord_map(projection = "lambert", lat0=30, lat1 = 15)+
    # coord_proj(paste0("+proj=lcc +lat_0=30 +lon_0=95 +lat_1=15 +lat_2=65 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))+
    #### Theme ####
    theme(
    plot.background = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.title = element_text(),
    legend.key = element_blank(),
    legend.justification = c("center"),               # left, top, right, bottom
    legend.text = element_text(size = 8),
    legend.position = "right",
    # legend.position = "none",
    panel.background = element_blank(),
    # panel.spacing = unit(0.7, "lines"),
    # plot.margin=unit(c(0,0,0,0),"cm")
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")
  )
  
  #### Export plot map ####
  library(patchwork)
  W = 850
  H = 500
  # dev.copy(png,'myplot.png')
  # pmap.occur
  # dev.off()
  ggsave(filename = "Figures/ACA/Trait/Map/Map_ACA_occurences.pdf", pmap.occur, width = W*0.026458333, height = H*0.026458333, units = "cm")
  # ggsave(filename = "Figures/ACA/Trait/Map/Map_ACA_DEM.pdf", pmap.dem, width = W*0.026458333, height = H*0.026458333, units = "cm")
  # ggsave(filename = "Figures/ACA/Trait/Map/Map_ACA_VP_SPo.pdf", pmap, width = W*0.026458333, height = H*0.026458333, units = "cm")
  # ptot <- pmap / (pmap.occur + pmap.dem) + plot_layout(guides="collect")
  # W = 2000
  # H = 2000
  # ggsave(filename = "Figures/ACA/Trait/Map/Map_ACA_occurences_V2.pdf", ptot, width = W*0.026458333, height = H*0.026458333, units = "cm", dpi = 20)
}

#### Test Maroc ####
Maroc.Eric = F
if(Maroc.Eric == T){
  EuroMed <- as_tibble(read.csv("Import/Maroc/Traits/EUROSL-2021.csv", sep = "\t", dec=".", header = T, stringsAsFactors = F))
  APD <- as_tibble(read.csv("Import/Maroc/Traits/APD_NorthAfrica_avril2021.csv", sep = "\t", dec=".", header = T, stringsAsFactors = F))
  
  Common.taxa <- intersect(EuroMed$ScientificName, APD$NOM_SIMPLE) 
  Keep.taxa.APD <- setdiff(APD$NOM_SIMPLE, EuroMed$ScientificName)
  APD.nodup <- APD[which(APD$NOM_SIMPLE %in% Keep.taxa.APD),]
  names(APD.nodup)[1] <- "ID"
  APD.nodup$ID <- paste("APD", APD.nodup$ID, sep = "_")
  
  Keep.taxa.EuroM <- setdiff(EuroMed$ScientificName, APD$NOM_SIMPLE)
  EuroMed.nodup <- EuroMed[which(EuroMed$ScientificName %in% Keep.taxa.EuroM),]
  EuroMed.nodup$ID <- paste("EM", EuroMed.nodup$ID, sep = "_")
  EuroMed.nodup <- EuroMed.nodup[-c(2:5, 7, 11, 12, 15)]
  
  #### Plus simple ####
  All.taxa <- union(EuroMed$ScientificName, APD$NOM_SIMPLE) 
  All.taxa <- All.taxa[which(lengths(gregexpr("\\W+", All.taxa)) + 1 >= 2)] # Remove genus only
  #A <- data.frame(X1 = All.taxa, X2 = lengths(gregexpr("\\W", All.taxa)))
  
  source("Scripts/TNRS.R")
  TL.check.med <- data.frame(Name_submitted = NA, Accepted_name = NA, match.score = NA, Unmatched_terms = NA,  Accepted_family = NA)
  for(i in 1:ceiling(length(All.taxa)/5000)){
    Min <- 1 + (i - 1)*5000
    Max <- i*5000
    print(paste("From", Min, "to", Max))
    Add.list <- TNRS.taxa.accepted(Taxa.vector = All.taxa[Min:Max])
    TL.check.med <- rbind(TL.check.med, Add.list)
    }
  TL.check.med <- TL.check.med[-c(1),]
  TL.check.med$EuroMed <- F
  TL.check.med[which(TL.check.med$Name_submitted %in% EuroMed$ScientificName), "EuroMed"] <- T
  TL.check.med$APD <- F
  TL.check.med[which(TL.check.med$Name_submitted %in% APD$NOM_SIMPLE), "APD"] <- T
  
  write.csv(TL.check.med, file = "Resultats/Maroc/Traits/APD_EuromeD_TNRS_taxa_list_raw.csv")
  
  # TL.check <- as_tibble(TL.check[TL.check$match.score >= 0.95,]) 
  # names(TL.merge)[names(TL.merge) == "species"] <- "Name_submitted"
   }

#### Treatment ####
Taxon.list.ACA = F
if(Taxon.list.ACA == T){
  #### Import data GBIF / BIEN / splot ####
  library(BIEN)
  
  if(exists("GBIF.oc.ACA") == F){GBIF.oc.ACA <- readRDS(paste(DB.path, "Vegetation/Occurences/GBIF_v2.1.0/GBIF_ACA_occ_V2.Rds", sep = ""))}
  if(exists("BIEN.oc.ACA") == F){BIEN.oc.ACA <- readRDS(paste(DB.path, "Vegetation/Occurences/BIEN/BIEN_ACA_plantlist.Rds", sep = ""))}
  if(exists("MV_TL") == F){MV_TL <- readRDS("Resultats/ACA/Vegetation/MV_ACA_taxonlist.Rds")}
  Table.conv.taxa_sl <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/Table.conv.taxa_sl.Rds")
  Table.conv.taxa_ss <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/Table.conv.taxa_ss.Rds")
  
  #### Merge taxon list from GBIF and BIEN #### 
  TL.GBIF <- GBIF.oc.ACA[GBIF.oc.ACA$taxonRank == "SPECIES",] 
  # TL.GBIF <- TL.GBIF[TL.GBIF$taxonomicStatus == "ACCEPTED", c(11,15,17,19)]
  TL.GBIF <- TL.GBIF[c("class", "order", "family", "genus", "species")]
  TL.GBIF <- TL.GBIF[order(TL.GBIF$species),]
  TL.GBIF <- TL.GBIF[!duplicated(TL.GBIF),]
  TL.BIEN <- sort(unique(BIEN.oc.ACA$scrubbed_species_binomial[!is.na(BIEN.oc.ACA$scrubbed_species_binomial)]))
  TL.BIEN.taxonomy <- BIEN_taxonomy_species(species = TL.BIEN)
  TL.BIEN.taxonomy.clean <- unique(TL.BIEN.taxonomy[TL.BIEN.taxonomy$scrubbed_taxonomic_status == "Accepted", c("class", "order", "scrubbed_family", "scrubbed_genus", "scrubbed_species_binomial")])
  names(TL.BIEN.taxonomy.clean) <- names(TL.GBIF)
  Add.to.GBIF <- sort(setdiff(TL.BIEN.taxonomy.clean$species, TL.GBIF$species))
  TL.BIEN.taxonomy.clean <- TL.BIEN.taxonomy.clean[match(Add.to.GBIF, TL.BIEN.taxonomy.clean$species),]
  TL.merge <- rbind(TL.GBIF, TL.BIEN.taxonomy.clean)
  TL.merge <- TL.merge[order(TL.merge$species),]
  
  #### Check TNRS accepted especies ####
  TNRS = F
  if(TNRS == T){
    source("Scripts/TNRS.R")
    TL.check <- data.frame(Name_submitted = NA, Accepted_name = NA, match.score = NA, Unmatched_terms = NA,  Accepted_family = NA, Accepted_name_TXSD = NA)
    for(i in 1:ceiling(nrow(TL.merge)/5000)){
      Min <- 1 + (i - 1)*5000
      Max <- i*5000
      print(paste("From", Min, "to", Max))
      Add.list <- TNRS.taxa.accepted(Taxa.vector = TL.merge[Min:Max,"species"], Taxonstand.test = T)
      TL.check <- rbind(TL.check, Add.list)
      }
    TL.check <- TL.check[-c(1),]
    # TL.check <- as_tibble(TL.check[c(1,5,6,3)])
    # names(TL.check)[3] <- "Accepted_name"
    saveRDS(TL.check, paste(DB.path, "Vegetation/Occurences/Merge_DB/Taxon_list_ACA_V0_tnrs.Rds", sep = ""))
    }
  if(TNRS == F){TL.check <- readRDS(paste(DB.path, "Vegetation/Occurences/Merge_DB/Taxon_list_ACA_V0_tnrs.Rds", sep = ""))}
  
  #### Merge with pollen taxa if missing ####
  TP.add <- unique(c(Table.conv.taxa_sl$New.name, Table.conv.taxa_ss$New.name))
  TP.add <- TP.add[!(TP.add %in% TL.check$Accepted_name)]
  TP.add <- TP.add[!(TP.add %in% unique(gsub("\\s.*", "", TL.check$Accepted_name)))]
  TP.add <- TP.add[!(TP.add %in% TL.check$Accepted_family)]
  Pollen.check = F
  if(Pollen.check == T){
    Tax.TNRS.pol <- TNRS.taxa.accepted(Taxa.vector = TP.add, Taxonstand.test = T)
    saveRDS(Tax.TNRS.pol, "Resultats/ACA/Traits/Corresp_tax/Pollen_recheck_taxa_TNRS.Rds")
    }
  else{Tax.TNRS.pol <- readRDS("Resultats/ACA/Traits/Corresp_tax/Pollen_recheck_taxa_TNRS.Rds")}
  Tax.TNRS.pol <- Tax.TNRS.pol[Tax.TNRS.pol$match.score >= 0.99,]
  Tax.TNRS.pol.spe <- Tax.TNRS.pol[grepl(" ", Tax.TNRS.pol$Accepted_name),]
  Tax.TNRS.pol.oth <- Tax.TNRS.pol[!grepl(" ", Tax.TNRS.pol$Accepted_name),]
  
  TL.check <- rbind(TL.check, Tax.TNRS.pol.spe)
  
  #### Merge with vegetation taxa list ####
  TL.check <- TL.check[TL.check$match.score >= 0.95,]
  TL.check$Accepted_name <- TL.check$Accepted_name_TXSD
  TL.check <- TL.check[-ncol(TL.check)]
  
  names(TL.merge)[names(TL.merge) == "species"] <- "Name_submitted"
  TL.clean <- merge(TL.check, TL.merge, by = "Name_submitted", all = T)
  TL.clean$Accepted_family[is.na(TL.clean$Accepted_family)] <- TL.clean$family[is.na(TL.clean$Accepted_family)]
  TL.clean$Accepted_family[TL.clean$Accepted_family %in% ""] <- TL.clean$family[TL.clean$Accepted_family %in% ""]
  
  MV_TL.to.merge <- MV_TL[MV_TL$Taxa_level != 3,]
  MV_TL.to.merge <- MV_TL.to.merge[1:3]
  
  TL.clean <- full_join(TL.clean, MV_TL.to.merge, by = c("Accepted_name", "Accepted_family", "genus"))
  TL.clean <- as_tibble(TL.clean[TL.clean$Accepted_name != "",])
  TL.clean <- TL.clean[-c(1,3,4,8)]
  
  TL.clean[TL.clean == ""] <- NA
  TL.clean$Accepted_family[is.na(TL.clean$Accepted_family)] <- TL.clean$Accepted_family[match(TL.clean$genus[is.na(TL.clean$Accepted_family)], TL.clean$genus)]
  TL.clean$class[is.na(TL.clean$class)] <- TL.clean$class[match(TL.clean$Accepted_family[is.na(TL.clean$class)], TL.clean$Accepted_family)]
  TL.clean$order[is.na(TL.clean$order)] <- TL.clean$order[match(TL.clean$Accepted_family[is.na(TL.clean$order)], TL.clean$Accepted_family)]
  
  #### Clean genus, familly, variety... ####
  TL.clean$Accepted_name[which(substring(TL.clean$Accepted_name,1,2)=="x ")] <- substring(TL.clean$Accepted_name[which(substring(TL.clean$Accepted_name,1,2)=="x ")], 3)
  
  TL.clean$genus <- gsub("\\s.*", "", TL.clean$Accepted_name)
  # TL.clean$Accepted_family[which(TL.clean$Accepted_family == "")] <- TL.clean$family[which(TL.clean$Accepted_family == "")]
  #TL.clean <- TL.clean[which(lengths(gregexpr("\\W+", TL.clean$Accepted_name)) + 1 > 1),] # Remove genus only
  TL.clean <- TL.clean[!grepl("var\\.", TL.clean$Accepted_name),] # Remove variety 
  TL.clean <- TL.clean[!grepl("subsp\\.", TL.clean$Accepted_name),] # Remove variety 
  TL.clean <- TL.clean[!grepl("ssp\\.", TL.clean$Accepted_name),] # Remove variety 
  TL.clean <- TL.clean[-c(4)]
  
  Species.by.x <- gsub(".*\\s", "", TL.clean$Accepted_name)
  Species.by.x <- Species.by.x[which(substring(Species.by.x,1,1)=="x")]
  Species.by.x <- Species.by.x[grep("[^aeiouy]", substring(Species.by.x,2,2))]
  
  if(length(Species.by.x > 0)){
    for(i in 1:length(Species.by.x)){TL.clean$Accepted_name[grep(Species.by.x[i], TL.clean$Accepted_name)] <- paste(TL.clean$genus[grep(Species.by.x[i], TL.clean$Accepted_name)], substring(Species.by.x[i], 1, 1), substring(Species.by.x[i], 2), sep = " ")}
    }
  TL.clean <- TL.clean[c(3,2,4,1)]
  names(TL.clean) <- c("class", "family", "genus", "species")
  
  TL.clean <- TL.clean[order(TL.clean$genus),]
  TL.clean <- TL.clean[!is.na(TL.clean$family),] # Remove empty family 
  Not.duplicated.genus <- TL.clean$genus[!duplicated(TL.clean$genus) & !duplicated(TL.clean$genus, fromLast = T)]
  Genus.in.species <- TL.clean$species[TL.clean$species %in% TL.clean$genus]
  
  TL.clean <- TL.clean[!TL.clean$species %in% Genus.in.species[!Genus.in.species %in% Not.duplicated.genus],] # Remove duplicate genus
  TL.clean <- TL.clean[!duplicated(TL.clean[2:4]),]
  
  
  #### Remove unless familly (aquatic...) ####
  Remove.class = c("", NA, "Zygnematophyceae", "Ulvophyceae", "Ulvophyceae", "Bryopsida", "Bangiophyceae", "Chlorophyceae",
                   "Trebouxiophyceae", "Sphagnopsida", "Rhodophyceae", "Prasinophyceae", "Porphyridiophyceae", "Charophyceae",
                   "Polypodiopsida", "Pedinophyceae", "Marchantiopsida", "Lycopodiopsida", "Andreaeopsida", "Haplomitriopsida", "Anthocerotopsida",
                   "Florideophyceae", "Jungermanniopsida")
  Remove.famille <- c("", NA, "Aspleniaceae", "Athyriaceae", "Blechnaceae", "Osmundaceae", "Equisetaceae", 
                      "Pteridaceae", "Polypodiaceae", "Woodsiaceae", "Unknown", "Selaginellaceae",
                      "Lycopodiaceae", "Dryopteridaceae", "Dennstaedtiaceae", "Cystopteridaceae",
                      "Davalliaceae", "Gleicheniaceae", "Hymenophyllaceae", "Lindsaeaceae", "Ophioglossaceae",
                      "Marsileaceae", "Salviniaceae", "Thelypteridaceae", "Tectariaceae", "Lygodiaceae", "Onocleaceae"
                      )
  TL.clean <- TL.clean[setdiff(seq(1, nrow(TL.clean)), which(TL.clean$class %in% Remove.class)),]
  TL.clean <- TL.clean[setdiff(seq(1, nrow(TL.clean)), which(TL.clean$family %in% Remove.famille)),]
  TL.clean$class <- gsub("Equisetopsida", NA, TL.clean$class)
  TL.clean <- TL.clean[order(TL.clean$family),]
  TL.clean$class <- zoo::na.locf(TL.clean$class)
  
  #### Check family fake names ####
  TL.clean$family[which(TL.clean$family == "Chenopodiaceae")] <- "Amaranthaceae"
  TL.clean$family[which(TL.clean$family == "Viburnaceae")] <- "Adoxaceae"
  TL.clean$family[which(TL.clean$family == "Petiveriaceae")] <- "Phytolaccaceae"
  TL.clean$family[which(TL.clean$family == "Corbichoniaceae")] <- "Aizoaceae"
  TL.clean$family[which(TL.clean$family == "Compositae")] <- "Asteraceae"
  TL.clean$family[which(TL.clean$family == "Leguminosae")] <- "Fabaceae"
  
  #### Check if all genus are in the same family ####
  TL.clean$family[which(TL.clean$genus == "Ilex")] <- "Aquifoliaceae"
  TL.clean$family[which(TL.clean$genus == "Carex")] <- "Cyperaceae"
  TL.clean$family[which(TL.clean$genus == "Scilla")] <- "Asparagaceae"
  TL.clean$family[which(TL.clean$genus == "Ornithogalum")] <- "Asparagaceae"
  TL.clean$family[which(TL.clean$genus == "Rohdea")] <- "Asparagaceae"
  TL.clean$family[which(TL.clean$genus == "Bassia")] <- "Amaranthaceae"
  TL.clean$family[which(TL.clean$genus == "Saussurea")] <- "Asteraceae"
  TL.clean$family[which(TL.clean$genus == "Calligonum")] <- "Polygonaceae"
  TL.clean$family[which(TL.clean$genus == "Stellaria")] <- "Caryophyllaceae"
  TL.clean$family[which(TL.clean$genus == "Lasianthus")] <- "Rubiaceae"
  TL.clean$family[which(TL.clean$genus == "Cynoglossum")] <- "Boraginaceae"
  TL.clean$family[which(TL.clean$genus == "Symplocos")] <- "Symplocaceae"
  TL.clean$family[which(TL.clean$genus == "Dodartia")] <- "Mazaceae"
  TL.clean$family[which(TL.clean$genus == "Stauntonia")] <- "Lardizabalaceae"
  TL.clean$family[which(TL.clean$genus == "Holboellia")] <- "Lardizabalaceae"
  
  Pb.gen <- c()
  for(gen in unique(TL.clean$genus)){
    A <- TL.clean[TL.clean$genus == gen,]
    Fam <- unique(A$family)
    if(length(Fam) > 1){
      Pb.gen <- append(Pb.gen, gen)
    }
  }
  
  if(length(Pb.gen) > 0){stop("Some genus are in more than one family !!!")}
  
  #### Quercus case, phenology from TRY ####
  TRY.quercus.extract = F
  if(TRY.quercus.extract == T){
    if(exists("TRY") == F){TRY <- vroom::vroom(file = paste(DB.path, "Traits/TRY_dec_2019/8066.csv", sep = ""))}
    TRY.quercus <- TRY[grep("Quercus", TRY$AccSpeciesName),]
    TRY.quercus <- TRY.quercus[TRY.quercus$DataID == 42,]
    TRY.quercus <- TRY.quercus[,c(7, 15)] 
    Quercus.type <- unique(data.frame(TRY.quercus))
    Quercus <- read.csv(file = "Import/World_DB/Taxonomie/Quercus_phenology.csv", sep = ",", header = T) 
    Quercus <- unique(Quercus)
    Quercus$Phenology <- gsub("evergreen", "Quercus evergreen", Quercus$Phenology)
    Quercus$Phenology <- gsub("deciduous", "Quercus deciduous", Quercus$Phenology)
    write.csv(Quercus, file = "Import/World_DB/Taxonomie/Quercus_phenology.csv")
    }
  if(TRY.quercus.extract == F){Quercus <- read.csv(file = "Import/World_DB/Taxonomie/Quercus_phenology.csv", sep = ",", header = T)}
  
  #### Cereales case ####
  Cereal.clean = F
  if(Cereal.clean == T){
    library(stringr)
    Cerealia <- read.csv(file = "Import/World_DB/Taxonomie/Cerealia.csv", sep = "\t", header = F) # cropsreview.com
    Cerealia_usda <- read.csv(file = "Import/World_DB/Taxonomie/Cerealia-USDA.csv", sep = "\t", header = F) # cropsreview.com
    Cerealia_usda <- str_split(Cerealia_usda$V2, " ")
    Cerealia_usda <-  unlist(lapply(Cerealia_usda, function(x) paste(x[1], x[2])))
    C2 <- gsub("\\(([^()]*)\\)|.", "\\1", Cerealia, perl=T)
    C2 <- gsub("\\(([^()]*)\\)|.", "\\1", paste(gsub("\\)", ",\\)", Cerealia$V1), collapse = " "), perl = F)
    C2 <- str_split(C2, ",")
    C2[[1]][which(substring(C2[[1]],1,1)==" ")] <- substring(C2[[1]][which(substring(C2[[1]],1,1)==" ")], 2)
    C2[[1]][which(substring(C2[[1]],1,1)=="x")] <- substring(C2[[1]][which(substring(C2[[1]],1,1)=="x")], 2)
    C2 <- C2[[1]]
    C2 <- C2[C2 != ""]
    C2 <- gsub("mainly ", "", C2)
    C2 <- gsub(" spp\\.", "", C2)
    C2 <- gsub(" sp\\.", "", C2)
    C2 <- unique(c(C2, Cerealia_usda))
    C3 <- str_split(C2, " ")
    Cg <- C2[lapply(C3, function(x) length(x)) == 1]
    Cs <- C2[lapply(C3, function(x) length(x)) == 2]
    C4 <- data.frame(Genre = c(Cg, rep(NA, (length(C2) - length(Cg)))))
    C4 <- cbind(C4, Espece = c(rep(NA, (length(C2) - length(Cs))), Cs))
    C4$Subgenus <- "Cerealia"
    Cerealia <- C4
    write.csv(C4, file = "Import/World_DB/Taxonomie/Cerealia-type.csv")
  }
  if(Cereal.clean == F){Cerealia <- read.csv(file = "Import/World_DB/Taxonomie/Cerealia-type.csv", sep = ",", header = T)}
  
  #### Rosaceae case ####
  Rosaceae.cal = T
  if(Rosaceae.cal == T){
    GrowthForm.ACA <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_GrowthForm_ACA.Rds", sep = ""))
    Rosaceae <- GrowthForm.ACA[GrowthForm.ACA$family == "Rosaceae", c(4,1,5)]
    if(exists("Herbier") == F){Herbier <- data.frame(read_ods("/home/lucas.dugerdil/Documents/Recherche/Herbier/Herbier_numerique/Index.ods"))}
    Rosaceae2 <- Herbier[Herbier$Family == "Rosaceae", c(4,5,18)]
    Rosaceae2 <- Rosaceae2[!is.na(Rosaceae2$species),]
    names(Rosaceae2) <- c("genus", "species", "GrowthForm")
    Rosaceae <- rbind(Rosaceae, Rosaceae2)
    Rosaceae$GrowthForm[Rosaceae$GrowthForm == "Unknown"] <- NA
    Rosaceae <- Rosaceae[!is.na(Rosaceae$GrowthForm),]
    Rosaceae <- data.table::setDT(Rosaceae)[, GrowthForm := zoo::na.locf(zoo::na.locf(GrowthForm, na.rm = FALSE), fromLast = TRUE), genus]
    Rosaceae <- as_tibble(Rosaceae)
    Rosaceae$Other.clade <- paste("Rosaceae (", Rosaceae$GrowthForm, ")", sep = "")
    Rosaceae$Other.clade[Rosaceae$Other.clade == "Rosaceae (NA)"] <- NA
    saveRDS(Rosaceae, "Import/World_DB/Taxonomie/Rosaceae_type.Rds")
  }
  
  
  #### Ephedra case ####
  Ephedra.cal = T
  if(Ephedra.cal == T){
    GrowthForm.ACA <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_GrowthForm_ACA.Rds", sep = ""))
    Ephedra <- GrowthForm.ACA[GrowthForm.ACA$family == "Ephedra", c(4,1,5)]
    if(exists("Herbier") == F){Herbier <- data.frame(read_ods("/home/lucas.dugerdil/Documents/Recherche/Herbier/Herbier_numerique/Index.ods"))}
    Ephedra2 <- Herbier[Herbier$Family == "Ephedra", c(4,5,18)]
    Ephedra2 <- Ephedra2[!is.na(Ephedra2$species),]
    names(Ephedra2) <- c("genus", "species", "GrowthForm")
    Ephedra <- rbind(Ephedra, Ephedra2)
    Ephedra$GrowthForm[Ephedra$GrowthForm == "Unknown"] <- NA
    Ephedra <- Ephedra[!is.na(Ephedra$GrowthForm),]
    Ephedra <- data.table::setDT(Ephedra)[, GrowthForm := zoo::na.locf(zoo::na.locf(GrowthForm, na.rm = FALSE), fromLast = TRUE), genus]
    Ephedra <- as_tibble(Ephedra)
    Ephedra$Other.clade <- paste("Ephedra (", Ephedra$GrowthForm, ")", sep = "")
    Ephedra$Other.clade[Ephedra$Other.clade == "Ephedra (NA)"] <- NA
    saveRDS(Ephedra, "Import/World_DB/Taxonomie/Ephedra_type.Rds")
  }
  
  #### Add the subfamily / subgenera ####
  Aster <- read.csv(file = "Import/World_DB/Taxonomie/Asteraceae.csv", sep = "\t", header = F)               # from NCBI Nat. Cent. Biotech. Info
  Pinus <- read.csv(file = "Import/World_DB/Taxonomie/Pinus_subgenre.csv", sep = ",", header = T)            # from NCBI Nat. Cent. Biotech. Info
  A = Aster[match(TL.clean$genus, Aster$V1),2]
  B = Pinus[match(TL.clean$species, Pinus$Espece),2]
  D = Quercus[match(TL.clean$species, Quercus$Espece),3]
  E = Cerealia[match(TL.clean$genus, Cerealia$Genre),4]
  G = Cerealia[match(TL.clean$species, Cerealia$Espece),4]
  C = data.frame(A, B, D, E, G)
  C <- mutate(C, Other.clade = coalesce(A,B,D,E,G))
  TL.clean$Other.clade <- C$Other.clade
  
  #### APG IV clean #### 
  APGIV <- as_tibble(read.table(file = "Import/World_DB/Taxonomie/APGIV/apg4.txt", sep = "\t", header = T))  # from APG III (wikipedia)
  APGIV <- APGIV[-c(3,6)]
  APGIV <- rbind(APGIV, c("Asphodelaceae", "family", "Liliales", ""))
  # APGIV <- rbind(APGIV, c("Xanthorrhoeaceae", "family", "Liliales", ""))
  APGIV <- rbind(APGIV, c("Cordiaceae", "family", "Boraginales", ""))
  APGIV <- rbind(APGIV, c("Heliotropiaceae", "family", "Boraginales", ""))
  APGIV <- rbind(APGIV, c("Namaceae", "family", "Boraginales", ""))
  APGIV <- rbind(APGIV, c("Hydrophyllaceae", "family", "Boraginales", ""))
  APGIV <- rbind(APGIV, c("Ehretiaceae", "family", "Boraginales", ""))
  APGIV <- rbind(APGIV, c("Ixiolirionaceae", "family", "Asparagales", ""))
  
  APGIV.clean <- as_tibble(setNames(data.frame(APGIV[APGIV$taxonRank == "family", c(3,1)]), c("order", "family"))) 
  APGIV.clade <-  as_tibble(setNames(data.frame(APGIV[APGIV$taxonRank == "order", c(3,1)]), c("kingdom", "order"))) 
  APGIV.clean <- merge(APGIV.clean, APGIV.clade, by = "order")
  APGIV.king <-  as_tibble(setNames(data.frame(APGIV[APGIV$taxonRank == "kingdom", c(3,1)]), c("kingdom", "clade"))) 
  APGIV.clean <- full_join(APGIV.clean, APGIV.king, by = "kingdom")
  
  APGIV.usage <- APGIV[!APGIV$acceptedNameUsage == "",]

  TL.APG <- left_join(TL.clean, APGIV.clean, by = "family")
  TL.APG <- TL.APG[c(7,8,6,5,2,3,4)]
  TL.APG$Subreign <- NA 
  TL.APG$Subreign[!is.na(TL.APG$order)] <- "Angiospermes"
  
  TL.APG$order[TL.APG$family == "Cycadaceae"] <- "Cycadales"
  TL.APG$order[TL.APG$family == "Ginkgoaceae"] <- "Ginkgoales"
  TL.APG$order[TL.APG$family == "Araucariaceae"] <- "Pinales"
  TL.APG$order[TL.APG$family == "Cupressaceae"] <- "Pinales"
  TL.APG$order[TL.APG$family == "Pinaceae"] <- "Pinales"
  TL.APG$order[TL.APG$family == "Podocarpaceae"] <- "Pinales"
  TL.APG$order[TL.APG$family == "Ephedraceae"] <- "Pinales"
  TL.APG$order[TL.APG$family == "Gnetaceae"] <- "Pinales"
  TL.APG$order[TL.APG$family == "Taxaceae"] <- "Pinales"
  TL.APG$order[TL.APG$family == "Sciadopityaceae"] <- "Pinales"
  
  TL.APG$Subreign[is.na(TL.APG$kingdom)] <- "Gymnospermes"
  TL.APG$kingdom[is.na(TL.APG$kingdom)] <- TL.APG$order[is.na(TL.APG$kingdom)]
  TL.APG <- TL.APG[!duplicated(TL.APG),]
  
  TL.APG <- TL.APG[c(8,1,3:7)]
  TL.ACA <- TL.APG[c(2,5,4,6,7)]
  
  
  #### Export ####
  saveRDS(TL.ACA,  paste(DB.path, "Vegetation/Occurences/Merge_DB/Taxon_list_ACA_V1.Rds", sep = ""))
  saveRDS(TL.APG,  paste(DB.path, "Vegetation/Occurences/Merge_DB/Taxon_list_ACA_full_taxo.Rds", sep = ""))
  write.csv(TL.ACA,  file = paste(DB.path, "Vegetation/Occurences/Merge_DB/Taxon_list_ACA_V1.csv", sep = ""))
  write.csv(TL.APG,  file = paste(DB.path, "Vegetation/Occurences/Merge_DB/Taxon_list_ACA_full_taxo.csv", sep = ""))
  write.csv(APGIV.clean,  file = "Resultats/World_DB/Taxonomie/APG_IV_clean.csv")
  
  #### Stats ####
  print(paste("ACA Species number from GBIF:", length(TL.GBIF$species)))
  print(paste("ACA Species number from BIEN:", length(TL.BIEN)))
  print(paste("ACA Species number from sPlot:", length(MV_TL$Accepted_name)))
  print(paste("ACA Species number not accepted:", length(TL.merge$Name_submitted)))
  print(paste("ACA Species number TNRS/Taxonstand-cheked:", length(TL.check$Accepted_name)))
  print(paste("ACA Species number monocot + dycot:", length(TL.clean$species)))
  print(paste("ACA Genus number:", length(unique(TL.clean$genus))))
  print(paste("ACA Family number:", length(unique(TL.clean$family))))
  }

Taxa.convertion = F
if(Taxa.convertion == T){
  #### Import ####
  Table.conv.taxa_sl <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/Table.conv.taxa_sl.Rds")
  Table.conv.taxa_ss <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/Table.conv.taxa_ss.Rds")
  
  if(exists("TL.ACA") == F){
    TL.ACA <- readRDS(paste(DB.path, "Vegetation/Occurences/Merge_DB/Taxon_list_ACA_V1.Rds", sep = ""))
    }
  
  #### TNRS check table corresp pollen #### 
  TNRS.TPT = F
  if(TNRS.TPT == T){
    Table.Taxon <- read.csv(file="Import/World_DB/Pollen/Odile_DB/Corresp_name_full_V10.csv", sep=",",dec=".", row.names = 1,  header=T, stringsAsFactors = F)
    TPT <- Table.Taxon
    Espece.cor <- TNRS.taxa.accepted(TPT$Espece, Taxonstand.test = T)
    
    TPT$Espece[which(TPT$Espece %in% Espece.cor$Name_submitted)] <- Espece.cor$Accepted_name_TXSD[match(TPT$Espece[which(TPT$Espece %in% Espece.cor$Name_submitted)], Espece.cor$Name_submitted)]
    Genre.cor <- TNRS.taxa.accepted(unique(TPT$Genre), Taxonstand.test = F)
    TPT$Genre[which(TPT$Genre %in% Genre.cor$Name_submitted)] <- Genre.cor$Accepted_name[match(TPT$Genre[which(TPT$Genre %in% Genre.cor$Name_submitted)], Genre.cor$Name_submitted)]
    
    saveRDS(TPT, "Import/World_DB/Pollen/Odile_DB/Corresp_name_full_V10_TNRS_check.csv")
    
    }
  else{
    TPT <- readRDS("Import/World_DB/Pollen/Odile_DB/Corresp_name_full_V10_TNRS_check.csv")
    }
  
  #### Convertion type pollen #### 
  
  Table.conv.taxa_sl <- unique(Table.conv.taxa_sl$New.name)[order(unique(Table.conv.taxa_sl$New.name))]
  Table.conv.taxa_ss <- unique(Table.conv.taxa_ss$New.name)[order(unique(Table.conv.taxa_ss$New.name))]
  
  TP.sl <- setNames(data.frame(matrix(NA, ncol = length(Table.conv.taxa_sl), nrow = nrow(TL.ACA))), Table.conv.taxa_sl)
  TP.sl <- cbind(TL.ACA, TP.sl)
  for(i in 6:ncol(TP.sl)){
    Type.p <- names(TP.sl)[i]
    if(Type.p %in% TPT$Espece){TP.sl[which(TP.sl$species == Type.p),Type.p] <- T}
    if(Type.p %in% TPT$Subgenus){TP.sl[which(TP.sl$Other.clade == Type.p),Type.p] <- T}
    if(Type.p %in% TPT$Genre){TP.sl[which(TP.sl$genus == Type.p),Type.p] <- T}
    if(Type.p %in% TPT$Subfamille){TP.sl[which(TP.sl$Other.clade == Type.p),Type.p] <- T}
    if(Type.p %in% TPT$Famille){TP.sl[which(TP.sl$family == Type.p),Type.p] <- T}
    }
  
  TP.sl[6:ncol(TP.sl)][is.na(TP.sl[6:ncol(TP.sl)])] <- F
  
  TP.ss <- setNames(data.frame(matrix(NA, ncol = length(Table.conv.taxa_ss), nrow = nrow(TL.ACA))), Table.conv.taxa_ss)
  TP.ss <- cbind(TL.ACA, TP.ss)
  for(i in 6:ncol(TP.ss)){
    Type.p <- names(TP.ss)[i]
    if(Type.p %in% TPT$Espece){TP.ss[which(TP.ss$species == Type.p),Type.p] <- T}
    if(Type.p %in% TPT$Subgenus){TP.ss[which(TP.ss$Other.clade == Type.p),Type.p] <- T}
    if(Type.p %in% TPT$Genre){TP.ss[which(TP.ss$genus == Type.p),Type.p] <- T}
    if(Type.p %in% TPT$Subfamille){TP.ss[which(TP.ss$Other.clade == Type.p),Type.p] <- T}
    if(Type.p %in% TPT$Famille){TP.ss[which(TP.ss$family == Type.p),Type.p] <- T}
  }
  
  TP.ss[6:ncol(TP.ss)][is.na(TP.ss[6:ncol(TP.ss)])] <- F
  
  saveRDS(TP.ss, file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_plant_type_bool_ss.Rds")
  saveRDS(TP.sl,file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_plant_type_bool_sl.Rds")
  
  #### Pollen type AP / NAP ####
  TPT$PT_sl <- NA
  TPT$PT_sl[which(TPT$Level == 5)] <- TPT$Espece[which(TPT$Level == 5)] 
  TPT$PT_sl[which(TPT$Level == 4)] <- TPT$Subgenus[which(TPT$Level == 4)] 
  TPT$PT_sl[which(TPT$Level == 3)] <- TPT$Genre[which(TPT$Level == 3)] 
  TPT$PT_sl[which(TPT$Level == 2)] <- TPT$Subfamille[which(TPT$Level == 2)] 
  TPT$PT_sl[which(TPT$Level == 1)] <- TPT$Famille[which(TPT$Level == 1)] 
  TPT$PT_ss <- NA
  TPT$PT_ss[which(!TPT$Espece == "")] <- TPT$Espece[which(!TPT$Espece == "")] 
  TPT$PT_ss[intersect(which(!TPT$Subgenus == ""), which(is.na(TPT$PT_ss)))] <- TPT$Subgenus[intersect(which(!TPT$Subgenus == ""), which(is.na(TPT$PT_ss)))] 
  TPT$PT_ss[intersect(which(!TPT$Genre == ""), which(is.na(TPT$PT_ss)))] <- TPT$Genre[intersect(which(!TPT$Genre == ""), which(is.na(TPT$PT_ss)))] 
  TPT$PT_ss[intersect(which(!TPT$Subfamille == ""), which(is.na(TPT$PT_ss)))] <- TPT$Subfamille[intersect(which(!TPT$Subfamille == ""), which(is.na(TPT$PT_ss)))] 
  TPT$PT_ss[intersect(which(!TPT$Famille == ""), which(is.na(TPT$PT_ss)))] <- TPT$Famille[intersect(which(!TPT$Famille == ""), which(is.na(TPT$PT_ss)))] 
  
  #### Old aggregation method ####
  Old.convert = F
  if(Old.convert == T){

    # TPT$PT_sl <- factor(TPT$PT_sl)
    # A <- levels(TPT$PT_sl)
    # A[which(A %in% Table.conv.taxa_sl$Old.name)] <- Table.conv.taxa_sl$New.name[which(Table.conv.taxa_sl$Old.name %in% A)]
    
    #### Convertion Taxa bota -> type pollen #### 
    TL.ACA$PT_sl <- NA
    TL.ACA$PT_sl[na.omit(match(TPT$Espece[which(TPT$Level == 5)], TL.ACA$species))] <- TPT$Espece[which(TPT$Level == 5)][!is.na(match(TPT$Espece[which(TPT$Level == 5)], TL.ACA$species))]
    TL.ACA$PT_sl[na.omit(intersect(which(TL.ACA$Other.clade %in% TPT$Subgenus[which(TPT$Level == 4)]), which(is.na(TL.ACA$PT_sl))))] <- TL.ACA$Other.clade[na.omit(intersect(which(TL.ACA$Other.clade %in% TPT$Subgenus[which(TPT$Level == 4)]), which(is.na(TL.ACA$PT_sl))))]
    TL.ACA$PT_sl[na.omit(intersect(which(TL.ACA$genus %in% TPT$Genre[which(TPT$Level == 3)]), which(is.na(TL.ACA$PT_sl))))] <- TL.ACA$genus[na.omit(intersect(which(TL.ACA$genus %in% TPT$Genre[which(TPT$Level == 3)]), which(is.na(TL.ACA$PT_sl))))]
    TL.ACA$PT_sl[na.omit(intersect(which(TL.ACA$Other.clade %in% TPT$Subfamille[which(TPT$Level == 2)]), which(is.na(TL.ACA$PT_sl))))] <- TL.ACA$Other.clade[na.omit(intersect(which(TL.ACA$Other.clade %in% TPT$Subfamille[which(TPT$Level == 2)]), which(is.na(TL.ACA$PT_sl))))]
    TL.ACA$PT_sl[na.omit(intersect(which(TL.ACA$family %in% TPT$Famille[which(TPT$Level == 1)]), which(is.na(TL.ACA$PT_sl))))] <- TL.ACA$family[na.omit(intersect(which(TL.ACA$family %in% TPT$Famille[which(TPT$Level == 1)]), which(is.na(TL.ACA$PT_sl))))]
    
    TL.ACA$PT_ss <- NA
    TL.ACA$PT_ss[which(TL.ACA$species %in% TPT$Espece[which(!TPT$Espece == "")])] <- TL.ACA$species[which(TL.ACA$species %in% TPT$Espece[which(!TPT$Espece == "")])]  
    TL.ACA$PT_ss[intersect(which(TL.ACA$Other.clade %in% TPT$Subgenus[which(!TPT$Subgenus == "")]), which(is.na(TL.ACA$PT_ss)))] <- TL.ACA$Other.clade[intersect(which(TL.ACA$Other.clade %in% TPT$Subgenus[which(!TPT$Subgenus == "")]), which(is.na(TL.ACA$PT_ss)))]  
    TL.ACA$PT_ss[intersect(which(TL.ACA$genus %in% TPT$Genre[which(!TPT$Genre == "")]), which(is.na(TL.ACA$PT_ss)))] <- TL.ACA$genus[intersect(which(TL.ACA$genus %in% TPT$Genre[which(!TPT$Genre == "")]), which(is.na(TL.ACA$PT_ss)))]  
    TL.ACA$PT_ss[intersect(which(TL.ACA$Other.clade %in% TPT$Subfamille[which(!TPT$Subfamille == "")]), which(is.na(TL.ACA$PT_ss)))] <- TL.ACA$Other.clade[intersect(which(TL.ACA$Other.clade %in% TPT$Subfamille[which(!TPT$Subfamille == "")]), which(is.na(TL.ACA$PT_ss)))]  
    TL.ACA$PT_ss[intersect(which(TL.ACA$family %in% TPT$Famille[which(!TPT$Famille == "")]), which(is.na(TL.ACA$PT_ss)))] <- TL.ACA$family[intersect(which(TL.ACA$family %in% TPT$Famille[which(!TPT$Famille == "")]), which(is.na(TL.ACA$PT_ss)))]  
    
    #### Check for new pollen names (TNRS) #### 
    for(i in 1:nrow(TL.ACA)){
      if(is.na(Table.conv.taxa_sl$New.name[match(TL.ACA$PT_sl[i], Table.conv.taxa_sl$Old.name)]) == F){
        TL.ACA$PT_sl[i] <- Table.conv.taxa_sl$New.name[match(TL.ACA$PT_sl[i], Table.conv.taxa_sl$Old.name)]}
      
      if(is.na(Table.conv.taxa_ss$New.name[match(TL.ACA$PT_ss[i], Table.conv.taxa_ss$Old.name)]) == F){
        TL.ACA$PT_ss[i] <- Table.conv.taxa_ss$New.name[match(TL.ACA$PT_ss[i], Table.conv.taxa_ss$Old.name)]}
      }
    
    pout <- TL.ACA[TL.ACA$PT_ss %in% Table.conv.taxa_ss$New.name,]
    
    missing.taxa.zeub <- Table.conv.taxa_ss[!Table.conv.taxa_ss$New.name %in% pout$PT_ss,]
    #### Convertion totale /!\ SLOW ####
    # Bool.cores <- F
    # if(Bool.cores == T){
    #   Naam <- paste(TL.ACA$family, TL.ACA$genus, TL.ACA$species, TL.ACA$Other.clade)
    #   Naam_ss <- sapply(unique(TPT$PT_ss)[!is.na(unique(TPT$PT_ss))], function (y) sapply(Naam, function (x) grepl(y,x)))
    #   row.names(Naam_ss) <- row.names(TL.ACA)
    #   TL.ACA2_PTss <- as_tibble(cbind(TL.ACA[,1:4], Naam_ss))
    #   
    #   Naam_sl <- sapply(unique(TPT$PT_sl)[!is.na(unique(TPT$PT_sl))], function (y) sapply(Naam, function (x) grepl(y,x)))
    #   row.names(Naam_sl) <- row.names(TL.ACA)
    #   TL.ACA2_PTsl <- as_tibble(cbind(TL.ACA[,1:4], Naam_sl))
    # }
    
    
    #### Convertion MV avec taxon list ACA valide /!\ SLOW ####
    # MV.conv.T <- F
    # if(MV.conv.T == T){
    #   MV_TL <- MV_TL[MV_TL$Accepted_name %in% TL.ACA$species | MV_TL$Accepted_name %in% TL.ACA$genus | MV_TL$Accepted_name %in% TL.ACA$family,]
    #   MV_TL <- MV_TL[!duplicated(MV_TL),]
    #   Total.list.of.accepted.name <- unique(c(MV_TL$Accepted_name, MV_TL$Accepted_family, MV_TL$genus))
    #   
    #   MV.releve <- MV.releve[MV.releve$Species %in% Total.list.of.accepted.name,]
    #   
    #   all_abund = aggregate(FA ~ ID, MV.releve, sum)
    #   colnames(all_abund)[2] = "FA_tot"
    #   MV.releve = merge(MV.releve, all_abund, by = "ID")
    #   MV.releve$FA = MV.releve$FA / MV.releve$FA_tot
    #   MV.releve <- MV.releve[-c(4)]
    #   
    #   MV.releve = aggregate(FA ~ ID + Species, MV.releve, sum)
    #   MV.releve = tidyr::spread(MV.releve, Species, FA, drop = T, fill = 0.00)
    #   row.names(MV.releve) <- MV.releve$ID
    #   MV.releve <- MV.releve[-c(1)]
    #   
    #   saveRDS(MV.releve, file = "Resultats/ACA/Vegetation/MV_ACA_plot_only_spermatophyte.Rds")
    #   saveRDS(MV_TL, file = "Resultats/ACA/Vegetation/MV_ACA_taxonlist.Rds")
    # }
    # else{
    #   MV.releve <- readRDS("Resultats/ACA/Vegetation/MV_ACA_plot_only_spermatophyte.Rds")
    #   MV_TL <- readRDS("Resultats/ACA/Vegetation/MV_ACA_taxonlist.Rds")
    # }
    
    
    }
  
  #### Export new taxon ####
  saveRDS(TPT, file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_pollen_type.Rds")
  # saveRDS(TL.ACA,file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_plant_type.Rds")
  write.csv(TPT, file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_pollen_type.csv")
  # write.csv(TL.ACA,file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_plant_type.csv")
  
  # if(Bool.cores == T){
    # saveRDS(TL.ACA2_PTss, file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_plant_type_bool_ss.Rds")
    # saveRDS(TL.ACA2_PTsl,file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_plant_type_bool_sl.Rds")
    # write.csv(TL.ACA2_PTss, file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_plant_type_bool_ss.csv")
    # write.csv(TL.ACA2_PTsl,file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_plant_type_bool_sl.csv")
  # }
  
  #### Stats report ####
  # length(which(is.na(TL.ACA$PT_sl)))
  # length(which(is.na(TL.ACA$PT_ss)))
  length(unique(TL.ACA$family))
  length(unique(TL.ACA$genus))
  length(unique(TL.ACA$species))
  # length(unique(TL.ACA$PT_sl))
  # length(unique(TL.ACA$PT_ss))
  # PT_sl.sans.bota <- TPT[setdiff(seq(1, nrow(TPT)), which(TPT$PT_sl %in% TL.ACA$PT_sl)),]
  # PT_ss.sans.bota <- TPT[setdiff(seq(1, nrow(TPT)), which(TPT$PT_ss %in% TL.ACA$PT_ss)),]
}

Climat.extract = F
if(Climat.extract == T){
  #### Import ####
  source("Scripts/Climat_extract.R")
  MPS.ACA.Cord <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Coord.Rds")
  MV.ACA.Cord <- readRDS("Resultats/ACA/Vegetation/MV_ACA_cord.Rds")
  
  #### Extract sites pollen ####
  MPS.ACA.Cord.noAlt <- MPS.ACA.Cord[which(is.na(MPS.ACA.Cord$Alt)), -c(5:7)]
  MPS.ACA.Cord.noAlt <- Clim.param.extraction(M = MPS.ACA.Cord.noAlt, Clim.cal = F, Altitude = T)
  MPS.ACA.Cord$Alt[which(is.na(MPS.ACA.Cord$Alt))] <- MPS.ACA.Cord.noAlt$Altitude
  
  MPS.ACA.Clim_chel <- Clim.param.extraction(M = MPS.ACA.Cord[,-c(5:7)],
                                              All.param = F, Chelsa = T,
                                              Altitude = T, Aridity = T)
  
  MPS.ACA.Clim_wc <- Clim.param.extraction(M = MPS.ACA.Cord[,-c(5:7)],
                                             All.param = F, Chelsa = F,
                                             Altitude = T, Aridity = T)
  
  MPS.ACA.Biom <- Clim.param.extraction(M = MPS.ACA.Cord[,-c(5:7)], Clim.cal = F, Biome = T,
                                             All.param = F, Chelsa = T,
                                             Altitude = F, Aridity = T)
  
  MPS.ACA.Land <- Clim.param.extraction(M = MPS.ACA.Cord[,-c(5:7)], Clim.cal = F, Biome = F,
                                        All.param = F, Chelsa = F, Land.cover = T,
                                        Altitude = F, Aridity = F)
  
  saveRDS(MPS.ACA.Clim_chel, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Clim_chel.Rds")
  saveRDS(MPS.ACA.Clim_wc, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Clim_wc.Rds")
  saveRDS(MPS.ACA.Land, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_LandCover.Rds")
  saveRDS(MPS.ACA.Biom, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Biom.Rds")
  
  #### Extract sites veget ####
  MV.ACA.Cord <- data.frame(MV.ACA.Cord)
  row.names(MV.ACA.Cord) <- MV.ACA.Cord$ID
  MV.ACA.Clim_chel <- Clim.param.extraction(M = MV.ACA.Cord,
                                             All.param = F, Chelsa = T,
                                             Altitude = T, Aridity = T)
  
  MV.ACA.Clim_wc <- Clim.param.extraction(M = MV.ACA.Cord,
                                           All.param = F, Chelsa = F,
                                           Altitude = T, Aridity = T)
  
  MV.ACA.Biom <- Clim.param.extraction(M = MV.ACA.Cord, Clim.cal = F, Biome = T,
                                        All.param = F, Chelsa = T,
                                        Altitude = F, Aridity = T)
  
  MV.ACA.Land <- Clim.param.extraction(M = MV.ACA.Cord, Clim.cal = F, Biome = F,
                                        All.param = F, Chelsa = T, Land.cover = T,
                                        Altitude = F, Aridity = F)
  
  saveRDS(MV.ACA.Land, "Resultats/ACA/Vegetation/DBACA_LandCover.Rds")
  saveRDS(MV.ACA.Clim_chel, "Resultats/ACA/Vegetation/DBACA_Clim_chel.Rds")
  saveRDS(MV.ACA.Clim_wc, "Resultats/ACA/Vegetation/DBACA_Clim_wc.Rds")
  saveRDS(MV.ACA.Biom, "Resultats/ACA/Vegetation/DBACA_Biom.Rds")
  #### Stats ####
  Ref.num <- table(MPS.ACA.Cord$Ref)
  Other.sites <- sum(Ref.num[which(Ref.num < 10)])
  Ref.num <- data.frame(Ref.num[which(Ref.num >= 10)])
  Ref.num <- rbind(Ref.num, data.frame(Var1 = "Other studies", Freq = as.numeric(Other.sites)))
  Ref.num$Var1 <- gsub(" ", "_", Ref.num$Var1)
  Ref.num$Var1 <- gsub("et_al_", "", Ref.num$Var1)
  Ref.num <- data.frame(Ref.num[order(Ref.num$Freq, decreasing = T),])
  names(Ref.num) <- c("Reference", "N° sites")
  library(stringr)
  Date <- purrr::map(str_extract_all(Ref.num$Reference, "[0-9]+"), 1)
  Date[sapply(Date, function(x) is.null(x))] <- 0
  Ref.num$Date <- unlist(Date)
  Ref.num$Name <- sub("_.*","",Ref.num$Reference)
  Ref.num$Latex <- paste("\\citet{", Ref.num$Name, "_", Ref.num$Date, "}", sep = "")
  write.csv(Ref.num[c(5,2)], "Resultats/ACA/Traits/Pollen_biblio/Biblio_latex_pollen.csv")
  }

TRY.BROT.BIEN.extract.ACA = F
if(TRY.BROT.BIEN.extract.ACA == T){
  #### Import taxa list ####
  TL.ACA <- readRDS(paste(DB.path, "Vegetation/Occurences/Merge_DB/Taxon_list_ACA_V1.Rds", sep = ""))
  
  #### Import BIEN  ####
  BIEN.trait.extract = T
  if(BIEN.trait.extract == T){
    library(BIEN)
    Trait.list.Barboni2004 <- c("leaf thickness", "leaf area", "leaf life span")
    TraitID.list.Barboni2004 <- c(46, 3115, 37, 343, 22) # TRY (thick, area /!\ plein de différents, leaf phenol, habit, pathway)
    TraitID.list.Frenette2013 <- c(155, 4, 3115, 47, 89, 329) # TRY (onset of flowering, stem density, SLA, LDMC, leaf C13, clonality )
    BIEN.tr.ACA <- BIEN_trait_species(species = TL.ACA$species)
    saveRDS(BIEN.tr.ACA,  paste(DB.path, "Vegetation/Occurences/BIEN/BIEN_ACA_trait.Rds", sep = ""))
    }
  
  #### Import list taxa ####
  library("vroom")
  if(exists("TL.TRY") == F){TL.TRY <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_plantlist.Rds", sep = ""))}
  if(exists("TRY") == F){TRY <- vroom::vroom(file = paste(DB.path, "Traits/TRY_dec_2019/8066.csv", sep = ""))}
 
  #### Extract TRY for taxa from ACA ####
  TL.TRY.ACA <- intersect(TL.TRY$AccSpeciesName, TL.ACA$species)
  TRY.ACA <- TRY[which(TRY$AccSpeciesName %in% TL.TRY.ACA),]

  #### Import BROT 2.0 + homog ####
  brot <- read.csv("Import/World_DB/Traits/BROT2.0/BROT2_dat.csv", row.names = 1, stringsAsFactors = F)
  brot.sou <- read.csv("Import/World_DB/Traits/BROT2.0/BROT2_sou.csv", row.names = 1, stringsAsFactors = F, encoding = "UTF-8")
  brot.tax <- read.csv("Import/World_DB/Traits/BROT2.0/BROT2_tax.csv", row.names = 1, stringsAsFactors = F, encoding = "UTF-8")
  brot.syn <- read.csv("Import/World_DB/Traits/BROT2.0/BROT2_syn.csv", stringsAsFactors = F, encoding = "UTF-8")
  
  TL.BROT.ACA <- intersect(brot$Taxon, TL.ACA$species)
  brot.ACA <- brot[which(brot$Taxon %in% TL.BROT.ACA),]
  
  #### Import Chinese trait DB + homog ####
  CPT.trait.conti <- as_tibble(read.csv("Import/China/Traits/China_Plant_TraitsDB_csv/Hard Traits.csv", row.names = 1, stringsAsFactors = F))
  CPT.trait.unconti <- as_tibble(read.csv("Import/China/Traits/China_Plant_TraitsDB_csv/Morphometric traits.csv", row.names = 1, stringsAsFactors = F))
  CPT.trait.pathway <- as_tibble(read.csv("Import/China/Traits/China_Plant_TraitsDB_csv/Photo Pathway.csv"))
  CPT.PFT <- as_tibble(read.csv("Import/China/Traits/China_Plant_TraitsDB_csv/PFT data.csv"))
  CPT.species.samp <- as_tibble(read.csv("Import/China/Traits/China_Plant_TraitsDB_csv/Species translations.csv"))
  
  CPT.species.samp$Species <- paste(CPT.species.samp$ACCEPTED.GENUS, CPT.species.samp$ACCEPTED.SPECIES, sep = " ")
  TL.CPT.ACA <- intersect(CPT.species.samp$Species, TL.ACA$species)
  CPT.species.samp <- CPT.species.samp[which(CPT.species.samp$Species %in% TL.CPT.ACA),]
    
  CPT.aggreg <- function(M, Aggeg.by, Aggreg){
    M <- merge(CPT.species.samp[c(Aggeg.by, "Species")], M[which(M[[Aggeg.by]] %in% CPT.species.samp[[Aggeg.by]]),], by = Aggeg.by)
    M <- M[-1]
    if(Aggreg == T){
      M <- aggregate(M, list(M$Species), FUN = mean, na.action = na.pass, na.rm = T)
      names(M)[1] <- "Species"
      M <- M[-2]}
    return(as_tibble(M))
    }
  
  CPT.trait.conti.ACA <- CPT.aggreg(CPT.trait.conti, "SAMPLE.ID", Aggreg = T)
  CPT.trait.unconti.ACA <- CPT.aggreg(CPT.trait.unconti, "SAMPLE.ID", Aggreg = F)
  CPT.PFT.ACA <- CPT.aggreg(CPT.PFT, "SAMPLE.ID", Aggreg = F)
  CPT.pathway.ACA <- CPT.aggreg(CPT.trait.pathway, "SPECIES.ID", Aggreg = F)
  CPT.ACA <- merge(CPT.trait.conti.ACA, CPT.trait.unconti.ACA, by = "Species", all = T)
  CPT.ACA <- merge(CPT.ACA, CPT.PFT.ACA, by = "Species", all = T)
  CPT.ACA <- merge(CPT.ACA, CPT.pathway.ACA, by = "Species", all = T)
  CPT.ACA <- CPT.ACA[!duplicated(CPT.ACA),]
  
  #### Save ####
  saveRDS(TRY.ACA,  paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_taxa_ACA.Rds", sep = ""))
  saveRDS(brot.ACA, "Import/World_DB/Traits/BROT2.0/BROT_ACA.Rds")
  saveRDS(CPT.ACA, "Import/China/Traits/China_Plant_TraitsDB_csv/CPT_ACA.Rds")
  
  }

Trait.DB.pulishing = F
if(Trait.DB.pulishing == T){
  #### Import DB ####
  BIEN.tr.ACA <- readRDS(paste(DB.path, "Vegetation/Occurences/BIEN/BIEN_ACA_trait.Rds", sep = ""))
  TRY.ACA <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_taxa_ACA.Rds", sep = ""))
  brot.ACA <- as_tibble(readRDS("Import/World_DB/Traits/BROT2.0/BROT_ACA.Rds"))
  
  #### TRY clean ####
  TRY.ACA <- TRY.ACA[!(grepl("Guy, A. L., J. M. Mischkolz, and E. G. Lamb.", TRY.ACA$Reference) & TRY.ACA$TraitID %in% 46),]
  TRY.ACA <- TRY.ACA[!(grepl("Wright JP, Sutton", TRY.ACA$Reference) & TRY.ACA$TraitID %in% 46),]
  TRY.ACA[TRY.ACA$TraitID %in% 46 & TRY.ACA$AccSpeciesID %in% 387342,21][3,] <- TRY.ACA[TRY.ACA$TraitID %in% 46 & TRY.ACA$AccSpeciesID %in% 387342,21][3,]/10
  TRY.ACA[TRY.ACA$TraitID %in% 46 & TRY.ACA$AccSpeciesID %in% 31479,21] <- TRY.ACA[TRY.ACA$TraitID %in% 46 & TRY.ACA$AccSpeciesID %in% 31479,21]/10
  TRY.ACA$StdValue[(grepl("Li, Y. and Shipley", TRY.ACA$Reference) & TRY.ACA$TraitID %in% 4)] <- TRY.ACA$StdValue[(grepl("Li, Y. and Shipley", TRY.ACA$Reference) & TRY.ACA$TraitID %in% 4)]*10
  TRY.ACA$StdValue[(TRY.ACA$DatasetID == 415 & TRY.ACA$TraitID %in% 146)] <- TRY.ACA$StdValue[(TRY.ACA$DatasetID == 415 & TRY.ACA$TraitID %in% 146)]/10
  TRY.ACA <- TRY.ACA[setdiff(seq(1,nrow(TRY.ACA)), which(TRY.ACA$TraitID %in% 15 & TRY.ACA$StdValue <= 0.02)),] 
  TRY.ACA <- TRY.ACA[setdiff(seq(1,nrow(TRY.ACA)), which(TRY.ACA$TraitID %in% 14 & TRY.ACA$StdValue == 0)),] 
  TRY.ACA[which(TRY.ACA$TraitID %in% 14 & TRY.ACA$StdValue > 100),"StdValue"] <- TRY.ACA[which(TRY.ACA$TraitID %in% 14 & TRY.ACA$StdValue > 100), "StdValue"]/10
  
  #### TRY references ####
  TRY.ref <- unique(TRY.ACA[c(3,4,26)])
  TRY.ref <- TRY.ref[!duplicated(TRY.ref$DatasetID),]
  TRY.ref <- TRY.ref[!TRY.ref$Reference == "unpub.",]
  TRY.ref <- TRY.ref[order(TRY.ref$DatasetID),]
  
  #### BIEN clean ####
  BIEN.tr.ACA <- BIEN.tr.ACA[!grepl("Albrectsen BR", BIEN.tr.ACA$project_pi),]
  BIEN.tr.ACA$trait_value[grepl("We sampled leaves from 5 to 20 random", BIEN.tr.ACA$method) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("We sampled leaves from 5 to 20 random", BIEN.tr.ACA$method) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)])/10
  BIEN.tr.ACA$trait_value[grepl("Fu B", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("Fu B", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)])*1000
  BIEN.tr.ACA$trait_value[grepl("Stevens JT", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("Stevens JT", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)])*1000
  BIEN.tr.ACA$trait_value[grepl("Ameztegui", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("Ameztegui", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)])*10^6
  BIEN.tr.ACA$trait_value[grepl("Kurokawa", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("Kurokawa", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)])*10^6
  BIEN.tr.ACA$trait_value[grepl("Edwards EJ", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("Edwards EJ", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)])*10^6
  BIEN.tr.ACA$trait_value[grepl("Mason CM", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("Mason CM", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)])*10^6
  BIEN.tr.ACA$trait_value[grepl("Muir CD", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("Muir CD", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)])*10^6
  BIEN.tr.ACA$trait_value[grepl("Price CA", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("Price CA", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)])*10^6
  BIEN.tr.ACA$trait_value[grepl("Storkey J", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("Storkey J", BIEN.tr.ACA$project_pi) & grepl("leaf area per leaf dry", BIEN.tr.ACA$trait_name)])*10^6
  BIEN.tr.ACA$trait_value[grepl("Muir CD", BIEN.tr.ACA$project_pi) & grepl("leaf dry mass per leaf fresh mass", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("Muir CD", BIEN.tr.ACA$project_pi) & grepl("leaf dry mass per leaf fresh mass", BIEN.tr.ACA$trait_name)])*10^2
  BIEN.tr.ACA$trait_value[grepl("Milla R", BIEN.tr.ACA$project_pi) & grepl("leaf dry mass per leaf fresh mass", BIEN.tr.ACA$trait_name)] <- as.numeric(BIEN.tr.ACA$trait_value[grepl("Milla R", BIEN.tr.ACA$project_pi) & grepl("leaf dry mass per leaf fresh mass", BIEN.tr.ACA$trait_name)])*10^2
  
  #### Save DB cleaned ####
  write.table(TRY.ref, "Resultats/ACA/Traits/TRY_ACA_references.csv", sep = ",", row.names = F)
  saveRDS(BIEN.tr.ACA, paste(DB.path, "Vegetation/Occurences/BIEN/BIEN_ACA_trait_clean.Rds", sep = ""))
  saveRDS(TRY.ACA, paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_taxa_ACA_clean.Rds", sep = ""))
  saveRDS(brot.ACA, "Import/World_DB/Traits/BROT2.0/BROT_ACA_clean.Rds")
}

Trait.extract = F
if(Trait.extract == T){
  #### Calculate.ext.trait ####
  Calculate.ext.trait = T
  if(Calculate.ext.trait == T){
    #### Import ####
    BIEN.tr.ACA <- readRDS(paste(DB.path, "Vegetation/Occurences/BIEN/BIEN_ACA_trait_clean.Rds", sep = ""))
    TRY.ACA <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_taxa_ACA_clean.Rds", sep = ""))
    TL.ACA <- readRDS(paste(DB.path, "Vegetation/Occurences/Merge_DB/Taxon_list_ACA_V1.Rds", sep = ""))
    
    brot.ACA <- as_tibble(readRDS("Import/World_DB/Traits/BROT2.0/BROT_ACA_clean.Rds"))
    CPT.ACA <- as_tibble(readRDS("Import/China/Traits/China_Plant_TraitsDB_csv/CPT_ACA.Rds"))
    LTr.TRY.ACA <- unique(TRY.ACA[,c("TraitName", "TraitID", "Comment")])
    LTr.TRY.ACA <- LTr.TRY.ACA[!is.na(LTr.TRY.ACA$TraitID),]
    TaTr.ACA <- as_tibble(TL.ACA[c(-1)])
    
    GIFT.sp.ACA <- readRDS("Import/ACA/Traits/GIFT/20210617_central_asia_species.rds")
    GIFT.tr.ACA <- readRDS("Import/ACA/Traits/GIFT/20210617_central_asia_traits.rds")
    #GIFT.tr.cov.ACA <- readRDS("Import/ACA/Traits/GIFT/20210617_central_asia_traits_coverage.rds")
    Table.trait.GIFT   <- data.frame(read.csv(file="Import/ACA/Traits/GIFT/GIFT_table_trait_ACA.csv",sep=",",dec=".",header=T, row.names = 1, stringsAsFactors = F))
    
    #### Trait names from List.trait.TRY ####
    Trait.labels <- c(
                     "146" = "CNRatio",           # leaf carbon content per leaf nitrogen content (g/g)
                     "47" = "LDMC",               # Leaf dry mass per leaf fresh mass (g.g-1)
                     "3106" = "Height",           # Plant height vegetative
                       "3108" = "LeafArea",  # Leaf area
                     "46" = "LeafThick",          # Leaf thickness 
                     "4" = "SSD",                # Stem specitif density (mg.mm −3)
                     "26" = "SeedMass",           # Seed dry mass 
                     "3115" = "SLA",        # Specific Leaf Area (mm2 mg-1)
                     "59" = "Lifespan",           # Plant lifespan (longevity)
                     # "3108" = "LeafArea_nopet1",  # Leaf area
                     # "3109" = "LeafArea_nopet2",  # Leaf area
                     # "3110" = "LeafArea_pet1",    # Leaf area
                     # "3111" = "LeafArea_pet2",    # Leaf area
                     # "3112" = "LeafArea_undef1",  # Leaf area
                     # "3113" = "LeafArea_undef2",  # Leaf area
                     # "3114" = "LeafArea_undef3",  # Leaf area
                     # "3115" = "SLA_nopet",        # Specific Leaf Area
                     # "3116" = "SLA_pet",          # Specific Leaf Area
                     # "3117" = "SLA_undef",        # Specific Leaf Area
                     #"42" = "PlantForm",          # Plant Growth Form (QUALI) 
                     "37" = "LeafPhenology",      # Leaf phenology type (QUALI)
                     "38" = "Wood",               # Plant Woodiness (QUALI)
                     "22" = "LPP",                 # Leaf photosynthesis pathway (QUALI)
                     "14" = "LeafN",               # Leaf N mass (mg.g-1)
                     "15" = "LeafP"               # Leaf P mass (mg.g-1)
                     #"155" = "Flower_start"    #"plant flowering begin",
                     #"335" = "plant flowering duration" 
    )
    Trait.log.dist <- c("Height", "SeedMass", "LeafArea", "SLA", "Lifespan", "LeafP")
    Trait.to.log.trans <- c("SSD", "LeafN", "LeafP", "SeedMass", "LeafThick", "LDMC", "CNRatio", "Height", "LeafArea", "SLA")
    # Trait.to.log.trans <- c("SSD", "LeafN", "LeafP", "LeafThick", "LDMC", "CNRatio",  "SLA")
    
    saveRDS(Trait.labels, "Resultats/ACA/Traits/Corresp_tax/Trait_ID_names.Rds")
  
    #### BIEN data prep to merge with TRY ####
    Table.trait.TRY.BIEN <- c("4" = "stem wood density",
                              "46" = "leaf thickness",
                              "47" = "leaf dry mass per leaf fresh mass",
                              "48" = "whole plant growth form",
                              "59" = "longest whole plant longevity",
                              "59" = "maximum whole plant longevity",
                              "3115" = "leaf area per leaf dry mass",
                              "26" = "seed mass",
                              "3108" = "leaf area",
                              "146" = "leaf carbon content per leaf nitrogen content" ,
                              "3106" = "maximum whole plant height",
                              "38" = "whole plant woodiness",
                              "37" = "whole plant vegetative phenology",
                              "155" = "plant flowering begin",
                              "335" = "plant flowering duration",
                              "14" = "leaf nitrogen content per leaf dry mass",
                              "15" = "leaf phosphorus content per leaf dry mass"
    )
  
  
    # BIEN convert to TRY datashape
    BIEN.tr.ACA <- as_tibble(BIEN.tr.ACA[which(BIEN.tr.ACA$trait_name %in% data.frame(Table.trait.TRY.BIEN)[[1]]), c(1:3)]) 
    for(i in 1:length(Table.trait.TRY.BIEN)){BIEN.tr.ACA$trait_name[BIEN.tr.ACA$trait_name == Table.trait.TRY.BIEN[[i]]] <- names(Table.trait.TRY.BIEN)[[i]]}
    names(BIEN.tr.ACA) <- c("AccSpeciesName", "TraitID", "StdValue")
    
    #### BROT 2.0 data prep to merge with TRY ####
    Table.trait.TRY.BROT2.0 <- c("4" = "StemDensity",
                              #"46" = "leaf thickness",
                              "47" = "LDMC", #(mg.g-1)
                              "59" = "Lifespan",
                              "3115" = "SLA",
                              "26" = "SeedMass",
                              "3108" = "LeafArea",
                              #"146" = "leaf carbon content per leaf nitrogen content" ,
                              "3106" = "Height",
                              #"38" = "whole plant woodiness",
                              "14" = "LNCm",
                              "37" = "LeafPhenology"#,
                              #"155" = "plant flowering begin",
                              #"335" = "plant flowering duration" 
    )
    
    brot.ACA <- as_tibble(brot.ACA[which(brot.ACA$Trait %in% data.frame(Table.trait.TRY.BROT2.0)[[1]]), c(2:4)]) 
    for(i in 1:length(Table.trait.TRY.BROT2.0)){brot.ACA$Trait[brot.ACA$Trait == Table.trait.TRY.BROT2.0[[i]]] <- names(Table.trait.TRY.BROT2.0)[[i]]}
    names(brot.ACA) <- c("AccSpeciesName", "TraitID", "StdValue")
    
    #### Traits continus selection (TRY + BIEN) ####
    Trait.list.conti.V1 <- c(4, 46, 47, 3115, 3116, 3117, 26, 3108, 3109, 3110, 3111, 3112, 3113, 3114, 146, 3106, 14, 15) #59,
    MC <- Clean.trait.con(TRY.ACA, Trait.ID = Trait.list.conti.V1, Average = T)
    MC.BIEN <- Clean.trait.con(BIEN.tr.ACA, Trait.ID = Trait.list.conti.V1)
    MC.BIEN$`47` <- MC.BIEN$`47`/1000
    MC.BROT <- Clean.trait.con(brot.ACA, Trait.ID = Trait.list.conti.V1)
    MC.BROT$`47` <- MC.BROT$`47`/1000
    
    #### GIFT merge TRY ####
    Trait.GIFT.TRY <- Table.trait.GIFT[which(Table.trait.GIFT$TRY_ID %in% names(MC)),]
    row.names(GIFT.tr.ACA) <- GIFT.tr.ACA$species
    GIFT.tr.TRY <- GIFT.tr.ACA[,c(which(names(GIFT.tr.ACA) %in% Trait.GIFT.TRY$trait))]
    GIFT.tr.TRY <- lapply(GIFT.tr.TRY, function(x) as.double(as.character(x)))   # convert data in numeric
    GIFT.tr.TRY <- as_tibble(c(GIFT.tr.ACA[c(length(GIFT.tr.ACA))], GIFT.tr.TRY))
    names(GIFT.tr.TRY) <- c("species", Trait.GIFT.TRY$TRY_ID)
    GIFT.tr.TRY$`26` <- GIFT.tr.TRY$`26`*100
    GIFT.tr.TRY$`3115` <- GIFT.tr.TRY$`3115`/10
    GIFT.tr.TRY$`4` <- GIFT.tr.TRY$`4`/1000
    GIFT.tr.TRY$`47` <- GIFT.tr.TRY$`47`/1000
    
    
    
    #### All continuous trait merge ####
    MC.ACA <- full_join(TaTr.ACA, MC, by = "species")
    MC.ACA.BIEN <- full_join(TaTr.ACA, MC.BIEN, by = "species")
    MC.ACA.BROT <- full_join(TaTr.ACA, MC.BROT, by = "species")
    MC.ACA.GIFT <- full_join(TaTr.ACA, GIFT.tr.TRY, by = "species")
    MC.ACA.GIFT <- MC.ACA.GIFT[-which(MC.ACA.GIFT$species %in% setdiff(MC.ACA.GIFT$species, TaTr.ACA$species)),]
    
    #### Traits discontinus selection (TRY + BIEN) ####
    MD <- Clean.trait.discon(TRY.ACA, Trait.ID = c(37, 22, 38)) # , 155, 335 /!\ work in progress pour les deux derniers !)
    MT.ACA <- full_join(MC.ACA, MD, by = "species")
     
    names(BIEN.tr.ACA) <- c("AccSpeciesName", "TraitID", "OrigValueStr")
    MD.BIEN <- Clean.trait.discon(BIEN.tr.ACA, Trait.ID = c(37, 22, 38, 48), Keep.string = c(48)) # , 155, 335 /!\ work in progress pour les deux derniers !)
    MT.ACA.BIEN <- full_join(MC.ACA.BIEN, MD.BIEN, by = "species")
    
    names(brot.ACA) <- c("AccSpeciesName", "TraitID", "OrigValueStr")
    MD.BROT <- Clean.trait.discon(brot.ACA, Trait.ID = c(37, 22, 38)) # , 155, 335 /!\ work in progress pour les deux derniers !)
    MT.ACA.BROT <- full_join(MC.ACA.BROT, MD.BROT, by = "species")
    
    #### Traits discontinuous GIFT ####
    Trait.GIFT.TRY <- Table.trait.GIFT[which(Table.trait.GIFT$TRY_ID %in% c(37, 22, 38, 48)),]
    row.names(GIFT.tr.ACA) <- GIFT.tr.ACA$species
    GIFT.tr.TRY <- GIFT.tr.ACA[,c(which(names(GIFT.tr.ACA) %in% Trait.GIFT.TRY$trait))]
    GIFT.tr.TRY <- as_tibble(c(GIFT.tr.ACA[c(length(GIFT.tr.ACA))], GIFT.tr.TRY))
    names(GIFT.tr.TRY) <- c("species", Trait.GIFT.TRY$TRY_ID)
    GIFT.tr.TRY[["38"]][which(GIFT.tr.TRY[["38"]] == "variable")] <- 0.5
    GIFT.tr.TRY[["38"]][which(GIFT.tr.TRY[["38"]] == "non-woody")] <- 0
    GIFT.tr.TRY[["38"]][which(GIFT.tr.TRY[["38"]] == "woody")] <- 1
    GIFT.tr.TRY[["37"]][which(GIFT.tr.TRY[["37"]] == "evergreen")] <- 1
    GIFT.tr.TRY[["37"]][which(GIFT.tr.TRY[["37"]] == "deciduous")] <- 0
    GIFT.tr.TRY[["37"]][which(GIFT.tr.TRY[["37"]] == "variable")] <- 0.5
    GIFT.tr.TRY[["22"]] <- gsub("C3", 0, GIFT.tr.TRY[["22"]])
    GIFT.tr.TRY[["22"]] <- gsub("C4", 1, GIFT.tr.TRY[["22"]])
    GIFT.tr.TRY[["22"]] <- gsub("CAM", -1, GIFT.tr.TRY[["22"]])
    GIFT.tr.TRY[["48"]][which(GIFT.tr.TRY[["48"]] == "subshrub")] <- "Shrub"
    GIFT.tr.TRY[["48"]][which(GIFT.tr.TRY[["48"]] == "shrub")] <- "Shrub"
    GIFT.tr.TRY[["48"]][which(GIFT.tr.TRY[["48"]] == "herb")] <- "Herb"
    GIFT.tr.TRY[["48"]][which(GIFT.tr.TRY[["48"]] == "tree")] <- "Tree"
    GIFT.tr.TRY[["48"]][which(GIFT.tr.TRY[["48"]] == "forb")] <- "Herb"
    GIFT.tr.TRY[["48"]][which(GIFT.tr.TRY[["48"]] == "sedge")] <- "Herb"
    GIFT.tr.TRY[["48"]][which(GIFT.tr.TRY[["48"]] == "graminoid")] <- "Herb"
    GIFT.tr.TRY[["48"]][which(GIFT.tr.TRY[["48"]] == "other")] <- "Other"
    GIFT.tr.TRY[["48"]][which(GIFT.tr.TRY[["48"]] == "parasite")] <- "Other"
    GIFT.tr.TRY[["48"]][which(GIFT.tr.TRY[["48"]] == "climber")] <- "Other"
    GIFT.tr.TRY[["48"]][which(GIFT.tr.TRY[["48"]] == "aquatic")] <- "Other"
    
    stop()
    GrowthForm.ACA <- GIFT.tr.TRY[c("species", "48")]
    names(GrowthForm.ACA)[2] <- "GrowthForm"
    GIFT.tr.TRY <- as_tibble(data.frame(GIFT.tr.TRY[c(1)], as_tibble(lapply(GIFT.tr.TRY[-c(1)], function(x) as.double(as.character(x))))))   # convert data in numeric
    names(GIFT.tr.TRY) <- c("species", Trait.GIFT.TRY$TRY_ID)
    GIFT.tr.TRY <- GIFT.tr.TRY[-c(3)]
    MT.ACA.GIFT <- full_join(MC.ACA.GIFT, GIFT.tr.TRY, by = "species")
    MT.ACA.GIFT <- MT.ACA.GIFT[-which(MT.ACA.GIFT$species %in% setdiff(MT.ACA.GIFT$species, TaTr.ACA$species)),]
    
    GrowthForm.ACA <- full_join(TaTr.ACA, GrowthForm.ACA, by = "species")
    GrowthForm.ACA <- full_join(MT.ACA.BIEN[c(4,ncol(MT.ACA.BIEN))], GrowthForm.ACA, by = "species")
    GrowthForm.ACA$GrowthForm[is.na(GrowthForm.ACA$GrowthForm) & !is.na(GrowthForm.ACA$`48`)] <- GrowthForm.ACA$`48`[is.na(GrowthForm.ACA$GrowthForm) & !is.na(GrowthForm.ACA$`48`)]
    GrowthForm.ACA <- GrowthForm.ACA[-c(2)]
    
    GrowthForm.ACA$GrowthForm[is.na(GrowthForm.ACA$GrowthForm)] <- "Unknown"
    MT.ACA.BIEN <- MT.ACA.BIEN[-ncol(MT.ACA.BIEN)]
    
    #### Merge TRY and BIEN traits ####
    MT.ACA.BROT <- melt(data.frame(MT.ACA.BROT))
    MT.ACA.BROT$DB <- "BROT"
    MT.ACA <- melt(data.frame(MT.ACA))
    MT.ACA$DB <- "TRY"
    MT.ACA.BIEN <- melt(data.frame(MT.ACA.BIEN))
    MT.ACA.BIEN$DB <- "BIEN"
    MT.ACA.GIFT <- melt(data.frame(MT.ACA.GIFT))
    MT.ACA.GIFT$DB <- "GIFT"
    
    A = rbind(MT.ACA.BIEN, MT.ACA, MT.ACA.BROT, MT.ACA.GIFT)
    B = A[c("species", "variable", "value")]
    B <- stats::aggregate(B, list(B$species, B$variable), FUN = mean, na.rm = T, na.action = na.pass)
    
    MT.ACA = reshape2::dcast(B[-c(3,4)], Group.1 ~ Group.2)
    names(MT.ACA)[1] <- "species"
    MT.ACA <- full_join(TaTr.ACA, MT.ACA, by = "species")
    
    # names(MT.ACA.no.scale)[grep("X", names(MT.ACA.no.scale))] <- as.vector(Trait.labels[match(gsub("X", "", names(MT.ACA.no.scale)[grep("X", names(MT.ACA.no.scale))]), names(Trait.labels))])
    MT.ACA$X26[which(MT.ACA$X26 == 0)] <- NA
    names(MT.ACA)[grep("X", names(MT.ACA))] <- as.vector(Trait.labels[match(gsub("X", "", names(MT.ACA)[grep("X", names(MT.ACA))]), names(Trait.labels))])
    
    #### Linear interpolation SSD from LDMC ####
    SSD.interpol = T
    if(SSD.interpol == T){
      TL.APG <- readRDS(paste(DB.path, "Vegetation/Occurences/Merge_DB/Taxon_list_ACA_full_taxo.Rds", sep = ""))
      List.faba <- TL.APG[TL.APG$family == "Fabaceae","species"][[1]]
      List.dicot <- TL.APG[TL.APG$Subreign == "Angiospermes","species"][[1]] 
      List.monoc <- TL.APG[TL.APG$kingdom == "monocots","species"][[1]] 
      List.dicot <- setdiff(List.dicot, c(List.monoc,List.faba))
      
      List.herb <- GrowthForm.ACA[GrowthForm.ACA$GrowthForm == "Herb", "species"][[1]]
      List.faba <- intersect(List.faba, List.herb)
      List.monoc <- intersect(List.monoc, List.herb)
      List.dicot <- intersect(List.dicot, List.herb)
          
      Interp.faba <- 0.692*MT.ACA[is.na(MT.ACA$SSD) & MT.ACA$species %in% List.faba, "LDMC"] + 0.048
      Interp.mono <- 0.888*MT.ACA[is.na(MT.ACA$SSD) & MT.ACA$species %in% List.monoc, "LDMC"] + 0.027
      Interp.dico <- 0.524*MT.ACA[is.na(MT.ACA$SSD) & MT.ACA$species %in% List.dicot, "LDMC"] + 0.096
      
      MT.ACA[is.na(MT.ACA$SSD) & MT.ACA$species %in% List.faba, "SSD"] <- Interp.faba
      MT.ACA[is.na(MT.ACA$SSD) & MT.ACA$species %in% List.monoc, "SSD"] <- Interp.mono
      MT.ACA[is.na(MT.ACA$SSD) & MT.ACA$species %in% List.dicot, "SSD"] <- Interp.dico
        
      MT.hist5 <- setNames(rbind(Interp.faba, Interp.mono, Interp.dico), "value")  
      MT.hist5$PT_sl <- NA 
      MT.hist5$PT_ss <- NA 
      MT.hist5$variable <- "SSD" 
      MT.hist5 <- MT.hist5[c(2,3,4,1)] 
      
      MT.ACA$SLA <- MT.ACA$SLA*0.001
      MT.hist5$Data.Base <- "Imputation of SSD"
      
      
      }
    
    #### Plot Histo Distribution of traits ####
    Histo.trait = F
    if(Histo.trait == T){
      MT.hist1 <- reshape2::melt(MC.ACA[5:ncol(MC.ACA)])
      MT.hist1$Data.Base <- "TRY 4.0"
      MT.hist2 <- reshape2::melt(MC.ACA.BIEN[5:ncol(MC.ACA.BIEN)])
      MT.hist2$Data.Base <- "BIEN 3.4"
      MT.hist3 <- reshape2::melt(MC.ACA.BROT[5:ncol(MC.ACA.BROT)])
      MT.hist3$Data.Base <- "BROT 2.0"
      MT.hist4 <- reshape2::melt(MC.ACA.GIFT[5:ncol(MC.ACA.GIFT)])
      MT.hist4$Data.Base <- "GIFT"
      
      MT.hist <- rbind(MT.hist1, MT.hist2, MT.hist3, MT.hist4) 
      MT.hist$variable <- as.vector(Trait.labels[match(MT.hist$variable, names(Trait.labels))])
      if(SSD.interpol == T){MT.hist <- rbind(MT.hist, MT.hist5)}
      
      Keep.trait <- c("SSD", "LeafN", "SeedMass", "Height", "LeafArea", "SLA")
      MT.hist <- MT.hist[MT.hist$variable %in% Keep.trait,]
      
      MT.hist.log <- MT.hist[MT.hist$variable %in% Trait.log.dist,]
      MT.hist <- MT.hist[!MT.hist$variable %in% Trait.log.dist,]
      p1 <- ggplot(data = MT.hist.log, aes(x = value)) + 
        geom_histogram(aes(fill = Data.Base), bins = 30, na.rm = T, alpha = 1, position="stack") + # dodge, identity
        # facet_wrap(vars(variable), scales = "free", ncol = 3)+
        facet_wrap(vars(variable), scales = "free", nrow = 2)+
        scale_x_continuous(trans = 'log10',
                           breaks = trans_breaks('log10', function(x) 10^x),
                           labels = trans_format('log10', math_format(10^.x)))+
        theme_par()+
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))
      
      p2 <- ggplot(data = MT.hist, aes(x = value)) + 
        geom_histogram(aes(fill = Data.Base), bins = 30, na.rm = T, alpha = 1, position="stack") +
        facet_wrap(vars(variable), scales = "free", nrow = 2)+
        theme_par()+
        ylab("Number of observations")+
        theme(axis.title.x = element_blank(), plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))
      
      # pf <- p2 + p1 + plot_layout(widths = c(2/5, 3/5), guides = "collect")
      pf <- p2 + p1 + plot_layout(widths = c(1/3, 2/3), guides = "collect")
      W = 1700
      H = 1000
      ggsave(filename = "Figures/ACA/Trait/Distribution/Trait_dist_hist_V3.pdf", pf, width = W*0.026458333, height = H*0.026458333, units = "cm")
      }
    
    #### Stats sur les traits ####
    MT.ACA.no.scale <- MT.ACA 
    
    Table.trait.stats <- data.frame(read.csv("Import/ACA/Traits/Table_traits.csv", sep = "\t", header = T))
    Stats.traits <- data.frame(lapply(MT.ACA.no.scale, function(x) signif(mean(x, na.rm = T), digits = 3)))
    Stats.traits <- rbind(Stats.traits, data.frame(lapply(MT.ACA.no.scale, function(x) signif(sd(x, na.rm = T), digits = 3))))
    Stats.traits <- rbind(Stats.traits, data.frame(lapply(MT.ACA.no.scale, function(x) min(x, na.rm = T))))
    Stats.traits <- rbind(Stats.traits, data.frame(lapply(MT.ACA.no.scale, function(x) max(x, na.rm = T))))
    Stats.traits <- rbind(Stats.traits, data.frame(lapply(MT.ACA.no.scale, function(x) round(length(which(is.na(x)==F)), digits = 0))))
    row.names(Stats.traits) <- c("Mean", "SD", "Min", "Max", "mathbb{N_{species}}")
    # Stats.traits <- t(Stats.traits[c(14,7,15,16,8,10)])
    Stats.traits <- t(Stats.traits[c(12,13,6,8,14,5)])
    row.names(Table.trait.stats) <- row.names(Stats.traits)
    Table.trait.stats <- cbind(Table.trait.stats, Stats.traits)
    Table.trait.stats$Min <- signif(Table.trait.stats$Min, digits = 3)
    Table.trait.stats$Max <- signif(Table.trait.stats$Max, digits = 3)
    Table.trait.stats <- Table.trait.stats[c(1,3,5,6,4,2),]
    
    print(Table.trait.stats)
    write.csv(Table.trait.stats, "Resultats/ACA/Traits/Corresp_tax/Table_traits_stats.csv", row.names = F)
    
    #### Scaling (z-scores) ####
    print("For the RESCALE hip hip hiiiip")
    
    # MT.ACA <- MT.ACA[-c(5,6,17,19)]
    # MT.ACA <- MT.ACA[-c(17,19)]
    MT.ACA <- MT.ACA[-c(15,17)]
    MT.ACA[,names(MT.ACA) %in% Trait.to.log.trans] <- log10(MT.ACA[,names(MT.ACA) %in% Trait.to.log.trans])
    MT.ACA.scaled <- MT.ACA
    # MT.ACA.scaled[7:ncol(MT.ACA.scaled)] <- scale(MT.ACA.scaled[7:ncol(MT.ACA.scaled)])
    MT.ACA.scaled[5:ncol(MT.ACA.scaled)] <- scale(MT.ACA.scaled[5:ncol(MT.ACA.scaled)])
    #### PLOT check z-score distrib ####
    Check.distri = T
    if(Check.distri == T){
      print("Plot gaussian scale distribution.")
      MT.ACA.hist <- reshape2::melt(MT.ACA[5:ncol(MT.ACA)-1])
      hist.ACA <- ggplot(data = MT.ACA.hist, aes(x = value)) + 
        geom_histogram(bins = 30, na.rm = T, alpha = 0.3, position="identity") +
        facet_wrap(vars(variable), scales = "free", nrow = 2)+
        theme_par()+
        ylab("Number of observations")+
        theme(axis.title.x = element_blank(), plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))
      
      MT.ACA.hist.scaled <- reshape2::melt(MT.ACA.scaled[5:ncol(MT.ACA.scaled)-1])
      hist.ACA.scaled <- ggplot(data = MT.ACA.hist.scaled, aes(x = value)) + 
        geom_histogram(bins = 30, na.rm = T, alpha = 0.3, position="identity") +
        facet_wrap(vars(variable), scales = "free", nrow = 2)+
        theme_par()+
        ylab("Number of observations")+
        theme(axis.title.x = element_blank(), plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))
      
      W = 2100
      H = 800
      ggsave(filename = "Figures/ACA/Trait/Distribution/Trait_dist_hist_ACA_V2.pdf", hist.ACA, width = W*0.026458333, height = H*0.026458333, units = "cm")
      ggsave(filename = "Figures/ACA/Trait/Distribution/Trait_dist_hist_ACA_scaled_V2.pdf", hist.ACA.scaled, width = W*0.026458333, height = H*0.026458333, units = "cm")
    }
    
    MT.ACA <- MT.ACA.scaled
    MT.ACA <- MT.ACA[!duplicated(MT.ACA),]
    
    saveRDS(MT.ACA.no.scale, paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_trait_ACA_no_scale.Rds", sep = ""))
    saveRDS(MT.ACA, paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_trait_ACA.Rds", sep = ""))
    saveRDS(GrowthForm.ACA, paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_GrowthForm_ACA.Rds", sep = ""))
    
    #### Gap-filled the matrix ####
    Gap.filled = F
    if(Gap.filled == T){
      #### Import matrix ####
      print("Gap-filling beginning ! Yipaaa")
      library("BHPMF") # en cas de probleme réinstaller le package library(devtools) install_github("fisw10/BHPMF")
      MT.ACA <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_trait_ACA.Rds", sep = ""))
      TL.APG <- readRDS(paste(DB.path, "Vegetation/Occurences/Merge_DB/Taxon_list_ACA_full_taxo.Rds", sep = ""))
      
      MT.ACA <- left_join(MT.ACA, TL.APG, by = c("family","genus", "species", "Other.clade"))
      MT.ACA <- MT.ACA[setdiff(1:ncol(MT.ACA), which(names(MT.ACA) %in% c("LeafPhenology", "Wood", "LPP", "PT_sl", "PT_ss")))]
      MT.ACA <- MT.ACA[c(15,16,17,1:14)]
      
      MT.ACA[MT.ACA == "NaN"] <- NA
      MT.ACA <- MT.ACA[which(!is.nan(rowMeans(MT.ACA[8:ncol(MT.ACA)], na.rm = T))),]
      
      #### Cut hierarchie / trait DB ####
      # Only.trait <- as.matrix(MT.ACA[,7:ncol(MT.ACA)])
      Only.trait <- as.matrix(MT.ACA[,8:ncol(MT.ACA)])
      Only.hiera <- MT.ACA[,c(7,6,4,3,2,1)]
      Only.hiera <- data.frame(lapply(Only.hiera, as.factor))
      Only.hiera <- cbind(plant_id = row.names(MT.ACA), Only.hiera)
      
      new.folder <- "/home/lucas.dugerdil/Documents/Recherche/R_stats/Resultats/World_DB/Traits/"
      saveRDS(Only.hiera, paste(new.folder, "Hierarchie_gap-filling.Rds", sep = ""))
      saveRDS(Only.trait, paste(new.folder, "TM_before_gap-filling.Rds", sep = ""))
      
      #### Gapfilling model ####
      GapFilling(Only.trait, Only.hiera, verbose = T,
                 # prediction.level = 4,
                 # used.num.hierarchy.levels = 2,
                 mean.gap.filled.output.path = "/tmp/MTrait_ACA_gf_V1.txt",
                 std.gap.filled.output.path = "/tmp/MTrait_ACA_gf_sd_V1.txt")

      list.of.files <- list.files("/tmp", "*txt", full.names = T)
      file.copy(list.of.files, new.folder, overwrite = T)
      
      MT.ACA.sd <- read.table(file = "Resultats/World_DB/Traits/MTrait_ACA_gf_sd_V1.txt", sep="\t",dec=".", header=T, stringsAsFactors = F)
      MT.ACA.gf <- read.table(file = "Resultats/World_DB/Traits/MTrait_ACA_gf_V1.txt", sep = "\t", dec=".", header=T)
      MT.ACA <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_trait_ACA.Rds", sep = ""))
      
      # Calculate average RMSE with the default values
      # out1 <- CalculateCvRmse(Only.trait, Only.hiera)
      # avg.rmse <- out1$avg.rmse
      # std.rmse <- out1$std.rmse
      }
    else{
      MT.ACA <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_trait_ACA.Rds", sep = ""))
      new.folder <- "/home/lucas.dugerdil/Documents/Recherche/R_stats/Resultats/World_DB/Traits/"
      Only.hiera <- readRDS(paste(new.folder, "Hierarchie_gap-filling.Rds", sep = ""))
      Only.trait <- readRDS(paste(new.folder, "TM_before_gap-filling.Rds", sep = ""))
      MT.ACA.sd <- read.table(file = "Resultats/World_DB/Traits/MTrait_ACA_gf_sd_V1.txt", sep="\t",dec=".", header=T, stringsAsFactors = F)
      MT.ACA.gf <- read.table(file = "Resultats/World_DB/Traits/MTrait_ACA_gf_V1.txt", sep = "\t", dec=".", header=T)
      }
    
    #### Gap-filling vegetation trait matrice ####
    MT.ACA.cv <- MT.ACA.sd/MT.ACA.gf
    MT.ACA.gf[MT.ACA.cv > 1 & is.na(Only.trait)] <- NA
    MT.ACA.gf <- as_tibble(cbind(Only.hiera[c(2)], MT.ACA.gf))
    # Keep.trait.conti <- MT.ACA[c("family","genus", "species", "PT_sl", "PT_ss")]
    Keep.trait.conti <- MT.ACA[c("family","genus", "species")]
    
    MT.ACA.gf <- full_join(MT.ACA.gf, MT.ACA, by = intersect(names(MT.ACA.gf), names(MT.ACA)))
    MT.ACA.gf <- MT.ACA.gf[match(names(MT.ACA), names(MT.ACA.gf))]
    MT.ACA.gf <- aggregate(MT.ACA.gf, list(MT.ACA.gf$species), FUN = mean, na.action = na.pass, na.rm = T)
    names(MT.ACA.gf)[1] <- "species"
    # MT.ACA.gf <- MT.ACA.gf[-c(2,3,4,5)]
    MT.ACA.gf <- MT.ACA.gf[-c(2:5)]
    MT.ACA.gf <- full_join(MT.ACA.gf, Keep.trait.conti, by = "species")
    MT.ACA.gf <- MT.ACA.gf[match(names(MT.ACA), names(MT.ACA.gf))[!is.na(match(names(MT.ACA), names(MT.ACA.gf)))]]
    saveRDS(MT.ACA.gf, "Resultats/World_DB/Traits/MT_ACA_gapfilled.Rds")
    
    #### Aggregate trait family / genus ####
    print("Let's aggregate pollen-type, bro !")
    Trait.aggregate.by.class <- function(M){
      M <- M[which((!is.na(M[[1]]))),]
      # return(M)
      M <- aggregate(M, list(M[[1]]), FUN = mean, na.action = na.pass, na.rm = T)
      names(M)[1] <- names(M)[2]
      M <- M[-2]
      return(M)
      }
    
    Tr.PT_fam <- Trait.aggregate.by.class(MT.ACA[-c(2,3,4)]) # Familly merge
    Tr.PT_gen <- Trait.aggregate.by.class(MT.ACA[-c(1,2,4)]) # Genus merge
    
    Tr.PT_fam_gf <- Trait.aggregate.by.class(MT.ACA.gf[-c(2,3)])
    Tr.PT_gen_gf <- Trait.aggregate.by.class(MT.ACA.gf[-c(1,3)])
    
    #### Aggregation pollen-type ####
    MT.bool.ss <- readRDS(file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_plant_type_bool_ss.Rds")
    MT.bool.sl <- readRDS(file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_plant_type_bool_sl.Rds")
    
    names(MT.bool.sl)[names(MT.bool.sl) == "Capparidaceae"] <- "Capparaceae"
    MT.bool.sl[MT.bool.sl$family == "Capparaceae", names(MT.bool.sl) == "Capparaceae"] <- T
    MT.bool.sl[MT.bool.sl$family == "Xanthorrhoeaceae", names(MT.bool.sl) == "Xanthorrhoeaceae"] <- T
    # MT.bool.sl <- MT.bool.sl[!names(MT.bool.sl) == "Hololachna"]
    MT.bool.ss[MT.bool.ss$family == "Capparaceae", names(MT.bool.ss) == "Capparaceae"] <- T
    names(MT.bool.ss)[names(MT.bool.ss) == "Capparidaceae"] <- "Capparaceae"
    # MT.bool.ss[MT.bool.ss$genus == "Dianthus petrorhagia", names(MT.bool.ss) == "Dianthus"] <- T
    # names(MT.bool.ss)[names(MT.bool.ss) == "Dianthus petrorhagia"] <- "Dianthus"
    
    MT.ACA <- MT.ACA[-c(2)]
    
    Trait.aggregate.by.type <- function(TP, MT, name.var){
      Fuck.off <- c("species", "family", "genus", "PT_ss", "PT_sl", "Subreign", "kingdom", "order", "Other.clade")
      Col.to.keep <- setdiff(names(MT), Fuck.off)
      Row.to.keep <- setdiff(names(TP), Fuck.off)
      
      TP.work <- setNames(data.frame(matrix(NA, ncol = length(Col.to.keep)+1, nrow = length(Row.to.keep))), c(name.var,Col.to.keep))
      TP.work[[name.var]] <- Row.to.keep
      names(TP.work)[1] <- "species"
      
      T.m <- TP.work
      T.sd <- TP.work
      TP <- TP[c("species",Row.to.keep)]
      MT <- data.frame(MT[c("species",Col.to.keep)])
      pb = txtProgressBar(min = 0, max = ncol(TP), initial = 0) 
      print(paste("Aggregation with", name.var))
      
      for(i in 2:ncol(TP)){
        setTxtProgressBar(pb,i)
        S.to.pick <- TP$species[which(TP[i] == T)]
        for(j in 2:ncol(MT)){
          Val <- MT[which(MT$species %in% S.to.pick), j]
          Trait.mean <- mean(Val, na.rm = T)
          Trait.sd <- sd(Val, na.rm = T)
          T.m[i-1,j] <- Trait.mean
          T.sd[i-1,j] <- Trait.sd
        }
      }
      close(pb)
      # return(list(Mean = T.m, SD = T.sd, Nval = T.n))
      return(list(Mean = T.m, SD = T.sd))
    }
    
    
    Tr.PT_sl.full <- Trait.aggregate.by.type(MT.bool.sl, MT.ACA, "PT_sl") # Normal sl
    Tr.PT_ss.full <- Trait.aggregate.by.type(MT.bool.ss, MT.ACA, "PT_ss") # Normal ss
    Tr.PT_sl.g.full <- Trait.aggregate.by.type(MT.bool.sl, MT.ACA.gf, "PT_sl") # Normal sl
    Tr.PT_ss.g.full <- Trait.aggregate.by.type(MT.bool.ss, MT.ACA.gf, "PT_ss") # Normal sl
    
    Tr.PT_sl.sd <- Tr.PT_sl.full$SD
    Tr.PT_ss.sd <- Tr.PT_ss.full$SD
    Tr.PT_sl_gf.sd <- Tr.PT_sl.g.full$SD
    Tr.PT_ss_gf.sd <- Tr.PT_ss.g.full$SD
    
    Tr.PT_sl <- Tr.PT_sl.full$Mean
    Tr.PT_ss <- Tr.PT_ss.full$Mean
    Tr.PT_sl_gf <- Tr.PT_sl.g.full$Mean
    Tr.PT_ss_gf <- Tr.PT_ss.g.full$Mean
    
  
    #### PLOT SD by traits ####
    MSD <- cbind(melt(Tr.PT_sl.sd, id = "species"), type = "sl")
    MSD <- rbind(MSD, cbind(melt(Tr.PT_sl_gf.sd, id = "species"), type = "sl_gf"))
    MSD2 <- cbind(melt(Tr.PT_ss.sd, id = "species"), type = "ss")
    MSD2 <- rbind(MSD2, cbind(melt(Tr.PT_ss_gf.sd, id = "species"), type = "ss_gf"))
    names(MSD)[1] <- "species"
    names(MSD2)[1] <- "species"
    MSD <- rbind(MSD, MSD2)
    Trait.selected <- c("Height", "LeafArea", "LeafN", "SeedMass", "SLA", "SSD")
    MSD <- MSD[MSD$variable %in% Trait.selected,]
    
    # levels(MSD$variable) <- Trait.selected
    MSD$variable <- factor(MSD$variable, levels = Trait.selected, ordered = T)
    
    MSD <- MSD[MSD$type %in% c("ss_gf", "sl_gf"),]
    # MSD <- MSD[MSD$type %in% c("sl_gf"),]
    
    Main.taxa.4 <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/4Main_taxa_PCA.Rds")
    Main.taxa.22 <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/22Main_taxa_PCA.Rds")
    MSD.maj4 <- MSD[MSD$species %in% Main.taxa.4,]
    MSD.maj4$type <- paste(MSD.maj4$type, "maj4", sep = "_")
    MSD.maj22 <- MSD[MSD$species %in% Main.taxa.22,]
    MSD.maj22$type <- paste(MSD.maj22$type, "maj22", sep = "_")
    # MSD <- rbind(MSD, MSD.maj22)
    MSD2 <- rbind(MSD, MSD.maj22, MSD.maj4)
    
    Col.vec <- c("sl_gf" = "grey10",
                 "sl_gf_maj22" = "grey30",
                 "sl_gf_maj4" = "grey60", "ss_gf" = "#335286ff",
                 "ss_gf_maj22" = "#3d77a1ff",
                 "ss_gf_maj4" = "#5fbfc2ff")
    
    p1 <- ggplot(MSD,  aes(x = variable, y = value, fill = type))+
      ylim(0,1.3)+ xlab("Traits")+  ylab("z-score SD")+
      scale_fill_manual(values = Col.vec, name = "ACASP aggregation scheme", drop = T)+
      geom_boxplot(outlier.colour = "red", outlier.shape = NA, alpha = 0.7, notch = F, notchwidth = 0.7, varwidth = F, na.rm = T, show.legend = T)+
      theme(legend.position = c(0.4,0.8), panel.background = element_blank(),
            legend.key = element_blank(),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 10),
            plot.margin = unit(x = c(1, 2, 2, 0),units="mm"), # Bas / Droite / Haut / Gauche
            axis.line = element_line(colour = "black"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 10),
            plot.background = element_blank())
      
    p2 <- ggplot(MSD2,  aes(x = variable, y = value, fill = type))+
      ylim(0,1.3)+ ylab(NULL)+ xlab("Traits")+
      # geom_boxplot()
      scale_fill_manual(values = Col.vec, drop = T)+
      geom_boxplot(outlier.colour = "red", outlier.shape = NA, alpha = 0.7, notch = F, notchwidth = 0.7, varwidth = F, na.rm = T, show.legend = T)+
      theme(legend.position = "none", panel.background = element_blank(),
            legend.key = element_blank(),
            axis.title = element_text(size = 14),
            axis.text.y = element_blank(),
            axis.text = element_text(size = 10),
            plot.margin = unit(x = c(1, 2, 2, 0),units="mm"), # Bas / Droite / Haut / Gauche
            axis.line = element_line(colour = "black"),
            plot.background = element_blank())
      
      
    
    p <- (p1 + p2) + plot_layout(widths = c(1/3, 2/3))
    W = 1100
    H = 500
    ggsave(filename = "Figures/ACA/Trait/Dispersion/SD_trait_main.pdf", p, width = W*0.026458333, height = H*0.026458333, units = "cm")

    
    
    #### Fullfill the fam / gen values in PT ####
    # Fullfill.fam.gen.PT <- function(M.to.fill, M.to.add){
    #   for(i in 1:ncol(M.to.fill)){
    #     for(j in 1:nrow(M.to.fill)){
    #       if(is.na(M.to.fill[j,i])){
    #         Taxa <- M.to.fill[j,1]
    #         
    #         if(grepl("eae", Taxa) == T){
    #           # print(M.to.fill[j,i])
    #           Trait <- names(M.to.fill)[i]
    #           print(Trait)
    #           print(Taxa)
    #         # print(M.to.fill[j,i])
    #         # New.val <- M.to.add[M.to.add$family == Taxa, i]
    #         # print(New.val)
    #         # if(!is.na(New.val)){M.to.fill[j,i] <- New.val}
    #         }
    #       }
    #     }
    #   }
    #   return(M.to.fill)}
    # 
    # 
    # Tr.PT_ss.c <- Fullfill.fam.gen.PT(Tr.PT_ss, Tr.PT_fam)
    # Tr.PT_ss_gf <- Fullfill.fam.gen.PT(Tr.PT_ss_gf, Tr.PT_fam)
    # Tr.PT_ss_gf <- Fullfill.fam.gen.PT(Tr.PT_ss_gf, Tr.PT_fam_gf)
    # Tr.PT_sl <- Fullfill.fam.gen.PT(Tr.PT_sl, Tr.PT_fam)
    # Tr.PT_sl_gf <- Fullfill.fam.gen.PT(Tr.PT_sl_gf, Tr.PT_fam)
    # Tr.PT_sl_gf <- Fullfill.fam.gen.PT(Tr.PT_sl_gf, Tr.PT_fam_gf)
    # Tr.PT_ss <- Fullfill.fam.gen.PT(Tr.PT_ss, Tr.PT_gen)
    # Tr.PT_ss_gf <- Fullfill.fam.gen.PT(Tr.PT_ss_gf, Tr.PT_gen)
    # Tr.PT_ss_gf <- Fullfill.fam.gen.PT(Tr.PT_ss_gf, Tr.PT_gen_gf)
    # Tr.PT_sl <- Fullfill.fam.gen.PT(Tr.PT_sl, Tr.PT_gen)
    # Tr.PT_sl_gf <- Fullfill.fam.gen.PT(Tr.PT_sl_gf, Tr.PT_gen)
    # Tr.PT_sl_gf <- Fullfill.fam.gen.PT(Tr.PT_sl_gf, Tr.PT_gen_gf)
    
    #### Replace par value non gf quand c'est possible (qui est normalement plus fiable que la gf) ####
    Keep.no.gf <- function(M.to.fill, M.to.add, id){
      Miss.type <- which(M.to.fill[[id]] %in% M.to.add[[1]])
      To.add <- which(M.to.add[[1]] %in% M.to.fill[[id]] )
      for(i in 2:ncol(M.to.fill)){
        for(j in 1:length(Miss.type)){
          if(is.na(M.to.add[To.add[j],i]) == F){
            M.to.fill[Miss.type[j], i] <- M.to.add[To.add[j],i]
          }
        }
      }
      return(M.to.fill)}
    
    Tr.PT_ss_gf <- Keep.no.gf(Tr.PT_ss_gf, Tr.PT_fam, 1)
    Tr.PT_ss_gf <- Keep.no.gf(Tr.PT_ss_gf, Tr.PT_gen, 1)
    Tr.PT_sl_gf <- Keep.no.gf(Tr.PT_sl_gf, Tr.PT_fam, 1)
    Tr.PT_sl_gf <- Keep.no.gf(Tr.PT_sl_gf, Tr.PT_gen, 1)
    Tr.PT_fam_gf <- Keep.no.gf(Tr.PT_fam_gf, Tr.PT_fam, 1)
    Tr.PT_gen_gf <- Keep.no.gf(Tr.PT_gen_gf, Tr.PT_gen, 1)
    
    #### Extract / aggregate trait pour MV_TL VEGETATION PLOT #### 
    MV_TL <- readRDS("Resultats/ACA/Vegetation/MV_ACA_taxonlist.Rds")
    MV_TL[which(MV_TL$Accepted_family == "Fabaceae"), "Accepted_family"] <- "Fabaceae"
    MV.releve <- readRDS("Resultats/ACA/Vegetation/MV_ACA_plot_only_spermatophyte.Rds")
    names(MV.releve)[which(names(MV.releve) == "Leguminosae")] <- "Fabaceae"
    MT.MV.gf <- MT.ACA.gf[which(MT.ACA.gf$species %in% MV_TL$Accepted_name),]
    MV_TL.full <- names(MV.releve)
    
    MT.MV.gen.gf <- Tr.PT_gen_gf[which(Tr.PT_gen_gf$genus %in% MV_TL.full),]
    MT.MV.fam.gf <- Tr.PT_fam_gf[which(Tr.PT_fam_gf$family %in% MV_TL.full),]
    names(MT.MV.gen.gf)[1] <- "species"
    names(MT.MV.fam.gf)[1] <- "species"
    Tr.MV_gf <- full_join(MT.MV.gen.gf, MT.MV.gf, by = intersect(names(MT.MV.gf), names(MT.MV.gen.gf)))
    Tr.MV_gf <- full_join(MT.MV.fam.gf, Tr.MV_gf, by = intersect(names(Tr.MV_gf), names(MT.MV.fam.gf)))
    Tr.MV_gf <- Tr.MV_gf[match(names(MT.ACA.gf), names(Tr.MV_gf))]
    
    MT.MV <- MT.ACA[which(MT.ACA$species %in% MV_TL$Accepted_name),]
    
    MT.MV.gen <- Tr.PT_gen[which(Tr.PT_gen$genus %in% MV_TL$genus),]
    MT.MV.fam <- Tr.PT_fam[which(Tr.PT_fam$family %in% MV_TL$Accepted_name),]
    names(MT.MV.gen)[1] <- "species"
    names(MT.MV.fam)[1] <- "species"
    Tr.MV <- full_join(MT.MV.gen, MT.MV, by = intersect(names(MT.MV), names(MT.MV.gen)))
    Tr.MV <- full_join(MT.MV.fam, Tr.MV, by = intersect(names(Tr.MV), names(MT.MV.fam)))
    Tr.MV <- Tr.MV[match(names(MT.ACA), names(Tr.MV))]
  
    Tr.MV_gf <- Tr.MV_gf[-c(1,2)]
    Tr.MV <- Tr.MV[-c(1,2)]
    
    #### Check if missing taxa in front of splot and surf pol DB ####
    if(exists("MV.releve") == F){MV.releve <- readRDS("Resultats/ACA/Vegetation/MV_ACA_plot_only_spermatophyte.Rds")}
    MP_sl <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_sl_clean.Rds")
    MP_ss <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_ss_clean.Rds")
    
    
    #### Check problems ####
    Check.missing.trait.1 <- function(MP, MT, type){
      if(type == "non-gf"){
        Gen <- Tr.PT_gen
        Fam <- Tr.PT_fam
      }
      if(type == "gf"){
        Gen <- Tr.PT_gen_gf
        Fam <- Tr.PT_fam_gf
      }
      Missing.trait <- setdiff(names(MP),MT[[1]])
      # print(Missing.trait)
      Add.fam <- Fam[which(Fam$family %in% Missing.trait),]
      names(Add.fam)[1] <- names(MT)[1]
      MT <- rbind(MT, Add.fam)

      Add.gen <- Gen[which(Gen$genus %in% Missing.trait),]
      names(Add.gen)[1] <- names(MT)[1]
      MT <- rbind(MT, Add.gen)
      Missing.trait <- setdiff(names(MP),MT[[1]])
      Found <- c(Add.fam[[1]], Add.gen[[1]])
      print("*** Found taxa ***")
      print(Found)
      print("*** Still missing ***")
      print(Missing.trait)

      return(MT)
      }

    Check.missing.trait <- function(MP, MT, TNRS.check, return.missing){
      #### Check missing ####
      if(missing(TNRS.check)){TNRS.check = F}
      if(missing(return.missing)){return.missing = F}
      APG_diff <- data.frame(read.csv("Import/World_DB/Taxonomie/APG_II_APG_IV_diff.csv" ,sep=",",dec=".",header=T, stringsAsFactors = F))
      
      Missing.trait <- setdiff(names(MP), MT[[1]])
      
      print("*** Previously missing ***")
      # print(length(Missing.trait))
      print(Missing.trait)
      
      #### Transpose ####
      t_df <- as.data.frame(t(MP[-1]),row.names = NA)
      t_df <- cbind(names(MP)[-1],t_df)
      names(t_df)[1] <- "ID"
      
      #### TNRS missing ####
      if(TNRS.check == T){
        # Missing.trait.tnrs <- TNRS.taxa.accepted(Taxa.vector = Missing.trait)
        library("TNRS")
        Missing.trait.tnrs <- TNRS(Missing.trait)
        # print(Missing.trait.tnrs)
        Missing.trait.tnrs <- Missing.trait.tnrs[Missing.trait.tnrs$match.score >= 0.80,]
        t_df$ID[which(t_df$ID %in% Missing.trait.tnrs$Name_submitted)] <- Missing.trait.tnrs$Accepted_name[which(Missing.trait.tnrs$Name_submitted %in% t_df$ID)]
      }
      
      #### Check diff APG II et APG IV    
      t_df$ID[which(t_df$ID %in% APG_diff$Old_name)] <- APG_diff$New_name[which(APG_diff$Old_name %in% t_df$ID)]
      
      #### Aggregate ####
      new_df=aggregate(t_df[-1], by=list(t_df[["ID"]]), sum)
      t_df2 <- as.data.frame(t(new_df[-1]),row.names = NA)
      names(t_df2) <- new_df$Group.1
      row.names(t_df2) <- names(new_df)[-1]
      MP <- t_df2
      
      Missing.trait <- setdiff(names(MP),MT[[1]])
      
      print("*** Still missing ***")
      # print(length(Missing.trait))
      print(Missing.trait)
      
      if(return.missing == T){return(Missing.trait)}
      else{return(MP)}
    }
    
    MP_sl <- Check.missing.trait(MP_sl, Tr.PT_sl, TNRS.check = T)
    MP_ss <- Check.missing.trait(MP_ss, Tr.PT_ss, TNRS.check = T)
    # Tr.PT_sl_gf <- Check.missing.trait(MP_sl, Tr.PT_sl_gf, "gf")
    # Tr.PT_ss_gf <- Check.missing.trait(MP_ss, Tr.PT_ss_gf, "gf")
    
    #### Test ####  
    Trait.to.test <- c("SSD", "SLA", "Height")
    # Taxa.to.test <- "Carpinus betulus" # "Astragalus", "Fabaceae"
    # print(paste("*** TEST the ", Trait.to.test, " on the ", Taxa.to.test, " taxa. ***"))
    # print(Tr.PT_fam[Tr.PT_fam$family == Taxa.to.test,Trait.to.test])
    # print(Tr.PT_fam_gf[Tr.PT_fam_gf$family == Taxa.to.test,Trait.to.test])
    # print(Tr.PT_ss[Tr.PT_ss$PT_ss == Taxa.to.test,Trait.to.test])
    # print(Tr.PT_sl[Tr.PT_sl$PT_sl == Taxa.to.test,Trait.to.test])
    # print(Tr.MV[Tr.MV$species == Taxa.to.test,Trait.to.test])
    # print(Tr.PT_ss_gf[Tr.PT_ss_gf$PT_ss == Taxa.to.test,Trait.to.test])
    # print(Tr.PT_sl_gf[Tr.PT_sl_gf$PT_sl == Taxa.to.test,Trait.to.test])
    # print(Tr.MV_gf[Tr.MV_gf$species == Taxa.to.test,Trait.to.test])
    
    #### PLOT Hist traits variance by pollen-type ####
    Hist.trait.by.spe = T
    if(Hist.trait.by.spe == T){
      Trait.distribution(MT = MT.ACA.no.scale, MP1 = MT.bool.ss, Show.var = T, Diaz.compa = T,
                         Selec.tax.to.shown = c("Amaranthaceae", "Cyperaceae", "Ericaceae", "Fabaceae", "Poaceae", "Artemisia", "Nitraria", "Ephedra", "Pinus diploxylon", "Pinus haploxylon", "Betula pubescens", "Carpinus betulus"),
                         Select.trait = c("SSD", "LeafN", "SeedMass", "Height", "LeafArea", "SLA"),
                         Plot.gaussian = F, Log.scale = T, Boundaries = c(6.5,9.5,11.5), Scale.var = F,
                         W = 1200, H = 750, Save.plot = "Figures/ACA/Trait/Dispersion/Hist_disp_tot_no_scale.pdf")
      
      # stop("here")
      # 
      # Trait.distribution(MT = MT.ACA.no.scale, MP1 = MT.bool.ss, MP2 = MT.bool.sl, Show.var = T,
      #                    # Selec.tax.to.shown = c("Amaranthaceae", "Ericaceae", "Fabaceae", "Poaceae", "Artemisia", "Ephedra", "Pinus diploxylon", "Pinus haploxylon", "Betula pubescens", "Carpinus betulus"),
      #                    # Selec.tax.to.shown = c("Amaranthaceae", "Ericaceae", "Fabaceae", "Cistaceae", "Poaceae", "Artemisia", "Nitraria", "Ephedra", "Pinus diploxylon", "Pinus haploxylon", "Betula pubescens", "Carpinus betulus"),
      #                    Selec.tax.to.shown = c("Amaranthaceae", "Cyperaceae", "Ericaceae", "Fabaceae", "Nitrariaceae", "Poaceae", "Artemisia", "Nitraria", "Ephedra", "Pinus diploxylon", "Pinus haploxylon", "Betula pubescens", "Carpinus betulus"),
      #                    Select.trait = c("SSD", "LeafN", "SeedMass", "Height", "LeafArea", "SLA"),
      #                    Plot.gaussian = F, Log.scale = T, Boundaries = c(6.5,9.5,11.5),
      #                    W = 1200, H = 750, Save.plot = "Figures/ACA/Trait/Dispersion/Hist_disp_tot_no_scale.pdf")
      # 
      # Trait.distribution(MT = MT.ACA.no.scale, MP1 = MT.bool.ss, MP2 = MT.bool.sl,
      #                    Selec.tax.to.shown = c("Poaceae", "Amaranthaceae", "Artemisia", "Pinus haploxylon", "Pinus diploxylon", "Fabaceae", "Ericaceae"),
      #                    # Selec.tax.to.shown = c("Pinus haploxylon", "Pinus diploxylon"),
      #                    # Selec.tax.to.shown = c("Poaceae", "Amaranthaceae", "Artemisia", "Plantago", "Asteroideae"),
      #                    # Select.trait = c("Height", "LeafArea"),
      #                    Select.trait = c("SSD", "LeafN", "SeedMass", "Height", "LeafArea", "SLA"),
      #                    # Plot.gaussian = F, Type.taxa = c("PT_ss", "PT_sl"),
      #                    Plot.gaussian = F, 
      #                    W = 1000, H = 650, Save.plot = "Figures/ACA/Trait/Dispersion/Hist_disp_ss.pdf")
      # 
      # Trait.distribution(MT = MT.ACA.no.scale, MP1 = MT.bool.ss, MP2 = MT.bool.sl,
      #                    Selec.tax.to.shown = c("Poaceae", "Amaranthaceae", "Artemisia", "Pinus haploxylon", "Pinus diploxylon", "Fabaceae", "Ericaceae"),
      #                    Select.trait = c("SSD", "LeafN", "SeedMass", "Height", "LeafArea", "SLA"),
      #                    Plot.gaussian = T, 
      #                    W = 900, H = 800, Save.plot = "Figures/ACA/Trait/Dispersion/Hist_disp_gauss_ss.pdf")
      
      # Trait.distribution(MT = MT.ACA.no.scale, MP1 = MT.bool.ss, MP2 = MT.bool.sl,
      #                    Selec.tax.to.shown = c("Poaceae", "Amaranthaceae", "Artemisia", "Pinus haploxylon", "Pinus diploxylon", "Fabaceae", "Ericaceae"),
      #                    Select.trait = c("SSD", "LeafN", "SeedMass", "Height", "LeafArea", "SLA"),
      #                    Plot.gaussian = F, 
      #                    W = 800, H = 500, Save.plot = "Figures/ACA/Trait/Dispersion/Hist_disp_sl.pdf")
      
      # Trait.distribution(MT = MT.ACA.no.scale, MP1 = MT.bool.ss, MP2 = MT.bool.sl,
      #                    Selec.tax.to.shown = c("Poaceae", "Amaranthaceae", "Artemisia", "Pinus haploxylon", "Pinus diploxylon", "Fabaceae", "Ericaceae"),
      #                    Select.trait = c("SSD", "LeafN", "SeedMass", "Height", "LeafArea", "SLA"),
      #                    Plot.gaussian = T, 
      #                    W = 900, H = 800, Save.plot = "Figures/ACA/Trait/Dispersion/Hist_disp_gauss_sl.pdf")      
      
    }
    
    #### Save data traits ####
    print("Export data")
    saveRDS(MP_sl, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_sl_clean2.Rds")
    saveRDS(MP_ss, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_ss_clean2.Rds")
    
    saveRDS(Tr.PT_fam, "Resultats/ACA/Traits/Tr.PT_fam.Rds")
    saveRDS(Tr.PT_gen, "Resultats/ACA/Traits/Tr.PT_gen.Rds")
    saveRDS(Tr.PT_ss, "Resultats/ACA/Traits/Tr.PT_ss.Rds")
    saveRDS(Tr.PT_sl, "Resultats/ACA/Traits/Tr.PT_sl.Rds")
    saveRDS(Tr.MV, "Resultats/ACA/Traits/MT_MV_full.Rds")
    saveRDS(Tr.PT_fam_gf, "Resultats/ACA/Traits/Tr.PT_fam_gf.Rds")
    saveRDS(Tr.PT_gen_gf, "Resultats/ACA/Traits/Tr.PT_gen_gf.Rds")
    saveRDS(Tr.PT_ss_gf, "Resultats/ACA/Traits/Tr.PT_ss_gf.Rds")
    saveRDS(Tr.PT_sl_gf, "Resultats/ACA/Traits/Tr.PT_sl_gf.Rds")
    saveRDS(Tr.MV_gf, "Resultats/ACA/Traits/MT_MV_full_gf.Rds")
    }
  else{
    GrowthForm.ACA <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_GrowthForm_ACA.Rds", sep = ""))
    MT.ACA <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_trait_ACA.Rds", sep = ""))
    MT.ACA.gf <- readRDS("Resultats/World_DB/Traits/MT_ACA_gapfilled.Rds")
    Tr.PT_ss <- readRDS("Resultats/ACA/Traits/Tr.PT_ss.Rds")
    Tr.PT_sl <- readRDS("Resultats/ACA/Traits/Tr.PT_sl.Rds")
    Tr.MV_gf <- readRDS("Resultats/ACA/Traits/MT_MV_full_gf.Rds")
    Tr.MV <- readRDS("Resultats/ACA/Traits/MT_MV_full.Rds")
    Tr.PT_ss_gf <- readRDS("Resultats/ACA/Traits/Tr.PT_ss_gf.Rds")
    Tr.PT_sl_gf <- readRDS("Resultats/ACA/Traits/Tr.PT_sl_gf.Rds")
    Tr.PT_fam_gf <- readRDS("Resultats/ACA/Traits/Tr.PT_fam_gf.Rds")
    Tr.PT_fam <- readRDS("Resultats/ACA/Traits/Tr.PT_fam.Rds")
    Tr.PT_gen_gf <- readRDS("Resultats/ACA/Traits/Tr.PT_gen_gf.Rds")
    Tr.PT_gen <- readRDS("Resultats/ACA/Traits/Tr.PT_gen.Rds")
  }
  #### Correlation traits a traits ####
  Trait.cor = F
  if(Trait.cor == T){
    H = 600
    W = 700
    Save.matcor.trait = "Figures/ACA/Trait/Correlation/Mat_correlation/Traits/Matcor_traits_4schemes.pdf"
    pdf(file = Save.matcor.trait, width = W*0.01041666666667, height = H*0.01041666666667)
    par(mfrow = c(2,2))

    MT.full.cor <- Mat.corel.CWT.clim(MT.ACA.gf[c(11,12,5,7,13,4)], MT.ACA.gf[c(11,12,5,7,13,4)],
                                  I.confiance = 0.95,
                                  Display.pval = "pch",
                                  Disp.R = "number", Display = "lower", Average = F, Label = F, Bar.pos = "n",
                                  Save.pval = "Resultats/ACA/Traits/Permutations/Traits_cheklist_perm.csv",
                                  Permutation.test = T, Nb.permutations = 10000, 
                                  
                                  # Title = "Correlation between CWT-ACASP-ss and bioclimate parameters.",
                                  # Save.path = "Resultats/ACA/Traits/Matcorr_CWT_MT_6t.csv",
                                  # Save.plot = "Figures/ACA/Trait/Correlation/Mat_correlation/Traits/Matcor_6traits_MT.pdf",
                                  H = 550,
                                  W = 400)

    MV.gf.cor <- Mat.corel.CWT.clim(Tr.MV_gf[c(9,10,3,5,11,2)], Tr.MV_gf[c(9,10,3,5,11,2)],
                                  I.confiance = 0.95,
                                  Display.pval = "pch",
                                  Disp.R = "number", Display = "lower", Average = F, Label = F, Bar.pos = "n",
                                  Save.pval = "Resultats/ACA/Traits/Permutations/Traits_ACAV_perm.csv",
                                  Permutation.test = T, Nb.permutations = 10000, 
                                  
                                  # Title = "Correlation between CWT-ACASP-ss and bioclimate parameters.",
                                  # Save.path = "Resultats/ACA/Traits/Matcorr_CWT_MV_gf_6t.csv",
                                  # Save.plot = "Figures/ACA/Trait/Correlation/Mat_correlation/Traits/Matcor_6traits_MV_gf.pdf",
                                  H = 550,
                                  W = 400)

    MC.ss.gf.cor <- Mat.corel.CWT.clim(Tr.PT_ss_gf[c(9,10,3,5,11,2)], Tr.PT_ss_gf[c(9,10,3,5,11,2)],
                                  I.confiance = 0.95,
                                  Display.pval = "pch",
                                  Disp.R = "number", Display = "lower", Average = F, Label = F, Bar.pos = "n",
                                  Save.pval = "Resultats/ACA/Traits/Permutations/Traits_ss_gf_perm.csv",
                                  Permutation.test = T, Nb.permutations = 10000, 
                                  
                                  # Title = "Correlation between CWT-ACASP-ss and bioclimate parameters.",
                                  # Save.path = "Resultats/ACA/Traits/Matcorr_CWT_ss_gf_6t.csv",
                                  # Save.plot = "Figures/ACA/Trait/Correlation/Mat_correlation/Traits/Matcor_6traits_ss_gf.pdf",
                                  H = 550,
                                  W = 400)

    MC.sl.gf.cor <- Mat.corel.CWT.clim(Tr.PT_sl_gf[c(9,10,3,5,11,2)], Tr.PT_sl_gf[c(9,10,3,5,11,2)],
                                  I.confiance = 0.95,
                                  Display.pval = "pch",
                                  Disp.R = "number", Display = "lower", Average = F, Label = F, Bar.pos = "n",
                                  Save.pval = "Resultats/ACA/Traits/Permutations/Traits_sl_gf_perm.csv",
                                  Permutation.test = T, Nb.permutations = 10000, 
                                  
                                  # Title = "Correlation between CWT-ACASP-ss and bioclimate parameters.",
                                  # Save.path = "Resultats/ACA/Traits/Matcorr_CWT_sl_gf_6t.csv",
                                  # Save.plot = "Figures/ACA/Trait/Correlation/Mat_correlation/Traits/Matcor_6traits_sl_gf.pdf",
                                  H = 550,
                                  W = 400)

    dev.off()
     
    MC.ss.cor <- Mat.corel.CWT.clim(Tr.PT_ss[-c(1)], Tr.PT_ss[-c(1)],
                                  I.confiance = 0.95, 
                                  Display.pval = "pch", 
                                  Disp.R = "number", Display = "lower", Average = F, 
                                  Save.pval = "Resultats/ACA/Traits/Permutations/Traits_ss_perm.csv",
                                  Permutation.test = T, Nb.permutations = 10000, 
                                  # Title = "Correlation between CWT-ACASP-ss and bioclimate parameters.",
                                  Save.path = "Resultats/ACA/Traits/Matcorr_CWT_ss.csv",
                                  Save.plot = "Figures/ACA/Trait/Correlation/Mat_correlation/Traits/Matcor_traits.pdf",
                                  H = 900,
                                  W = 800)
      
    H = 400
    W = 700
    Save.matcor.trait = "Figures/ACA/Trait/Correlation/R2_compar/Traits/R2_comp_all.pdf"
    #### R2 compar traits ####
    R2.compar.t = F
    if(R2.compar.t == T){
        
      p1 <- R2.compar(MV = MV.gf.cor, MP = MC.ss.gf.cor, Repel.outliers = F,
                Xlab = "r (ACASP-fi)", Ylab = "r (ACAV)", Nb.labs = 4, Panel.lab = "(A)",
                # H = 500, W = 500, Save.plot = "Figures/ACA/Trait/Correlation/R2_compar/Traits/R2_comp.pdf"
                )
      p2 <- R2.compar(MV = MV.gf.cor, MP = MC.sl.gf.cor,
                Xlab = "r (ACASP-co)", Ylab = "r (ACAV)", Nb.labs = 4, Panel.lab = "(B)",
                # H = 500, W = 500, Save.plot = "Figures/ACA/Trait/Correlation/R2_compar/Traits/R2_comp.pdf"
                )
      p3 <- R2.compar(MV = MC.ss.gf.cor, MP = MC.sl.gf.cor,
                Xlab = "r (ACASP-co)", Ylab = "r (ACASP-fi)", Nb.labs = 4,
                # H = 500, W = 500, Save.plot = "Figures/ACA/Trait/Correlation/R2_compar/Traits/R2_comp.pdf"
                )
      pall <- p1 + p2 #+ p3 #+ plot_layout(guides = "collect")
      ggsave(file = Save.matcor.trait, pall, width = W*0.01041666666667, height = H*0.01041666666667)}
  }
 
  #### Comparaison p-values (permutation) ####
  Comp.pval = F
  if(Comp.pval == T){
    pval.Checklist <- readRDS("Resultats/ACA/Traits/Permutations/Traits_cheklist_perm.Rds")
    pval.ACAV <- readRDS("Resultats/ACA/Traits/Permutations/Traits_ACAV_perm.Rds")
    pval.ACASP_fi <- readRDS("Resultats/ACA/Traits/Permutations/Traits_ss_gf_perm.Rds")
    pval.ACASP_co <- readRDS("Resultats/ACA/Traits/Permutations/Traits_sl_gf_perm.Rds")
    
    pval.Checklist <- cbind(suppressWarnings(melt(pval.Checklist)), Sampling = "Cheklist")
    pval.ACAV <- cbind(suppressWarnings(melt(pval.ACAV)), Sampling = "ACAV")
    pval.ACASP_fi <- cbind(suppressWarnings(melt(pval.ACASP_fi)), Sampling = "ACASP-fine")
    pval.ACASP_co <- cbind(suppressWarnings(melt(pval.ACASP_co)), Sampling = "ACASP-coarse")
    Mpval <- rbind(pval.Checklist, pval.ACAV, pval.ACASP_fi, pval.ACASP_co)
    
    # ggplot(Mpval, aes(value, fill = Sampling)) + geom_histogram(binwidth = 0.1)
    
    }
  #### PCA ####
  PCA.trait = F
  if(PCA.trait == T){
    #### Prep data ####
    TPT <- readRDS(file = "Resultats/ACA/Traits/Corresp_tax/Table_corresp_pollen_type.Rds")
    Tr.PT_ss_gf.pca <- cbind(AP_NAP = TPT$AP_NAP[match(Tr.PT_ss_gf$species, TPT$PT_ss)], Tr.PT_ss_gf)
    Tr.PT_sl_gf.pca <- cbind(AP_NAP = TPT$AP_NAP[match(Tr.PT_sl_gf$species, TPT$PT_sl)], Tr.PT_sl_gf)
    Tr.PT_ss.pca <- cbind(AP_NAP = TPT$AP_NAP[match(Tr.PT_ss$species, TPT$PT_ss)], Tr.PT_ss)
    Tr.PT_sl.pca <- cbind(AP_NAP = TPT$AP_NAP[match(Tr.PT_sl$species, TPT$PT_sl)], Tr.PT_sl)
    
    MT.ACA.gf.pca <- cbind(GrowthForm = GrowthForm.ACA$GrowthForm[match(MT.ACA.gf$species, GrowthForm.ACA$species)], MT.ACA.gf)
    Tr.MV_gf.pca <- cbind(GrowthForm = GrowthForm.ACA$GrowthForm[match(Tr.MV_gf$species, GrowthForm.ACA$species)], Tr.MV_gf)
    MT.ACA.pca <- cbind(GrowthForm = GrowthForm.ACA$GrowthForm[match(MT.ACA$species, GrowthForm.ACA$species)], MT.ACA)
    Tr.MV.pca <- cbind(GrowthForm = GrowthForm.ACA$GrowthForm[match(Tr.MV$species, GrowthForm.ACA$species)], Tr.MV)
    
    Trait.to.keep.pca1 <- c(1,5,6,8,12,13,14)
    Trait.to.keep.pca2 <- c(1,3,4,6,10,11,12)
    
    #### PCA plots ####
    BioCl.pca.full.MT <- PCA.bioclim(MT.ACA.gf.pca[Trait.to.keep.pca1], Dot.opac = 0.65, Dot.size = 1, Cluster.core = "GrowthForm",
                                     Ellipse = F, Cluster.core.lab = "Growth Forms", 
                                     Density.contour = F, Opa.range = c(0.75,1), Density.type = "contour",
                                     transp_OK = T, Scale.PCA = 10, return.pick = T, Num.facet = "(A)",
                                     Site.name = "All ACA plant checklist", Legend.position = "none",
                                     # Save.plot = "Figures/ACA/Trait/Correlation/PCA/Traits/PCA_trait_ACA_fullplantlist.pdf", H = 650, W = 650
                                     )

    BioCl.pca.MV.gf <- PCA.bioclim(Tr.MV_gf.pca[Trait.to.keep.pca2], Dot.opac = 0.65, Dot.size = 1, Cluster.core = "GrowthForm",
                                     Ellipse = F, Cluster.core.lab = "Growth Forms", Legend.size = 10, 
                                     Density.contour = F, Opa.range = c(0,1), Density.type = "contour",
                                     transp_OK = T, Scale.PCA = 10, return.pick = T, Num.facet = "(B)",
                                     Site.name = "ACAV plant checklist (gap-filled)", #Reverse.dim = T, 
                                     # Manu.lim.x = c(7.5,-5),
                                     # Manu.lim.y = c(7.5,-5)
                                     # Save.plot = "Figures/ACA/Trait/Correlation/PCA/Traits/PCA_trait_ACA_vegetation_plot_taxa_gf.pdf", H = 650, W = 650
                                   )

    BioCl.pca.ss.gf <- PCA.bioclim(Tr.PT_ss_gf.pca[Trait.to.keep.pca2], Cluster.core = "AP_NAP", Dot.size = 1.5, Dot.opac = 0.8,
                                  return.pick = T, Ellipse = F, Cluster.core.lab = "Pollen types", Num.facet = "(C)",
                                  Legend.position = "none", Density.contour = F, Opa.range = c(0.05,0.3),
                                  Site.name = "ACASP-fine (gap-filled)", transp_OK = T, Scale.PCA = 5,
                                  # Save.plot = "Figures/ACA/Trait/Correlation/PCA/Traits/PCA_trait_ACA_ss_gf.pdf", H = 650, W = 650
                                  )


    BioCl.pca.sl.gf <- PCA.bioclim(Tr.PT_sl_gf.pca[Trait.to.keep.pca2], Cluster.core = "AP_NAP", Dot.size = 1.5, Dot.opac = 0.8,
                                  return.pick = T, Ellipse = F, Cluster.core.lab = "Pollen types", Legend.size = 10,
                                  Density.contour = F, Opa.range = c(0.03,0.15), Num.facet = "(D)",
                                  Site.name = "ACASP-coarse (gap-filled)", transp_OK = T, Scale.PCA = 5,
                                  # Save.plot = "Figures/ACA/Trait/Correlation/PCA/Traits/PCA_trait_ACA_sl_gf.pdf", H = 650, W = 1000
                                  )

    # pca.full <- (BioCl.pca.full.MT + BioCl.pca.MV.gf) /(BioCl.pca.ss.gf + BioCl.pca.sl.gf) / plot_layout(guides = "keep")
    pca.full <- BioCl.pca.full.MT | BioCl.pca.MV.gf | BioCl.pca.ss.gf | BioCl.pca.sl.gf | plot_layout(guides = "keep")
    W = 1700
    H = 380
    # W = 850
    # H = 720
    ggsave(filename = "Figures/ACA/Trait/Correlation/PCA/Traits/PCA_trait_full_gf.pdf", pca.full, width = W*0.026458333, height = H*0.026458333, units = "cm")
  
    
   
      }
  
  #### Test Box's M test ####
  Box.M = F
  if(Box.M == T){
    library(heplots)
    
    Trait.to.keep.pca1 <- c(5,6,8,12,13,14)
    Trait.to.keep.pca2 <- c(3,4,6,10,11,12)
    
    M1 = MT.ACA.gf.pca[Trait.to.keep.pca1]
    M2 = Tr.MV_gf.pca[Trait.to.keep.pca2]
    M3 = Tr.PT_ss_gf.pca[Trait.to.keep.pca2]
    M4 = Tr.PT_sl_gf.pca[Trait.to.keep.pca2]
    
    M1$Sampling <- "Checklist"
    M2$Sampling <- "ACAV"
    M3$Sampling <- "ACASP-fine"
    M4$Sampling <- "ACASP-coarse"
    M <- rbind(M1, M2, M3, M4)
    M <- na.omit(M)
    
    res <- boxM(M[!names(M) %in% "Sampling"], M[, "Sampling"])
    
    A = summary(res)
    
    #### M calculation ####
    Tab.n = count(M, c("Sampling"))
    Tab.n$nn <- Tab.n$freq -1
    Tab.n <- rbind(Tab.n, c("Pooled", sum(Tab.n$freq), sum(Tab.n$nn)))
    Tab.n$ln.det <- c(res$logDet)
    Tab.n$nn.ln.det <- Tab.n$ln.det * as.numeric(Tab.n$nn)
    M.statistic = Tab.n$nn.ln.det[Tab.n$Sampling == "Pooled"] - sum(Tab.n$nn.ln.det[Tab.n$Sampling != "Pooled"])
    
    # visualize (what is done in the plot method) 
    dets <- res$logDet
    dets <- dets[order(dets)]
    ng <- length(res$logDet)-1
    dotchart(dets, xlab = "log determinant")
    points(dets , 1:length(dets),  
           cex=c(rep(1.5, ng), 2.5), 
           pch=c(rep(16, ng), 15),
           col= c(rep("blue", ng), "red"))
  }
  
  #### Test permutation generale ####
  Test.permut = F
  if(Test.permut == T){
    #### import data ####
    Trait.to.keep.pca1 <- c(4,5,7,11,12,13)
    Trait.to.keep.pca2 <- c(2,3,5,9,10,11)
    
    M1 = MT.ACA.gf[Trait.to.keep.pca1]
    M2 = Tr.MV_gf[Trait.to.keep.pca2]
    M3 = Tr.PT_ss_gf[Trait.to.keep.pca2]
    M4 = Tr.PT_sl_gf[Trait.to.keep.pca2]
    
    M1$Sampling <- "Checklist"
    M2$Sampling <- "ACAV"
    M3$Sampling <- "ACASP-fine"
    M4$Sampling <- "ACASP-coarse"
    M <- rbind(M1, M2, M3, M4)
    M <- na.omit(M)
    
    # Relation.trait.sampling <- function(M,Groups)
    Groups = "Sampling"
    Nb.simul = 2
    Nb.permu = 1
    alpha = 0.1
    alpha.r = 0.15
    names(M)[names(M) == Groups] <- "My_Group"
    Traits <- names(M)[names(M) != "My_Group"]
    M$My_Group <- factor(M$My_Group, levels = unique(M$My_Group), ordered = T)
    Group_type <- M[!duplicated(M$My_Group), !names(M) %in% c(Traits)]

    print(M)
    #### Measuring trait-based functional indices ####
# 
#     for (i in 1:length(Group_type)) {
#       presence <- subset(M, My_Group == Group_type[i]) %>%
#         mutate(present = 1) %>%
#         pivot_wider(names_from = TreeID, values_from = present)  # Presence/absence matrix for a tree
# 
#       vector <- data.frame(t(presence[, 17]))
#       names(vector) <- presence$LeafID
#       leaves.distance <- gowdis(presence[, 8:16]) # Leaf distance matrix
#       leaves.distance <- as.dist(as.matrix(leaves.distance))
#       attr(leaves.distance, "Labels") <- names(vector)
# 
#       functional <- dbFD(leaves.distance, vector,
#                          calc.FGR = FALSE, calc.CWM = FALSE, calc.FDiv = FALSE,
#                          corr = "lingoes") # We assess the functional indices
# 
#       FRic[i] <- functional$FRic # We store the results for FRic
#       FDis[i] <- functional$FDis # We store the results for FDis
#     }

    
    #### SD for each traits ####
    Mat.sd <- data.frame(Group_type = Group_type)
    for(i in 1:length(Traits)){
      sd.trait <- aggregate(M[Traits[i]], by = list(M$My_Group), FUN = sd)
      Mat.sd <- cbind(Mat.sd, sd.trait[c(2)])
      }
    
    #### Network properties between traits ####
    trait.properties <- matrix(NA, length(Group_type), 1+length(Traits))
    for (i in 1:length(Group_type)) {
      significant.correlations <- matrix(NA, length(Traits), length(Traits))
      colnames(significant.correlations) <- Traits
      rownames(significant.correlations) <- Traits
      
      for (k in 1:length(Traits)) {
        for (j in k:length(Traits)) {
          set.seed(0)
          significant.correlations[k, j] <- perm.cor.test(subset(M, My_Group == Group_type[i])[, k],
                                                          subset(M, My_Group == Group_type[i])[, j], num.sim = Nb.permu)$p.value
        }}
      
      correlations <- cor(subset(M, My_Group == Group_type[i])[, 1:length(Traits)])
      
      p.val <- c()
      r <- c()
      Var1 <- c() 
      Var2 <- c()
      for (k in 1:length(Traits)) {
        p.val <- c(p.val, significant.correlations[k, k:length(Traits)]) 
        r <- c(r, correlations[k, k:length(Traits)]) 
        Var1 <- c(Var1, rep(rownames(correlations)[k], (length(Traits)+1-k)))
        Var2 <- c(Var2, rownames(correlations)[k:length(Traits)])
        }
      
      cor <- data.frame(cbind(Var1, Var2, p.val, r))
      cor$p.val <- as.numeric(cor$p.val)
      cor$r <- as.numeric(cor$r)
       
      # Based on our selected correlations (the significant ones) we build our network
      my_adj_list <- cor %>% filter(p.val < alpha)
      my_adj_list <- my_adj_list %>% filter(abs(r) > alpha.r)
      names(my_adj_list) <- c("from", "to", "weight", "r") # We define the origin and end of our nodes
      net <- graph.data.frame(my_adj_list, directed = FALSE) # We build the network
      HOP <- c(edge_density(net, loops = FALSE), degree(net)) # We store edge density and degree for every tree in our matrix
      trait.properties[i,] <- c(edge_density(net, loops = FALSE), degree(net)) # We store edge density and degree for every tree in our matrix
      }
    
    trait.properties <- as.data.frame(trait.properties)
    names(trait.properties) <- c("edge.density", paste("deg", Traits, sep = "."))
    
    #### Observed correlations between metrics ####
    Cor <- c()
    for(i in 2:(length(Traits)+1)){Cor <- c(Cor, suppressWarnings(cor(trait.properties[[i]], Mat.sd[[i]])))}
    # Cor[is.na(Cor)] <- 1
    names(Cor) <- Traits
    Cor <- melt(Cor)
    Cor$variable <- row.names(Cor)
    
    #### On fait la même chose mais randomisé ####
    random.cor.traits <- matrix(NA, Nb.simul, length(Traits)) # Empty vector to store the values of pearson correlations
    for(p in 1:Nb.simul) {
      #### Creation d'une matrice randomisée ####
      Mat.rand <- M["My_Group"]
      for(i in 1:length(Traits)){
        set.seed(p + sample.int(1000, 1))
        Mat.rand <- cbind(Mat.rand, sample(M[[i]]))
      }
      names(Mat.rand) <- c("My_Group", Traits)
      
      ##### TRUC A COMPRENDRE PLUS TARD ####
      # random.FRic <- rep(NA, Nb.truc) # Empty vector to store the FRic data
      # random.FDis <- rep(NA, Nb.truc) # Empty vector to store the FDis data
      # 
      # for (i in 1:Nb.truc) {
      #   presence <- subset(data, TreeID == trees[i, ]$TreeID) %>%
      #     mutate(present = 1) %>%
      #     pivot_wider(names_from = TreeID, values_from = present)  # Presence/absence matrix for a tree
      #   
      #   vector <- data.frame(t(presence[, 17]))
      #   names(vector) <- presence$LeafID              
      #   leaves.distance <- gowdis(presence[, 17:25]) # Leaf distance matrix 
      #   leaves.distance <- as.dist(as.matrix(leaves.distance))
      #   attr(leaves.distance, "Labels") <- names(vector)
      #   
      #   functional <- dbFD(leaves.distance, vector, 
      #                      calc.FGR = FALSE, calc.CWM = FALSE, calc.FDiv = FALSE,
      #                      corr = "lingoes") # We assess the functional indices
      #   
      #   random.FRic[i] <- functional$FRic # We store the results for FRic
      #   random.FDis[i] <- functional$FDis # We store the results for FDis
      # }
      
      
      #### SD for each traits ####
      Mat.sd.random <- data.frame(Group_type = Group_type)
      for(i in 1:length(Traits)){
        sd.trait <- aggregate(Mat.rand[Traits[i]], by = list(Mat.rand$My_Group), FUN = sd)
        Mat.sd.random <- cbind(Mat.sd.random, sd.trait[c(2)])
      }
    
      Mat.rand <- Mat.rand[c(2:ncol(Mat.rand),1)]
      # print(head(M))
      # print(head(Mat.sd.random))
      
      #### Network properties between traits ####
      random.trait.properties <- matrix(NA, length(Group_type), 1+length(Traits))
      for (i in 1:length(Group_type)){
        significant.correlations <- matrix(NA, length(Traits), length(Traits))
        colnames(significant.correlations) <- Traits
        rownames(significant.correlations) <- Traits
        for (k in 1:length(Traits)){
          for (j in k:length(Traits)){
            set.seed(0)
            significant.correlations[k, j] <- perm.cor.test(subset(Mat.rand, My_Group == Group_type[i])[, k],
                                                            subset(Mat.rand, My_Group == Group_type[i])[, j], num.sim = Nb.permu)$p.value
          }}
        correlations <- cor(subset(Mat.rand, My_Group == Group_type[i])[, 1:length(Traits)])
        p.val <- c()
        r <- c()
        Var1 <- c()
        Var2 <- c()
        for (k in 1:length(Traits)) {
          p.val <- c(p.val, significant.correlations[k, k:length(Traits)])
          r <- c(r, correlations[k, k:length(Traits)])
          Var1 <- c(Var1, rep(rownames(correlations)[k], (length(Traits)+1-k)))
          Var2 <- c(Var2, rownames(correlations)[k:length(Traits)])
        }
        cor <- data.frame(cbind(Var1, Var2, p.val, r))
        cor$p.val <- as.numeric(cor$p.val)
        cor$r <- as.numeric(cor$r)
        # print(cor)

        # Based on our selected correlations (the significant ones) we build our network
        # my_adj_list <- cor %>% filter(p.val < 0.05)
        my_adj_list <- cor %>% filter(p.val < alpha)
        # my_adj_list <- my_adj_list %>% filter(abs(r) > 0.6)
        my_adj_list <- my_adj_list %>% filter(abs(r) > alpha.r)
        names(my_adj_list) <- c("from", "to", "weight", "r") # We define the origin and end of our nodes
        net <- graph.data.frame(my_adj_list, directed = FALSE) # We build the network
        HOP <- c(edge_density(net, loops = FALSE), degree(net)) # We store edge density and degree for every tree in our matrix
        random.trait.properties[i,] <- c(edge_density(net, loops = FALSE), degree(net)) # We store edge density and degree for every tree in our matrix
        }

      random.trait.properties <- as.data.frame(random.trait.properties)
      names(random.trait.properties) <- c("edge.density", paste("deg", Traits, sep = "."))
      # print(random.trait.properties)
      #### Observed correlations between metrics ####
      Random.Cor <- c()
      for(i in 2:(length(Traits)+1)){Random.Cor <- c(Random.Cor, suppressWarnings(cor(random.trait.properties[[i]], Mat.sd.random[[i]])))}
      # Random.Cor[is.na(Random.Cor)] <- 1
      random.cor.traits[p,] <- Random.Cor
      # print(Random.Cor)
      print(c("p", p))
      }
    
    random.cor.traits <- setNames(as.data.frame(random.cor.traits), Traits) 
    random.cor.traits <- melt(random.cor.traits)
    
    #### Plots ####
    p <- ggplot(random.cor.traits, aes(x = value, fill = variable)) +
      # annotate("rect", xmin = quantile(random.cor.fdis, 0.025), xmax = quantile(random.cor.fdis, 0.975), 
      #          ymin = 0, ymax = 40, alpha = 0.2) +
      # geom_vline(aes(xintercept = quantile(random.cor.fdis, 0.025)), color = "grey80", size = 1) +
      # geom_vline(aes(xintercept = quantile(random.cor.fdis, 0.975)), color = "grey80", size = 1) +
      geom_histogram(color = "black", fill = "white", binwidth = 0.1) +
      geom_vline(aes(xintercept = 0), linetype = "dotted", size = 1.1) +
      facet_wrap(vars(variable))+
      geom_vline(data = Cor, aes(xintercept = value), linetype = "dashed", size = 1.1, color = "red") +
      theme_bw() +
      # coord_cartesian(ylim = c(1.65, 35), xlim = c(-0.5, 0.5)) +
      labs(y = "Count", x = "")
    print(p)
   
    
    }
  #### Stats cover and export ####
  Stats.cover = F
  if(Stats.cover == T){
    #### Functions calcul de cover ####
    Cal.matrice.NA.pourc <- function(M){
      To.remove <- c("species", "Other.clade", "genus", "family", "PT_sl", "PT_ss", "CNRatio", "LPP", "Wood", "LeafPhenology", "LeafP")
      M <- M[,!names(M) %in% To.remove]
      Tot.trait.val <- prod(dim(M)[1], dim(M)[2])
      Tot.empty <- length(which(is.na(M)))  
      Pour.remplissage <- round((Tot.trait.val - Tot.empty)*100/Tot.trait.val, digits = 0) 
      # print(Pour.remplissage)
      return(Pour.remplissage)}
    
    Cal.trait.NA.pourc <- function(M){
      To.remove <- c("species", "Other.clade", "genus", "family", "PT_sl", "PT_ss", "CNRatio", "LPP", "Wood", "LeafPhenology", "LeafP")
      M <- M[,!names(M) %in% To.remove]
      Save.val <- c()
      for (i in 1:ncol(M)){
        # print(paste("Calculation for the trait", names(M)[i]))
        Tot.trait.val <- length(M[[i]])
        Tot.empty <- length(which(is.na(M[[i]])))
        Pour.remplissage <- round((Tot.trait.val - Tot.empty)*100/Tot.trait.val, digits = 0) 
        Save.val[i] <- Pour.remplissage
        }
      Save.val <- data.frame(Save.val)
      row.names(Save.val) <- names(M)
      return(Save.val)}
   
    #### Calcule par matrice type ####
    PR.MT.ACA <- Cal.matrice.NA.pourc(MT.ACA)
    PR.MT.ACA.gf <- Cal.matrice.NA.pourc(MT.ACA.gf)
    PR.MV.ACA <- Cal.matrice.NA.pourc(Tr.MV)
    PR.MV.ACA.gf <- Cal.matrice.NA.pourc(Tr.MV_gf)
    PR.Tr.PT_ss <- Cal.matrice.NA.pourc(Tr.PT_ss)
    PR.Tr.PT_ss_gf <- Cal.matrice.NA.pourc(Tr.PT_ss_gf)
    PR.Tr.PT_sl <- Cal.matrice.NA.pourc(Tr.PT_sl)
    PR.Tr.PT_sl_gf <- Cal.matrice.NA.pourc(Tr.PT_sl_gf)
    PR.Tr.PT_gen <- Cal.matrice.NA.pourc(Tr.PT_gen)
    PR.Tr.PT_gen_gf <- Cal.matrice.NA.pourc(Tr.PT_gen_gf)
    PR.Tr.PT_fam <- Cal.matrice.NA.pourc(Tr.PT_fam)
    PR.Tr.PT_fam_gf <- Cal.matrice.NA.pourc(Tr.PT_fam_gf)
    
    Stats.pour.remp <- data.frame(Matrix = c("MT.ACA", "MV.ACA", "Tr.PT_ss", "Tr.PT_sl", "Tr.PT_fam", "Tr.PT_gen"),
                                  Pour.remp = c(PR.MT.ACA, PR.MV.ACA, PR.Tr.PT_ss, PR.Tr.PT_sl, PR.Tr.PT_fam, PR.Tr.PT_gen),
                                  Pour.remp.gf = c(PR.MT.ACA.gf, PR.MV.ACA.gf, PR.Tr.PT_ss_gf, PR.Tr.PT_sl_gf, PR.Tr.PT_fam_gf, PR.Tr.PT_gen_gf)
                                  )
    
    Line.tot <- data.frame(PR.MT.ACA, PR.MT.ACA.gf, PR.MV.ACA, PR.MV.ACA.gf, PR.Tr.PT_ss, PR.Tr.PT_ss_gf, PR.Tr.PT_sl, PR.Tr.PT_sl_gf, PR.Tr.PT_fam, PR.Tr.PT_fam_gf, PR.Tr.PT_gen, PR.Tr.PT_gen_gf)
    row.names(Line.tot) <- "Mean coverage"
    
    PR.MT.ACA <- Cal.trait.NA.pourc(MT.ACA)
    PR.MT.ACA.gf <- Cal.trait.NA.pourc(MT.ACA.gf)
    PR.MV.ACA <- Cal.trait.NA.pourc(Tr.MV)
    PR.MV.ACA.gf <- Cal.trait.NA.pourc(Tr.MV_gf)
    PR.Tr.PT_ss <- Cal.trait.NA.pourc(Tr.PT_ss)
    PR.Tr.PT_ss_gf <- Cal.trait.NA.pourc(Tr.PT_ss_gf)
    PR.Tr.PT_sl <- Cal.trait.NA.pourc(Tr.PT_sl)
    PR.Tr.PT_sl_gf <- Cal.trait.NA.pourc(Tr.PT_sl_gf)
    PR.Tr.PT_gen <- Cal.trait.NA.pourc(Tr.PT_gen)
    PR.Tr.PT_gen_gf <- Cal.trait.NA.pourc(Tr.PT_gen_gf)
    PR.Tr.PT_fam <- Cal.trait.NA.pourc(Tr.PT_fam)
    PR.Tr.PT_fam_gf <- Cal.trait.NA.pourc(Tr.PT_fam_gf)
    
    Pour.remp.tr = data.frame(PR.MT.ACA, PR.MT.ACA.gf, PR.MV.ACA, PR.MV.ACA.gf, PR.Tr.PT_ss, PR.Tr.PT_ss_gf, PR.Tr.PT_sl, PR.Tr.PT_sl_gf, PR.Tr.PT_fam, PR.Tr.PT_fam_gf,  PR.Tr.PT_gen, PR.Tr.PT_gen_gf)
    Stats.pour.remp.tr <- c("Traits ACA","Traits ACA (gf)", "ACAV","ACAV (gf)", "ACASP-ss", "ACASP-ss (gf)", "ACASP-sl", "ACASP-sl (gf)", "ACASP-family", "ACASP-family (gf)", "ACASP-genus", "ACASP-genus (gf)")
    names(Pour.remp.tr) <- Stats.pour.remp.tr
    names(Line.tot) <- Stats.pour.remp.tr
    Pour.remp.tr <- Pour.remp.tr[c(1:3,6:8),]
    # Pour.remp.tr <- Pour.remp.tr[order(Pour.remp.tr$`Traits ACA`),]
    Pour.remp.tr <- Pour.remp.tr[order(row.names(Pour.remp.tr)),]
    Pour.remp.tr <- rbind(Pour.remp.tr, Line.tot)
    Aver.by.tr <- data.frame(round(rowMeans(Pour.remp.tr), digits = 0))
    names(Aver.by.tr) <- "Mean coverage"
    Pour.remp.tr <- cbind(Pour.remp.tr, Aver.by.tr)
    Pour.remp.tr <- Pour.remp.tr[,c(1:8,13)]
    
    #### Export ####
    P.fam.sl <- round(length(Tr.PT_sl$species[grepl("aceae", Tr.PT_sl$species)])/length(Tr.PT_sl$species)*100, digits = 0)
    P.fam.ss <- round(length(Tr.PT_ss$species[grepl("aceae", Tr.PT_ss$species)])/length(Tr.PT_ss$species)*100, digits = 0)
    print(paste("Percentage of pollen-types coarse has family", P.fam.sl))
    print(paste("Percentage of pollen-types fine has family", P.fam.ss))
    print(Stats.pour.remp)
    print(Pour.remp.tr)
    write.csv(Pour.remp.tr, "Resultats/ACA/Traits/Statistic_coverage.csv")
    }

  }

Clean.DB.pol = F
if(Clean.DB.pol == T){
  #### Import data base ####
  if(exists("MV.releve") == F){MV.releve <- readRDS("Resultats/ACA/Vegetation/MV_ACA_plot_only_spermatophyte.Rds")}
  MP_sl <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_sl_clean.Rds")
  MP_ss <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_ss_clean.Rds")
  Tr.PT_ss <- readRDS("Resultats/ACA/Traits/Tr.PT_ss.Rds")
  Tr.PT_sl <- readRDS("Resultats/ACA/Traits/Tr.PT_sl.Rds")
  Tr.MV_gf <- readRDS("Resultats/ACA/Traits/MT_MV_full_gf.Rds")
  Tr.MV <- readRDS("Resultats/ACA/Traits/MT_MV_full.Rds")
  Tr.PT_ss_gf <- readRDS("Resultats/ACA/Traits/Tr.PT_ss_gf.Rds")
  Tr.PT_sl_gf <- readRDS("Resultats/ACA/Traits/Tr.PT_sl_gf.Rds")
  Tr.PT_fam_gf <- readRDS("Resultats/ACA/Traits/Tr.PT_fam_gf.Rds")
  Tr.PT_fam <- readRDS("Resultats/ACA/Traits/Tr.PT_fam.Rds")
  Tr.PT_gen_gf<- readRDS("Resultats/ACA/Traits/Tr.PT_gen_gf.Rds")
  Tr.PT_gen<- readRDS("Resultats/ACA/Traits/Tr.PT_gen.Rds")
  Tr.PT_fam_gf <- readRDS("Resultats/ACA/Traits/Tr.PT_fam_gf.Rds")
  Tr.PT_fam <- readRDS("Resultats/ACA/Traits/Tr.PT_fam.Rds")
  Tr.PT_gen_gf<- readRDS("Resultats/ACA/Traits/Tr.PT_gen_gf.Rds")
  Tr.PT_gen <- readRDS("Resultats/ACA/Traits/Tr.PT_gen.Rds")
  APG_diff <- data.frame(read.csv("Import/World_DB/Taxonomie/APG_II_APG_IV_diff.csv" ,sep=",",dec=".",header=T, stringsAsFactors = F))

  #### Functions ####
  Check.missing.trait <- function(MP, MT, TNRS.check, return.missing){
    #### Check missing ####
    if(missing(TNRS.check)){TNRS.check = F}
    if(missing(return.missing)){return.missing = F}
    Missing.trait <- setdiff(names(MP),MT[[1]])
    print("*** Previously missing ***")
    # print(length(Missing.trait))
    print(Missing.trait)
    
    #### Transpose ####
    t_df <- as.data.frame(t(MP),row.names = NA)
    t_df <- cbind(names(MP),t_df)
    names(t_df)[1] <- "ID"
    
    #### TNRS missing ####
    if(TNRS.check == T){
      Missing.trait.tnrs <- TNRS.taxa.accepted(Taxa.vector = Missing.trait)
      print(Missing.trait.tnrs)
      Missing.trait.tnrs <- Missing.trait.tnrs[Missing.trait.tnrs$match.score >= 0.80,]
      t_df$ID[which(t_df$ID %in% Missing.trait.tnrs$Name_submitted)] <- Missing.trait.tnrs$Accepted_name[which(Missing.trait.tnrs$Name_submitted %in% t_df$ID)]
      }
    
    #### Check diff APG II et APG IV    
    t_df$ID[which(t_df$ID %in% APG_diff$Old_name)] <- APG_diff$New_name[which(APG_diff$Old_name %in% t_df$ID)]
  
    #### Aggregate ####
    new_df <- aggregate(t_df[-1], by=list(t_df[["ID"]]), sum)
    t_df2 <- as.data.frame(t(new_df[-1]),row.names = NA)
    names(t_df2) <- new_df$Group.1
    row.names(t_df2) <- names(new_df)[-1]
    MP <- t_df2
    
    Missing.trait <- setdiff(names(MP),MT[[1]])
    
    print("*** Still missing ***")
    # print(length(Missing.trait))
    print(Missing.trait)
    
    if(return.missing == T){return(Missing.trait)}
    else{return(MP)}
  }
  
  #### Applications ####
  MP_sl.c <- Check.missing.trait(MP_sl, Tr.PT_sl, TNRS.check = T)
  MP_ss.c <- Check.missing.trait(MP_ss, Tr.PT_ss, TNRS.check = T)
  Missing_taxa_trait_sl <- Check.missing.trait(MP_sl, Tr.PT_sl, TNRS.check = T, return.missing = T)
  Missing_taxa_trait_ss <- Check.missing.trait(MP_ss, Tr.PT_ss, TNRS.check = T, return.missing = T)
  
  MP_ss.c$Polygonum <- MP_ss.c$Polygonum + MP_ss.c$Knorringia
  MP_ss.c <- subset(MP_ss.c, select = -c(Knorringia))

  names(MV.releve)[names(MV.releve) == "Leguminosae"] <- "Fabaceae" 
  
  #### Export ####
  saveRDS(MP_sl.c, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_sl_clean2.Rds")
  saveRDS(MP_ss.c, "Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_ss_clean2.Rds")
  saveRDS(Missing_taxa_trait_sl, "Resultats/World_DB/Traits/Missing_taxa_trait_sl.Rds")
  saveRDS(Missing_taxa_trait_ss, "Resultats/World_DB/Traits/Missing_taxa_trait_ss.Rds")
  saveRDS(MV.releve, "Resultats/ACA/Vegetation/MV_ACA_plot_only_spermatophyte.Rds")
  }

CWT.calcul = F
if(CWT.calcul == T){
  #### Settings ####
  if(exists("MV_climtot") == F){
    #### Import data base ####
    MP_sl <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_sl_clean2.Rds")
    MP_ss <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_ss_clean2.Rds")
    
    if(exists("MV.releve") == F){MV.releve <- readRDS("Resultats/ACA/Vegetation/MV_ACA_plot_only_spermatophyte.Rds")}
    MV.ACA.Clim_chel <- readRDS("Resultats/ACA/Vegetation/DBACA_Clim_chel.Rds")
    MV.ACA.Clim_wc <- readRDS("Resultats/ACA/Vegetation/DBACA_Clim_wc.Rds")
    MV.ACA.Biom <- readRDS("Resultats/ACA/Vegetation/DBACA_Biom.Rds")
    MV.ACA.Land <- readRDS("Resultats/ACA/Vegetation/DBACA_LandCover.Rds")
    MV.ACA.Biom <- merge(MV.ACA.Land, MV.ACA.Biom, by = 0)
    row.names(MV.ACA.Biom) <- MV.ACA.Biom$Row.names
    MV.ACA.Biom <- MV.ACA.Biom[-c(1,6,7,8)]
    names(MV.ACA.Biom)[1:2] <- c("Longitude", "Latitude")
    
    MPS.ACA.Clim_wc <- data.frame(readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Clim_wc.Rds"))
    MPS.ACA.Clim_chel <- data.frame(readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Clim_chel.Rds"))
    MPS.ACA.Biom <- data.frame(readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Biom.Rds"))
    MPS.ACA.Land <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_LandCover.Rds")
    MPS.ACA.Biom <- merge(MPS.ACA.Land, MPS.ACA.Biom, by = 0)
    row.names(MPS.ACA.Biom) <- MPS.ACA.Biom$Row.names
    MPS.ACA.Biom <- MPS.ACA.Biom[-c(1,6,7,8)]
    names(MPS.ACA.Biom)[1:2] <- c("Longitude", "Latitude")
    
    Tr.PT_ss <- readRDS("Resultats/ACA/Traits/Tr.PT_ss.Rds")
    Tr.PT_sl <- readRDS("Resultats/ACA/Traits/Tr.PT_sl.Rds")
    Tr.PT_ss_gf <- readRDS("Resultats/ACA/Traits/Tr.PT_ss_gf.Rds")
    Tr.PT_sl_gf <- readRDS("Resultats/ACA/Traits/Tr.PT_sl_gf.Rds")
    Tr.MV_gf <- readRDS("Resultats/ACA/Traits/MT_MV_full_gf.Rds")
    Tr.MV <- readRDS("Resultats/ACA/Traits/MT_MV_full.Rds")
    MT.ACA.gf <- readRDS("Resultats/World_DB/Traits/MT_ACA_gapfilled.Rds")
    MT.ACA <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_trait_ACA.Rds", sep = ""))
  
    #### Clim merge ####
    MP_climtot <- MPS.ACA.Clim_chel
    names(MP_climtot)[4:9] <- paste(names(MP_climtot)[4:9], "chel", sep = "_")
    MP_climtot <- cbind(MP_climtot, MPS.ACA.Clim_wc[4:9])
    names(MP_climtot)[11:16] <- paste(names(MP_climtot)[11:16], "wc", sep = "_")
    MP_climtot <- MP_climtot[order(names(MP_climtot))]
    
    MV_climtot <- MV.ACA.Clim_chel
    names(MV_climtot)[4:9] <- paste(names(MV_climtot)[4:9], "chel", sep = "_")
    MV_climtot <- cbind(MV_climtot, MV.ACA.Clim_wc[4:9])
    names(MV_climtot)[11:16] <- paste(names(MV_climtot)[11:16], "wc", sep = "_")
    MV_climtot <- MV_climtot[order(names(MV_climtot))]
    saveRDS(MV_climtot, "Resultats/ACA/Vegetation/DBACA_Clim_tot.Rds")
    saveRDS(MP_climtot, "Resultats/World_DB/Pollen/Merging_DB/DB13549/MP_climtot.Rds")
    }
  
  #### CWT calculation ####
  Seuil.pcover = 0.2 # ?
  # Seuil.pcover = 0.5 # ?
  # Seuil.pcover = 0.6 # Borgy et al 2014
  # Seuil.pcover = 0.7 # Garnier personal. com.
  # Seuil.pcover = 0.8 # Garnier personal. com.
  MCWT.clim.PT_ss <- CWT.calculation(MT = Tr.PT_ss$species, MP = MP_ss, Mclim = MP_climtot, MPS.ACA.Biom = MPS.ACA.Biom, Accep.seuil = Seuil.pcover, Remove.biom = c("N/A", "Tundra", NA))
  MCWT.clim.PT_ss_gf <- CWT.calculation(MT = Tr.PT_ss_gf, MP = MP_ss, Mclim = MP_climtot, MPS.ACA.Biom = MPS.ACA.Biom, Accep.seuil = Seuil.pcover, Remove.biom = c("N/A", "Tundra", NA))
  MCWT.clim.PT_sl <- CWT.calculation(MT = Tr.PT_sl, MP = MP_sl, Mclim = MP_climtot, MPS.ACA.Biom = MPS.ACA.Biom, Accep.seuil = Seuil.pcover, Remove.biom = c("N/A", "Tundra", NA))
  MCWT.clim.PT_sl_gf <- CWT.calculation(MT = Tr.PT_sl_gf, MP = MP_sl, Mclim = MP_climtot, MPS.ACA.Biom = MPS.ACA.Biom, Accep.seuil = Seuil.pcover, Remove.biom = c("N/A", "Tundra", NA))
  # MCWT.clim.MV <- CWT.calculation(MT = Tr.MV, MP = MV.releve, Mclim = MV_climtot, MPS.ACA.Biom = MV.ACA.Biom, Accep.seuil = Seuil.pcover, Remove.biom = c("N/A", "Tundra", NA))
  MCWT.clim.MV_gf <- CWT.calculation(MT = Tr.MV_gf, MP = MV.releve, Mclim = MV_climtot, MPS.ACA.Biom = MV.ACA.Biom, Accep.seuil = Seuil.pcover, Remove.biom = c("N/A", "Tundra", NA))

  #### Export ####
  if(Seuil.pcover == 0.2){
    saveRDS(MCWT.clim.PT_ss, "Resultats/ACA/Traits/CWT_ACA_PT_ss_20p.Rds")
    saveRDS(MCWT.clim.PT_ss_gf, "Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_20p.Rds")
    saveRDS(MCWT.clim.PT_sl, "Resultats/ACA/Traits/CWT_ACA_PT_sl_20p.Rds")
    saveRDS(MCWT.clim.PT_sl_gf, "Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_20p.Rds")
    saveRDS(MCWT.clim.MV, "Resultats/ACA/Traits/CWT_ACA_MV_20p.Rds")
    saveRDS(MCWT.clim.MV_gf, "Resultats/ACA/Traits/CWT_ACA_MV_gf_20p.Rds")}
  if(Seuil.pcover == 0.3){
    saveRDS(MCWT.clim.PT_ss, "Resultats/ACA/Traits/CWT_ACA_PT_ss_30p.Rds")
    saveRDS(MCWT.clim.PT_ss_gf, "Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_30p.Rds")
    saveRDS(MCWT.clim.PT_sl, "Resultats/ACA/Traits/CWT_ACA_PT_sl_30p.Rds")
    saveRDS(MCWT.clim.PT_sl_gf, "Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_30p.Rds")
    saveRDS(MCWT.clim.MV, "Resultats/ACA/Traits/CWT_ACA_MV_30p.Rds")
    saveRDS(MCWT.clim.MV_gf, "Resultats/ACA/Traits/CWT_ACA_MV_gf_30p.Rds")}
  if(Seuil.pcover == 0.4){
    saveRDS(MCWT.clim.PT_ss, "Resultats/ACA/Traits/CWT_ACA_PT_ss_40p.Rds")
    saveRDS(MCWT.clim.PT_ss_gf, "Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_40p.Rds")
    saveRDS(MCWT.clim.PT_sl, "Resultats/ACA/Traits/CWT_ACA_PT_sl_40p.Rds")
    saveRDS(MCWT.clim.PT_sl_gf, "Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_40p.Rds")
    saveRDS(MCWT.clim.MV, "Resultats/ACA/Traits/CWT_ACA_MV_40p.Rds")
    saveRDS(MCWT.clim.MV_gf, "Resultats/ACA/Traits/CWT_ACA_MV_gf_40p.Rds")}
  if(Seuil.pcover == 0.5){
    saveRDS(MCWT.clim.PT_ss, "Resultats/ACA/Traits/CWT_ACA_PT_ss_50p.Rds")
    saveRDS(MCWT.clim.PT_ss_gf, "Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_50p.Rds")
    saveRDS(MCWT.clim.PT_sl, "Resultats/ACA/Traits/CWT_ACA_PT_sl_50p.Rds")
    saveRDS(MCWT.clim.PT_sl_gf, "Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_50p.Rds")
    saveRDS(MCWT.clim.MV, "Resultats/ACA/Traits/CWT_ACA_MV_50p.Rds")
    saveRDS(MCWT.clim.MV_gf, "Resultats/ACA/Traits/CWT_ACA_MV_gf_50p.Rds")}
  if(Seuil.pcover == 0.6){
    saveRDS(MCWT.clim.PT_ss, "Resultats/ACA/Traits/CWT_ACA_PT_ss_60p.Rds")
    saveRDS(MCWT.clim.PT_ss_gf, "Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_60p.Rds")
    saveRDS(MCWT.clim.PT_sl, "Resultats/ACA/Traits/CWT_ACA_PT_sl_60p.Rds")
    saveRDS(MCWT.clim.PT_sl_gf, "Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_60p.Rds")
    saveRDS(MCWT.clim.MV, "Resultats/ACA/Traits/CWT_ACA_MV_60p.Rds")
    saveRDS(MCWT.clim.MV_gf, "Resultats/ACA/Traits/CWT_ACA_MV_gf_60p.Rds")}
  if(Seuil.pcover == 0.7){
    saveRDS(MCWT.clim.PT_ss, "Resultats/ACA/Traits/CWT_ACA_PT_ss_70p.Rds")
    saveRDS(MCWT.clim.PT_ss_gf, "Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_70p.Rds")
    saveRDS(MCWT.clim.PT_sl, "Resultats/ACA/Traits/CWT_ACA_PT_sl_70p.Rds")
    saveRDS(MCWT.clim.PT_sl_gf, "Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_70p.Rds")
    saveRDS(MCWT.clim.MV, "Resultats/ACA/Traits/CWT_ACA_MV_70p.Rds")
    saveRDS(MCWT.clim.MV_gf, "Resultats/ACA/Traits/CWT_ACA_MV_gf_70p.Rds")}
  if(Seuil.pcover == 0.8){
    saveRDS(MCWT.clim.PT_ss, "Resultats/ACA/Traits/CWT_ACA_PT_ss_80p.Rds")
    saveRDS(MCWT.clim.PT_ss_gf, "Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_80p.Rds")
    saveRDS(MCWT.clim.PT_sl, "Resultats/ACA/Traits/CWT_ACA_PT_sl_80p.Rds")
    saveRDS(MCWT.clim.PT_sl_gf, "Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_80p.Rds")
    saveRDS(MCWT.clim.MV, "Resultats/ACA/Traits/CWT_ACA_MV_80p.Rds")
    saveRDS(MCWT.clim.MV_gf, "Resultats/ACA/Traits/CWT_ACA_MV_gf_80p.Rds")}
  if(Seuil.pcover == 0.9){
    saveRDS(MCWT.clim.PT_ss, "Resultats/ACA/Traits/CWT_ACA_PT_ss_90p.Rds")
    saveRDS(MCWT.clim.PT_ss_gf, "Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_90p.Rds")
    saveRDS(MCWT.clim.PT_sl, "Resultats/ACA/Traits/CWT_ACA_PT_sl_90p.Rds")
    saveRDS(MCWT.clim.PT_sl_gf, "Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_90p.Rds")
    saveRDS(MCWT.clim.MV, "Resultats/ACA/Traits/CWT_ACA_MV_90p.Rds")
    saveRDS(MCWT.clim.MV_gf, "Resultats/ACA/Traits/CWT_ACA_MV_gf_90p.Rds")}
  if(Seuil.pcover == 1){
    saveRDS(MCWT.clim.PT_ss, "Resultats/ACA/Traits/CWT_ACA_PT_ss_100p.Rds")
    saveRDS(MCWT.clim.PT_ss_gf, "Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_100p.Rds")
    saveRDS(MCWT.clim.PT_sl, "Resultats/ACA/Traits/CWT_ACA_PT_sl_100p.Rds")
    saveRDS(MCWT.clim.PT_sl_gf, "Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_100p.Rds")
    saveRDS(MCWT.clim.MV, "Resultats/ACA/Traits/CWT_ACA_MV_100p.Rds")
    saveRDS(MCWT.clim.MV_gf, "Resultats/ACA/Traits/CWT_ACA_MV_gf_100p.Rds")}
  }

#### Plot for biomes reconstruction ####
Biom.modelling = F
if(Biom.modelling == T){
  #### Import ####
  if(exists("MV.releve") == F){
    header.ACA <- readRDS("Import/ACA/Vegetation/sPlot.ACA.metadata.Rds")
    MV.releve <- readRDS("Resultats/ACA/Vegetation/MV_ACA_plot_only_spermatophyte.Rds")
    }
    
  CWM.selectcov = "PCover60"
  if(CWM.selectcov == "PCover20"){
    MCWT.clim.PT_ss <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_20p.Rds")
    MCWT.clim.PT_ss_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_20p.Rds")
    MCWT.clim.PT_sl <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_20p.Rds")
    MCWT.clim.PT_sl_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_20p.Rds")
    # MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_20p.Rds")
    MCWT.clim.MV_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_gf_20p.Rds")}  
  if(CWM.selectcov == "PCover30"){
    MCWT.clim.PT_ss <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_30p.Rds")
    MCWT.clim.PT_ss_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_30p.Rds")
    MCWT.clim.PT_sl <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_30p.Rds")
    MCWT.clim.PT_sl_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_30p.Rds")
    # MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_30p.Rds")
    MCWT.clim.MV_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_gf_30p.Rds")}  
  if(CWM.selectcov == "PCover40"){
    MCWT.clim.PT_ss <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_40p.Rds")
    MCWT.clim.PT_ss_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_40p.Rds")
    MCWT.clim.PT_sl <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_40p.Rds")
    MCWT.clim.PT_sl_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_40p.Rds")
    # MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_40p.Rds")
    MCWT.clim.MV_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_gf_40p.Rds")}  
  if(CWM.selectcov == "PCover50"){
    MCWT.clim.PT_ss <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_50p.Rds")
    MCWT.clim.PT_ss_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_50p.Rds")
    MCWT.clim.PT_sl <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_50p.Rds")
    MCWT.clim.PT_sl_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_50p.Rds")
    MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_50p.Rds")
    MCWT.clim.MV_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_gf_50p.Rds")}
  if(CWM.selectcov == "PCover60"){
    MCWT.clim.PT_ss <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_60p.Rds")
    MCWT.clim.PT_ss_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_60p.Rds")
    MCWT.clim.PT_sl <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_60p.Rds")
    MCWT.clim.PT_sl_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_60p.Rds")
    MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_60p.Rds")
    MCWT.clim.MV_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_gf_60p.Rds")}
  if(CWM.selectcov == "PCover70"){
    MCWT.clim.PT_ss <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_70p.Rds")
    MCWT.clim.PT_ss_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_70p.Rds")
    MCWT.clim.PT_sl <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_70p.Rds")
    MCWT.clim.PT_sl_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_70p.Rds")
    MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_70p.Rds")
    MCWT.clim.MV_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_gf_70p.Rds")}
  if(CWM.selectcov == "PCover80"){
    MCWT.clim.PT_ss <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_80p.Rds")
    MCWT.clim.PT_ss_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_80p.Rds")
    MCWT.clim.PT_sl <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_80p.Rds")
    MCWT.clim.PT_sl_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_80p.Rds")
    MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_80p.Rds")
    MCWT.clim.MV_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_gf_80p.Rds")}
  if(CWM.selectcov == "PCover90"){
    MCWT.clim.PT_ss <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_90p.Rds")
    MCWT.clim.PT_ss_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_90p.Rds")
    MCWT.clim.PT_sl <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_90p.Rds")
    MCWT.clim.PT_sl_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_90p.Rds")
    # MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_90p.Rds")
    MCWT.clim.MV_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_gf_90p.Rds")}
  if(CWM.selectcov == "PCover100"){
    MCWT.clim.PT_ss <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_100p.Rds")
    MCWT.clim.PT_ss_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_100p.Rds")
    MCWT.clim.PT_sl <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_100p.Rds")
    MCWT.clim.PT_sl_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_100p.Rds")
    # MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_100p.Rds")
    MCWT.clim.MV_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_gf_100p.Rds")}
  MCWT.clim.splot <- data.frame(readRDS("Import/ACA/Vegetation/sPlot_ACA_CWM_clim.Rds"))
  ### !!!! enlever le commentaire en dessous si CWM.selectcov est diff de 60
  # MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_60p.Rds")
  
  #### Table cover CWM ####
  Tab.cover.cwm = F
  if(Tab.cover.cwm == T){
    Stat.PT_ss <- MCWT.clim.PT_ss$MCWT.stat
    Stat.PT_ss_gf <- MCWT.clim.PT_ss_gf$MCWT.stat
    Stat.PT_sl <- MCWT.clim.PT_sl$MCWT.stat
    Stat.PT_sl_gf <- MCWT.clim.PT_sl_gf$MCWT.stat
    Stat.MV <- MCWT.clim.MV$MCWT.stat
    Stat.MV_gf <- MCWT.clim.MV_gf$MCWT.stat
    
    Table.cover.CWM <- data.frame(cbind(mv = as.numeric(Stat.MV$Pour.sites.kept)*100, 
                             mv_gf = as.numeric(Stat.MV_gf$Pour.sites.kept)*100,
                             ss = as.numeric(Stat.PT_ss$Pour.sites.kept)*100, 
                             ss_gf = as.numeric(Stat.PT_ss_gf$Pour.sites.kept)*100, 
                             sl = as.numeric(Stat.PT_sl$Pour.sites.kept)*100, 
                             sl_gf = as.numeric(Stat.PT_sl_gf$Pour.sites.kept)*100 
                             ))
    row.names(Table.cover.CWM) <- gsub("TRY_", "", Stat.PT_ss$Trait)
    colnames(Table.cover.CWM) <- c("ACAV","ACAV (gf)", "ACASP-ss", "ACASP-ss (gf)", "ACASP-sl", "ACASP-sl (gf)")
    Table.cover.CWM <- Table.cover.CWM[c(8,9,2,4,10,1),]
    Table.cover.CWM[["Mean coverage"]] <- round(rowMeans(Table.cover.CWM), digits = 0)
    Table.cover.CWM["Mean coverage",] <- round(colMeans(Table.cover.CWM), digits = 0)
    print(Table.cover.CWM)
    write.csv(Table.cover.CWM, "Resultats/ACA/Traits/Corresp_tax/Table_CWM_cover.csv", row.names = T)}
  
  #### Clean import ####
  MCWT.clim.PT_ss <- MCWT.clim.PT_ss$MCWT
  MCWT.clim.PT_ss_gf <- MCWT.clim.PT_ss_gf$MCWT
  MCWT.clim.PT_sl <- MCWT.clim.PT_sl$MCWT
  MCWT.clim.PT_sl_gf <- MCWT.clim.PT_sl_gf$MCWT
  MCWT.clim.MV <- MCWT.clim.MV$MCWT
  MCWT.clim.MV_gf <- MCWT.clim.MV_gf$MCWT
  
  #### Test Taiga ####
  Test.taiga = F
  if(Test.taiga == T){
    if(exists("MV.releve") == F){MV.releve <- readRDS("Resultats/ACA/Vegetation/MV_ACA_plot_only_spermatophyte.Rds")}
    MV.ACA.Biom <- readRDS("Resultats/ACA/Vegetation/DBACA_Biom.Rds")
    GrowthForm.ACA <- readRDS(paste(DB.path, "Traits/TRY_dec_2019/Extraction/TRY_GrowthForm_ACA.Rds", sep = ""))
    
    Taiga.site <- MCWT.clim.MV[which(MCWT.clim.MV$Biome == "Boreal Forests/Taiga"),]$Site
    MV.releve.taig <- MV.releve[which(row.names(MV.releve) %in% Taiga.site),]
    MV.releve.taig <- MV.releve.taig[colMeans(MV.releve.taig) > 0]
    
    Mean.taiga <- data.frame(Mean = colMeans(MV.releve.taig))
    Mean.taiga$Nsup1 <- NA
    Mean.taiga$Nsup2 <- NA
    Mean.taiga$Nsup3 <- NA
    Mean.taiga$Nsup3 <- NA
    Mean.taiga$species <- row.names(Mean.taiga)
    for(i in 1:ncol(MV.releve.taig)){
      A = length(which(MV.releve.taig[i] >= 1))
      B = length(which(MV.releve.taig[i] >= 0.5))
      C = length(which(MV.releve.taig[i] >= 0.2))
      D = length(which(MV.releve.taig[i] >= 0.05))
      Mean.taiga$Nsup1[i] <- A
      Mean.taiga$Nsup2[i] <- B
      Mean.taiga$Nsup3[i] <- C
      Mean.taiga$Nsup4[i] <- D
      print(A)
    }
    MV.releve.taig.tree <- left_join(Mean.taiga, GrowthForm.ACA, id = "species")
    # MV.releve.taig.tree <- MV.releve.taig.tree[MV.releve.taig.tree$GrowthForm == "Tree",]
    MV.releve.taig.tree <- MV.releve.taig.tree[MV.releve.taig.tree$GrowthForm != "Herb",]
    MV.releve.taig.tree <- MV.releve.taig.tree[MV.releve.taig.tree$GrowthForm != "Unknown",]
    MV.releve.taig.tree <- MV.releve.taig.tree[MV.releve.taig.tree$GrowthForm != "Other",]
    
    Mean.taiga.keep <- MV.releve.taig[MV.releve.taig.tree$species[!is.na(MV.releve.taig.tree$species)]]
    Mean.taiga.keep$havetree <- rowSums(Mean.taiga.keep)
    Mean.taiga.keep <- Mean.taiga.keep[Mean.taiga.keep$havetree>0,]
    Keep.taig <- row.names(Mean.taiga.keep)
    Remove.taig <- setdiff(row.names(MV.releve.taig), row.names(Mean.taiga.keep))
    
    saveRDS(Remove.taig, "Import/ACA/Vegetation/sPlot_remove_fake_taiga.Rds")
    Mean.taiga.rep <- Mean.taiga[Mean.taiga$Mean >= 0.015 & Mean.taiga$Nsup2 >= 30,]
    Mean.taiga <- left_join(Mean.taiga, GrowthForm.ACA, id = "species")
    Mean.taiga <- Mean.taiga[order(Mean.taiga$GrowthForm),]
    
    values.gf = c("Herb" = "#b5ab32ff",
                  "Shrub" = "#aa373aff",
                  "Other" = "grey90",
                  "Unknown" = "grey90",
                  "Tree" = "#0f6b31ff")
    
    library(ggrepel)
    ptaig <- ggplot()+ 
      # geom_point(data = Mean.taiga, aes(x = Mean, y = Nsup1), alpha = 0.4, color = "grey10")+
      # geom_point(data = Mean.taiga, aes(x = Mean, y = Nsup2), color = "darkblue", alpha = 0.4)+
      geom_point(data = Mean.taiga, aes(x = Mean, y = Nsup2, color = GrowthForm), alpha = 0.8)+
      xlab("Mean taxa abundance")+ ylab("Nb of plot with taxa abundance > 50%")+
      scale_color_manual(values = values.gf, name = "Growth forms", drop = T)+
      annotate("text", x = 0.03, y = 0.5, label = paste("n_site = ", length(Taiga.site), "\nn_taxa = ", length(unique(Mean.taiga$species))), hjust = 0)+
      geom_text_repel(data = Mean.taiga.rep, aes(x = Mean, y = Nsup2, label = species), size = 3)+#, nudge_y = 3, nudge_x = 0.005) +
      theme(axis.line = element_line(colour = "grey30"), 
            plot.background = element_blank(), panel.background = element_blank(),
            panel.grid = element_blank(), 
            legend.key = element_blank(),
            panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 1.5)
            
            )
    # geom_point(aes(x = Mean, y = Nsup3), color = "royalblue", alpha = 0.6)#+
    # geom_point(aes(x = Mean, y = Nsup4), color = "lightblue", alpha = 0.8)
    W = 700
    H = 600
    ggsave(filename = "Figures/ACA/Trait/Dispersion/Taiga_CWM_mean.pdf", ptaig, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    MV.coord.taig <- MV.ACA.Biom[which(row.names(MV.ACA.Biom) %in% Taiga.site),]
    MV.releve.taig.maj <- MV.releve.taig[colSums(MV.releve.taig) >1]
    MV.releve.taig.maj <- t(MV.releve.taig.maj)
    names(MV.releve.taig.maj) <- row.names(MV.coord.taig)
    MV.releve.taig.maj$gender <- row.names(MV.releve.taig.maj)
    MV.releve.taig.maj$gender <- gsub("\\s.*", "", MV.releve.taig.maj$gender)
    
    #### Map points taiga ####
    Save.plot = "Figures/ACA/Trait/Map/Taiga_samples.pdf"
    pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)
    map("world", regions = c("Russia"))
    points(MV.coord.taig$Longitude, MV.coord.taig$Latitude, col = "darkred")
    # sort(colSums(MV.releve.taig), decreasing = T)[1:100]
    dev.off()
    # Faba.hei <- MT.ACA[MT.ACA$family == "Fabaceae", c("species", "Height")]
  }
  
  Test.Gobi = F
  if(Test.Gobi == T){
    set.seed(12345)       
    Remove.gobi <- read.csv("Import/ACA/Vegetation/sPlot_Riparian_Gobi.csv")
    Remove.gobi <- Remove.gobi$Site
    MV.releve.gobi <- MV.releve[match(Remove.gobi,row.names(MV.releve)),]  
    MV.releve.gobi <- MV.releve.gobi[,match(names(which(colSums(MV.releve.gobi)>0)), names(MV.releve.gobi))]
    header.Gobi <- header.ACA[match(Remove.gobi,header.ACA$ID),]  
    }
  
  #### Remove points pourris ####
  Remove.taiga.wetland = T
  if(Remove.taiga.wetland == T){
    Remove.taig <- readRDS("Import/ACA/Vegetation/sPlot_remove_fake_taiga.Rds")
    MCWT.clim.MV <- MCWT.clim.MV[!MCWT.clim.MV$Site %in% Remove.taig,]
    MCWT.clim.MV_gf <- MCWT.clim.MV_gf[!MCWT.clim.MV_gf$Site %in% Remove.taig,]
    MCWT.clim.splot <- MCWT.clim.splot[!MCWT.clim.splot$Site %in% Remove.taig,]
    }
  
  Remove.ladakh = F
  if(Remove.ladakh == T){
    set.seed(12345)       
    Remove.lada <- read.csv("Import/ACA/Vegetation/sPlot_Ladakh.csv")
    
    Remove.lada <- Remove.lada$Site
    Keep.lada.10p <- sample(Remove.lada, 100)
    Remove.lada <- setdiff(Remove.lada, Keep.lada.10p)
    
    MCWT.clim.MV <- MCWT.clim.MV[!MCWT.clim.MV$Site %in% Remove.lada,]
    MCWT.clim.MV_gf <- MCWT.clim.MV_gf[!MCWT.clim.MV_gf$Site %in% Remove.lada,]
    MCWT.clim.splot <- MCWT.clim.splot[!MCWT.clim.splot$Site %in% Remove.lada,]
  }
    
  Remove.local.forest = T
  if(Remove.local.forest == T){
    header.ACA.fo <- header.ACA[c(1,36)]
    names(header.ACA.fo)[1] <- "Site"
    MCWT.clim.MV.clean <- left_join(MCWT.clim.MV, header.ACA.fo, by = "Site")
    Row.to.remove <- which(MCWT.clim.MV.clean$is.forest == T & MCWT.clim.MV.clean$Biome.no %in% c(8, 10, 13))
    Row.to.remove2 <- which(MCWT.clim.MV.clean$is.forest == F & MCWT.clim.MV.clean$Biome.no %in% c(4, 5, 6))
    # Row.to.remove <- NULL
    Row.to.keep <- setdiff(seq(1:nrow(MCWT.clim.MV.clean)), c(Row.to.remove, Row.to.remove2))
    Site.to.keep <- MCWT.clim.MV.clean[Row.to.keep,"Site"]
    
    MCWT.clim.MV <- MCWT.clim.MV[MCWT.clim.MV$Site %in% Site.to.keep,]
    MCWT.clim.MV_gf <- MCWT.clim.MV_gf[MCWT.clim.MV_gf$Site %in% Site.to.keep,]
    MCWT.clim.splot <- MCWT.clim.splot[MCWT.clim.splot$Site %in% Site.to.keep,]
    }
  
  #### Test taille DB #### 
  Test.biome.rep = F
  if(Test.biome.rep == T){
    # MCWT.clim.MV_gf <- na.omit(MCWT.clim.MV_gf)
    Nbiome.MV <- round(table(MCWT.clim.PT_ss_gf$Biome)/length(MCWT.clim.PT_ss_gf$Biome), digits = 2)
    Nbiome.MP <-round(table(MCWT.clim.MV_gf$Biome)/length(MCWT.clim.MV_gf$Biome), digits = 2)
    
    # Nbiome.MV <- table(MCWT.clim.MV_gf$Biome)
    # Nbiome.MP <- table(MCWT.clim.PT_ss_gf$Biome)
    
    set.seed(12345)       
    Keep.MV <- data.frame()
    for(i in names(Nbiome.MV)){
      Data.biome <- MCWT.clim.MV_gf[MCWT.clim.MV_gf$Biome == i,]
      # print(Data.biome)
      print(Nbiome.MP[[i]]*nrow(MCWT.clim.MV_gf))
      # data_s1 <- Data.biome[sample(1:nrow(Data.biome), Nbiome.MP[[i]]*nrow(MCWT.clim.MV_gf)), ]
      
      # print(Nbiome.MP[[i]])
      # print(round(Nbiome.MP[[i]]*nrow(MCWT.clim.MV_gf), digits = 0))
      # Keep.MV <- rbind(Keep.MV, data_s1)
    }
    # Nbiome.MV2 <- table(Keep.MV$Biome)
    Nbiome.MV2 <- round(table(Keep.MV$Biome)/length(Keep.MV$Biome), digits = 2)
    
    MCWT.clim.MV_gf <- Keep.MV
  }
  
  #### Export clean data pour plot ####
  saveRDS(MCWT.clim.PT_ss, "Resultats/ACA/Traits/CWT_ACA_PT_ss_clean.Rds")
  saveRDS(MCWT.clim.PT_ss_gf, "Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_clean.Rds")
  saveRDS(MCWT.clim.PT_sl, "Resultats/ACA/Traits/CWT_ACA_PT_sl_clean.Rds")
  saveRDS(MCWT.clim.PT_sl_gf, "Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_clean.Rds")
  saveRDS(MCWT.clim.MV, "Resultats/ACA/Traits/CWT_ACA_MV_clean.Rds")
  saveRDS(MCWT.clim.MV_gf, "Resultats/ACA/Traits/CWT_ACA_MV_gf_clean.Rds")
  saveRDS(MCWT.clim.splot, "Resultats/ACA/Traits/sPlot_ACA_CWM_clim_clean.Rds")
  write.csv(MCWT.clim.MV_gf, "Resultats/ACA/Traits/CWT_ACA_MV_gf_clean.csv")
  
  #### Show number of sites by sampling ####
  Number.site.by.biom = F
  if(Number.site.by.biom == T){
    Area.biomes <- read.table("Import/ACA/Site/ACA_area/ACA_stats.csv", sep = ",", header = T)
    Area.biomes <- Area.biomes[which(Area.biomes$Biome.no %in% c(4,5,6,8,10,13)),]
    Area.biomes$Area.pour <- round(Area.biomes$Area / sum(Area.biomes$Area) * 100, digits = 2)
    print(Area.biomes)
    
    Site.clust <- function(Mveg, Mpol, Category, Save.plot, H, W){
      #### Settings ####
      if(missing(Save.plot)){Save.plot = NULL}
      if(missing(W)){W = NULL}
      if(missing(H)){H = NULL}
      if(missing(Mveg)){Mveg = NULL}
      if(missing(Mpol)){Mpol = NULL}
      
      #### Save plots ####
      if(is.null(Save.plot) == F){
        Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
        dir.create(file.path(Path.to.create), showWarnings = FALSE)
        if(is.null(W) == F & is.null(H) == F){
          pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
        else{pdf(file = Save.plot)}}
      
      #### Calculation ####
      if(is.null(Mveg) == F){
        Mveg$Type <- "Vegetation"
        M <- Mveg}
      if(is.null(Mpol) == F){
        Mpol$Type <- "Pollen"
        M <- Mpol}
      if(is.null(Mpol) == T & is.null(Mveg) == T){
        print("**** WARNING: a matrix is needed ****")
        return(NULL)
        }
      if(is.null(Mpol) == F & is.null(Mveg) == F){
        Mveg <- Mveg[Category]
        Mveg$Type <- "Vegetation"
        
        Mpol <- Mpol[Category]
        Mpol$Type <- "Pollen"
        
        M <- rbind(Mveg, Mpol)
        # print(Mveg)
      
        }
      
      #### Color and fill settings ####
      Value.bi <- c(
        "Boreal Forests/Taiga" = "#8FB8E6",
        "Deserts & Xeric Shrublands" = "#C88282",
        "Montane Grasslands & Shrublands" = "#D0C3A7",
        "Temperate Broadleaf & Mixed Forests" = "#3E8A70",
        "Temperate Conifer Forests" = "#6B9A88",
        "Temperate Grasslands, Savannas & Shrublands" = "#ECED8A",
        "N/A" = "#FFEAAF",
        "Tundra" = "#A9D1C2",
        "Tropical & Subtropical Coniferous Forests" = "#99CA81",
        "Mangroves" = "#FE01C4",
        "Flooded Grasslands & Savannas" = "#BEE7FF")
       
      Value.bi <- Value.bi[which(names(Value.bi) %in% unique(M$Biome))]
      My_fill <-  scale_fill_manual(values = Value.bi, name = "") # "Biomes (Dinerstein et al., 2017)"#, labels = levels(ACA.biom.proj$BIOME_NAME)
      
      #### Plot ####
      # p <- ggplot(M, aes(x = Biome, fill = Biome, color = Type))+ 
      p <- ggplot(M, aes(x = Biome, fill = Type))+ 
        geom_histogram(stat = "count", position = "dodge", alpha = 1)+
        # My_fill +
        scale_fill_manual(values = c(Vegetation = "#3f6b3dff", Pollen = "#d9a15eff"), "Sampling type")+
        stat_count(aes(y=..count..,label=..count..), geom="text", vjust=-1, position = "dodge") +
        ylab("Number of sites")+ xlab("Biomes")+
        #### Theme ####
      theme(
        panel.background = element_blank(),
        legend.position = "top", panel.grid = element_blank(),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1),
        # legend.text = element_text(size = Legend.size), 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.border = element_rect(fill = NA),
        strip.text = element_text(size = 12), 
        strip.background = element_blank()
      )
      
      #### Export ####
      print(p)
      if(is.null(Save.plot) == F){dev.off()}
      
    }
    
    Site.clust(Mveg = MCWT.clim.MV_gf, Mpol = MCWT.clim.PT_sl_gf, Category = "Biome",
               H = 500, W = 700, Save.plot = "Figures/ACA/Trait/Distribution/Site_by_biomes.pdf")
  }
  #### Correl matrice 10 traits ####
  Correl.mat.10t = F
  if(Correl.mat.10t == T){
    # Trait <- c(4,5,7,11,12,13)
    Trait <- c(4:16)
    Clim.chel <- c(2,3,15,16,17,19,21,23,25,27)
    Clim.wc <- c(2,3,15,16,17,18,20,22,24,26,28)
    Clim.both <- c(Clim.chel, Clim.wc)
    Clim <- Clim.chel
    # Clim <- Clim.wc
    
    H = 700
    W = 800
    Save.matcor.CWM = "Figures/ACA/Trait/Correlation/Mat_correlation/CWM/Matcor_clim_CWT-10t-merge_wc.pdf"
    pdf(file = Save.matcor.CWM, width = W*0.01041666666667, height = H*0.01041666666667)
    par(mfrow = c(2,3))
    
    MC.ss <- Mat.corel.CWT.clim(MCWT.clim.PT_ss[Clim], MCWT.clim.PT_ss[Trait],
                            I.confiance = 0.95, Display.pval = "blank", Disp.R = "number",
                            Title = "Correlation between CWT-ACASP-ss and bioclimate parameters.",
                            Save.path = "Resultats/ACA/Traits/Matcorr_CWT_ss.csv",
                            )
    
    MC.ss.gf <- Mat.corel.CWT.clim(MCWT.clim.PT_ss_gf[Clim], MCWT.clim.PT_ss_gf[Trait],
                                I.confiance = 0.95, Display.pval = "blank", Disp.R = "number",
                                Title = "Correlation between CWT-ACASP-ss-gf and bioclimate parameters.",
                                Save.path = "Resultats/ACA/Traits/Matcorr_CWT_ss_gf.csv",
                                )
    
    MC.sl <- Mat.corel.CWT.clim(MCWT.clim.PT_sl[Clim], MCWT.clim.PT_sl[Trait],
                                I.confiance = 0.95, Display.pval = "blank", Disp.R = "number",
                                Title = "Correlation between CWT-ACASP-sl and bioclimate parameters.",
                                Save.path = "Resultats/ACA/Traits/Matcorr_CWT_sl.csv",
                                )
    
    MC.sl.gf <- Mat.corel.CWT.clim(MCWT.clim.PT_sl_gf[Clim], MCWT.clim.PT_sl_gf[Trait],
                                   I.confiance = 0.95, Display.pval = "blank", Disp.R = "number",
                                   Title = "Correlation between CWT-ACASP-sl-gf and bioclimate parameters.",
                                   Save.path = "Resultats/ACA/Traits/Matcorr_CWT_sl_gf.csv",
                                   )
    
    MC.mv <- Mat.corel.CWT.clim(MCWT.clim.MV[Clim], MCWT.clim.MV[Trait],
                                I.confiance = 0.95, Display.pval = "blank", Disp.R = "number",
                                Title = "Correlation between CWT-ACAV and bioclimate parameters.",
                                Save.path = "Resultats/ACA/Traits/Matcorr_CWT_mv.csv")
    
    MC.mv.gf <- Mat.corel.CWT.clim(MCWT.clim.MV_gf[Clim], MCWT.clim.MV_gf[Trait],
                                   I.confiance = 0.95, Display.pval = "blank", Disp.R = "number",
                                   Title = "Correlation between CWT-ACAV-gf and bioclimate parameters.",
                                   Save.path = "Resultats/ACA/Traits/Matcorr_CWT_mv_gf.csv")
    
    # MC.sP.gf <- Mat.corel.CWT.clim(MCWT.clim.splot[c(2,15:27)], MCWT.clim.splot[c(4:13)],
    #                                I.confiance = 0.95,
    #                                Display.pval = "blank",
    #                                Disp.R = "number",
    #                                Title = "Correlation between CWT-ACAV-gf and bioclimate parameters.",
    #                                Save.path = "Resultats/ACA/Traits/Matcorr_CWT_mv_gf.csv",
    #                                Save.plot = "Figures/ACA/Trait/Correlation/Mat_correlation/CWM/Matcor_clim_CWT_splot_gf-10t.pdf",
    #                                H = 1200,
    #                                W = 800)
    dev.off()
    
    
    }
  
  #### Correl matrice all data ####
  Correl.mat.full = F
  if(Correl.mat.full == T){
    M1 <- MCWT.clim.PT_ss
    names(M1)[grep("TRY", names(M1))] <- paste(names(M1)[grep("TRY", names(M1))], "ss", sep = "_")
    M2 <- MCWT.clim.PT_ss_gf[,grep("TRY", names(MCWT.clim.PT_ss_gf))]
    names(M2) <- paste(names(M2), "ss_gf", sep = "_")
    M3 <- MCWT.clim.PT_sl[,grep("TRY", names(MCWT.clim.PT_sl))]
    names(M3) <- paste(names(M3), "sl", sep = "_")
    M4 <- MCWT.clim.PT_sl_gf[,grep("TRY", names(MCWT.clim.PT_sl_gf))]
    names(M4) <- paste(names(M4), "sl_gf", sep = "_")
    
    M5 <- MCWT.clim.MV
    # M5 <- MCWT.clim.MV[,grep("TRY", names(MCWT.clim.MV))]
    names(M5)[grep("TRY", names(M5))] <- paste(names(M5)[grep("TRY", names(M5))], "ACAV", sep = "_")
    # names(M5) <- paste(names(M5), "ACAV", sep = "_")
    M6 <- MCWT.clim.MV_gf[,grep("TRY", names(MCWT.clim.MV_gf))]
    names(M6) <- paste(names(M6), "ACAV_gf", sep = "_")
    
    M.cor.full <- cbind(M1, M2, M3, M4)
    M.cor.full <- M.cor.full[order(names(M.cor.full))]
    MCWT.clim.ACAP <- M.cor.full[c(22, 1:21, 23:ncol(M.cor.full))]
    
    M7 <- MCWT.clim.splot[,grep("TRY", names(MCWT.clim.splot))]
    names(M7) <- paste(names(M7), "splot", sep = "_")
    
    MCWT.clim.ACAP.full <- cbind(M5, M6)
    MCWT.clim.ACAP.full <- MCWT.clim.ACAP.full[order(names(MCWT.clim.ACAP.full))]
    MCWT.clim.ACAP.full <- MCWT.clim.ACAP.full[c(22, 1:21, 23:ncol(MCWT.clim.ACAP.full))]
    
    saveRDS(MCWT.clim.ACAP, "Resultats/ACA/Traits/CWT_ACA_PT_all.Rds")
    saveRDS(MCWT.clim.ACAP.full, "Resultats/ACA/Traits/CWT_ACA_ACAV_all.Rds")
    
    # Clim.wc <- c(2,3,9,10,12,14,16,18,20,22)
    # Clim.full <- c(2,3,9:22)
    Clim.full <- c(2,3,9:22)
    
    MC.full <- Mat.corel.CWT.clim(MCWT.clim.ACAP[c(27:30,35:42,51:62)], MCWT.clim.ACAP[Clim.full], # full M.cor.full[c(22:65)], M.cor.full[c(1,7:20)], 
                                I.confiance = 0.95, 
                                Display.pval = "blank",
                                Disp.R = "number",
                                Title = "Correlation between CWT-ACASP (ss / sl / gf) and bioclimate parameters.",
                                Save.path = "Resultats/ACA/Traits/Matcorr_CWT_full.csv",
                                Save.plot = "Figures/ACA/Trait/Correlation/Mat_correlation/CWM/Matcor_clim_CWT_full_ACAP.pdf",
                                H = 1900,
                                W = 1000
                                )
    
    # Clim.wc <- c(2,3,9,10,12,14,16,18,20,22)
    # Clim.full <- c(2,3,9:22)
    
    MC.ACAV <- Mat.corel.CWT.clim(MCWT.clim.ACAP.full[c(25,26,29,30,31,32,37,38,39,40,41,42)], MCWT.clim.ACAP.full[Clim.full],
                                  I.confiance = 0.95,
                                  Display.pval = "blank",
                                  Disp.R = "number",
                                  Title = "Correlation between CWT-ACAV and bioclimate parameters.",
                                  Save.path = "Resultats/ACA/Traits/Matcorr_CWT_ACAV.csv",
                                  Save.plot = "Figures/ACA/Trait/Correlation/Mat_correlation/CWM/Matcor_clim_CWT_ACAV.pdf",
                                  H = 1100,
                                  W = 800
    )
    # print(MC.ACAV)
    }
    
  #### Correl matrice 6 traits + wc ####
  Correl.mat = F
  if(Correl.mat == T){
    #### Select traits ####
    Trait <- c(11,12,5,7,13,4)
    Clim.chel <- c(2,3,15,16,17,19,21,23,25,27)
    Clim.wc <- c(2,3,15,16,18,20,22,24,26,28)
    Clim.both <- c(Clim.chel, Clim.wc)
    # Clim <- Clim.chel
    Clim <- Clim.wc
    
    #### Matrice correl ####
    H = 700
    W = 800
    Save.matcor.CWM = "Figures/ACA/Trait/Correlation/Mat_correlation/CWM/Matcor_clim_CWT-6t-merge_wc.pdf"
    pdf(file = Save.matcor.CWM, width = W*0.01041666666667, height = H*0.01041666666667)
    par(mfrow = c(2,3))
    # names(MCWT.clim.MV_gf) <- gsub("_chel", "", names(MCWT.clim.MV_gf))
    # names(MCWT.clim.PT_ss_gf) <- gsub("_chel", "", names(MCWT.clim.PT_ss_gf))
    # names(MCWT.clim.PT_sl_gf) <- gsub("_chel", "", names(MCWT.clim.PT_sl_gf))
    names(MCWT.clim.MV_gf) <- gsub("_wc", "", names(MCWT.clim.MV_gf))
    names(MCWT.clim.PT_ss_gf) <- gsub("_wc", "", names(MCWT.clim.PT_ss_gf))
    names(MCWT.clim.PT_sl_gf) <- gsub("_wc", "", names(MCWT.clim.PT_sl_gf))
    my_color <- colorRampPalette(c("royalblue", "royalblue", "grey95", "grey95", "darkorange", "darkorange"))(20)
    print(names(MCWT.clim.MV_gf))
    
    ##### "everything", "all.obs", "complete.obs", "na.or.complete",  "pairwise.complete.obs"
    My_use.cor <- "pairwise.complete.obs"
    #### "pearson", "kendall", "spearman"
    My_method.cor <- "pearson"
    
    MC.mv.gf <- Mat.corel.CWT.clim(MCWT.clim.MV_gf[Clim], MCWT.clim.MV_gf[Trait], 
                                   I.confiance = 0.95,
                                   Display.pval = "pch", my_color = my_color, 
                                   Disp.R = "number", Label = F, Average = F,  Bar.pos = "n",
                                   Use.cor = My_use.cor, Method.cor = My_method.cor, #Xlab.pos = "d",
                                   # Title = "Correlation between CWT-ACAV-gf and bioclimate parameters.",
                                   # Save.path = "Resultats/ACA/Traits/Matcorr_CWT_mv_gf.csv",
                                   # Save.plot = "Figures/ACA/Trait/Correlation/Mat_correlation/CWM/Matcor_clim_CWT_MV_gf-6t.pdf",
                                   H = 800,
                                   W = 500)
    
    MC.ss.gf <- Mat.corel.CWT.clim(MCWT.clim.PT_ss_gf[Clim], MCWT.clim.PT_ss_gf[Trait],
                                   I.confiance = 0.95, return.pick = F,
                                   my_color = my_color, Display.pval = "pch", #blank
                                   Disp.R = "number", Label = F, Average = F,  Bar.pos = "n",
                                   Use.cor = My_use.cor, Method.cor = My_method.cor,# Xlab.pos = "l",
                                   # Title = "Correlation between CWT-ACASP-ss-gf and bioclimate parameters.",
                                   # Save.path = "Resultats/ACA/Traits/Matcorr_CWT_ss_gf.csv",
                                   # Save.plot = "Figures/ACA/Trait/Correlation/Mat_correlation/CWM/Matcor_clim_CWT_ss_gf-6t.pdf",
                                   H = 800,
                                   W = 500)
    
    MC.sl.gf <- Mat.corel.CWT.clim(MCWT.clim.PT_sl_gf[Clim], MCWT.clim.PT_sl_gf[Trait], 
                                   I.confiance = 0.95,
                                   Display.pval = "pch", my_color = my_color,
                                   Disp.R = "number", Label = F, Average = F,  Bar.pos = "n",
                                   Use.cor = My_use.cor, Method.cor = My_method.cor,# Xlab.pos = "l",
                                   # Title = "Correlation between CWT-ACASP-sl-gf and bioclimate parameters.",
                                   # Save.path = "Resultats/ACA/Traits/Matcorr_CWT_sl_gf.csv",
                                   # Save.plot = "Figures/ACA/Trait/Correlation/Mat_correlation/CWM/Matcor_clim_CWT_sl_gf-6t.pdf",
                                   H = 800,
                                   W = 500)

    plot.new()
    Diff.mat <- round(MC.mv.gf$R2 - MC.ss.gf$R2, digits = 2)
    CP1 <- corrplot(Diff.mat, col = my_color, #tl.pos = "n",
                   tl.col="black", tl.srt=45, tl.cex = .7, method = "number", cl.align.text = "l", cl.pos = "n",
                   sig.level = 0.05, insig = "pch", pch.cex = 2)

    Diff.mat <- round(MC.mv.gf$R2 - MC.sl.gf$R2, digits = 2) #tl.pos = "l",
    CP2 <- corrplot(Diff.mat, col = my_color, tl.col="black", tl.srt=45,
                    tl.cex = .7, method = "number", cl.align.text = "l", cl.pos = "n",
                   sig.level = 0.05, insig = "pch", pch.cex = 2)
    dev.off()

    #### r compar ####
    H = 400
    W = 700
    Save.matcor.trait = "Figures/ACA/Trait/Correlation/R2_compar/CWM/R2_comp_all_clim.pdf"
    
    p1 <- R2.compar(MV = MC.mv.gf, MP = MC.ss.gf, Repel.outliers = F, Show.Plotly = T, Bisectrice.area = 10, Leg.pos = "none",
                    Xlab = "r (CWM-ACASP-fi)", Ylab = "r (CWM-ACAV)", Nb.labs = 4, Panel.lab = "(A)",
                    H = 500, W = 500, Save.plot = "Figures/ACA/Trait/Correlation/R2_compar/CWM/R2_comp_ACAV_clim.pdf"
                    )
    p2 <- R2.compar(MV = MC.mv.gf, MP = MC.sl.gf, Bisectrice.area = 10, Repel.outliers = T,  #Show.Plotly = T,
                    Xlab = "r (CWM-ACASP-co)", Ylab = "r (CWM-ACAV)", Nb.labs = 4, Panel.lab = "(B)")
    # p3 <- R2.compar(MV = MC.ss.gf, MP = MC.sl.gf, Bisectrice.area = 10,
    #                 Xlab = "r (CWM-ACASP-co)", Ylab = "r (CWM-ACASP-fi)", Nb.labs = 4,  Show.Plotly = T,
    #                 H = 500, W = 500, Save.plot = "Figures/ACA/Trait/Correlation/R2_compar/CWM/R2_comp_ACASP_clim.pdf")
    pall <- p1 + p2 #+ p3 #+ plot_layout(guides = "collect")
    # pall <- p1 + p2 + plot_layout(guides = "collect")
    ggsave(file = Save.matcor.trait, pall, width = W*0.01041666666667, height = H*0.01041666666667)

    Test.Pcover = F
    if(Test.Pcover == T){
      H = W = 500
      Save.R2compPcov = "Figures/ACA/Trait/Correlation/R2_compar/CWM/R2_evol_Pcover.pdf"
      Mat.R2.Pcov <- data.frame(read.csv("Import/ACA/Traits/R2_comp_Pcover.csv"))
      Mat.R2.Pcov <- melt(Mat.R2.Pcov, id = c("P_cover", "P_site_ACAV"))
      Mat.R2.Pcov$variable <- as.character(Mat.R2.Pcov$variable)
      Mat.R2.Pcov$variable[grep("_co", Mat.R2.Pcov$variable)] <- "ACAV vs. ACASP-co"
      Mat.R2.Pcov$variable[grep("_fi", Mat.R2.Pcov$variable)] <- "ACAV vs. ACASP-fi"
      coeff = 2
      
      p <- ggplot(Mat.R2.Pcov, aes(x = P_cover))+ 
        geom_vline(xintercept = 50, linewidth = .8, color = "grey60", linetype = "dashed")+
        xlab("P[cover] threshold (%)") +
        geom_line(aes(y = value, shape = variable), color = "darkorange") +
        geom_point(aes(y = value, shape = variable), color = "darkorange", size = 3) +
        geom_line(aes(y = P_site_ACAV / coeff), color = "royalblue") +
        geom_point(aes(y = P_site_ACAV / coeff), color = "royalblue") + 
        scale_x_continuous(n.breaks = 10)+
        scale_y_continuous(
          name = "R²",
          sec.axis = sec_axis(~.*coeff, name="Average %[sPlot-sites used in CWM]")) +
        
        theme(axis.line = element_line(colour = "grey30"), 
              plot.background = element_blank(), panel.background = element_blank(),
              panel.grid = element_blank(), legend.position = "top",
              axis.title.y = element_text(color = "darkorange"),
              legend.key = element_blank(), axis.title.y.right = element_text(color = "royalblue"),
              panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 1.5))

      ggsave(file = Save.R2compPcov, p, width = W*0.01041666666667, height = H*0.01041666666667)
    }
  }
  
  #### Correl CLIM / CLIM ####
  Correl.mat.t.t = F
  if(Correl.mat.t.t == T){
    #### Trait select ####
    Trait <- c(11,12,5,7,13,4)
    Clim.chel <- c(2,3,15,16,17,19,21,23,25,27)
    Clim.wc <- c(2,3,15,16,17,18,20,22,24,26,28)
    Clim.both <- c(Clim.chel, Clim.wc)
    Clim <- Clim.chel
    Trait <- Clim.chel
    # Clim <- Clim.wc
    
    H = 850
    W = 950
    Save.matcor.CWM = "Figures/ACA/Trait/Correlation/Mat_correlation/CWM/Matcor_clim_CWT-6t-merge_clim_clim.pdf"
    pdf(file = Save.matcor.CWM, width = W*0.01041666666667, height = H*0.01041666666667)
    par(mfrow = c(2,2))
    names(MCWT.clim.MV_gf) <- gsub("_chel", "", names(MCWT.clim.MV_gf))
    names(MCWT.clim.PT_ss_gf) <- gsub("_chel", "", names(MCWT.clim.PT_ss_gf))
    names(MCWT.clim.PT_sl_gf) <- gsub("_chel", "", names(MCWT.clim.PT_sl_gf))
    names(MCWT.clim.MV_gf) <- gsub("_wc", "", names(MCWT.clim.MV_gf))
    names(MCWT.clim.PT_ss_gf) <- gsub("_wc", "", names(MCWT.clim.PT_ss_gf))
    names(MCWT.clim.PT_sl_gf) <- gsub("_wc", "", names(MCWT.clim.PT_sl_gf))
    my_color <- colorRampPalette(c("royalblue", "royalblue", "grey95", "grey95", "darkorange", "darkorange"))(20)
    
    names(MCWT.clim.MV_gf) <- gsub("TRY_", "", names(MCWT.clim.MV_gf))
    names(MCWT.clim.PT_ss_gf) <- gsub("TRY_", "", names(MCWT.clim.PT_ss_gf))
    names(MCWT.clim.PT_sl_gf) <- gsub("TRY_", "", names(MCWT.clim.PT_sl_gf))
    
    #### Mat correl ####
    MC.mv.gf <- Mat.corel.CWT.clim(MCWT.clim.MV_gf[Trait], MCWT.clim.MV_gf[Trait], 
                                   I.confiance = 0.95,
                                   Display.pval = "pch", my_color = my_color, Display = "lower",
                                   Disp.R = "number", Label = F, Average = F,  Bar.pos = "n")
    
    MC.ss.gf <- Mat.corel.CWT.clim(MCWT.clim.PT_ss_gf[Trait], MCWT.clim.PT_ss_gf[Trait],
                                   I.confiance = 0.95, return.pick = F,
                                   my_color = my_color, Display.pval = "pch", Display = "lower",
                                   Disp.R = "number", Label = F, Average = F,  Bar.pos = "n")
    
    MC.sl.gf <- Mat.corel.CWT.clim(MCWT.clim.PT_sl_gf[Trait], MCWT.clim.PT_sl_gf[Trait], 
                                   I.confiance = 0.95,
                                   Display.pval = "pch", my_color = my_color, Display = "lower",
                                   Disp.R = "number", Label = F, Average = F,  Bar.pos = "n")

    # Diff.mat <- round(abs(MC.mv.gf$R2 - MC.ss.gf$R2), digits = 2)
    Diff.mat <- round(MC.mv.gf$R2 - MC.ss.gf$R2, digits = 2)
    CP <- corrplot(Diff.mat, type = "lower", col = my_color,
                   tl.col="black", tl.srt=45, tl.cex = .7, method = "number", cl.align.text = "l", cl.pos = "n", 
                   sig.level = 0.05, insig = "pch", pch.cex = 2)
    dev.off()
    
    #### R2 compar ####
    H = 400
    W = 700
    Save.matcor.trait = "Figures/ACA/Trait/Correlation/R2_compar/CWM/R2_comp_onlyclim.pdf"
    
    p1 <- R2.compar(MV = MC.mv.gf, MP = MC.ss.gf, Repel.outliers = F, shade.areas = T,
                    Xlab = "r (CWM-ACASP-fi)", Ylab = "r (CWM-ACAV)", Nb.labs = 0, Panel.lab = "(A)")
    p2 <- R2.compar(MV = MC.mv.gf, MP = MC.sl.gf, Leg.pos = "none", shade.areas = T,
                    Xlab = "r (CWM-ACASP-co)", Ylab = "r (CWM-ACAV)", Nb.labs = 6, Panel.lab = "(B)")
    # p3 <- R2.compar(MV = MC.ss.gf, MP = MC.sl.gf,
    #                 Xlab = "r (CWM-ACASP-co)", Ylab = "r (CWM-ACASP-fi)", Nb.labs = 4,
    #                 H = 500, W = 500, Save.plot = "Figures/ACA/Trait/Correlation/R2_compar/CWM/R2_comp_ACASP.pdf")
    pall <- p1 + p2 #+ p3 #+ plot_layout(guides = "collect")
    ggsave(file = Save.matcor.trait, pall, width = W*0.01041666666667, height = H*0.01041666666667)
    
    

  }
  
  #### Correl matrice 6 traits / traits ####
  Correl.mat.t.t = F
  if(Correl.mat.t.t == T){
    #### Trait select ####
    Trait <- c(11,12,5,7,13,4)
    Clim.chel <- c(2,3,15,16,17,19,21,23,25,27)
    Clim.wc <- c(2,3,15,16,17,18,20,22,24,26,28)
    Clim.both <- c(Clim.chel, Clim.wc)
    Clim <- Trait
    Trait <- Trait
    # Clim <- Clim.wc
    
    #### Mat correl ####
    H = 650
    W = 700
    Save.matcor.CWM = "Figures/ACA/Trait/Correlation/Mat_correlation/CWM/Matcor_clim_CWT-6t-trait-trait.pdf"
    pdf(file = Save.matcor.CWM, width = W*0.01041666666667, height = H*0.01041666666667)
    par(mfrow = c(2,2))
    names(MCWT.clim.MV_gf) <- gsub("_chel", "", names(MCWT.clim.MV_gf))
    names(MCWT.clim.PT_ss_gf) <- gsub("_chel", "", names(MCWT.clim.PT_ss_gf))
    names(MCWT.clim.PT_sl_gf) <- gsub("_chel", "", names(MCWT.clim.PT_sl_gf))
    names(MCWT.clim.MV_gf) <- gsub("_wc", "", names(MCWT.clim.MV_gf))
    names(MCWT.clim.PT_ss_gf) <- gsub("_wc", "", names(MCWT.clim.PT_ss_gf))
    names(MCWT.clim.PT_sl_gf) <- gsub("_wc", "", names(MCWT.clim.PT_sl_gf))
    my_color <- colorRampPalette(c("royalblue", "royalblue", "grey95", "grey95", "darkorange", "darkorange"))(20)
    
    names(MCWT.clim.MV_gf) <- gsub("TRY_", "", names(MCWT.clim.MV_gf))
    names(MCWT.clim.PT_ss_gf) <- gsub("TRY_", "", names(MCWT.clim.PT_ss_gf))
    names(MCWT.clim.PT_sl_gf) <- gsub("TRY_", "", names(MCWT.clim.PT_sl_gf))
    
    MC.mv.gf <- Mat.corel.CWT.clim(MCWT.clim.MV_gf[Trait], MCWT.clim.MV_gf[Trait], 
                                   I.confiance = 0.95,
                                   Display.pval = "pch", my_color = my_color, Display = "lower",
                                   Disp.R = "number", Label = F, Average = F,  Bar.pos = "n")
    
    MC.ss.gf <- Mat.corel.CWT.clim(MCWT.clim.PT_ss_gf[Trait], MCWT.clim.PT_ss_gf[Trait],
                                   I.confiance = 0.95, return.pick = F,
                                   my_color = my_color, Display.pval = "pch", Display = "lower",
                                   Disp.R = "number", Label = F, Average = F,  Bar.pos = "n")
    
    MC.sl.gf <- Mat.corel.CWT.clim(MCWT.clim.PT_sl_gf[Trait], MCWT.clim.PT_sl_gf[Trait], 
                                   I.confiance = 0.95,
                                   Display.pval = "pch", my_color = my_color, Display = "lower",
                                   Disp.R = "number", Label = F, Average = F,  Bar.pos = "n")
    
    # Diff.mat <- round(abs(MC.mv.gf$R2 - MC.ss.gf$R2), digits = 2)
    Diff.mat <- round(MC.mv.gf$R2 - MC.ss.gf$R2, digits = 2)
    CP <- corrplot(Diff.mat, type = "lower", col = my_color,
                   tl.col="black", tl.srt=45, tl.cex = .7, method = "number", cl.align.text = "l", cl.pos = "n", 
                   sig.level = 0.05, insig = "pch", pch.cex = 2)
    dev.off()
    
    #### R2 compar ####
    H = 400
    W = 700
    Save.matcor.trait = "Figures/ACA/Trait/Correlation/R2_compar/CWM/R2_comp_all.pdf"
    
    p1 <- R2.compar(MV = MC.mv.gf, MP = MC.ss.gf, Repel.outliers = F, p.val.show = T, Pourc.coerrent = F,
                    Xlab = "r (CWM-ACASP-fine)", Ylab = "r (CWM-ACAV)", Nb.labs = 4, Panel.lab = "(A)")
    p2 <- R2.compar(MV = MC.mv.gf, MP = MC.sl.gf, Leg.pos = "none", Pourc.coerrent = T, p.val.show = T,
                    Xlab = "r (CWM-ACASP-coarse)", Ylab = "r (CWM-ACAV)", Nb.labs = 4, Panel.lab = "(B)")
    # p3 <- R2.compar(MV = MC.ss.gf, MP = MC.sl.gf,
    #                 Xlab = "r (CWM-ACASP-co)", Ylab = "r (CWM-ACASP-fi)", Nb.labs = 4,
    #                 H = 500, W = 500, Save.plot = "Figures/ACA/Trait/Correlation/R2_compar/CWM/R2_comp_ACASP.pdf")
    pall <- p1 + p2 #+ p3 #+ plot_layout(guides = "collect")
    ggsave(file = Save.matcor.trait, pall, width = W*0.01041666666667, height = H*0.01041666666667)
    
  }
  
  #### Graph network ####
  Correl.graph = F
  if(Correl.graph == T){
    # Type.graph "linear", "kk", circular = T/F 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl', 'lgl'.
    DB.full <- c(4,5,7,11,12,13,18,20,22,24,26,28) # world clim
    T1gr = "stress"
    
    Gf.PT.ss.gf <- Mat.corel.network(MCWT.clim.PT_ss_gf[,DB.full], Color.scale.bar = T, Type.graph = T1gr)
    Gf.PT.sl.gf <- Mat.corel.network(MCWT.clim.PT_sl_gf[,DB.full], Color.scale.bar = F, Type.graph = T1gr)
    Gf.MV.gf <- Mat.corel.network(MCWT.clim.MV_gf[,DB.full], Color.scale.bar = F, Type.graph = T1gr)
    
    pgraph <- Gf.MV.gf + Gf.PT.ss.gf + Gf.PT.sl.gf
    W = 2700
    H = 1000
    ggsave(filename = "Figures/ACA/Trait/Correlation/Graphs/MCor_network_full_gf_wc.pdf", pgraph, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    DB.full <- c(4,5,7,11,12,13,16,17,19,21,23,25,27) # chelsa
    
    Gf.PT.ss.gf <- Mat.corel.network(MCWT.clim.PT_ss_gf[,DB.full], Color.scale.bar = T, Type.graph = T1gr)
    Gf.PT.sl.gf <- Mat.corel.network(MCWT.clim.PT_sl_gf[,DB.full], Color.scale.bar = F, Type.graph = T1gr)
    Gf.MV.gf <- Mat.corel.network(MCWT.clim.MV_gf[,DB.full], Color.scale.bar = F, Type.graph = T1gr)
    
    pgraph <- Gf.MV.gf + Gf.PT.ss.gf + Gf.PT.sl.gf
    ggsave(filename = "Figures/ACA/Trait/Correlation/Graphs/MCor_network_full_gf_chel.pdf", pgraph, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    }
  
  #### SSM test ####
  SSM.clim = F
  if(SSM.clim == T){
    Keep.trait <- c("TRY_38", "TRY_59", "TRY_3108", "TRY_3114", "TRY_22", "TRY_4", "TRY_46", "TRY_3116")
    MT1 <- MCWT.clim[which(names(MCWT.clim) %in% Keep.trait)]
    #Keep.trait <- c(grepl("TRY", names(MCWT.clim)))
    # MT1 <- MCWT.clim[Keep.trait]
    MC1 <- subset(MCWT.clim, select = c(MAP, MAAT, MTWAQ, MPWAQ, MPCOQ, MTCOQ))
    SG <- SSM.CWT.clim(MT1, MC1, Full.stats = F, Nb.max = 2, 
                       Save.path = "Resultats/ACA/Traits/MR/SSM_trait_clim_v1.csv")
    
    print(SG)
    }
  
  #### Prediction categorical values ####
  Biom.classifier <- function(M, Trait, Biome, Titre, Method, Only.stat,
                              Scaling, Normalization, Save.plot, W, H){
    #### Libraries & settings ####
    #library(e1071)
    #library(caTools)
    library(randomForest)
    library(gbm)
    library(dismo) # BRT
    library(caret)
    # library(neuralnet) # Neural Network
    library(nnet) # Neural Network
    library(NeuralNetTools)
    if(missing(Normalization)){Normalization = F}
    if(missing(Only.stat)){Only.stat = F}
    if(missing(Scaling)){Scaling = T}
    if(missing(Save.plot)){Save.plot = NULL}
    if(missing(W)){W = NULL}
    if(missing(H)){H = NULL}
    if(missing(Method)){Method = "NBC"}
    
    
    M <- M[c(Biome, Trait)]
    names(M)[1] <- "Species"
    Keep.species <- names(which(table(M$Species)>5))
    M <- M[which(M$Species %in% Keep.species),]
    M$Species = as.factor(M$Species)
    M$Species <- droplevels(M$Species)
    
    #### Splitting train and test data ####
    set.seed(1234)
    data1 <- sample(2, nrow(M), 
                    replace = T, 
                    prob = c(0.8, 0.2))
    train_cl <- M[data1 == 1,]
    test_cl <- M[data1 == 2,]
    
    #### Feature Scaling / standardization #### 
    if(Scaling == T){
      train_cl[, 2:ncol(train_cl)] <- scale(train_cl[, 2:ncol(train_cl)])
      test_cl[, 2:ncol(test_cl)] <- scale(test_cl[, 2:ncol(test_cl)])
      }
    
    if(Normalization == T){
    normalize <- function(x){return((x - min(x)) / (max(x) - min(x)))}
      train_cl[, 2:ncol(train_cl)] <- as.data.frame(lapply(train_cl[, 2:ncol(train_cl)], normalize))
      test_cl[, 2:ncol(test_cl)] <- as.data.frame(lapply(test_cl[, 2:ncol(test_cl)], normalize))
      }
    
    #### Fitting on train #### 
    set.seed(120)  # Setting Seed
    if(Method == "NBC"){classifier_cl <- naiveBayes(Species ~ ., data = train_cl)} # Naive Bayes Model
    if(Method == "BRT"){classifier_cl <- gbm.step(data = train_cl, gbm.x = 2:ncol(train_cl), gbm.y = 1, family = "poisson", 
                                                  #verbose = F, call = F,
                                                  tree.complexity = 5, #tolerance.method= "fixed", tolerance = 0.1, 
                                                  learning.rate = 0.01, bag.fraction = 0.5#, max.trees = 7000
                                                  )} # Boosted Regression Trees
    if(Method == "GBM"){classifier_cl <- gbm(Species ~.,
                                             data = train_cl,
                                             distribution = "gaussian",
                                             cv.folds = 5,
                                             shrinkage = 0.1,
                                             n.minobsinnode = 10,
                                             n.trees = 100
                                             )} # Gradient Boosting Model
    if(Method == "RF"){classifier_cl <- randomForest(Species ~ ., data = train_cl, ntree = 1000, na.action = na.omit)} # Random Forest
    if(Method == "NN"){#classifier_cl <- neuralnet(Species ~ .,data = train_cl, hidden = c(10,6), act.fct = "logistic",  # hidden=c(7,6,5) pour avoir plusieurs noeuds
                       #             linear.output = F, stepmax=10^5, threshold = 0.01)
                       classifier_cl <- nnet(Species ~ ., data = train_cl, size = 12,  decay = 0.0001, maxit = 2000, verbose = F)
                      X11()
                      plotnet(classifier_cl, bias = F, circle_cex = 6, line_stag = 0.025, alpha_val = 0.3, 
                              circle_col = "grey90", neg_col = "royalblue", pos_col = "darkorange", bord_col = "grey30",
                              max_sp = T, pad_x = 0.8)
                      } # Neural Network
    
    #### Prediction on test #### 
    if(Method == "GBM"){y_pred <- predict.gbm(object = classifier_cl,
                                              newdata = test_cl,
                                              n.trees = 100,
                                              type = "response")
                        y_pred <- colnames(y_pred)[apply(y_pred, 1, which.max)]}
    if(Method == "NN"){#y_pred <- compute(classifier_cl, test_cl)
    #                    y_pred <- data.frame(y_pred$net.result)
    #                    names(y_pred) <- unique(test_cl$Species)
    #                    y_pred <- as.factor(apply(y_pred, 1, function(x) return(names(x)[which(x == max(x))])))
    #                    levels(y_pred) <- levels(M$Species) 
    y_pred <- predict(classifier_cl, test_cl, type = "class")
                       }
    else{y_pred <- predict(classifier_cl, newdata = test_cl)}
    
    #### Confusion Matrix ####
    cm <- table(test_cl$Species, y_pred)
    df.cm <- confusionMatrix(cm)$table
    #print(cm)
    
    #convert confusion matrices to tables, and binding them together
    df.cm.col <- df.cm / rowSums(df.cm)
    df.table <- reshape2::melt(df.cm)
    df.table.col <- reshape2::melt(df.cm.col)
    df.table <- left_join(df.table, df.table.col, by =c("Var1", "y_pred"))
    
    #calculate accuracy and class accuracy
    acc.vector <- c(diag(df.cm)) / c(rowSums(df.cm))
    recall.vector <- c(diag(df.cm)) / c(colSums(df.cm))
    
    F1.vector <- round(2 * acc.vector * recall.vector / (acc.vector + recall.vector), 2)
    class.acc <- data.frame(y_pred = "Biome Acc.", Var1 = names(acc.vector), value = acc.vector)
    class.F1 <- data.frame(y_pred = "Biome F1", Var1 = names(F1.vector), value = F1.vector)
    acc <- sum(diag(df.cm)) / sum(df.cm)
    # recall.tot <- sum(diag(df.cm)) / c(colSums(df.cm))
    # F1.tot <- round(2 * acc * recall.tot / (acc + recall.tot), 2)
    F1.tot <- mean(2 * acc.vector * recall.vector / (acc.vector + recall.vector))
    
    #### Save plots ####
    if(is.null(Save.plot) == F){
      Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
      dir.create(file.path(Path.to.create), showWarnings = FALSE)
      if(is.null(W) == F & is.null(H) == F){
        pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
      else{pdf(file = Save.plot)}}
    
    #### Graphical settings ####
    true.lab = "True Biome"
    pred.lab ="Predicted Biome"
    high.col = "darkorange" 
    low.col = 'white'
    #### plot ####
    if(Only.stat == F){
      print(paste0("Accuracy: ", round(100*acc, 1), "%"))
      print(paste0("Accuracy: ", round(F1.tot, 2)))
      print(paste0("Nb test sites: ", nrow(test_cl)))
      p <- ggplot() +
        geom_tile(aes(x=y_pred, y=Var1, fill=value.y),
                  data=df.table, size=0.2, color=grey(0.5)) +
        geom_tile(aes(x=y_pred, y=Var1),
                  data=df.table[df.table$Var1==df.table$y_pred, ], size=1, color="black", fill = 'transparent') +
        scale_x_discrete(position = "top",  limits = c(levels(df.table$y_pred), "Biome Acc.", "Biome F1")) +
        scale_y_discrete(limits = rev(unique(levels(df.table$y_pred)))) +
        labs(x=pred.lab, y=true.lab, fill=NULL, 
             title = paste0(Titre, "\nAccuracy = ", round(100*acc, 1), "%", 
                            "\nF1 = ", round(F1.tot, 2),
                            "\nNb test sites: ", nrow(test_cl), "\nNb traits: ", length(Trait))) +
        geom_text(aes(x=y_pred, y=Var1, label=value.x),
                  data=df.table, size=4, colour="black") +
        geom_text(data = class.acc, aes(y_pred, Var1, label = paste0(round(100*value), "%"),  colour = value)) +
        geom_text(data = class.F1, aes(y_pred, Var1, label = value, colour = value)) +
        scale_fill_gradient(low=low.col, high=high.col, labels = scales::percent,
                            limits = c(0,1), breaks = c(0,0.5,1)) +
        scale_colour_gradient(low = "grey", high = "darkblue")+
        guides(size=F, colour = F) +
        theme_bw() +
        theme(panel.border = element_blank(), legend.position = "bottom",
              axis.text = element_text(color='black'), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.text.x.top = element_text(angle = 30, vjust = 0, hjust = 0)) +
        coord_fixed()
      print(p)
      
      if(is.null(Save.plot) == F){dev.off()}
      return(df.table)
      }
    #### Only.stats ####
    else{
      Return.stat <- t(data.frame(Accuracy = round(acc, digits = 2),
                      F1 = round(F1.tot, digits = 2),
                      Nb.test = nrow(test_cl),
                      Nb.trait = length(Trait)
                      ))
      # names(Return.stat) <- eval(parse(text = M))
      print(Return.stat)
      
    }

  }
  
  # NBC.CWT.pol <- Biom.classifier(MCWT.clim.PT_ss_gf, Trait = c("TRY_Height", "TRY_LeafN", "TRY_SLA", "TRY_SeedMass", "TRY_SSD", "TRY_LeafArea"), Biome = "Biome", Method = "RF",
  #                                Normalization = F, Scaling = T, Only.stat = T)
  # 
  
  Biom.model.matrix.confusion = F
  if(Biom.model.matrix.confusion == T){
    NBC.CWT.pol <- Biom.classifier(MCWT.clim.PT_ss_gf, Trait = c("TRY_Height", "TRY_LeafN", "TRY_SLA", "TRY_SeedMass", "TRY_SSD", "TRY_LeafArea"), Biome = "Biome", Method = "RF",
                                          Titre = "Biome reconstruction based on ACASP-ss (gf) Random Forest", Normalization = F, Scaling = T,
                                          H = 700, W = 900, Save.plot = "Figures/ACA/Trait/Predictions/Biome_pred_NBC_6traits.pdf")
  
    NBC.CWT.pol.sl.6t <- Biom.classifier(MCWT.clim.PT_sl_gf, Trait = c("TRY_Height", "TRY_LeafN", "TRY_SLA", "TRY_SeedMass", "TRY_SSD", "TRY_LeafArea"), Biome = "Biome", Method = "RF",
                                          Titre = "Biome reconstruction based on ACASP-sl (gf) Random Forest", Normalization = F, Scaling = T,
                                          H = 700, W = 900, Save.plot = "Figures/ACA/Trait/Predictions/Biome_pred_NBC_6traits_sl.pdf")
  
  
    NBC.CWT.pol <- Biom.classifier(MCWT.clim.PT_ss_gf, Trait = names(MCWT.clim.PT_ss)[grep("TRY", names(MCWT.clim.PT_ss))], Biome = "Biome", Method = "RF",
                                          Titre = "Biome reconstruction based on pollen-trait Random Forest", Normalization = F, Scaling = F,
                                          H = 700, W = 900, Save.plot = "Figures/ACA/Trait/Predictions/Biome_pred_NBC_alltraits_ss.pdf")
    
    NBC.CWT.pol.SL <- Biom.classifier(MCWT.clim.PT_sl_gf, Trait = names(MCWT.clim.PT_sl)[grep("TRY", names(MCWT.clim.PT_sl))], Biome = "Biome", Method = "RF",
                                          Titre = "Biome reconstruction based on pollen-trait Random Forest", Normalization = F, Scaling = F,
                                          H = 700, W = 900, Save.plot = "Figures/ACA/Trait/Predictions/Biome_pred_NBC_alltraits_sl.pdf")
     
    NBC.CWT.veg <- Biom.classifier(MCWT.clim.MV_gf, Trait = names(MCWT.clim.MV_gf)[grep("TRY", names(MCWT.clim.MV_gf))], Biome = "Biome", Method = "RF",
                                          Titre = "Biome reconstruction based on pollen-trait Random Forest", Normalization = F, Scaling = F,
                                          H = 700, W = 900, Save.plot = "Figures/ACA/Trait/Predictions/Biome_pred_NBC_alltraits_ACAV.pdf")
    
    NBC.CWT.veg.6t <- Biom.classifier(MCWT.clim.MV_gf, Trait = c("TRY_Height", "TRY_LeafN", "TRY_SLA", "TRY_SeedMass", "TRY_SSD", "TRY_LeafArea"), Biome = "Biome", Method = "RF",
                                   Titre = "Biome reconstruction based on ACAV (gf) Random Forest", Normalization = F, Scaling = T,
                                   H = 700, W = 900, Save.plot = "Figures/ACA/Trait/Predictions/Biome_pred_NBC_ACAV_6traits.pdf")
    
    All.models <- list(MCWT.clim.MV_gf, MCWT.clim.MV, MCWT.clim.PT_ss_gf)
    # for(i in )
    
    
  }
  }
  
#### Plot bioclimat param. ####
Plot.clim = F
if(Plot.clim == T){
  #### Import data ####
  MPS.ACA.Clim_chel <- as_tibble(readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Clim_chel.Rds"))
  MPS.ACA.Clim_wc <- as_tibble(readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Clim_wc.Rds"))
  MPS.ACA.Biom <- as_tibble(readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Biom.Rds"))
  
  MV.ACA.Clim_chel <- as_tibble(readRDS("Resultats/ACA/Vegetation/DBACA_Clim_chel.Rds"))
  MV.ACA.Clim_wc <- as_tibble(readRDS("Resultats/ACA/Vegetation/DBACA_Clim_wc.Rds"))
  MV.ACA.Biom <- as_tibble(readRDS("Resultats/ACA/Vegetation/DBACA_Biom.Rds"))
  
  #### Extract clim ACA cells ####
  Extract.bioclim.ACA = F
  if(Extract.bioclim.ACA == T){
    source("Scripts/Climat_extract.R")
    Xsamp <- seq(extent(ACA.bo)[1], extent(ACA.bo)[2], length.out = 200)
    Ysamp <- seq(extent(ACA.bo)[3], extent(ACA.bo)[4], length.out = 200)
    Psamp <- expand.grid(Xsamp, Ysamp)
    names(Psamp) <- c("Long", "Lat")
    Psamp <- SpatialPointsDataFrame(coords = Psamp, data = Psamp, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    Psamp <- Psamp[ACA.bo,]
    Psamp <- data.frame(Psamp)
    Psamp.var.chel <- Clim.param.extraction(M = Psamp[c(1,2)], Clim.cal = T, Biome = T,
                                       All.param = F, Chelsa = T,
                                       Altitude = T, Aridity = T)
    
    Psamp.var.wc <- Clim.param.extraction(M = Psamp[c(1,2)], Clim.cal = T, Biome = T,
                                            All.param = F, Chelsa = F,
                                            Altitude = T, Aridity = T)
    
    saveRDS(Psamp.var.chel, file = "Resultats/ACA/Climat/Chelsa_extract_ACA_df.Rds")
    saveRDS(Psamp.var.wc, file = "Resultats/ACA/Climat/WC_extract_ACA_df.Rds")
    }
  else{
    Psamp.var.chel <- as_tibble(readRDS("Resultats/ACA/Climat/Chelsa_extract_ACA_df.Rds"))  
    Psamp.var.wc <- as_tibble(readRDS("Resultats/ACA/Climat/WC_extract_ACA_df.Rds"))
    }  
  
  #### Data prep ####
  MC.extract.pca.prep <- function(MSample.Clim, MSample.Biom, Pclim){
                          MC <- cbind(MSample.Clim, MSample.Biom)
                          MC <- cbind(MC, Type = "MPS")
                          row.names(MC) <- paste("MPS", row.names(MC), sep = "-")
                          PV <- data.frame(cbind(Pclim, Type = "PV"))
                          row.names(PV) <- paste("PV", row.names(PV), sep = "-")
                          PV <- PV[names(MC)]
                          MCPV <- rbind(PV, MC)
                          MCPV <- MCPV[MCPV$Biome.no %in% c(13,8,10,5,4,6),]
                          return(MCPV)
                          }
      
  MCPV.MPS <- MC.extract.pca.prep(MPS.ACA.Clim_chel, MPS.ACA.Biom[c(5,6)], Psamp.var.chel)
  MCPV.MV <- MC.extract.pca.prep(MV.ACA.Clim_chel, MV.ACA.Biom[c(5,6)], Psamp.var.chel)
  MCPV.MPS.wc <- MC.extract.pca.prep(MPS.ACA.Clim_wc, MPS.ACA.Biom[c(5,6)], Psamp.var.wc)
  MCPV.MV.wc <- MC.extract.pca.prep(MV.ACA.Clim_wc, MV.ACA.Biom[c(5,6)], Psamp.var.wc)
  
  
  #### PCA Bioclim ####
  Chelsa = F
  if(Chelsa == T){
    BioCl.pca.mps <- PCA.bioclim(MCPV.MPS, transp_OK = T, Scale.PCA = 7,
                        Cluster.core = "Biome", Shape = 21, Legend.size = 8,
                        Groupes = list(Water = c("AI", "MAP", "MPWAQ", "MPCOQ"),
                                       Temperature = c("MAAT", "MTWAQ", "MTCOQ"), Altitude = "Altitude"),
                        # Site.name = "Bioclim CHELSA+CGAR for ACA (res: 0.5°)",
                        Site.name = "", Num.facet = "(A) Pollen", return.pick = T, Legend.position = "bottom",
                        # Save.plot = "Figures/ACA/Trait/Surface/PCA_bioclim/PCA_bioclim_ACA_mps_ch.pdf", H = 500, W = 800
                        )
    
    BioCl.pca.mv <- PCA.bioclim(MCPV.MV, transp_OK = T, Scale.PCA = 7,
                             Cluster.core = "Biome", Shape = 24,
                             Groupes = list(Water = c("AI", "MAP", "MPWAQ", "MPCOQ"),
                                            Temperature = c("MAAT", "MTWAQ", "MTCOQ"), Altitude = "Altitude"),
                             # Site.name = "Bioclim CHELSA+CGAR for ACA (res: 0.5°)",
                             Site.name = "", Num.facet = "(B) Vegetation", return.pick = T, Legend.position = "none",
                             # Save.plot = "Figures/ACA/Trait/Surface/PCA_bioclim/PCA_bioclim_ACA_vm_ch.pdf", H = 500, W = 800
                             )
    
    full.pca.clim <- BioCl.pca.mps + BioCl.pca.mv / plot_layout(guides = "collect")
    W = 1000
    H = 600
    ggsave(filename = "Figures/ACA/Trait/Surface/PCA_bioclim/PCA_bioclim_ACA_tot_chel.pdf", full.pca.clim, width = W*0.026458333, height = H*0.026458333, units = "cm")
    }
    
  Worldclim = F
  if(Worldclim == T){
    BioCl.pca.mps <- PCA.bioclim(MCPV.MPS.wc, transp_OK = T, Scale.PCA = 7, 
                        Cluster.core = "Biome", Shape = 21, 
                        Legend.size = 8, Dot.size = 2,
                        Groupes = list(Water = c("AI", "MAP", "MPWAQ", "MPCOQ"),
                                       Temperature = c("MAAT", "MTWAQ", "MTCOQ"),
                                       Altitude = "Altitude"),
                        Site.name = "", Num.facet = "(A) Pollen", return.pick = T, Legend.position = "bottom", 
                        # Save.plot = "Figures/ACA/Trait/Surface/PCA_bioclim/PCA_bioclim_ACA_mps_wc.pdf", H = 500, W = 800
                        )
    
    BioCl.pca.mv <- PCA.bioclim(MCPV.MV.wc, transp_OK = T, Scale.PCA = 7,
                             Cluster.core = "Biome", Shape = 24,
                             Groupes = list(Water = c("AI", "MAP", "MPWAQ", "MPCOQ"),
                                            Temperature = c("MAAT", "MTWAQ", "MTCOQ"),
                                            Altitude = "Altitude"),
                             Site.name = "", Num.facet = "(B) Vegetation", return.pick = T, Legend.position = "none",
                             # Save.plot = "Figures/ACA/Trait/Surface/PCA_bioclim/PCA_bioclim_ACA_vm_wc.pdf", H = 500, W = 800
                             )
    
    full.pca.clim <- BioCl.pca.mps + BioCl.pca.mv / plot_layout(guides = "collect")
    W = 1000
    H = 600
    ggsave(filename = "Figures/ACA/Trait/Surface/PCA_bioclim/PCA_bioclim_ACA_tot_wc.pdf", full.pca.clim, width = W*0.026458333, height = H*0.026458333, units = "cm")
    }
  
  #### Biplot Climat ACA vs. pollen samp ####
  Biplot.clim = T
  if(Biplot.clim == T){
    BioCl.summer.plot <- Biplot.bioclim(
      MC.area = MCPV.MPS[MCPV.MPS$Type == "PV",],
      Same.mat = F, Limites = c(0, 1000, 26, 58),
      MC.samp = MCPV.MPS[MCPV.MPS$Type == "MPS",], Add.reg = T, R2.pos = "topright",
      Leg.pos = "none", Emprise.size = 60, Emprise.alpha = .5, Emprise.bin = T,
      PClim1 = "MPWAQ", PClim2 = "Latitude", return.pick = T, Show.density = T, 
      Save.plot = "Figures/ACA/Trait/Surface/PCA_bioclim/Biplot_bioclim_ACA_PollenSamples.pdf", H = 600, W = 500
      )
    
    BioCl.summer.plot.veg <- Biplot.bioclim(
      MC.area = MCPV.MV[MCPV.MV$Type == "PV",],
      Same.mat = F, Limites = c(0, 1000, 26, 58),
      MC.samp = MCPV.MV[MCPV.MV$Type == "MPS",], Add.reg = T, R2.pos = "topright",
      Leg.pos = "none", Emprise.size = 60, Emprise.alpha = .5, Emprise.bin = T,
      PClim1 = "MPWAQ", PClim2 = "Latitude", return.pick = T, Show.density = T,
      Save.plot = "Figures/ACA/Trait/Surface/PCA_bioclim/Biplot_bioclim_ACA_VegetPlot.pdf", H = 600, W = 500
      )
  
    
    # BioCl.tot <- (BioCl.summer.plot) + (BioCl.summer.plot.veg) + plot_layout(widths = c(2,1))

    # H = 600
    # W = 1000
    # Save.plot =  "Figures/ACA/Trait/Surface/PCA_bioclim/Biplot_bioclim_ACA_tot.pdf"
    # ggsave(filename = Save.plot, BioCl.tot, width = W*0.026458333, height = H*0.026458333, units = "cm")
    }
  }

#### Plot DB surface (pollen + veget) ####
Plot.surf.DB = F
if(Plot.surf.DB == T){
  #### Import data ####
  MP_sl <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_sl_clean.Rds")
  MP_ss <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Pol_PT_ss_clean.Rds")
  MPS.ACA.Biom <- data.frame(readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Biom.Rds"))
  
  #### Diag pol tot ####
  DP.t = F
  if(DP.t == T){
    # MPS.ACA.Biom <- MPS.ACA.Biom[order(MPS.ACA.Biom$Latitude),]
    MPS.ACA.Biom <- MPS.ACA.Biom[order(MPS.ACA.Biom$Biome.no),]
    MPS.ACA.Biom$Biom.ordin <- seq(1:nrow(MPS.ACA.Biom)) 
    MPS.ACA.Biom <- cbind(name = row.names(MPS.ACA.Biom), MPS.ACA.Biom)
    write.csv(MPS.ACA.Biom, "Resultats/World_DB/Pollen/Merging_DB/DB13549/Bioclim_ordin.csv")
    MP_sl_DP <- data.frame(t(MP_sl))
    # names(MP_sl_DP) <- gsub("\\.", "-", names(MP_sl_DP))
    
    source("Scripts/Pollen_surf.R")
    Plot.Mong.pollen.surf <- Diag.pol.surf(#MP = as.data.frame(MP_sl_DP[,c(1:1000)]),
                                           MP = MP_sl_DP,
                                           Ordin.path = "Resultats/World_DB/Pollen/Merging_DB/DB13549/Bioclim_ordin.csv",
                                           # Index = "Import/Mongolia/Pollen/Indexes/Index_pollen_Mongolia.csv",
                                           # Name.zone = "Mongolia - Bioclimate ordinated",
                                           # Y.legend = "Sample sites",
                                           # Sort.taxon = "Auto", # Auto, Manual
                                           # Manual.sort = c("Pinus.sylvestris", "Pinus.sibirica", "Betula","Alnus", "Picea", "Abies",  "Larix", "Salix",
                                           #                 "Artemisia", "Poaceae",
                                           #                 "Cyperaceae", "Brassicaceae", "Convolvulus", 
                                           #                 "Rumex",  "Amaranthaceae", "Caryophyllaceae","Thalictrum",
                                           #                 "Other"),
                                           # TaxonLabel = c("AP", "NAP", "Pinus sylvestris", "Pinus sibirica", "Betula sp.", "Alnus sp.", "Picea obovata", "Abies sibirica", 
                                           #                "Larix sibirica", "Salix sp.", "Artemisia sp.", "Poaceae", "Cyperaceae", "Brassicaceae", "Convolvulus-type",
                                           #                "Rumex sp.", "Amaranthaceae", "Caryophyllaceae", "Thalictrum sp.", "Other taxa"),
                                           Csv.sep = ",",
                                           # Max_seuil = 0.8, 
                                           Max_seuil = 0.1, 
                                           # AP.NAP = T,
                                           CONISS = T,
                                           Nzone = 7,
                                           # Abiot.plot = c("MAP", "MAAT", "Altitude"), # Fonction à faire ! Transform MAAT en MAAT.ordin (1, 2, 3... et pas 1.0° 2.3° 3.4°...)
                                           Sort = "Biom.ordin",  # Altitude, MAP.ordin, MAAT, Latitude, Longitude, Ordination
                                           # Label.eco = c("Light taiga", "Dark taiga-","birch subtaiga", "Forest-steppe", "Steppe", "Alpine meadow", "Steppe-desert", "Desert"),
                                           # Sort.eco = Sort.eco.mong,
                                           H = 5000, W = 1900,
                                           Save.plot = "Figures/ACA/Trait/Surface/Pollen/DP_ACASP_sl.pdf")
    
    W = 1100
    H = 550
    Save.plot = "Figures/ACA/Trait/Surface/Pollen/PCA_pollen_surf_sl_2.pdf"
    pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)
    par(mfrow=c(1,1))

    PCA.pollen.Mong <- PCA.pollen.surf(data.frame(t(Plot.Mong.pollen.surf)),               # Import Matrice pour la PCA
                                       # Cluster.path = Meco,
                                       # Cluster.groups = "Vegetation.community",
                                       Csv.sep =",", Simple.title = T,
                                       transp_OK = T,                         # Log trans (T), (F) sinon
                                       Scale.PCA = 1,                         # 1 or 2
                                       # Save.path = "Resultats/Mongolia/Pollen/Surface/MP_mong.csv",
                                       #Site.name = "Mongolia surface sample",
                                       # Sort.eco = Sort.eco.mong,
                                       Type.samples = "Pollen"
                                       )     # Nom du site

    dev.off()
    }
  
  #### Determine most important taxa ####
  Dominant.table = T
  if(Dominant.table == T){

    
  Table.Taxon <- read.csv(file="Import/World_DB/Pollen/Odile_DB/Corresp_name_full_V10.csv", sep=",",dec=".", row.names = 1,  header=T, stringsAsFactors = F)
    
  MP_sl <- cbind(MP_sl, MPS.ACA.Biom[c(5,6)])
  MP_ss <- cbind(MP_ss, MPS.ACA.Biom[c(5,6)])

  Pollen.PCA.sl <- PCA.bioclim(MP_sl, transp_OK = F, Scale.PCA = 5, Contrib = T,
                               Cluster.core = "Biome", Shape = 21,
                               Legend.size = 8, Dot.size = 1, #Ellipse = T, Show.centroid = F, Ellipse.opa = 0.3, #Groupes = T,
                               # Groupes = list(Water = c("AI", "MAP", "MPWAQ", "MPCOQ"),
                               #                Temperature = c("MAAT", "MTWAQ", "MTCOQ"),
                               #                Altitude = "Altitude"),
                               Site.name = "", Num.facet = "(A) Pollen", return.pick = F, 
                               Save.plot = "Figures/ACA/Trait/Surface/Pollen/PCA_pollen_surf_sl.pdf", H = 500, W = 800)

   Pollen.PCA.ss <- PCA.bioclim(MP_ss, transp_OK = F, Scale.PCA = 4, Contrib = T,
                               Cluster.core = "Biome", Shape = 21,
                               Legend.size = 8, Dot.size = 1,
                               # Groupes = list(Water = c("AI", "MAP", "MPWAQ", "MPCOQ"),
                               #                Temperature = c("MAAT", "MTWAQ", "MTCOQ"),
                               #                Altitude = "Altitude"),
                               Site.name = "", Num.facet = "(A) Pollen", return.pick = F,
                               Save.plot = "Figures/ACA/Trait/Surface/Pollen/PCA_pollen_surf_ss.pdf", H = 500, W = 800)


  Pollen.PCA.sl.contrib <- Pollen.PCA.sl[,1]+Pollen.PCA.sl[,2]
  Pollen.PCA.sl.contrib <- Pollen.PCA.sl.contrib[order(Pollen.PCA.sl.contrib, decreasing = T)]
  Pollen.PCA.sl.contrib1 <- Pollen.PCA.sl.contrib[Pollen.PCA.sl.contrib > 1]
  Pollen.PCA.sl.contrib2 <- Pollen.PCA.sl.contrib[Pollen.PCA.sl.contrib > 0.01]
  Tab.to.exp <- data.frame(Rank = seq(1, nrow(data.frame(Pollen.PCA.sl.contrib2))), 
                           Pollen.type = names(Pollen.PCA.sl.contrib2), 
                           "PCA contribution" = round(Pollen.PCA.sl.contrib2, digits = 3))
  
  Tab.to.exp$Identification.level = Table.Taxon[match(Tab.to.exp$Pollen.type, Table.Taxon$Nom),12]
  Tab.to.exp$Identification.level[Tab.to.exp$Identification.level == 1] <- "Family"
  Tab.to.exp$Identification.level[Tab.to.exp$Identification.level == 2] <- "Sub-family"
  Tab.to.exp$Identification.level[Tab.to.exp$Identification.level == 3] <- "Genus"
  Tab.to.exp$Identification.level[Tab.to.exp$Identification.level == 4] <- "Sub-genus"
  Tab.to.exp$Identification.level[Tab.to.exp$Identification.level == 5] <- "Species"
  Pollen.PCA.sl.contrib1 <- names(Pollen.PCA.sl.contrib1)
  Pollen.PCA.sl.contrib2 <- names(Pollen.PCA.sl.contrib2)
  write.table(Tab.to.exp, "Resultats/ACA/Pollen/PCA/PCA_22_type_sl_ACA_contrib.csv", row.names = F)
  # Pollen.PCA.sl.contrib3 <- names(colSums(MP_sl)[order(colSums(MP_sl))][colSums(MP_sl)[order(colSums(MP_sl))] >100])
  saveRDS(Pollen.PCA.sl.contrib1, "Resultats/World_DB/Pollen/Merging_DB/DB13549/4Main_taxa_PCA.Rds")
  saveRDS(Pollen.PCA.sl.contrib2, "Resultats/World_DB/Pollen/Merging_DB/DB13549/22Main_taxa_PCA.Rds")}
  
}
#### Plot CWM trait vs. climat #### 
Plot.CWT = T
if(Plot.CWT == T){
  #### Import data ####
  MCWT.clim.PT_ss <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_clean.Rds")
  MCWT.clim.PT_ss_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_clean.Rds")
  MCWT.clim.PT_sl <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_clean.Rds")
  MCWT.clim.PT_sl_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_clean.Rds")
  MCWT.clim.ACAP.full <- readRDS( "Resultats/ACA/Traits/CWT_ACA_PT_all.Rds")
  MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_clean.Rds")
  MCWT.clim.MV_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_gf_clean.Rds")
  MCWT.clim.MV.tot <- readRDS("Resultats/ACA/Traits/CWT_ACA_ACAV_all.Rds")
  MCWT.clim.splot <- data.frame(readRDS("Resultats/ACA/Traits/sPlot_ACA_CWM_clim_clean.Rds"))
  
  
  LR.MCWT.clim.PT_ss_gf <- MCWT.clim.PT_ss_gf
  LR.MCWT.clim.PT_ss_gf$Type <- "Pollen"
  LR.MCWT.clim.MV_gf <- MCWT.clim.MV_gf
  LR.MCWT.clim.MV_gf$Type <- "Vegetation"
  LR.MCWT.clim.splot <- MCWT.clim.splot
  LR.MCWT.clim.splot$Type <- "Vegetation"
  
  #### Plots linéaire ####
  RL.clim.CWT = F
  if(RL.clim.CWT == T){
    RL.CWT.clim.PT <- LRelation.CWT.clim(CWT = LR.MCWT.clim.PT_ss_gf,
                              Select.Pclim = c("MTCOQ_chel", "MPCOQ_chel", "MAP_chel"),
                              Select.trait = c("TRY_SSD", "TRY_Height", "TRY_LeafArea"),
                              Strip.lab = F,
                              Select.eco = c("Type"),
                              Leg.pos = "none", Add.linear = T, Alpha = .2, Trait.lim = c(-3,2),
                              Bit.map = T, 
                              # H = 800, W = 1000, Save.plot = "Figures/ACA/Trait/Correlation/LR/LR_clim_CWT.pdf"
                              )
    
    RL.CWT.clim.V <- LRelation.CWT.clim(CWT = LR.MCWT.clim.MV_gf,
                                        Select.Pclim = c("MTCOQ_chel", "MPCOQ_chel", "MAP_chel"),
                                        Select.trait = c("TRY_SSD", "TRY_Height", "TRY_LeafArea"),
                                        Select.eco = c("Type"), Strip.lab = F, Bit.map = T,
                                        Leg.pos = "none", Add.linear = T, Alpha = .05#,
                                        # H = 800, W = 1000, Save.plot = "Figures/ACA/Trait/Correlation/LR/LR_clim_CWT_ACAV.pdf"
                                        )
    
    # RL.CWT.clim.V.splot <- LRelation.CWT.clim(CWT = LR.MCWT.clim.splot,
    #                                     Select.Pclim = c("MTWAQ_chel", "MPCOQ_chel", "MAP_chel"),
    #                                     Select.trait = c("TRY_SSD", "TRY_Height", "TRY_LeafArea"),
    #                                     Select.eco = c("Type"), Strip.lab = F, Bit.map = T,
    #                                     Leg.pos = "none", Add.linear = T, Alpha = .05, Trait.lim = c(-1,5),
    #                                     # H = 800, W = 1000, Save.plot = "Figures/ACA/Trait/Correlation/LR/LR_clim_CWT_ACAV_splot.pdf"
    #   )
    
    RL.full <- RL.CWT.clim.PT/RL.CWT.clim.V
    W = 800
    H = 1200
    Save.plot = "Figures/ACA/Trait/Correlation/LR/LR_clim_CWT_full.pdf"
    ggsave(RL.full, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm", useDingbats = TRUE)
    }
  
  RL.clim.CWT.full = F
  if(RL.clim.CWT.full == T){
    RL.CWT.clim.PT.f <- LRelation.CWT.clim(CWT = LR.MCWT.clim.PT_ss_gf,
                                         Select.Pclim = c("MPCOQ_chel", "MTCOQ_chel", "AI", "MAAT_chel", "MAP_chel"),
                                         Select.trait = c("TRY_SSD", "TRY_Height", "TRY_LeafArea", "TRY_SLA", "TRY_LeafN", "TRY_SeedMass"),
                                         # Select.eco = c("Biome"),
                                         Select.eco = c("Type"),
                                         Strip.lab = F, Bit.map = T, Add.n = T,
                                         Leg.pos = "none", Add.linear = T, Alpha = .1, Trait.lim = c(-4,4),
                                         # H = 1000, W = 1300, Save.plot = "Figures/ACA/Trait/Correlation/LR/LR_clim_CWT_ACASP_ss_full.pdf"
                                         )
  
    RL.CWT.clim.V.f <- LRelation.CWT.clim(CWT = LR.MCWT.clim.MV_gf,
                                        Select.Pclim = c("MPCOQ_chel", "MTCOQ_chel", "AI", "MAAT_chel", "MAP_chel"),
                                        Select.trait = c("TRY_SSD", "TRY_Height", "TRY_LeafArea", "TRY_SLA", "TRY_LeafN", "TRY_SeedMass"),
                                        Select.eco = c("Type"),  Add.n = T,
                                        Leg.pos = "none", Add.linear = T, Alpha = .05, Strip.lab = F, Bit.map = T,
                                        # H = 1000, W = 1300, Save.plot = "Figures/ACA/Trait/Correlation/LR/LR_clim_CWT_ACAV_full.pdf"
                                        )
    
    RL.CWT.clim.V.splot.f <- LRelation.CWT.clim(CWT = LR.MCWT.clim.splot,
                                              Select.Pclim = c("MPCOQ_chel", "MTCOQ_chel", "AI", "MAAT_chel", "MAP_chel"),
                                              Select.trait = c("TRY_SSD", "TRY_Height", "TRY_LeafArea", "TRY_SLA", "TRY_LeafN", "TRY_SeedMass"),
                                              Select.eco = c("Type"), Strip.lab = F, Bit.map = T,
                                              Leg.pos = "none", Add.linear = T, Alpha = .05, Trait.lim = c(-4,4),
                                              # H = 800, W = 1000, Save.plot = "Figures/ACA/Trait/Correlation/LR/LR_clim_CWT_ACAV_splot_full.pdf"
                                              )

    RL.full.merge <- RL.CWT.clim.PT.f/RL.CWT.clim.V.f
    W = 1300
    H = 2000
    Save.plot.2 = "Figures/ACA/Trait/Correlation/LR/LR_clim_CWT_full_merge.pdf"
    ggsave(RL.full.merge, file = Save.plot.2, width = W*0.026458333, height = H*0.026458333, units = "cm", useDingbats = TRUE)
  }
  
  RL.clim = F
  if(RL.clim == T){
    RL.CWT.clim.PT <- LRelation.CWT.clim(CWT = LR.MCWT.clim.PT_ss_gf,
                                         Select.Pclim = c("MPWAQ_chel", "MTWAQ_chel", "Latitude"),
                                         Select.trait = c("MPWAQ_chel", "MTWAQ_chel", "Latitude"),
                                         Strip.lab = F, Facet.scale = "free",
                                         Select.eco = c("Type"),
                                         Leg.pos = "none", Add.linear = T, Alpha = .2, 
                                         Bit.map = T,
                                         # H = 800, W = 1000, Save.plot = "Figures/ACA/Trait/Correlation/LR/LR_clim_CWT.pdf"
    )
    
    RL.CWT.clim.V <- LRelation.CWT.clim(CWT = LR.MCWT.clim.MV_gf,
                                        Select.Pclim = c("MPWAQ_chel", "MTWAQ_chel", "Latitude"),
                                        Select.trait = c("MPWAQ_chel", "MTWAQ_chel", "Latitude"),
                                        Select.eco = c("Type"), Strip.lab = F, Bit.map = T,
                                        Leg.pos = "none", Add.linear = T, Alpha = .05, Facet.scale = "free",
                                        # H = 800, W = 1000, Save.plot = "Figures/ACA/Trait/Correlation/LR/LR_clim_CWT_ACAV.pdf"
    )
    
    RL.full <- RL.CWT.clim.PT/RL.CWT.clim.V
    W = 800
    H = 1200
    Save.plot = "Figures/ACA/Trait/Correlation/LR/LR_clim.pdf"
    ggsave(RL.full, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm", useDingbats = TRUE)
  }
    
  #### Box plot biomes ####
  Boxplot.CWT = F
  if(Boxplot.CWT == T){
    # Box.plt.CWT.clim <- BPlot.CWT.clim(CWT = MCWT.clim.PT_ss,
    #                           Select.Pclim = c("MAAT_wc"),
    #                           Select.trait = c("TRY_LeafN", "TRY_SSD", "TRY_Height", "TRY_SeedMass", "TRY_LeafArea", "TRY_SLA"), 
    #                           Select.eco = c("Biome"),
    #                           # Select.eco = c("GLC2000"),
    #                           H = 800, W = 1250, Save.plot = "Figures/ACA/Trait/Correlation/Boxplot/Boxplot_clim_CWT.pdf"
    #                           )
    # 
    # Box.plt.CWT.clim <- BPlot.CWT.clim(CWT = MCWT.clim.MV.tot,
    #                           Select.Pclim = c("MAAT"),
    #                           Select.trait = c("TRY_LeafN", "TRY_SSD", "TRY_Height", "TRY_SeedMass", "TRY_LeafArea", "TRY_SLA"), 
    #                           Select.eco = c("Biome"),
    #                           Subgroup = c("_ACAV_gf", "_ACAV"),
    #                           # Select.eco = c("GLC2000"),
    #                           H = 800, W = 1250, Save.plot = "Figures/ACA/Trait/Correlation/Boxplot/Boxplot_clim_CWT_ACAV.pdf"
    #                           )
    
    # Box.plt.CWT.clim <- BPlot.CWT.clim(CWT = MCWT.clim.ACAP.full,
    #                                    Select.Pclim = c("MPWAQ_wc", "MAAT_wc"), Subgroup = c("_ss_gf", "_sl_gf", "_sl", "_ss"),
    #                                    Select.trait = c("TRY_LeafN", "TRY_SSD", "TRY_Height", "TRY_SeedMass", "TRY_LeafArea", "TRY_SLA"),
    #                                    # Select.trait = c("TRY_LeafThick", "TRY_SSD"),
    #                                    Select.eco = c("Biome"),
    #                                    H = 900, W = 1350, Save.plot = "Figures/ACA/Trait/Correlation/Boxplot/Boxplot_clim_CWT_ACASP.pdf")
    # 
    Total.fun <- full_join(MCWT.clim.MV.tot, MCWT.clim.ACAP.full, by = intersect(names(MCWT.clim.MV.tot), names(MCWT.clim.ACAP.full)))
    # Box.plt.CWT.clim <- BPlot.CWT.clim(CWT = Total.fun, 
    #                                    Select.Pclim = c("MPWAQ_wc", "MAAT_wc"), Subgroup = c("_ACAV_gf", "_ACAV","_ss_gf", "_sl_gf", "_sl", "_ss"),
    #                                    Select.trait = c("TRY_LeafN", "TRY_SSD", "TRY_Height", "TRY_SeedMass", "TRY_LeafArea", "TRY_SLA"),
    #                                    # Select.trait = c("TRY_LeafThick", "TRY_SSD"),
    #                                    Select.eco = c("Biome"), Legend.size = 10,
    #                                    H = 1100, W = 1500, Save.plot = "Figures/ACA/Trait/Correlation/Boxplot/Boxplot_clim_CWT_tot.pdf")
    
    # Box.plt.CWT.clim.splot <- BPlot.CWT.clim(CWT = MCWT.clim.splot,
    #                                    Select.Pclim = c("MAAT_wc"),
    #                                    Select.trait = c("TRY_LeafN", "TRY_SSD", "TRY_Height", "TRY_SeedMass", "TRY_LeafArea", "TRY_SLA"),
    #                                    Select.eco = c("Biome"), Legend.size = 10,
    #                                    H = 1100, W = 1500, Save.plot = "Figures/ACA/Trait/Correlation/Boxplot_clim_CWT_splot.pdf")
    
    Total.fun.gf <- Total.fun[c(1:22, grep("_gf", names(Total.fun)))] 
    names(Total.fun.gf) <- gsub("_ss_", "_2_fine_", names(Total.fun.gf))
    names(Total.fun.gf) <- gsub("_sl_", "_3_coarse_", names(Total.fun.gf))
    names(Total.fun.gf) <- gsub("_ACAV_", "_1_ACAV_", names(Total.fun.gf))
    
    Box.plt.CWT.clim <- BPlot.CWT.clim(CWT = Total.fun.gf, 
                                       # Subgroup = c("_ACAV_gf", "_ss_gf", "_sl_gf"),
                                       Subgroup = c("_1_ACAV_gf", "_2_fine_gf", "_3_coarse_gf"),
                                       # Select.trait = c("TRY_LeafN", "TRY_SSD", "TRY_Height", "TRY_SeedMass", "TRY_LeafArea", "TRY_SLA"),
                                       Select.trait = c("TRY_Height", "TRY_LeafArea", "TRY_LeafN","TRY_SeedMass", "TRY_SLA", "TRY_SSD"),
                                       # Select.trait = c("TRY_LeafN", "TRY_SSD"),
                                       Select.eco = c("Biome"), Legend.size = 8, Annotate = F, 
                                       Alpha.point = 0.018, Only.one.col = T,
                                       # Alpha.point = 0.1, 
                                       Nuage.point = T,
                                       Save.anova = "Resultats/ACA/Traits/ANOVA/ACA_anova.csv", ANOVA.letters = T,
                                       H = 650, W = 1000, Save.plot = "Figures/ACA/Trait/Correlation/Boxplot/Boxplot_clim_CWT_gf.pdf")
    
    
    }
  
  #### Density zone ####
  Density.plot = F
  if(Density.plot == T){
    Density.zone.CWT.clim1 <- Biplot.bioclim(MC.area = MCWT.clim.PT_ss_gf, MC.samp = MCWT.clim.PT_ss_gf, PClim1 = "TRY_Height", PClim2 = "TRY_SLA",
                                       Leg.pos = c(0.85,0.85), Same.mat = T, #Limites = c(0, 6000, 0.1, 0.5),
                                        Save.plot = "Figures/ACA/Trait/Correlation/Biplots/Biplot_SS_1.pdf", H = 900, W = 1200)
    
    Density.zone.CWT.clim2 <- Biplot.bioclim(MC.area = MCWT.clim.PT_ss_gf, MC.samp = MCWT.clim.PT_ss_gf, PClim1 = "TRY_LeafN", PClim2 = "TRY_Height",
                                       Leg.pos = c(0.85,0.85), Same.mat = T, #Limites = c(0, 6000, 0.1, 0.5),
                                        Save.plot = "Figures/ACA/Trait/Correlation/Biplots/Biplot_SS_2.pdf", H = 900, W = 1200)
    
    Density.zone.CWT.clim3 <- Biplot.bioclim(MC.area = MCWT.clim.MV_gf, MC.samp = MCWT.clim.MV_gf, PClim1 = "TRY_LeafN", PClim2 = "TRY_Height",
                                             Leg.pos = c(0.85,0.85), Same.mat = T, #Limites = c(0, 6000, 0.1, 0.5),
                                             Save.plot = "Figures/ACA/Trait/Correlation/Biplots/Biplot_MV_2.pdf", H = 900, W = 1200)
    
    Density.zone.CWT.clim2 <- Biplot.bioclim(MC.area = MCWT.clim.PT_ss_gf, MC.samp = MCWT.clim.PT_ss_gf, PClim1 = "MAAT_wc", PClim2 = "TRY_Height",
                                       Leg.pos = c(0.85,0.85), Same.mat = T, Limites = c(-15,20,-3,3),
                                        Save.plot = "Figures/ACA/Trait/Correlation/Biplots/Biplot_SS_3.pdf", H = 900, W = 1200)
    
    Density.zone.CWT.clim3 <- Biplot.bioclim(MC.area = MCWT.clim.MV_gf, MC.samp = MCWT.clim.MV_gf, 
                                             PClim1 = "MAAT_wc", PClim2 = "TRY_Height",
                                             # Leg.pos = c(0.85,0.85), 
                                             Leg.pos = c(0.1,0.85), 
                                             Same.mat = T, #Limites = c(0, 6000, 0.1, 0.5),
                                             Save.plot = "Figures/ACA/Trait/Correlation/Biplots/Biplot_MV_3.pdf", H = 900, W = 1200)
    
    Biplot.clim.MV <- Biplot.bioclim(MC.area = MCWT.clim.MV_gf, MC.samp = MCWT.clim.MV_gf, 
                                             PClim1 = "MPWAQ_wc", PClim2 = "Latitude",
                                             # Leg.pos = c(0.85,0.85), 
                                             Leg.pos = "none", Show.density = F,
                                             Same.mat = T, Limites = c(0, 500, 25, 55),
                                             # Save.plot = "Figures/ACA/Trait/Correlation/Biplots/Biplot_MV_clim.pdf", H = 900, W = 1200
                                             )
    
    Biplot.clim.PT <- Biplot.bioclim(MC.area = MCWT.clim.PT_ss_gf, MC.samp = MCWT.clim.PT_ss_gf, 
                                             PClim1 = "MPWAQ_wc", PClim2 = "Latitude",
                                             # Leg.pos = c(0.85,0.85), 
                                             Leg.pos = "none", Show.density = F,
                                             Same.mat = T, Limites = c(0, 500, 25, 55)
                                             # Save.plot = "Figures/ACA/Trait/Correlation/Biplots/Biplot_PT_clim.pdf", H = 900, W = 1200
                                             )
    
    H = 800
    W = 1300
    Biplot2 <- Biplot.clim.PT + Biplot.clim.MV 
    Save.plot = "Figures/ACA/Trait/Correlation/Biplots/Biplot_PT_MV_clim.pdf"
    ggsave(Biplot2, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm", useDingbats = TRUE)
    }
  
  #### PCA ####
  Trait.to.keep <- c(4,5,7,11,12,13,33)
  
  PCA.CWM.gf = F
  if(PCA.CWM.gf == T){
    BioCl.pca.pt.ss.gf <- PCA.bioclim(MCWT.clim.PT_ss_gf[,Trait.to.keep], transp_OK = T, Scale.PCA = 6,
                             Cluster.core = "Biome", Ellipse = F, return.pick = T,Legend.position = "none", Show.centroid = F, Dot.size = 1.2, # Ellipse.opa = 0.2,
                             Site.name = "CWM, ACASP-fine  (gap-filled)",  Num.facet = "(B)", Legend.size = 12, Marg.density.plot = T,
                             Save.plot = "Figures/ACA/Trait/Correlation/PCA/CWM/BioCl.pca.pt.ss.gf.pdf", H = 460, W = 440
                             )

    BioCl.pca.pt.sl.gf <- PCA.bioclim(MCWT.clim.PT_sl_gf[,Trait.to.keep], transp_OK = T, Scale.PCA = 6, Contrib = T,
                             Cluster.core = "Biome", Ellipse = F, return.pick = T, Show.centroid = F,  Dot.size = 1.2,
                             Legend.position = "bottom", Marg.density.plot = T,
                             Site.name = "CWM, ACASP-coarse (gap-filled)", Num.facet = "(C)", Legend.size = 12,
                             Save.plot = "Figures/ACA/Trait/Correlation/PCA/CWM/BioCl.pca.pt.sl.gf.pdf", H = 460, W = 440
                             )

    BioCl.pca.mv.gf <- PCA.bioclim(MCWT.clim.MV_gf[,Trait.to.keep], transp_OK = T, Scale.PCA = 6, Show.centroid = F, Dot.size = 1, #0.05
                                Cluster.core = "Biome", Ellipse = F, return.pick = T, Ellipse.opa = 0, #Manu.lim.x = c(6,-6),
                                Reverse.dim = F, Legend.position = "none", Marg.density.plot = T, # Manu.lim.y = c(6,-6),
                                Site.name = "CWM, ACAV (gap-filled)",  Num.facet = "(A)", Legend.size = 12,
                                Save.plot = "Figures/ACA/Trait/Correlation/PCA/CWM/BioCl.pca.mv.gf.pdf", H = 460, W = 440
                                )

    # BioCl.pca.mv.gf.mini <- PCA.bioclim(MCWT.clim.MV_gf[,Trait.to.keep], transp_OK = T, Scale.PCA = 6, Show.centroid = T, Dot.size = 0.03,
    #                                Cluster.core = "Biome", Ellipse = T, return.pick = T, Ellipse.opa = 0.25, Manu.lim.x = c(7,-7),
    #                                Reverse.dim = F, Legend.position = "none", Show.annot = F,
    #                                Site.name = "ACAV",  Num.facet = "(D)", Legend.size = 12, Show.arrow = F
    #                                )
    # 
    # BioCl.pca.pt.ss.gf.mini <- PCA.bioclim(MCWT.clim.PT_ss_gf[,Trait.to.keep], transp_OK = T, Scale.PCA = 6,
    #                                   Cluster.core = "Biome", Ellipse = T, return.pick = T,  Show.annot = F,
    #                                   Show.centroid = T,Dot.size = 0.03, Ellipse.opa = 0.25, Show.arrow = F,
    #                                   Site.name = "ACASP-ss",  Num.facet = "(E)", Legend.size = 12)
    
    # pca.CWT.full <- (BioCl.pca.mv.gf + BioCl.pca.pt.ss.gf) / (BioCl.pca.pt.sl.gf + ((BioCl.pca.mv.gf.mini+BioCl.pca.pt.ss.gf.mini)/guide_area() + plot_layout(guides = "collect")))
    pca.CWT.full <- (BioCl.pca.mv.gf + BioCl.pca.pt.ss.gf + BioCl.pca.pt.sl.gf) / 
      guide_area()  + plot_layout(guides = "collect", heights = c(9/10, 1/10))
                                                              
    W = 1400
    H = 550
    ggsave(filename = "Figures/ACA/Trait/Correlation/PCA/CWM/PCA_CWT_all_gf.pdf", pca.CWT.full, width = W*0.026458333, height = H*0.026458333, units = "cm")
    }
  
  PCA.CWT.normal = F # Run the PCA.CWM.gf before
  if(PCA.CWT.normal == T){
    BioCl.pca.pt.ss <- PCA.bioclim(MCWT.clim.PT_ss[,Trait.to.keep], transp_OK = T, Scale.PCA = 3,
                             Cluster.core = "Biome", Ellipse = F, return.pick = T,
                             # Groupes = list(Desert = c("TRY_LPP", "TRY_LeafThick", "TRY_SeedMass", "TRY_SSD"),
                             #                Taiga = c("TRY_Wood", "TRY_SLA", "TRY_LDMC", "TRY_LeafArea", "TRY_Lifespan", "TRY_Height", "TRY_CNRatio", "TRY_LeafPhenology")),
                             Site.name = "CWM-traits, ACASP-sensus stricto (non gap-filled)", Legend.position = "none",
                             # Save.plot = "Figures/ACA/Trait/Correlation/PCA/CWM/PCA_CWT_ACASP-ss.pdf", H = 650, W = 650
                             )


    BioCl.pca.pt.sl <- PCA.bioclim(MCWT.clim.PT_sl[,Trait.to.keep], transp_OK = T, Scale.PCA = 3,
                             Cluster.core = "Biome", Ellipse = F, return.pick = T, Legend.position = "none",
                             # Groupes = list(Desert = c("TRY_LPP", "TRY_LeafThick", "TRY_SeedMass", "TRY_SSD"),
                             #                Taiga = c("TRY_Wood", "TRY_SLA", "TRY_LDMC", "TRY_LeafArea", "TRY_Lifespan", "TRY_Height", "TRY_CNRatio", "TRY_LeafPhenology")),
                             Site.name = "CWM-traits, ACASP-sensus largo (non gap-filled)", #Num.facet = "(A)",
                             # Save.plot = "Figures/ACA/Trait/Correlation/PCA/CWM/PCA_CWT_ACA_sl.pdf", H = 650, W = 650
                             )


    BioCl.pca.mv <- PCA.bioclim(MCWT.clim.MV[,Trait.to.keep], transp_OK = T, Scale.PCA = 3, Manu.lim.x = c(8,-9),
                             Cluster.core = "Biome", Ellipse = F, return.pick = T, Legend.position = "none",
                             # Groupes = list(Desert = c("TRY_LPP", "TRY_LeafThick", "TRY_SeedMass", "TRY_SSD"),
                             #                Taiga = c("TRY_Wood", "TRY_SLA", "TRY_LDMC", "TRY_LeafArea", "TRY_Lifespan", "TRY_Height", "TRY_CNRatio", "TRY_LeafPhenology")),
                             Site.name = "CWM-traits, Vegetation plot (non gap-filled)",
                             # Save.plot = "Figures/ACA/Trait/Correlation/PCA/CWM/PCA_CWT_MV.pdf", H = 650, W = 1000
                             )


   BioCl.pca.mv.splot <- PCA.bioclim(MCWT.clim.splot[,c(4,5,7,11,12,13,32)], transp_OK = T, Scale.PCA = 1,
                             Cluster.core = "Biome", Ellipse = F, return.pick = F, Legend.position = "none",
                             # Groupes = list(Desert = c("TRY_LPP", "TRY_LeafThick", "TRY_SeedMass", "TRY_SSD"),
                             #                Taiga = c("TRY_Wood", "TRY_SLA", "TRY_LDMC", "TRY_LeafArea", "TRY_Lifespan", "TRY_Height", "TRY_CNRatio", "TRY_LeafPhenology")),
                             Site.name = "CWM-traits, sPlot traits", Reverse.dim = F,
                             Manu.lim.x = c(-20,20), #c(-2,5), 
                             Manu.lim.y = c(-20,20),#c(-3,5),
                             # Save.plot = "Figures/ACA/Trait/Correlation/PCA/CWM/PCA_CWT_MV_splot.pdf", H = 650, W = 1000
                             )

    pca.CWT.full <- ((BioCl.pca.mv + BioCl.pca.mv.gf + BioCl.pca.pt.ss) /(BioCl.pca.pt.ss.gf + BioCl.pca.pt.sl +  BioCl.pca.pt.sl.gf)) /
     guide_area() + plot_layout(guides = "collect", heights = c(1, 1, 1/10))
    W = 1600
    H = 1000
    ggsave(filename = "Figures/ACA/Trait/Correlation/PCA/CWM/PCA_CWT_full.pdf", pca.CWT.full, width = W*0.026458333, height = H*0.026458333, units = "cm")
    }
  
  #### Test Box's M test ####
  Box.M = F
  if(Box.M == T){
    library(heplots)
    
    Trait.to.keep.Box <- c(4,5,7,11,12,13)
    
    M1 = MCWT.clim.PT_ss_gf[,Trait.to.keep.Box]
    M2 = MCWT.clim.PT_sl_gf[Trait.to.keep.Box]
    M3 = MCWT.clim.MV_gf[Trait.to.keep.Box]
    
    M1$Sampling <- "ACASP-fine"
    M2$Sampling <- "ACASP-coarse"
    M3$Sampling <- "ACAV"
    M <- rbind(M1, M2, M3)
    M <- na.omit(M)
    res <- boxM(M[!names(M) %in% "Sampling"], M[, "Sampling"])
    
    A = summary(res)
    
    #### M calculation ####
    Tab.n = count(M, c("Sampling"))
    Tab.n$nn <- Tab.n$freq -1
    Tab.n <- rbind(Tab.n, c("Pooled", sum(Tab.n$freq), sum(Tab.n$nn)))
    Tab.n$ln.det <- c(res$logDet)
    Tab.n$nn.ln.det <- Tab.n$ln.det * as.numeric(Tab.n$nn)
    M.statistic = Tab.n$nn.ln.det[Tab.n$Sampling == "Pooled"] - sum(Tab.n$nn.ln.det[Tab.n$Sampling != "Pooled"])
    
    # visualize (what is done in the plot method) 
    dets <- res$logDet
    dets <- dets[order(dets)]
    ng <- length(res$logDet)-1
    dotchart(dets, xlab = "log determinant")
    points(dets , 1:length(dets),  
           cex=c(rep(1.5, ng), 2.5), 
           pch=c(rep(16, ng), 15),
           col= c(rep("blue", ng), "red"))
  }
  
  #### Procrustean Co-innertia Analysis ####
  PCoIA = F
  if(PCoIA == T){
    # PCoIA.CWM.PT.norm.gf <- PCoIA.plot(MCWT.clim.PT_ss_gf[,Trait.to.keep], MCWT.clim.PT_ss[,Trait.to.keep],
    #                                    H = 800, W = 1200, transp_OK = T, return.pick = F, Show.outliers = F, Ylim = c(0,1),
    #                                    Save.plot = "Figures/ACA/Trait/Correlation/PCoI/Procrustian_CWM_ACASP_ss_ssgf.pdf")
    
    PCoIA.CWM.PT <- PCoIA.plot(MCWT.clim.PT_sl_gf[,Trait.to.keep], MCWT.clim.PT_ss_gf[,Trait.to.keep],
                               H = 700, W = 1100, transp_OK = T, return.pick = F, Show.outliers = F, Ylim = c(0,0.01),
                               Save.plot = "Figures/ACA/Trait/Correlation/PCoI/Procrustian_CWM_ACASP_ss_sl.pdf"
                               )
    }
  }

#### CWM trait biogeography ####
Biogeo.CWM = F
if(Biogeo.CWM == T){
  #### Import data ####
  if(exists("MCWT.clim.PT_ss") == F){
    library(rgdal)
    library(raster)
    MCWT.clim.PT_ss <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_clean.Rds")
    MCWT.clim.PT_ss_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_ss_gf_clean.Rds")
    MCWT.clim.PT_sl <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_clean.Rds")
    MCWT.clim.PT_sl_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_PT_sl_gf_clean.Rds")
    MCWT.clim.ACAP.full <- readRDS( "Resultats/ACA/Traits/CWT_ACA_PT_all.Rds")
    MCWT.clim.MV <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_clean.Rds")
    MCWT.clim.MV_gf <- readRDS("Resultats/ACA/Traits/CWT_ACA_MV_gf_clean.Rds")
    MCWT.clim.MV.tot <- readRDS("Resultats/ACA/Traits/CWT_ACA_ACAV_all.Rds")
    MCWT.clim.splot <- data.frame(readRDS("Resultats/ACA/Traits/sPlot_ACA_CWM_clim_clean.Rds"))
    MPS.ACA.Cord <- readRDS("Resultats/World_DB/Pollen/Merging_DB/DB13549/DBACA_Coord.Rds")
    MV.Cord <- readRDS("Resultats/ACA/Vegetation/MV_ACA_cord.Rds")
  }
  
  #### Import data SIG ####
  if(exists("ACA.bo.co") == F){
    Path.ACA.border = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/Projets/ACA/Borders_ACA/Extern_border/Out_border_ACA.shp"
    # Path.ACA.border = "/media/lucas.dugerdil/Climacoptera/Documents/Recherche/SIG/Projets/ACA/Borders_ACA/Extern_border/Out_border_ACA.shp"
    ACA.bo = readOGR(Path.ACA.border)
    proj4string(ACA.bo) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    # proj4string(ACA.bo) <- CRS("+proj=lcc +lat_0=30 +lon_0=95 +lat_1=15 +lat_2=65 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    ACA.bo.proj = fortify(ACA.bo)
    
    Path.ACA.border.co = "/media/lucas.dugerdil/Extreme SSD/Documents/Recherche/SIG/Projets/ACA/Borders_ACA/Full/ACA_border.shp"
    ACA.bo.co = readOGR(Path.ACA.border.co)
    proj4string(ACA.bo.co) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    ACA.bo.co.proj = fortify(ACA.bo.co)
  }
  
  #### Applications #### 
  Map.fig6 = T
  if(Map.fig6 == T){
    Map2 <- Map.biogeo.CWM(MCWT = MCWT.clim.PT_ss_gf, Type1 = "ACASP-fine",
                          MCWT2 = MCWT.clim.MV_gf, Type2 = "ACA-vegetation",
                          Strip.lab = F, Hex.size = 40, Show.diff = F,
                          Select.trait = c("TRY_LeafArea", "TRY_SSD", "TRY_Height"),
                          W = 1200, H = 700, Save.plot = "Figures/ACA/Trait/Map/CWM_map/Map_ACA_3CWT_best.pdf"
                          )
    
    # Map2 <- Map.biogeo.CWM(MCWT = MCWT.clim.MV_gf, Type1 = "ACA-vegetation",
    #                         MCWT2 = MCWT.clim.PT_ss_gf, Type2 = "ACASP-fine",
    #                         Strip.lab = F, Hex.size = 40, Show.diff = T,
    #                         Select.trait = c("MPWAQ_chel"),
    #                         W = 600, H = 700, Save.plot = "Figures/ACA/Trait/Map/CWM_map/Map_ACA_clim_diff.pdf"
    #                         )
    }
  
  #### Test maps ####
  # Map1 <- Map.biogeo.CWM(MCWT = MCWT.clim.MV_gf, Type1 = "ACA-vegetation",
  #                         MCWT2 = MCWT.clim.PT_ss_gf, Type2 = "Pollen",
  #                         Strip.lab = T,
  #                         Select.trait = c("TRY_LeafArea", "TRY_SSD", "TRY_LeafN", "TRY_Height", "TRY_SLA", "TRY_SeedMass"),
  #                         W = 2000, H = 750, Save.plot = "Figures/ACA/Trait/Map/CWM_map/Map_ACA_6CWT.pdf"
  #                         )
  # 
  # Map3 <- Map.biogeo.CWM(MCWT = MCWT.clim.MV_gf, Type1 = "ACA-vegetation",
  #                         MCWT2 = MCWT.clim.PT_ss_gf, Type2 = "Pollen",
  #                         Strip.lab = F, Crop.zone = c("China","Mongolia"), Hex.size = 30,
  #                         # Strip.lab = F, Crop.zone = c("Tajikistan","Kyrgyzstan"), Hex.size = 30,
  #                         Select.trait = c("TRY_LeafArea", "TRY_SSD", "TRY_LeafN", "TRY_Height", "TRY_SLA", "TRY_SeedMass"),
  #                         W = 2000, H = 700, Save.plot = "Figures/ACA/Trait/Map/CWM_map/Map_ACA_6CWT_Mong.pdf"
  #                         )
  
  # Map4 <- Map.biogeo.CWM(MCWT = MCWT.clim.MV_gf, Type1 = "ACA-Vegetation",
  #                         MCWT2 = MCWT.clim.PT_ss_gf, Type2 = "ACASP-ss",
  #                         Strip.lab = T, Hex.size = 60,
  #                         # Select.trait = c("TRY_LeafArea", "TRY_Height", "TRY_SLA"),
  #                         Select.trait = c("TRY_LeafArea"),
  #                         # W = 1100, H = 750, Save.plot = "Figures/ACA/Trait/Map/Map_ACA_best_CWT.pdf"
  #                         W = 600, H = 1000, Save.plot = "Figures/ACA/Trait/Map/CWM_map/Map_ACA_LA_d.pdf"
  #                         )
  # 
  # 
  # Map5 <- Map.biogeo.CWM(MCWT = MCWT.clim.PT_sl_gf, Type1 = "ACASP-sl",
  #                        MCWT2 = MCWT.clim.PT_ss_gf, Type2 = "ACASP-ss",
  #                        Strip.lab = F, Hex.size = 60,
  #                        Select.trait = c("TRY_LeafArea", "TRY_Height", "TRY_SLA"),
  #                        W = 1200, H = 800, Save.plot = "Figures/ACA/Trait/Map/CWM_map/Map_ACA_pollen_CWT.pdf"
  # )
  #### SI ####
  Figure.SI = F
  if(Figure.SI == T){
    Map.SI1 <- Map.biogeo.CWM(MCWT = MCWT.clim.MV_gf, Type1 = "ACA-vegetation",
                           MCWT2 = MCWT.clim.PT_ss_gf, Type2 = "ACASP-fine", 
                           Strip.lab = F, Hex.size = 40, Show.diff = F, Vertical = T, Leg.pos = "bottom",
                           Select.trait = c("TRY_LeafArea", "TRY_LeafN", "TRY_SeedMass", "TRY_SLA", "TRY_SSD", "TRY_Height"),
                           W = 700, H = 1500, Save.plot = "Figures/ACA/Trait/Map/CWM_map/Map_ACA_all_CWM_ss.pdf"
                           ) 
    Map.SI2 <- Map.biogeo.CWM(MCWT = MCWT.clim.MV_gf, Type1 = "ACA-vegetation",
                           MCWT2 = MCWT.clim.PT_sl_gf, Type2 = "ACASP-coarse",
                           Strip.lab = F, Hex.size = 40, Show.diff = F, Vertical = T, Leg.pos = "bottom",
                           Select.trait = c("TRY_LeafArea", "TRY_LeafN", "TRY_SeedMass", "TRY_SLA", "TRY_SSD", "TRY_Height"),
                           W = 700, H = 1500, Save.plot = "Figures/ACA/Trait/Map/CWM_map/Map_ACA_all_CWM_sl.pdf"
                           )


    # MAP.CWM.SI <- Map.SI1 + Map.SI2 #+ plot_layout(widths = c(2,1))
    # 
    # W = 1400
    # H = 2400
    # Save.plot = "Figures/ACA/Trait/Map/CWM_map/Map_ACA_full.pdf"
    # ggsave(filename = Save.plot, MAP.CWM.SI, width = W*0.026458333, height = H*0.026458333, units = "cm")
  }
}
