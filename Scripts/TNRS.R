#### Function ####
TNRS.taxa.accepted <- function(Taxa.vector, Taxonstand.test){
  #### Settings ####
  if(missing(Taxonstand.test)){Taxonstand.test = F}
  #rm(list=ls())
  # Base URL for TNRS api
  url = "https://tnrsapi.xyz/tnrs_api.php"	# Production on paramo
  #url = "http://vegbiendev.nceas.ucsb.edu:8975/tnrs_api.php"  # Dev on vegbiendev 
  # library(RCurl) # API requests
  library(httr)		# API requests
  library(jsonlite) # JSON coding/decoding
  # sources <- "tropicos,tpl,usda"	        # Taxonomic sources
  sources <- "tropicos,wfo,wcvp,usda"	# Taxonomic sources
  class <- "tropicos"										# Family classification
  mode <- "resolve"										# Processing mode
  matches <- "best"											# Return best match only "all" for all the matches
  opts <- data.frame(c(sources),c(class), c(mode), c(matches))   # Convert the options to data frame and then JSON
  names(opts) <- c("sources", "class", "mode", "matches")
  opts_json <-  jsonlite::toJSON(opts)
  opts_json <- gsub('\\[','',opts_json)
  opts_json <- gsub('\\]','',opts_json)
  
  #### Vector preparation ####
  Taxa.vector <- unique(Taxa.vector)
  if(length(Taxa.vector[duplicated(Taxa.vector)]>0)){print(paste("Attention ! The following taxa are duplicated in the Taxa vector:", Taxa.vector[duplicated(Taxa.vector)]))}
  data <- data.frame(V1 = seq(1:length(Taxa.vector)), V2 = Taxa.vector)
  data_json <- jsonlite::toJSON(unname(data))
  input_json <- paste0('{"opts":', opts_json, ',"data":', data_json, '}' )# Combine the options and data into single JSON object
  headers <- list('Accept' = 'application/json', 'Content-Type' = 'application/json', 'charset' = 'UTF-8')# Construct the request
  
  #### Send the API request ####
  # results_json <- postForm(url, .opts=list(postfields= input_json, httpheader=headers))
  results_json <- POST(url = url,
                       add_headers('Content-Type' = 'application/json'),
                       add_headers('Accept' = 'application/json'),
                       add_headers('charset' = 'UTF-8'),
                       body = input_json,
                       encode = "json")
  results_raw <- fromJSON(rawToChar(results_json$content)) 
  results <- as.data.frame(results_raw)
  print(results)
  # results <-  jsonlite::fromJSON(results_json) # Convert JSON results to a data frame
  results$match.score <- format(round(as.numeric(results$Overall_score),2), nsmall=2)
  # print(results)
  results.clean <- results[c("Name_submitted","Accepted_name", "match.score",  "Unmatched_terms", "Accepted_family")]
  # results.clean <- results
  results.clean[results$Taxonomic_status == "No opinion", c(2,5)] <- results[results$Taxonomic_status == "No opinion",c("Name_matched","Name_matched_accepted_family")]    
  
  #### Test The Plant List taxo ref. (Cayuela et al., 2012) ####
  if(Taxonstand.test == T){
    library(Taxonstand)
    r1 <- TPL(results.clean$Name_submitted, corr = TRUE)
    print(r1)
    results.clean[["Accepted_name_TXSD"]] <- paste(r1$New.Genus, r1$New.Hybrid.marker, r1$New.Species, sep = " ")
    results.clean[["Accepted_name_TXSD"]] <- gsub("  ", " ", results.clean[["Accepted_name_TXSD"]])
    results.clean[["Accepted_family"]] <- r1$Family
    }
  return(results.clean)
  }

#### Applications ####
# Add.list <- TNRS.taxa.accepted(Taxa.vector = TL.merge[0:20,"species"], Taxonstand.test = T)

Index.pollen.Reille = F
if(Index.pollen.Reille == T){
  # Reille <- data.frame(read.csv(file="Import/World_DB/Taxonomie/Index_reille.csv",sep=",",dec=".",header = T, stringsAsFactors = FALSE))
  Reille <- data.frame(read.csv(file="Import/World_DB/Taxonomie/Index_reille_p2.csv",sep=",",dec=".",header = T, stringsAsFactors = FALSE))
  Reille$species <- paste(Reille$Genus, Reille$Species, Reille$Subsp.)
  Reille$species <- gsub(" $","",Reille$species, perl=T)
  
  Reille.check <- TNRS.taxa.accepted(Taxa.vector = Reille[,"species"], Taxonstand.test = F)
  names(Reille.check)[2] <- "species"
  Reille.full <- merge(Reille, Reille.check, by = "species") 
  Reille.full <- Reille.full[,c(2:9,11)]
  Reille.full$Acc_Genus <- gsub("\\s.*", "", Reille.full$Accepted_name)
  Reille.full$Acc_Spe <- sub(".*? ", "", Reille.full$Accepted_name)
  # Reille.full$Acc_Subsp <- gsub(".* subsp", "", Reille.full$Acc_Spe)
  # Reille.full$Acc_Spe <- gsub("\\s.*", "", Reille.full$Acc_Spe)
  
  write.table(Reille.full, file = "Resultats/World_DB/Taxonomie/Index_Reille_checked_p2.csv", row.names=T, col.names=NA, sep=",", dec = ".")
  }

Check.reille.p2 = F
if(Check.reille.p2 == T){
  Reille <- read.csv(file="Import/World_DB/Taxonomie/Index_Reille_setdiff.csv",sep=",",dec=".",header = F, stringsAsFactors = FALSE)
  Reille.check <- TNRS.taxa.accepted(Taxa.vector = Reille$V1, Taxonstand.test = F)
  write.table(Reille.check, file = "Resultats/World_DB/Taxonomie/Index_Reille_setdiff_check.csv", row.names=T, col.names=NA, sep=",", dec = ".")
  
  Reille <- read.csv(file="Import/World_DB/Taxonomie/Index_Reille_im_setdiff.csv",sep=",",dec=".",header = F, stringsAsFactors = FALSE)
  Reille$ V1 <- substr(Reille$V1, start = 2, nchar(Reille$V1))
  Reille.check <- TNRS.taxa.accepted(Taxa.vector = Reille$V1, Taxonstand.test = F)
  write.table(Reille.check, file = "Resultats/World_DB/Taxonomie/Index_Reille_im_setdiff_check.csv", row.names=T, col.names=NA, sep=",", dec = ".")
  
}

Checklist.taj = F
if(Checklist.taj == T){
  CLT <- data.frame(read.csv(file="Import/Uzbekistan/Vegetation/Checklist_Tajikistan.csv",sep=",",dec=".",header = T, stringsAsFactors = FALSE))
  CLT <- CLT[-nrow(CLT),]
  CLT$ID <- as.integer(CLT$ID)
  CLT$Species <- gsub("\\[.*","",CLT$Species.name.according.to.Cherepanov.1995.with.basic.synonym.s.)
  CLT.check <- TNRS.taxa.accepted(Taxa.vector = CLT[,"Species"], Taxonstand.test = F)
  names(CLT.check)[2] <- "Species"
  CLT.check <- merge(CLT, CLT.check, by = "Species") 
  write.table(CLT.check, file = "Resultats/World_DB/Taxonomie/Checklist_Tajikistan_TNRS.csv", row.names=T, col.names=NA, sep=",", dec = ".")
  
}