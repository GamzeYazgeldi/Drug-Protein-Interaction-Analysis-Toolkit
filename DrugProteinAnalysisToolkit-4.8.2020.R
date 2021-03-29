# TODO: Find host-based and pathogen-based drugs which are targeting host and pathogen proteins within interested 
# Pathogen-Host Interaction (PHI) data.

# INPUT FILES: DRUG_LINKS, TARGET_LINKS, HIV_PHI AND HIV_GA.

# Relevant drug and drug target informations were retrieved from DrugBank.
drug_links <- read.csv2("drug links.csv", header=TRUE, sep=",");
target_links <-read.csv2("target links.csv",header=TRUE,sep=",");
# Tidyverse library was required for seperation DrugIDs column by ';'.
library(tidyverse)
target_links <- target_links %>% separate_rows(Drug.IDs)

# About input files:
# Interested HIV-Human PHI data were downloaded from PHISTO. 
# Topological values information of host proteins were downloaded from PHISTO througly Graph Analysis option.

# Finding host-based drugs:
func_hostbased <- function(phi_data, ga_data){
  phi_data <- read.csv2(file.choose(),header=TRUE,sep=',')     # to choose interested phi data
  ga_data <- read.csv2(file.choose(),header=TRUE,sep=',')   # to choose interested degree & btw centrality data
  
  phi_data <- phi_data[!duplicated(phi_data[,"Host.Uniprot.ID"]),]
  
  pro_drug <- merge(phi_data, target_links, by.x =c("Host.Uniprot.ID"), by.y=c("UniProt.ID"), all.x = FALSE, all.y = FALSE)
  length(unique(pro_drug$Drug.IDs))
  
  pro_top <- merge(ga_data, pro_drug, by="Host.Uniprot.ID",all = FALSE)
  
  # Other drug IDs in different sources are retrieving.
  drug_other <- merge(pro_top,drug_links,by.x = "Drug.IDs",by.y="DrugBank.ID",all=FALSE)
  
  drug_analysis <- data.frame(host_protein_ids=drug_other$Host.Uniprot.ID, host_protein_name=drug_other$Name.x,
                              pdb_id=drug_other$PDB.ID, species=drug_other$Species, drug_ids=drug_other$Drug.IDs,
                              drug_name=drug_other$Name.y, drug_type=drug_other$Drug.Type, kegg_id=drug_other$KEGG.Drug.ID,
                              chebi_id=drug_other$ChEBI.ID, degree=drug_other$Degree, 
                              betweenness_centrality=drug_other$Betweenness.Centrality)
}
result_hostbased <- func_hostbased(phi_data,ga_data);
write.csv2(result_hostbased, file="Host Based Drug Analysis Result.csv")


# Finding host-based drugs:
func_pathogenbased <- function(phi_data){
  phi_data <- read.csv2(file.choose(),header=TRUE,sep=',')     # to choose interested phi data
  
  phi_data <- phi_data[!duplicated(phi_data[,"Uniprot.ID"]),]
  
  pro_drug <- merge(phi_data, target_links, by.x= "Uniprot.ID", by.y = "UniProt.ID", all=FALSE)
  length(unique(pro_drug$Drug.IDs))
  
  # Other drug IDs in different sources are retrieving.
  drug_other <- merge(pro_drug, drug_links,by.x = "Drug.IDs",by.y="DrugBank.ID",all=FALSE)
  
  drug_analysis <- data.frame(pathogen_protein_ids=drug_other$Uniprot.ID, pathogen_protein_name=drug_other$Name.x,
                              pdb_id=drug_other$PDB.ID, species=drug_other$Species, drug_ids=drug_other$Drug.IDs,
                              drug_name=drug_other$Name.y, drug_type=drug_other$Drug.Type, kegg_id=drug_other$KEGG.Drug.ID,
                              chebi_id=drug_other$ChEBI.ID)
}
result_pathogenbased <- func_pathogenbased(phi_data)
write.csv2(result_pathogenbased, file="Pathogen Based Drug Analysis Result.csv")


#####################################################################################################################
# STUDY PART II: DRUG-PROTEIN INTERACTIONS NETWORK ANALYSIS
# TODO: Compare drug targets and non-drug targets by topological parameters analysis.

# INPUT FILES: HIV_GA.

# About input files:
# Interested HIV-Human PHI data were downloaded from PHISTO. 
# Topological values information of host proteins were downloaded from PHISTO througly Graph Analysis option.

ga_data <- read.csv(file.choose(), header=TRUE, sep=",")
mismatch_data <- as.data.frame(ga_data$Host.Uniprot.ID[!(ga_data$Host.Uniprot.ID %in% result_hostbased$host_protein_ids)])
non_drug_targets <- merge(ga_data, mismatch_data, by.x = "Host.Uniprot.ID", by.y ="ga_data$Host.Uniprot.ID[!(ga_data$Host.Uniprot.ID %in% result_hostbased$host_protein_ids)]")
non_drug_targets <- non_drug_targets[, 1:3]
drug_targets <- unique(data.frame(Host.UniProt.ID=result_hostbased$host_protein_ids, Degree=result_hostbased$degree,
                                  Betweenness.Centrality=result_hostbased$betweenness_centrality))

write.csv2(non_drug_targets, file="Non-drug target host proteins.csv")
write.csv2(drug_targets, file="Drug target host proteins.csv")

mean(non_drug_targets$Degree); 
mean(drug_targets$Degree) #degree
mean(non_drug_targets$Betweenness.Centrality);
#mean(drug_targets$Betweenness.Centrality)
