setwd("C:/Users/Gamze Yazgeldi/Desktop/DrugAnalysisFiles/Drugbank.ca_AllDownloadableFiles/")
drug_links <- read.csv2("drug links.csv", header=TRUE, sep=",")
target_links <-read.csv2("target links.csv",header=TRUE,sep=",")
library(tidyverse)
target_links <- target_links %>% separate_rows(Drug.IDs)

setwd("C:/Users/Gamze Yazgeldi/Desktop//CanonicalCorrelationAnalysis")
func_hostbased <- function(phi_data, ga_data){
  phi_data <- read.csv2(file.choose(),header=TRUE,sep=';')     # to choose interested phi data
  ga_data <- read.csv2(file.choose(),header=TRUE,sep=',')   # to choose interested degree & btw centrality data
  
  phi_data <- phi_data[!duplicated(phi_data[,"Host.Uniprot.ID"]),]
  pro_drug <- merge(phi_data, target_links, by.x =c("Host.Uniprot.ID"), by.y=c("UniProt.ID"), all.x = FALSE, all.y = FALSE)
  pro_top <- merge(ga_data, pro_drug, by="Host.Uniprot.ID",all = FALSE)
  
  drug_other <- merge(pro_top,drug_links,by.x = "Drug.IDs",by.y="DrugBank.ID",all=FALSE)
  drug_analysis <- data.frame(host_protein_ids=drug_other$Host.Uniprot.ID, host_protein_name=drug_other$Name.x,
                              pdb_id=drug_other$PDB.ID, species=drug_other$Species, drug_ids=drug_other$Drug.IDs,
                              drug_name=drug_other$Name.y, drug_type=drug_other$Drug.Type, kegg_id=drug_other$KEGG.Drug.ID,
                              chebi_id=drug_other$ChEBI.ID, degree=drug_other$Degree, 
                              betweenness_centrality=drug_other$Betweenness.Centrality)
}
result_hostbased <- func_hostbased(phi_data,ga_data);
write.csv2(result_hostbased, file="Host Based Drug Analysis Result.csv")


# PATHOGEN BASED:
func_pathogenbased <- function(phi_data){
  phi_data <- read.csv2(file.choose(),header=TRUE,sep=';')     # to choose interested phi data
  phi_data <- phi_data[!duplicated(phi_data[,"Uniprot.ID"]),]
  pro_drug <- merge(phi_data, target_links, by.x= "Uniprot.ID", by.y = "UniProt.ID", all=FALSE)
  drug_other <- merge(pro_drug, drug_links,by.x = "Drug.IDs",by.y="DrugBank.ID",all=FALSE)
  drug_analysis <- data.frame(pathogen_protein_ids=drug_other$Uniprot.ID, pathogen_protein_name=drug_other$Name.x,
                              pdb_id=drug_other$PDB.ID, species=drug_other$Species, drug_ids=drug_other$Drug.IDs,
                              drug_name=drug_other$Name.y, drug_type=drug_other$Drug.Type, kegg_id=drug_other$KEGG.Drug.ID,
                              chebi_id=drug_other$ChEBI.ID)
}
result_pathogenbased <- func_pathogenbased(phi_data)
write.csv2(result_pathogenbased, file="Pathogen Based Drug Analysis Result.csv")


#############################################################################################################
# DRUG-PROTEIN INTERACTIONS NETWORK ANALYSIS
# Interested pathogen-host interaction data's host protein analysis result is retrieved by users from online.
# Inputs are t_analysis and interaction file.
ga_data <- read.table("t_analysis.txt",sep=",",header = FALSE)
phi_data <- read.csv2(file.choose(),header=TRUE,sep=';')     

# Finding drug targeted human proteins:
mismatch_data <- as.data.frame(ga_data$V1[!(ga_data$V1 %in% result_hostbased$host_protein_ids)])
colnames(mismatch_data) <- "V1"
match_data <- as.data.frame(ga_data$V1[(ga_data$V1 %in% result_hostbased$host_protein_ids)])
colnames(match_data) <- "V1"
drug_targets <- merge(ga_data,match_data,by.x = "V1",by.y="V1")

# Finding non-drug targeted human proteins:
mismatch_data <- phi_data$Host.Uniprot.ID[!(phi_data$Host.Uniprot.ID %in% result_hostbased$host_protein_ids)]
mismatch_data <- as.data.frame(unique(mismatch_data))
non_drug_targets <- merge(ga_data, mismatch_data, by.x="V1", by.y="unique(mismatch_data)")

# Host proteins in these two datasets can be analyse wrt gene ontology and pathway enrichment.
write.csv2(non_drug_targets, file="Non-drug target host proteins in the SARS-CoV-2.csv")
write.csv2(drug_targets, file="Drug target host proteins in the SARS-CoV-2.csv")

mean(non_drug_targets$V2)
mean(drug_targets$V2)    

mean(non_drug_targets$V3)
mean(drug_targets$V3)