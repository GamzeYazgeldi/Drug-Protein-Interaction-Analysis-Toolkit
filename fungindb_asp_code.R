# FungiNDB project for Aspergillus fuminatus's gene-gene interactions
# by Gamze Yazgeldi in June-July, 2019
library(data.table)

# PROCESSING OF OPTIMIZED DATA AND PRIOR KNOWLEDGE DATA
setwd("/home/gamze/Documents");
opt_data <- read.csv("NETWORKS_Data_Penalthy_0.25.tab",sep="\t"); dim(opt_data) #32431 row, 51 columns
opt_data <- opt_data[, c(1,2,39)] #selected min RSS size factor that is X80 
write.table(opt_data, file="optimized_data.csv", row.names = FALSE,col.names = FALSE)

pr_data <- read.csv("GS_v2.csv", sep="\t", header = FALSE, stringsAsFactors = FALSE)
pr_data[1:2] <- lapply(pr_data, trimws)

plot(as.numeric(unlist(opt_data[,3])),xlab="index",ylab="weight",main="optimized_data")

###############################################################################################################################################################################################################
# CONSTRUCTION DATA TABLES 
# GENE TABLES:
gene_names_opt_data <-as.data.frame(union(opt_data[,1],opt_data[,2]), stringsAsFactors = FALSE) # Unique genes between source and target genes. There are unique 246 genes.
gene_names_pr_data <- as.data.frame(union(pr_data[,1],pr_data[,2]), stringsAsFactors = FALSE) # Unique genes between source and target genes. There are unique 214 genes.
gene_names_both <- as.data.frame(union(gene_names_pr_data[,1],gene_names_opt_data[,1]))
# Count of unique gene names belongs to Candida albicans will be added into new unique gene names belongs to Aspergillus fuminagatus.
# So, seq((nrow(For_candida_unique_gene_names)+1), (nrow(for_aspergillus_unique_gene_names)+nrow(For_candida_unique_gene_names)), by=1)
# k <- nrow(gene_names_for_candida)
k <- 6214
unid <- seq((k+1), (nrow(gene_names_both)+k), by=1)
gene_both <- data.frame(unid=unid, gene_names=gene_names_both)     # Another columns will be retrieved from related genome database. Interested strain is ATCC 46645 (NCPF 2109) and genome is Af293 (NCPF 7367).
write.table(gene_both, file = "both_unidID_gene_table.csv", row.names = FALSE,col.names = FALSE)

indices_opt <- match(gene_names_opt_data[,1], gene_names_both[,1])
gene_opt <- gene_both[indices_opt,]
write.table(gene_opt, file="optimized_data_unidID_gene_table.csv",col.names = FALSE)

indices_pr <- match(gene_names_pr_data[,1], gene_names_both[,1])
gene_pr <- gene_both[indices_pr,]
write.table(gene_opt, file="priorkn_data_unidID_gene_table.csv", col.names = FALSE)

intersect_genes <- intersect(gene_names_pr_data[,1],gene_names_opt_data[,1]) #203 genes are common...
diff_genes <- as.data.frame(setdiff(gene_names_pr_data[,1], gene_names_opt_data[,1]), stringsAsFactors = FALSE) # 13 different genes btw opt_data and pr_data were assigned
#subs <- subset(gene_opt, !(gene_names %in% gene_pr$gene_names))

############################################################################################################################################################################################################### 
# ANNOTATION TABLE:
# This annotation table was retrieved from ftp://ftp.ensemblgenomes.org/pub/release-38/fungi/gff3/aspergillus_fumigatusa1163/
# Here is annotation data and columns are seperated by tab. (retrieved from ensemblegenomes)
annotation_data <- readLines(con = "Aspergillus_fumigatusa1163.CADRE.38.gff3")
annotation_data_t <- grep(pattern = "^\\#\\#|^\\#\\!", x = annotation_data, value = TRUE, invert = TRUE)
annotation_data_t <- paste0(annotation_data_t, collapse = "\n")
annotation_data <- fread(input = annotation_data_t, sep="\t",header=FALSE)

separated_anndata <- strsplit(as.character(annotation_data$V9),";")  #interested part is here contains CADAFUBG* and AFUB* present in 9th column
gene_id <- sapply(separated_anndata, grep, pattern ="ID=gene:CADAFUBG[0-9]+", value = TRUE)   #Value gives your returned selected pattern.
gene_name <- sapply(separated_anndata, grep, pattern ="Name=AFUB\\_[0-9]+", value = TRUE)

gene_id <- gsub("ID=gene:", "",gene_id)
gene_name <- gsub("Name=", "",gene_name)  
length(gene_name)     #79441 gene names contain that their geneid, exonid and transcriptid
length(unique(gene_name))    #10097 genes
# This table contains gene_id corresponding to gene_name of Asp genome annotation. 
matched_table <- data.frame(gene_id, gene_name, stringsAsFactors = TRUE) 

f <- unique(matched_table[,1]);     #there is unique 10161 gene_id 
f1 <- unique(matched_table[,2]);     #there is unique 10097 gene_name
check3 <- unique(gene_pr[,2]) %in% unique(matched_table[,1]); all(check3)

indices_matched_1 <- match(gene_both[,2], matched_table[,1])    #indices are matched interested 260 gene_ids in annotation table
gene_name_table <-matched_table[indices_matched_1, 2]    #gene names (CADAFUBG*) are corresponding to optimized data gene_ids (AFUB_*),

matched_CAD_AFU <- data.frame(gene_both[,2],gene_name_table)
write.table(matched_CAD_AFU, file = "matched_259_geneid_genename_table.csv", row.names = FALSE,col.names = FALSE) 

# This annotation table was retrieved from: 
# http://www.aspergillusgenome.org/download/chromosomal_feature_files/A_fumigatus_Af293/A_fumigatus_Af293_current_chromosomal_feature.tab
ann_data <- read.csv2("asp_genome_annotation",sep="\t",header=FALSE)     #first explanations were removed to seperate columns.
gene_name_annt <- lapply(ann_data$V3, grep, pattern ="AFUB\\_[0-9]+", value = TRUE)       #contain AFUB_* and other relevant row elements
gene_name_annt_2 <- sapply(as.character(gene_name_annt), strsplit, split= "|", fixed=TRUE)      #splitted by |
gene_name_annt_3 <- lapply(gene_name_annt_2, grep, pattern ="AFUB\\_[0-9]+", value = TRUE)      #only contain AFUB_*
class(gene_name_annt_3)
a <- sapply(gene_name_annt_3, "length<-", max(lengths(gene_name_annt_3)))        #now, a has same length for each list elements while converting to data frame
a <- as.data.frame(t(a))    #transposon of a
indices_matched_2 <- sapply(a, function(col) match(trimws(gene_name_table), as.character(col)))      #whole matching elements' indices
num_nonna <- apply(indices_matched_2, 1, function(row) sum(!is.na(row)))      #check which elements containts only NA. These indices are 31, 71, and 205. 
table(num_nonna) 
ind_match_vec <- rep(NA, length(num_nonna))      #create new vector and then, will assign our indices into new vector
ind_match_vec[num_nonna == 1] <- unlist(apply(indices_matched_2, 1, function(row) which(!is.na(row))))     #i have no idea what we did in here. loooook!!
ind_match_vec_merge <- sapply(1:length(ind_match_vec), function(n) indices_matched_2[n, ind_match_vec[n]])

gene_symbol <- ann_data[ind_match_vec_merge, 18]
gene_description <- ann_data[ind_match_vec_merge, 11]
gene_systematic <- ann_data[ind_match_vec_merge, 1]
gene_alias <- ann_data[ind_match_vec_merge, 9]
chromosome <- ann_data[ind_match_vec_merge, 5]
start_pos <- ann_data[ind_match_vec_merge, 6]   #transcript start and end position
end_pos <- ann_data[ind_match_vec_merge,7]  
strand_type <- ann_data[ind_match_vec_merge,8]  #strand type whether sense (+ or watson/w) or antisense (- or crick/c)
annotation_table_1 <- data.frame(unid=unid, gene_id=gene_both[,2], gene_symbol=gene_symbol, gene_systematic=gene_systematic, 
                                 gene_alias=gene_alias, gene_name=gene_name_table, chromosome=chromosome, strand=strand_type, start_pos=start_pos, end_pos=end_pos, description=gene_description)       
write.table(annotation_table_1, file = "annotation_table_257genes.csv", row.names = FALSE, col.names = FALSE)
a <- duplicated(annotation_table_1[,4]); all(a); annotation_table_1[a,1]  

###############################################################################################################################################################################################################
# INTERACTION TABLES:
int_type <- 5
nodeA_opt <- (match(opt_data[,1],gene_both[,2])+k)
nodeB_opt <- (match(opt_data[,2],gene_both[,2])+k)
nodeA_pr <- (match(pr_data[,1],gene_both[,2])+k)
nodeB_pr <- (match(pr_data[,2],gene_both[,2])+k)
int_opt <- data.frame(nodeA=nodeA_opt,nodeB=nodeB_opt,int_type)
int_pr <- data.frame(nodeA=nodeA_pr,nodeB=nodeB_pr,int_type)
int_opt_pr <- rbind(int_opt,int_pr)
int_opt_pr <- int_opt_pr[!duplicated(int_opt_pr[,c(1,2)]),]


# l is number of interaction betweeen candida albicans genes. 
# l <- nrow(interaction between candida albicans genes) or max(interaction between candida albicans genes). l is 91894
l <- 91894
unid2 <- seq(l+1, (nrow(int_opt_pr)+l), by=1)
int_opt_pr_df <- data.frame(unid=unid2, nodeA=int_opt_pr$nodeA, nodeB=int_opt_pr$nodeB, int_type)
write.table(int_opt_pr_df, file = "both_data_interaction_table.csv", row.names = FALSE,col.names = FALSE)

###############################################################################################################################################################################################################
# NETWORK TABLE:
networkid <- c("3","4")
network_name <- c("optimized data network","prior knowledge data network")
node_number <- c(nrow(gene_opt), nrow(gene_pr))
edge_number <- c(nrow(opt_data), nrow(pr_data))
is.prior <- c('0','1')
description <- c("in silico prediction, regression model, least-angle regression (LARS), no prior knowledge, CV, contains all inferred edges",
                 "prior knowledge, gold standard, text mining")
network <- data.frame(unid=networkid, network_name=network_name, network_name_description=description, species_id="9",provider_name='_', node_number=node_number, edge_number=edge_number, is.prior, type='NULL', predict_method_id='NULL',data_create="10.7.19")
write.table(network, file = "network_table.csv", row.names = FALSE, quote=FALSE,col.names = FALSE, sep = "\t")

###############################################################################################################################################################################################################
# NETWORK EDGE TABLE:
edge_type <- int_type  
interaction_id <- int_opt_pr_df$unid
nodeA <- int_opt_pr_df$nodeA
nodeB <- int_opt_pr_df$nodeB

# m is number of edge/interaction betweeen candida albicans genes. 
# m <- nrow(interaction/edge between candida albicans genes) or max(interaction/edge between candida albicans genes). m is 135985
m <- 135985 
unid4 <- seq(m+1, (length(interaction_id)+m), by=1)      # SPECIFIC FOR EACH INTERACTION - EDGE ID  

network_id <- c( (rep('3',nrow(int_opt_pr_df))) , (rep('4',nrow(int_opt_pr_df))))  
is_prior <- c( (rep('0',nrow(int_opt_pr_df))) , (rep('1',nrow(int_opt_pr_df)))) 
networkedge <- data.frame(unid=unid4, network_id, interaction_id, nodeA, nodeB, edge_type, is_prior, references='NULL', is_homo='0', homo_species='NULL',homo_nodeA='NULL', homo_nodeB='NULL')
write.table(networkedge, file = "networkedge_table.csv", row.names = FALSE, quote=FALSE, col.names = FALSE)

###############################################################################################################################################################################################################
# NETWORK EDGE ATTRIBUTE TABLE:
indices_att <- which(opt_data$X80!=0)     #non-zero weights' indices in the specified parameter column
int_opt_att <- opt_data[indices_att,c(1,2,3)]   # gene-gene interaction data which weights are not equal to 0
a <- match(int_opt_att$source, gene_both[,2])     # sifira esit olmayanlarin sourcelarinin gene unidleri
b <- match(int_opt_att$target,gene_both[,2])
df_att <- data.frame(source=gene_both[a,1],target=gene_both[b,1],weight=opt_data[indices_att,3])      #gene-gene interaction table (!=0) with unid ids
  
  
ind_int_opt <- rep(0, length(df_att$source))
for (n in 1:length(df_att$source)){
  ind_int_opt[n] <- which(df_att$source[n]==int_opt$nodeA & df_att$target[n]==int_opt$nodeB)
}

inter_id_att <- int_opt[ind_int_opt,1]
df_att <- data.frame(unid=inter_id_att, source=gene_both[a,1],target=gene_both[b,1], weight=opt_data[indices_att,3])      #gene-gene interaction table (!=0) with unid ids

# n shows weight/correlation btw each asp genes. n is 87081 due to there is 87081 edge attr id.
n <- 87081
unid5 <- seq(n+1, (length(inter_id_att)+n), by=1)  # Weights must not be zero. Thats why, weights are belongs to prior knowledge data.
network_edge_attr <- data.frame(unid=unid5, inter_id=inter_id_att, weight=df_att$weight, rss=NULL, pearson=NULL, mutualinfo=NULL)     # SPECIFIC FOR EACH INTERACTION - EDGE ID
write.table(network_edge_attr, file = "networkedgeattribute_table.csv", row.names = FALSE, col.names = FALSE, quote=FALSE)

###############################################################################################################################################################################################################
# FOR CYTOSCAPE NETWORK TABLE:
d <- match(int_opt$nodeA, annotation_table_1[,1])
d <- annotation_table_1[d,3]

e <- match(int_opt$nodeB, annotation_table_1[,1])
e <- annotation_table_1[e,3]

f <- data.frame(d,e)
f <- f[-which(f$d == "" | f$e== ""), ]
write.table(f, file = "cytoscape_interaction_table.csv", row.names = FALSE, col.names = FALSE)
