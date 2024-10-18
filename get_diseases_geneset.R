library(igraph)
library(RCy3)

#import data
setwd("~/Documents/MU/Y2P1/Net Bio/project")
human_disease_textmining_full <- read.delim("human_disease_textmining_full.tsv", header = FALSE)
colnames(human_disease_textmining_full) <- c("g.id", "g.name", "d.id","d.name","z-core","confidence","url")

#set filename, DOID, confidence cut-off
#Alzheimers
ad_id <- c("AD", "DOID:10652", 2.9)
#Parkinsons
pd_id <- c("PD", "DOID:14330",2.8) 
#ALS
als_id <-c("ALS", "DOID:332",2.5)
#Huntingtons
hd_id <- c("HD", "DOID:12858",2.5)
#Prion
p_id <- c("P", "DOID:649",2)
#Frontotemporal dementia
ftd_id <- c("FDT", "DOID:9255",2.3)
#Multiple Sclerosis
ms_id <- c("MS", "DOID:2377",2.7)
#Lewy body dementia	
lbd_id <- c ("LBD", "DOID:12217",2)
# Spinal muscular atrophy
sma_id <- c("SMA", "DOID:12377",2)
# Friedreich ataxia
fa_id <- c("FA", "DOID:12705",2)
#Progressive supranuclear palsy
psp_id <- c("PSP", "DOID:678",2)
#Gerstmann-Straussler-Scheinker syndrome
gsss_id <- c("GSSS", "DOID:4249",1)
#Multiple system atrophy
msa_id <- c("MSA", "DOID:4752",1.75)
#Olivopontocerebellar atrophy
opca_id <- c("OPCA", "DOID:14784",1)
#Vascular dementia
vd_id <- c("VD", "DOID:8725",1.9)
#Progressive Multifocal Leukoencephalopathy
pml_id <- c("PML", "DOID:643",1.5)
#Wernicke-Korsakoff syndrome
wks_id<- c("WKS", "DOID:10915",1.2)
#Neurodegeneration with brain iron accumulation
ndbia_id<- c("NBDIA", "DOID:0110734",1.9)

nd_list <- list(ad_id, pd_id, als_id, hd_id, p_id, ftd_id, ms_id, lbd_id, 
                sma_id,fa_id, psp_id, gsss_id,msa_id,opca_id,vd_id,pml_id,wks_id,ndbia_id)

# finding genes above cutoff for each disease
gene_list <- data.frame()
diseases <- data.frame()
for (i in 1:length(nd_list)){
  disease <- nd_list[[i]][1]
  filename <- paste(disease, "tsv", sep=".")
  id <- nd_list[[i]][2]
  conf <- nd_list[[i]][3]
  disease_genes <- subset(human_disease_textmining_full, d.id ==id & confidence > conf)
  #hist(disease_genes$confidence)
  #write to file
  write.table(disease_genes, file = filename, row.names=FALSE, sep="\t")
  gene_temp <- data.frame(id = disease, genes = disease_genes$g.id)
  print(length(gene_temp$genes))
  gene_list <- rbind(gene_list, gene_temp)
  diseases <- rbind(diseases, disease)
}
colnames(diseases) <- "d.name"
#get disease names

similarity = matrix(, nrow = length(nd_list), ncol = length(nd_list))
rownames(similarity) <- diseases$d.name
colnames(similarity) <- diseases$d.name
#calculate overlap
for (i in 1:length(nd_list)){
  for(j in 1:length(nd_list)){
    if(i ==j){
      similarity[i, j] <- 0
    }
    else{
      gene_list_1 <- subset(gene_list, id ==nd_list[[i]][1])
      gene_list_2 <- subset(gene_list, id ==nd_list[[j]][1])
      
      gene_list_1 <- gene_list_1$genes
      gene_list_2 <- gene_list_2$genes
      
      overlap <- length(intersect(gene_list_1, gene_list_2))/min(length(gene_list_1), length(gene_list_2))
      similarity[i, j] <- overlap
    }
  }
}

write.table(similarity, file="similarity.txt", row.names=FALSE, col.names=FALSE)

hist(similarity)
#set cut-off to .2
threshold = 0.1
similarity[similarity < threshold] <- 0

#plot in cytoscape
g <- graph_from_adjacency_matrix(similarity, mode = "lower", weighted = "weight")
plot(g, edge.width = E(g)$weight, edge.label = E(g)$weight)

RCy3::createNetworkFromIgraph(g, title="Association Network", collection="igraph")

