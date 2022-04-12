setwd("~/all_samples_scales_csv/")
setwd("~/AIM_cells_channel_csv/")
setwd("~/Downloads/")

load("Results.RData")
library(flowCore)
library(stringr)
library(tidyverse)
library(data.table)
library(iMUBAC)
library(SingleCellExperiment)
setwd("/Users/mishrp03/all_samples_scales_csv/")
setwd("/Users/mishrp03/AIM_cells_channel_csv/test")

#"S071221" was removed
batch_names =c("S0813021","S073021","S071421")

getwd()
file_names= list.files(pattern="_DROPPED.csv")
md <- data.frame()
md_final <- data.frame()
for (batch in batch_names){
  
  if (batch == "S0813021"){
    file_names = intersect(list.files(pattern="._DROPPED.fcs"), list.files(pattern="Hea"))
    md <- data.frame(file_name=file_names, batch=paste0("Data_",batch,"_Hea"),
                     panel="Panel3", 
                     group="HD", treatment="base")
    
    #Add ids to the healthy samples
    md=mutate(md, donor_id = paste0(group,"_",batch,"_",row_number()))
    md$full_path = paste(getwd(),md$file_name,sep="/")
    #md=mutate(md, donor_id = paste0(group,row_number()))
    
    md_final<- rbind(md_final, md)
  }else{
    
    file_names = intersect(intersect(list.files(pattern=".fcs"), list.files(pattern="Lym")), list.files(pattern=batch))
    md <- data.frame(file_name=file_names, batch=paste0("Data_",batch,"_Lymph"),
                     panel="Panel3", 
                     group="Lymp", treatment="base" )
    #Add ids to the Lymphoma samples
    md=mutate(md, donor_id = paste0(group,"_",batch,"_",row_number()))
    md$full_path = paste(getwd(),md$file_name,sep="/")
    #md=mutate(md, donor_id = paste0(group,row_number()))
    
    md_final<- rbind(md_final, md)
  }                     
}

#Convert it into factors
md_final$batch = as.factor(md_final$batch)
md_final$group = as.factor(md_final$group)
md_final$treatment = as.factor(md_final$treatment)
#Convert into data table
md_final = setDT(md_final)
#md_final[,batch :=factor(batch, levels=c("Data_S0813021_Hea","Data_S073021_Lymph","Data_S071421_Lymph"))]
#md_final[,group :=factor(group, levels=c("HD","Lymp"))]



# To first convert it into FCS
file_names = list.files(pattern=".csv") 
readfile(file_names)






#md$full_path = paste(getwd(),md$file_name,sep="/")
#md_final = setDT(md_final)




#Read the FCS file to get their columnnames
#Considering healthy, Lmph etc have the same headers.
#pd <- list(
#data.table::fread(system.file("PanelInfo_Data23_Panel3.csv", package="iMUBAC")),
#data.table::fread(system.file("PanelInfo_Data29_Panel3.csv", package="iMUBAC"))
#) %>% magrittr::set_names(c("Data23","Data29"))
#
tempdata <- exprs(read.FCS(md_final$file_name[[1]], transformation = FALSE))
tempdata <- tempdata[1:nrow(tempdata),1:ncol(tempdata)]
colnames(tempdata)
dim(tempdata)
#CHECK FOR CD40L
markers =c("CD278", "CD38", "CD137", "GzmB", "CD185","CD279","CD194", "CD19","CD25", "CD23","HLADR","IgD","CD21",
           "CD95", "CD24","CD138",  "IgM", "CD197", "CD69", "CD183","CD134",
           "FOXP3", "Tbet", "CD8", "CD183","CD40", "CD11c", "IGg" ,"CD4", "CD127", "CD3","CD45RA","CD27", "CCR6" )

markers =c("CD278", "CD38", "CD137", "GzmB", "CD185","CD279","CD194","CD25","CD23","HLADR","CD95", "CD197",
           "CD69", "CD183","CD134", "FOXP3", "Tbet", "CD8", "CD40L","CD40", "CD4", "CD127", "CD3","CD45RA","CCR6" )
g=cbind(as.data.frame(colnames(tempdata)), markers)
names(g)=c("channel", "antigen")
g=setDT(g)
list_batches = unique(sort(md_final$batch))
length(list_batches)

#Just appended 3 times. Can be done better (Needs to be adjusted)
pd<-list(g,g,g) %>% magrittr::set_names(list_batches)

md_final$group


md_final$group
#Laod into a object
sce_cytof <- iMUBAC::prepSCE(
  md=md_final,
  pd=pd,
  channel_length=NULL,      ## Doublet detection by event length is disabled.
  channel_DNA=NULL,
  channel_LD=NULL,
  type="CyTOF"
)

sce_flow <- iMUBAC::prepSCE(
  md=md_final,
  pd=pd,
  channel_length=NULL,      ## Doublet detection by event length is disabled.
  channel_DNA=NULL,
  channel_LD=NULL,
  type="Flow"
)

options(java.parameters="-Xmx80G")

colData(sce_cytof) <- colData(sce_cytof)[c("batch","file_name","panel","group","treatment","donor_id" )]
colData(sce_flow) <- colData(sce_flow)[c("batch","file_name","panel","group","treatment","donor_id" )]
colData(sce_all_data) <-  colData(sce_all_data)[c("batch","file_name","panel","group","treatment","donor_id" )]
saveRDS(sce_cytof, "CyTek_SCE_onlyHealthy_cytof.rds")
saveRDS(sce_down_snn, "CyTek_SCE_down_snn.rds")

saveRDS(sce_flow, "CyTek_FLOW_SCE_flow_DROPPED.rds")
sce1 <- readRDS("CyTek_FLOW_SCE_flow_DROPPED.rds")
saveRDS(sce_flow_down, "Cytek_sce_flow_down.rds")
#Downsample the healthy controls
sce_cytof <- readRDS("CyTek_SCE_onlyHealthy_cytof.rds")
sce_flow <- readRDS("CyTek_FLOW_SCE_onlyHealthy_flow.rds")

sce_cytof_down <- sce_cytof[,sce_cytof$"group"=="HD"]
sce_flow_down <- sce_flow[,sce_flow$"group"=="HD"]

sce_flow_down <- sce_flow
sce_flow_down <- sce_flow
sce_flow_down <- iMUBAC::batchCorrection(
  sce_flow_down,
  maxN=50000, ## A maximum of 50000 cells are randomly selected from each batch.
  seed=12345  ## a random seed
)
sce_cytof_down <- iMUBAC::batchCorrection(
  sce_cytof_down,
  maxN=50000, ## A maximum of 50000 cells are randomly selected from each batch.
  seed=12345  ## a random seed
)


sce_all_data_down<- iMUBAC::batchCorrection(
  sce_all_data_down,
  maxN=50000, ## A maximum of 50000 cells are randomly selected from each batch.
  seed=12345  ## a random seed
)
#Tale all data and subset it 
sce_all_data <- iMUBAC::prepSCE(
  md=md_final,
  pd=pd,
  channel_length=NULL,      ## Doublet detection by event length is disabled.
  channel_DNA=NULL,
  channel_LD=NULL,
  type="Flow"
)


#After batch correctio use Umap
sce_all_data_down <- iMUBAC::runUMAP(
  sce_all_data_down,
  by_exprs_values="exprs",
  name="UMAP",
  n_neighbors=25,
  min_dist=0.4,
  scale=T,
  n_threads=4, ## the number of threads for parallel computing
  seed=12345
)
sce_all_data_down <- iMUBAC::runUMAP(
  sce_all_data_down,
  by_exprs_values="normexprs",
  name="UMAPnorm",
  n_neighbors=25,
  min_dist=0.4,
  scale=T,
  n_threads=4,
  seed=12345
)



sce_flow_down <- iMUBAC::runUMAP(
  sce_flow_down,
  by_exprs_values="exprs",
  name="UMAP",
  n_neighbors=25,
  min_dist=0.4,
  scale=T,
  n_threads=4,
  seed=12345
)

sce_flow_down <- iMUBAC::runUMAP(
  sce_flow_down,
  by_exprs_values="normexprs",
  name="UMAPnorm",
  n_neighbors=25,
  min_dist=0.4,
  scale=T,
  n_threads=4,
  seed=12345
)


sce_all_data <- iMUBAC::runUMAP(
  sce_all_data,
  by_exprs_values="exprs",
  name="UMAP",
  n_neighbors=25,
  min_dist=0.4,
  scale=T,
  n_threads=4, ## the number of threads for parallel computing
  seed=12345
)



pdf("Umap_by_cluster.pdf")
ggpubr::ggarrange(
  iMUBAC::plotDR(
    sce_down_snn,
    dimred="UMAPnorm",
    colour_by="cluster_id",
    lave
  ) +
    scale_colour_brewer(palette="Dark2"))
  dev.off()



pdf("Plot_all_data_Umpas_normexprs.pdf")
dev.off()
plt_all_data_down <- ggpubr::ggarrange(
  iMUBAC::plotDR(
    sce_down_snn,
    by_exprs_values="normexprs",
    colour_by="batch",
  ) +
    scale_colour_brewer(palette="Dark2")),
  iMUBAC::plotDR(
    sce_all_data_down,
    dimred="UMAPnorm",
    colour_by="batch"
  ) +
    scale_colour_brewer(palette="Dark2"),
  ncol=2, nrow=1, common.legend=T, legend="right"
)

pdf("Flow_Umaps_afterbatch_50Kcells_ALL.pdf")
plt_all_data_down
dev.off()
#plotUtility::savePDF(plt1, o="Flow_Plt1.pdf", w=15, h=6)


sce_down_snn<- iMUBAC::clustering(
  sce_flow_down,
  features=rownames(sce_flow_down), ## Using all markers for clustering
  by_exprs_values="normexprs",
  method="SNNGraph",
  n_components=10, ## the number of reduced dimentions for constructing the SNN graph
  n_neighbors=25,  ## the parameters for UMAP-based dimention reduction
  min_dist=0.4,    ## the parameters for UMAP-based dimention reduction
  seed=12345
)




pdf("Heatmap_all_n30.pdf")
plt_snn_sce_all_data <- iMUBAC::plotClusterHeatmap(
  sce_down_snn,
  features=rownames(sce_down_snn),
  clusters=sce_down_snn$"cluster_id",
  by_exprs_values="normexprs",
  fun="median",
  scale=T,
  cluster_rows=T,
  cluster_anno=T,
  draw_dend=T,
  draw_freqs=T
)

dev.off()
pdf("Umap_heatmap_after_clsuter_all_n30.pdf")
plt_snn_sce_all_data_umap <- iMUBAC::plotDR(
  sce_flow_down,
  dimred="UMAP",
  colour_by="cluster_id",
  text_by="cluster_id"    ## to overlay cluster ids on each of the clusters
)



#Read the CSV files and convert them into FCS files. These files are the ones that have header 
#starting at line 341
readfile <- function(file_names){
  for (i in file_names) {
      x=read.csv(i, header=TRUE)
      colnames(x) = gsub("\\.","", gsub("Comp","",gsub("APCFire","",gsub("PECy5A","",colnames(x)))))
                  
      #Still replace a few 
      colnames(x)=str_replace_all(colnames(x), c("BB700A","BUV395A", "AlexaFluor647A", "AlexaFluor700A","SparkBlue550A"),"")
      metadata <-  data.frame(name=dimnames(x)[[2]], desc=paste('column',dimnames(x)[[2]], 'from dataset'))
      metadata$minRange <- apply(x,2,min)
      metadata$maxRange <- apply(x,2,max)

      #Make an object from the CSV
      #This class represents the data contained in a FCS file or similar data structure. There are three parts of the data:
      # a numeric matrix of the raw measurement values with rows=events and columns=parameters
      #annotation for the parameters (e.g., the measurement channels, stains, dynamic range)
      #additional annotation provided through keywords in the FCS file
      #https://www.rdocumentation.org/packages/flowCore/versions/1.38.2/topics/flowFrame-class
      fcs_object <- new("flowFrame",exprs=as.matrix(x), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
      replace_name =gsub(".csv","",i)
      write.FCS(fcs_object, paste0(replace_name, ".fcs"))
      
  }
}
sce_down_snn$group


only_hd = subset(sce_down_snn, , group=="HD")
CD8=only_hd["CD8"]
CD4=only_hd["CD4"]
assays(CD8)$normexprs
only_lymp = subset(sce_down_snn, , group=="Lymp")
CD8_lymp= only_lymp["CD8"]
CD4_lymp= only_lymp["CD4"]
getwd()
pdf("CD4_Healthy.pdf")
iMUBAC::plotDR(
  CD8_lymp,
  by_exprs_values="UMAPnorm",
  dimred="UMAPnorm",
  color_by="batch"

)

BiocManager::install("tnaake/MatrixQCvis")
featurePlot(CD8)
sce_down_snn$celltype

CD8$celltype
dev.off()
pdf("CD4_lymp.pdf")
iMUBAC::plotDR(
  CD4_lymp,
  by_exprs_values="normexprs",
  colour_by="UMAPnorm",
  
)
assay(sce, by_exprs_values) 
dev.off()
only_hd = subset(sce_down_snn, , group=="HD")
only_lymp = subset(sce_down_snn, , group=="Lymp")

pdf("CD8.pdf")

head(reducedDim(sce_down_snn["CD4",]))

lypm8<-iMUBAC::plotDR(
  only_hf["CD8",],
  dimred="UMAPnorm",
  colour_by="cluster_id",
  text_by="cluster_id",
  scale=T ## to overlay cluster ids on each of the clusters
) + scale_colour_brewer(palette="Dark2")
dev.off()

pdf("Umap_by_cluster.pdf")
plt3
dev.off()
plt3 <- iMUBAC::plotDR(
  assays(sce_down_snn)[2]$normexprs,
  dimred="UMAPnorm",
  colour_by=h,
  by_exprs_values="normexprs",
  text_by="cluster_id",
  ## to overlay cluster ids on each of the clusters
) +
  ggpubr::rremove("legend")


library(flowCore)

# Median expression heatmap
plt2 <- iMUBAC::plotClusterHeatmap(
  sce_flow,
  features=rownames(sce_down_snn),
  clusters=sce_down_snn$"cluster_id",
  by_exprs_values="normexprs",
  fun="median",
  scale=T,
  cluster_rows=T,
  cluster_anno=T,
  draw_dend=T,
  draw_freqs=T
)
pdf("Heatmap_AIM.pdf")
plt2
dev.off()

pdf("Heatmap_SNN_50K.pdf")

dev.off()

pdf("Clustering_SNN_500k.pdf")
plt2
dev.off()


iMUBAC::plotDR(sce_flow)

sce_flow <- iMUBAC::clusterPropagation(
  sce_flow,          
  sce_down_snn,      
  by_exprs_values="exprs", 
  maxN=100,
  numThreads=4, 
  seed=12345
)


#Use excel sheet
df_celltype <- readxl::read_excel("Markers_all.xlsx")
sce <- clusterAnnotation(sce, df_celltype)
sce_down_snn<- clusterAnnotation(sce_down_snn, df_celltype)

sce_down_snn <- iMUBAC::runUMAP(
  sce_down_snn,
  by_exprs_values="normexprs",
  name="UMAP",
  n_neighbors=25,
  min_dist=0.4,
  scale=T,
  n_threads=4,
  seed=12345
)



dt_da <- colData(sce) %>%
  as.data.frame() %>%
  dplyr::group_by(donor_id, batch, group, celltype) %>%
  dplyr::tally() %>%
  dplyr::mutate(percent=n/sum(n)*100) %>% data.table::as.data.table() 




#EDGER by celltype
#Removed a few things
dt_da <- colData(sce) %>%
  as.data.frame() %>%
  dplyr::group_by(donor_id, batch, group, celltype) %>%
  dplyr::tally() %>%
  dplyr::mutate(percent=n/sum(n)*100) %>%
  data.table::as.data.table()
write.table(dt_da, file="Howmany", quote=FALSE, row.names=FALSE)
mat <- mat %>% replace(is.na(.), 0)
d <- edgeR::DGEList(counts=mat, lib.size=colSums(mat), group=groups, remove.zeros=T)
design <- model.matrix(~0+groups+batches, data=d$samples)

d <- edgeR::estimateDisp(d, design)
fit <- edgeR::glmQLFit(d, design, robust=T)
dt_da_qlf <- edgeR::glmQLFTest(
  fit, 
  contrast=makeContrasts(NRvsR=groupsNR-groupsR, levels=design)
) %>%
  edgeR::topTags(n=Inf) %>%
  as.data.frame() %>%
  dplyr::transmute(
    celltype_detailed=rownames(.),
    log2FC=logFC, 
    PValue=PValue,
    AdjPValue=FDR
  ) %>%
  data.table::as.data.table()

pdf("Heatmap_Masato_AIM.pdf")
iMUBAC::plotClusterHeatmap(
  sce_down_snn,
  features=rownames(sce_down_snn),
  clusters=sce_down_snn$"celltype",
  by_exprs_values="normexprs",
  fun="median",
  scale=T,
  cluster_rows=F,
  cluster_anno=F,
  draw_dend=T,
  draw_freqs=T
)
dev.off()