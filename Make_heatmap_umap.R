library(dplyr)
devtools::install_github("immunogenomics/harmony")

devtools::install_github("casanova-lab/iMUBAC")


#Write the right FCS files
norm_exp = assay(sce_down_snn,"normexprs")
file_name = as.data.frame(t(colData(sce_down_snn)$file_name))
#Make the column names as file names for the normalized expr
colnames(norm_exp) = file_name
#n_exp <-rbind(file_name, norm_exp)
#names(n_exp) <- as.character(unlist(n_exp[1,])
#n_exp<- n_exp[-1,]

n_exp <- norm_exp
colnames_list <-unique(sort(colnames(n_exp)))

for (i in colnames_list){
    print (i)
    make_files(i, n_exp)
  
}

#Read the patien file
patient_data= read.table("Final_info_patients.txt", sep=" ", header=TRUE, check.names=FALSE)

make_files<- function(fname, df){
  #sp <- sapply(n_exp, function(fname){ grep (fname, colnames(n_exp))})
  #sp<- grep(i, colnames(n_exp))
  sp<-(grep(fname, colnames(df)))
  new_n_exp <- n_exp[,sp]
  new_n_exp_t <- as.data.frame(t(new_n_exp))
  #Find the date to match with file
  date = gsub(pattern="S0",replacement = "", x=stringr::str_split(fname, "_", simplify = TRUE)[3])
  exp_id = gsub(pattern="S0",replacement = "", x=stringr::str_split(fname, "_", simplify = TRUE)[4])
  exp_id = gsub(pattern = "^0",  replacement ="", exp_id)
  print (paste(date, exp_id))
  
  index = grep(date ,patient_data$Date) 
  index =  grep (paste("^",exp_id,"$", sep=""), patient_data[index,]$Experiment_ID)
  

  

  new_n_exp_t <- new_n_exp_t %>% mutate(patient_data[index,]$Response, patient_data[index,]$Cancer_type, patient_data[index,]$Study_ID, patient_data[index,]$Experiment_ID  )
 
  colnames(new_n_exp_t) <- stringi::stri_replace_all_regex(colnames(new_n_exp_t), 
                                                         pattern=c('\\[', '\\]', 'index', '\\,', '\\$', 'patient_data', ' '), 
                                                         replacement=c('','','','','','', ''),
                                                         vectorize=FALSE)
  
    csv_file = stringr::str_replace(fname,".fcs","_.csv")
    csv_file <-paste0("/Users/mishrp03/all_samples_scales_csv/For_anna_for_flowjo/",csv_file)
    write.csv(new_n_exp_t, csv_file, row.names = FALSE, quote=FALSE)
    index=""
}


#Find the healthy
G<-as.data.frame(stringr::str_split(string=colnames_list[1:10], pattern="_" , simplify = TRUE))[3:4]

fcs_files = unique(sort(colData(sce_down_snn)$file_name))


for (i in f1){
  
  print (i)
  cluster<-sce_down_snn$"cluster_id"
  rows_markers_obj <- sce_down_snn[i]
  marker_norm_exp = assay(rows_markers_obj,"normexprs")
  #Umap values from the selected rows
  Umap_norm <- reducedDim(rows_markers_obj, "UMAPnorm")
   #Cluster information and nornalized values
   clust_expt <- as.data.frame(cbind(cluster, t(as.data.frame(marker_norm_exp))))
   
   #combine the umap info clusters
   final <- cbind(clust_expt,Umap_norm)
   write.table(final, file=i, row.names = FALSE, quote = FALSE)
    #j<-clust_expt %>% group_by(cluster) %>% summarise(median=median(i))
    filename = paste0(i,"_Umapnorm_normalizedexpr.pdf")
    print (filename)
    pdf(filename)
    colnames(final) =c("cluster_id", "Marker_expression", "UMAPnorm1", "UMAPnorm2")
    ggplot(final, aes(x=UMAPnorm1,y=UMAPnorm2, colour=Marker_expression)) + geom_point( ) + 
      scale_colour_gradient(low="blue", high="red")  + ggtitle(paste0(i, " Umap")) +
      theme(plot.title = element_text(family="Arial", face="bold-italics", size=15)) + 
      theme_bw()
      
    dev.off()

   
}


#x=read.table("UMAP_norm_clusters_CD8_normalized", header=TRUE, sep=",")
#head(x)
#new_cd8<-x%>% left_join(j, by="cluster")
#dim(new_cd8)
#ggplot(new_cd8, aes(x=UMAPnorm1,y=UMAPnorm2, fill="red",  label=median)) + geom_point(aes(color=median) ) + scale_colour_gradient(low="grey", high="red") +  geom_text(check_overlap = TRUE)
#head(new_cd8)
#ggplot(new_cd8, aes(x=UMAPnorm1,y=UMAPnorm2,  scale=T)) + geom_point(aes(color=median) )
  