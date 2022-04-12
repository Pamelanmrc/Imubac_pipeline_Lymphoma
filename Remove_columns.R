
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library("flowCore")

setwd("~/all_samples_scales_csv/")
setwd("~/AIM_cells_channel_csv/")
#CD11c
#CD138
#CD19
#CD21
#CD24
#CD27
#IgD
#IgG
#IgM

x <-exprs(read.FCS("export_Lymph_S073021_14_CD3.csv", transformation = FALSE))
x<-x[1:nrow(x),1:ncol(x)]
DataList=list()

setwd("~/all_samples_scales_csv/For_anna_for_flowjo/")
f="export_Lymph_S073021_14_CD3_TEST.csv"
lfiles= list.files(pattern=".csv",)
length(lfiles)
for (f in lfiles){
  #x=read.csv(f, header=TRUE, skip=341)
  x=read.csv(f, header=TRUE,stringsAsFactors = F, check.names = FALSE)
  num_data <- data.frame(data.matrix(x))
x=num_data
length(names(x))
  #colnames(x) = gsub("\\.","", gsub("Comp","",gsub("APCFire","",gsub("PECy5A","",colnames(x)))))
  #names(x %>% select(contains("CD11")))
  #grep("CD11c$", names(x))
  #grep("CD24$", names(x))
  #grep ("CD138$", names(x))
  #grep ("CD19$", names(x))
  #grep("IGg", names(x))
 
  #drops<-c(names(x)[29], names(x)[15],names(x)[17], names(x)[35], names(x)[18],names(x)[10],names(x)[14], names(x)[19], names(x)[30] )
  #drops <- c(  names(x)[8], names(x)[27], names(x)[16], names(x)[13], names(x)[15], names(x)[33], names(x)[12], names(x)[17], names(x)[28] )
  #x=x[ , !(names(x) %in% drops)]
  #new=str_replace(f,".csv", "_Dropped.csv")

  
  
  metadata <-  data.frame(name=dimnames(x)[[2]], desc=paste('column',dimnames(x)[[2]], 'from dataset'))
  metadata$minRange <- apply(x,2,min)
  metadata$maxRange <- apply(x,2,max)
  fcs_object <- new("flowFrame",exprs=as.matrix(x), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
  
  #f=str_replace(f, ".csv", "_DROPPED.csv")
  replace_name =gsub(".csv","",f)
  write.FCS(fcs_object, paste0(replace_name,".fcs"))
  #write.csv(x, file=new, quote=FALSE, row.names=FALSE)
  metadata <- data.frame()
  
}