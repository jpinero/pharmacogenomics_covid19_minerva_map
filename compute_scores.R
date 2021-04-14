rm(list = ls())
library(data.table)
library(dplyr)
library(tidyr)
setwd("/home/janet/Documents/COVID_PATHWAY_MAPS/DATA/processed_data/final_version/pharmacogenomics_covid19_minerva_map/")
# LOAD DATA
# The file contains the information about drugs associated to the proteins in the COVID-19 map
drug2target_data <- fread("data/drug2target.tsv")
# The file contains the genomic variants in the proteins in the COVID-19 map with pharmacogenomics information in PharmGKB, and the drugs
pharmacogenic_data <-fread("data/pharmacogenomics_annotations.tsv", quote = "")

# The folder contains GNOMAD exome information for the genes in the COVID-19 map 
gnomad_files <- list.files("data/exome_data/")
gnomad_files <- gsub(".tsv", "", gnomad_files)
 

##################
##################
### CAP SCORE 
##################
###################

populations <- c("", paste0("_", c("afr", "amr", "fin", "asj", "nfe", "eas", "sas")))
atts <- c("AF")
gender <- c("", "_male", "_female")
miscols <- c()
for (pp in populations){
  for (at in atts){
    aa <- paste(at, pp, sep = "")
    for (gg in gender){
      cc <- paste0(aa, gg)
      miscols <- c(miscols, cc)
    }
  }
}

miscols2 <- c()
for (pp in populations){
  aa <- paste("CAP_PKB", pp, sep = "")
  for (gg in gender){
    cc <- paste0(aa, gg)
    miscols2 <- c(miscols2, cc)
  }
}

atts <- c("nvar")
miscols3 <- c()
for (pp in populations){
  for (at in atts){
    aa <- paste(at, pp, sep = "")
    for (gg in gender){
      cc <- paste0(aa, gg)
      miscols3 <- c(miscols3, cc)
    }
  }
}

populations <- c( "afr", "amr", "fin", "asj", "nfe", "eas", "sas")
ddata <- data.frame(matrix(ncol = 2+ length(miscols3)+ length(miscols2), nrow = 0))
colnames(ddata) <-  c("gene", "tvar",  miscols2, miscols3)
i <- 1
for (gene in  intersect(drug2target_data$Symbol, gnomad_files)){
   print(gene)  

  outfile1 <- paste0( "data/exome_data/", gene, ".tsv")
  if (file.exists(outfile1)) {
    df <- fread(outfile1,  header=T, stringsAsFactors = F, quote = "")
    if ( length(intersect(colnames(df), "AF")) > 0){
      df <- subset(df, AF < 0.5 & V3!=".")
      df$id <- paste(df$V1, df$V2, df$V3, df$V4, df$V5, sep = "_")
      mutaciones <- length(unique(df$id))
      df3 <- merge(df, pharmacogenic_data, by.x =c("V3", "V5"), by.y=c("Variant", "Alleles")  )
      mutpkg <- length(unique(df3$id))
      if (mutpkg >0 ){
        
        ddata[i,"gene"]  <- gene
        ddata[i,"tvar"]  <- mutaciones
        ddata[i,"nvar"]  <- mutpkg
        prod <- 1
        prodfemale <-1
        prodmale <-1
        for (rr in 1:nrow(df3)) {
          prod <- prod*(1- df3$AF[rr])^2
          prodmale <- prodmale*(1- df3$AF_male[rr])^2
          prodfemale <- prodfemale*(1- df3$AF_female[rr])^2
        }
        ddata[i, c("CAP_PKB","CAP_PKB_male","CAP_PKB_female")] <- c(  1-prod, 1-prodmale, 1-prodfemale)
        ddata[i,"nvar_male"]  <- length(which(df3$AF_male>0))
        ddata[i,"nvar_female"]  <- length(which(df3$AF_female>0))
        
        for (pp in populations) {
          af <- paste("AF", pp, sep = "_")
          #print(af)
          af_female <- paste("AF", pp, "female", sep = "_")
          af_male <- paste("AF", pp, "male", sep = "_")
          prod <- 1
          prodfemale <-1
          prodmale <-1
          for (rr in 1:nrow(df3)) {
            prod <- prod*(1- df3[[af]][rr])^2
            prodmale <- prodmale*(1- df3[[af_male]][rr])^2
            prodfemale <- prodfemale*(1- df3[[af_female]][rr])^2
            
          }
          ddata[i,c( paste0("CAP_PKB_", pp), paste0("CAP_PKB_", pp, "_male"), paste0("CAP_PKB_", pp, "_female"))] <- c( 1-prod, 1-prodmale, 1-prodfemale)
          ddata[i,paste("nvar", pp, sep = "_")]  <- length(which(df3[[af]]>0))
          ddata[i,paste("nvar", pp, "male", sep = "_")]  <- length(which(df3[[af_male]] >0))
          ddata[i,paste("nvar", pp, "female", sep = "_")]  <- length(which(df3[[af_female]]>0))
          
        }
        i <- i +1 
      }
      
    } else {
      print("no data") 
      print(gene)
    }
  }
  else {print("does not exists") 
    print(outfile1)}
}  
# write the results
write.table(ddata, "cap_score_genes.tsv", sep = "\t", quote = F, row.names = F)


##################
### DRP SCORE for drugs 
##################
###################

populations <- c("", paste0("_", c("afr", "amr", "fin", "asj", "nfe", "eas", "sas")))
atts <- c("AF")
gender <- c("", "_male", "_female")
miscols <- c()
for (pp in populations){
  for (at in atts){
    aa <- paste(at, pp, sep = "")
    for (gg in gender){
      cc <- paste0(aa, gg)
      miscols <- c(miscols, cc)
    }
  }
}

miscols2 <- c()
for (pp in populations){
  aa <- paste("CAP_PKB", pp, sep = "")
  for (gg in gender){
    cc <- paste0(aa, gg)
    miscols2 <- c(miscols2, cc)
  }
}


atts <- c("nvar")
miscols3 <- c()
for (pp in populations){
  for (at in atts){
    aa <- paste(at, pp, sep = "")
    for (gg in gender){
      cc <- paste0(aa, gg)
      miscols3 <- c(miscols3, cc)
    }
  }
}

populations <- c( "afr", "amr", "fin", "asj", "nfe", "eas", "sas")

ddata <- data.frame(matrix(ncol = 2+length(miscols2)+length(miscols3), nrow = 0))
colnames(ddata) <-  c("drug",  "tvar", miscols2,  miscols3)
cc <- c("id","V3", "V5", miscols)
i <- 1


for (dd in intersect(drug2target_data$drugbank, pharmacogenic_data$drugbank)) {
  print(dd)
  gnomaddata <- data.frame(matrix(ncol = 3+length(miscols), nrow = 0))
  colnames(gnomaddata) <-  c("id","V3", "V5",   miscols)
  gg<- unique( drug2target_data[ drug2target_data$drugbank == dd, ]$Symbol )
  if (length(gg)>0) {
    for (gene in  gg ){
      outfile1 <- paste0("data/exome_data/", gene, ".tsv")
      if (file.exists(outfile1)) {
        df <- fread(outfile1,  header=T, stringsAsFactors = F)
        if ( length(intersect(colnames(df), "AF")) > 0){
          df <- subset(df, AF < 0.5 & V3!=".")
          df$id <- paste(df$V1, df$V2, df$V3, df$V4, df$V5, sep = "_")
          df <- df[, ..cc]
          gnomaddata <- rbind(gnomaddata, df)  
        }
      }  
      else {
        print(paste(outfile1, "does not exists") )
      }  
    }  
    # print(dim(gnomaddata))
    if (dim(gnomaddata)[1]> 0){
      
      prod <- 1
      prodfemale <-1
      prodmale <-1
      
      df3 <- merge(gnomaddata, pharmacogenic_data[ drugbank ==dd], by.x =c("V3", "V5"), by.y=c("Variant", "Alleles")  )
      mutpkg <- length(unique(df3$id))
      # ddata[i,"nvars_pkg"]  <- dim(gnomaddata)[1]
      if (mutpkg >0 ){
        ddata[i,"tvar"]  <- dim(gnomaddata)[1]
        ddata[i,"drug"]  <- dd
        ddata[i,"nvar"]  <- mutpkg
        prod <- 1
        prodfemale <-1
        prodmale <-1
        for (rr in 1:nrow(df3)) {
          prod <- prod*(1- df3$AF[rr])^2
          # print(prod)
          prodmale <- prodmale*(1- df3$AF_male[rr])^2
          prodfemale <- prodfemale*(1- df3$AF_female[rr])^2
        }
        ddata[i, c("CAP_PKB","CAP_PKB_male","CAP_PKB_female")] <- c(  1-prod, 1-prodmale, 1-prodfemale)
        ddata[i,"nvar_male"]  <- length(which(df3$AF_male>0))
        ddata[i,"nvar_female"]  <- length(which(df3$AF_female>0))
        
        for (pp in populations) {
          af <- paste("AF", pp, sep = "_")
          # print(af)
          af_female <- paste("AF", pp, "female", sep = "_")
          af_male <- paste("AF", pp, "male", sep = "_")
          prod <- 1
          prodfemale <-1
          prodmale <-1
          for (rr in 1:nrow(df3)) {
            prod <- prod*(1- df3[[af]][rr])^2
            prodmale <- prodmale*(1- df3[[af_male]][rr])^2
            prodfemale <- prodfemale*(1- df3[[af_female]][rr])^2
          }
          
          ddata[i,c( paste0("CAP_PKB_", pp), paste0("CAP_PKB_", pp, "_male"), paste0("CAP_PKB_", pp, "_female"))] <- c( 1-prod, 1-prodmale, 1-prodfemale)
          ddata[i,paste("nvar", pp, sep = "_")]  <- length(which(df3[[af]]>0))
          ddata[i,paste("nvar", pp, "male", sep = "_")]  <- length(which(df3[[af_male]] >0))
          ddata[i,paste("nvar", pp, "female", sep = "_")]  <- length(which(df3[[af_female]]>0))
        }
        i <- i +1 
      }
    }
  } else{
    print( dd)
    print(subset(dat_fin, drugbank ==dd ))
  }
}


write.table(ddata, "drp_score_drugs.tsv", sep = "\t", quote = F, row.names = F)


