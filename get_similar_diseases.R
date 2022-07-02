#' Estimates the p value of the Jaccard Coefficient for pairs of diseases
#'
#' This function estimates the statistic significance of the Jaccard Coefficient with a bootstrap
#' @param disease receives a disease, or a list of diseases, that are used as query to get the associated diseases
#' @param entity type of entity to connect diseases: gene or variant
#' @param min_entities minimum number of genes/variants shared between diseases
#' @param nboot Number of iterations sued to compute the pvalue associted
#' to the calculated Jaccard Index. Default: 1000.
#' @param ncores Number of cores used to calculate the pvalue associated to
#' the computed Jaccard Index. Default: 1.
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @examples
#' res <- get_similar_diseases( c("C0002395","C0011570"), entity = "gene", nboot = 10, min_entities = 50, verbose = T)
#' @export get_similar_diseases



get_similar_diseases <- function(disease, entity = "gene", min_entities = 10, nboot = 10, ncores = 1, verbose = FALSE) {
  
  if(entity =="gene"){
    data <- fread("/home/janet/Documents/all_gene_disease_associations.tsv.gz")
    data <- unique(data[ , c("geneId", "diseaseId", "diseaseName")])
    colnames(data) <- tolower(colnames(data))
    universe <- unique(data$geneid)
    
    
    diseaseid1 <-  unique(data[ diseaseid %in% disease]$diseaseid) 
    gene_list <- unique(data[ diseaseid %in% diseaseid1]$geneid)
    diseaseid2 <- unique(data[ geneid %in% gene_list]$diseaseid)
    
    
    aa <- data[diseaseid  %in% diseaseid1, c("geneid", "diseaseid")]
    bb <- data[diseaseid  %in%  diseaseid2, c("geneid", "diseaseid")]
    
    aa = data.table(aa, key=c("geneid"))
    bb = data.table(bb, key=c("geneid"))
    disease.pairs <- merge(aa, bb,  by=.EACHI, allow.cartesian=TRUE  )
    disease.pairs <- subset(disease.pairs, diseaseid.x != diseaseid.y)
    colnames(disease.pairs) <- c("geneid", "diseaseid1", "diseaseid2")
    
    disease.pairs <- disease.pairs[, .( ngenes = length(geneid)), by = c("diseaseid1","diseaseid2")]
    dt.agg <- aa[, .( ngenes1 = length(geneid)), by = c("diseaseid")]
    
    disease.pairs <- merge(disease.pairs, dt.agg, by.x = c("diseaseid1" ), by.y = c("diseaseid" ))
    dt.agg <- bb[, .( ngenes2 = length(geneid)), by = c("diseaseid")]
    disease.pairs <- merge(disease.pairs, dt.agg, by.x = c("diseaseid2" ), by.y = c("diseaseid" ))
    
    disease.pairs$union <- disease.pairs$ngenes1+disease.pairs$ngenes2 - disease.pairs$ngenes
    disease.pairs$jaccard_genes <- disease.pairs$ngenes/disease.pairs$union 
    disease.pairs <- merge(disease.pairs, unique(data[, c("diseaseid", "diseasename")]), by.x = "diseaseid1", by.y = "diseaseid")
    colnames(disease.pairs)[which(colnames(disease.pairs) == "diseasename")] <- "disease1_name"
    disease.pairs <- merge(disease.pairs, unique(data[, c("diseaseid", "diseasename")]), by.x = "diseaseid2", by.y = "diseaseid")
    colnames(disease.pairs)[which(colnames(disease.pairs) == "diseasename")] <- "disease2_name"
    
    message(paste0("There is total of ", length(diseaseid2),  " diseases associated to the input disease(s) ",
                   "with a total number of ", nrow(disease.pairs), "  possible combinations\n")
    )
    
    disease.pairs <- disease.pairs[ngenes >= min_entities]
    
    message(paste0("After applying cut off of number of entities, ", min_entities, 
                   " there is total   number of ", nrow(disease.pairs), "  possible combinations\n")
    )
    
  } else if(entity =="variant"){
    data <- fread("/home/janet/Documents/all_variant_disease_associations.tsv.gz")
    data <- unique(data[ , c("snpId", "diseaseId", "diseaseName")])
    colnames(data) <- tolower(colnames(data))
    universe <- data$snpid
    
    diseaseid1 <-  unique(data[ diseaseid %in% disease]$diseaseid) 
    variant_list <- unique(data[ diseaseid %in% diseaseid1]$snpid)
    diseaseid2 <- unique(data[ snpid %in% variant_list]$diseaseid)
    
    
    aa <- data[diseaseid  %in% diseaseid1, c("snpid", "diseaseid")]
    bb <- data[diseaseid  %in%  diseaseid2, c("snpid", "diseaseid")]
    
    aa = data.table(aa, key=c("snpid"))
    bb = data.table(bb, key=c("snpid"))
    disease.pairs <- merge(aa, bb,  by=.EACHI, allow.cartesian=TRUE  )
    disease.pairs <- subset(disease.pairs, diseaseid.x != diseaseid.y)
    colnames(disease.pairs) <- c("snpid", "diseaseid1", "diseaseid2")
    
    disease.pairs <- disease.pairs[, .( nvariants = length(snpid)), by = c("diseaseid1","diseaseid2")]
    dt.agg <- aa[, .( nvariants1 = length(snpid)), by = c("diseaseid")]
    
    disease.pairs <- merge(disease.pairs, dt.agg, by.x = c("diseaseid1" ), by.y = c("diseaseid" ))
    dt.agg <- bb[, .( nvariants2 = length(snpid)), by = c("diseaseid")]
    disease.pairs <- merge(disease.pairs, dt.agg, by.x = c("diseaseid2" ), by.y = c("diseaseid" ))
    
    disease.pairs$union <- disease.pairs$nvariants1+disease.pairs$nvariants2 - disease.pairs$nvariants
    disease.pairs$jaccard_variants <- disease.pairs$nvariants/disease.pairs$union 
    disease.pairs <- merge(disease.pairs, unique(data[, c("diseaseid", "diseasename")]), by.x = "diseaseid1", by.y = "diseaseid")
    colnames(disease.pairs)[which(colnames(disease.pairs) == "diseasename")] <- "disease1_name"
    disease.pairs <- merge(disease.pairs, unique(data[, c("diseaseid", "diseasename")]), by.x = "diseaseid2", by.y = "diseaseid")
    colnames(disease.pairs)[which(colnames(disease.pairs) == "diseasename")] <- "disease2_name"
    
    message(paste0("There is total of ", length(diseaseid2), " diseases associated to the input disease(s) ",
                   "with a total number of ", nrow(disease.pairs), "  possible combinations\n")
    )
    
    disease.pairs <- disease.pairs[nvariants >= min_entities]
    
    message(paste0("After applying cut off of number of entities, ", min_entities, 
                   " there is total   number of ", nrow(disease.pairs), "  possible combinations\n")
    )
  } else {
    stop(paste0(entity, " is a wrong entity! Allowed entities are: \n", 
                paste(c("gene","variant"), collapse = "\n")))
  }
  
  
  message("Pvalue estimation")
  
  for(i in 1:nrow( disease.pairs )){
    if( verbose ) {
      message(paste0("\n\t->Disease pair ", i, " of ", nrow(disease.pairs), " total diseases' pairs."))
    }
    if(entity  == "gene"){
      bb <- ji.internal(disease.pairs$ngenes1[i], disease.pairs$ngenes2[i], universe, nboot, ncores)
      pvalue <- (1 + sum(bb > disease.pairs$jaccard_genes[i]) ) / (nboot +1)
    }  else if(entity =="variant"){
      bb <- ji.internal(disease.pairs$nvariants1[i], disease.pairs$nvariants2[i], universe, nboot, ncores)
      pvalue <- (1 + sum(bb > disease.pairs$jaccard_variants[i]) ) / (nboot +1)
    }  

    disease.pairs$pvalue[i] <- round(pvalue, 5)
  }
  return(disease.pairs)
}


ji.internal <- function(len1, len2, universe, nboot, ncores) {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    pfun <- lapply
  } else {
    pfun <- parallel::mclapply
  }
  
  unlist(pfun(1:nboot, function(ii) {
    g1 <- sample( universe, len1 )
    g2 <- sample( universe, len2 )
    ja.coefr <- length(intersect(g1, g2)) / length(union(g1, g2))
  }, mc.cores = ncores))
}
