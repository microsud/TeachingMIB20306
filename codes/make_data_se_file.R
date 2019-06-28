
make_data_sefile <- function(proteomeData = NULL, SampleData = NULL){
  
  data <- read.csv(proteomeData, 
                   header = T, 
                   row.names = 1)
  # Read sample information
  sam_dat <- read.csv(SampleData, 
                      header = T,
                      row.names = 1, 
                      check.names = F)
  
  # Filter out contaminants
  data <- filter(data, Reverse != "+", Potential.contaminant != "+")
  data$Gene.names %>% duplicated() %>% any()
  #head(data)
  LFQ_columns <- grep("LFQ.", colnames(data)) 
  #rownames(sam_dat)
  #data$id
  # Focus here on Anaerobutyricum soehngenii strain L2-7
  # Make data compatible with DEP package
  sam_dat.l2 <- subset(sam_dat, Organism == "L2-7")
  experimental_design <- sam_dat.l2
  
  data$Gene.names <- gsub("EHL_","EHLA_", data$Gene.names)
  data <- tidyr::separate(data, Gene.names, c("Gene.names", "Product name"), sep = " ", remove = TRUE)
  
  data_unique <- make_unique(data, "Gene.names", "Fasta.headers", delim = ";")
  
  sam_dat.l2$condition <- sam_dat.l2$Condition
  sam_dat.l2$label <- rownames(sam_dat.l2)
  sam_dat.l2$replicate <- c("1","2", "3", "1","2", "3","1","2", "3")
  LFQ_columns <- grep("LFQ.", colnames(data_unique))
  
  data_se <- make_se(data_unique, LFQ_columns, sam_dat.l2)
  return(data_se)
}
