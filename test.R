library(DEP)
data_prot <- readRDS("data/data_se_intestinimonas.rds")

data_filt <- filter_missval(data_prot, thr = 1)
data_norm <- normalize_vsn(data_filt)
set.seed(2156)
# All possible imputation methods are printed in an error, if an invalid function name is given.
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")

dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = 1.5)

data_results <- get_results(dep)
data_results.sub <- data_results[,-7]
colnames(data_results.sub) <- c("Locus Tag",	"Name",	"Glucose plus Acetate vs Lysine p.val",	"Glucose plus Acetate vs Lysine p.adj",	"significant",	"Glucose.plus Acetate vs Lysine significant",	"Glucose plus Acetate centered",	"Lysine centered")
head(data_results.sub)

library(reshape2)  
library(tidyverse)
data_results.sub2 <- data_results.sub

#data_results.sub2$Name <- str_replace(data_results.sub2$Name, "^[>AF_)", "")
#head(data_results.sub2)

  
data_results.sub2$Name = colsplit(data_results.sub2$Name," ", names = c("name1","Name"))[,2]
data_results.sub2$Name = colsplit(data_results.sub2$Name," ", names = c("name1","Name"))[,2]
head(data_results.sub2)
data_results.sub2$Name <- gsub( " reverse ", " ", data_results.sub2$Name)
data_results.sub2$Name <- gsub( " forward ", " ", data_results.sub2$Name)
head(data_results.sub2)

data_results.sub2$Name <- gsub( "^ [0-9]\\:", "", data_results.sub2$Name)


######
dir.create("kegg_maps")
list.ibu <- c("00010",	"00020",	"00030",	"00040",	"00051",	"00052",	"00061",	"00071",	"00072",	"00190",	"00220",	"00230",	"00240",	"00250",	"00260",	"00261",	"00270",	"00280",	"00281",	"00290",	"00300",	"00310",	"00311",	"00330",	"00332",	"00340",	"00350",	"00360",	"00362",	"00380",	"00400",	"00401",	"00410",	"00430",	"00440",	"00450",	"00460",	"00471",	"00473",	"00480",	"00500",	"00511",	"00514",	"00515",	"00520",	"00521",	"00524",	"00550",	"00561",	"00564",	"00592",	"00600",	"00620",	"00625",	"00630",	"00633",	"00640",	"00650",	"00660",	"00670",	"00680",	"00730",	"00740",	"00750",	"00760",	"00770",	"00780",	"00790",	"00860",	"00900",	"00903",	"00910",	"00920",	"00970",	"01054",	"01100",	"01110",	"01120",	"01130",	"01200",	"01210",	"01212",	"01220",	"01230",	"01501",	"01502",	"01503",	"02010",	"02020",	"02024",	"02030",	"03010",	"03018",	"03020",	"03030",	"03060",	"03070",	"03410",	"03420",	"03430",	"03440",	"04122",	"04122")

download.kegg(pathway.id = list.ibu, species = "ibu", kegg.dir = "kegg_maps/",
              file.type=c("xml", "png"))


#df$value <- gsub("^\\D.*", "", df$value)
library(clusterProfiler)
data_results.sub2$Name <- gsub( "^>AF_*^", "", data_results.sub2$Name)
head(data_results.sub2)

data_results <- get_results(dep)
df_wide <- get_df_wide(dep)

foldchanges.view = data_results$Lysine_centered
names(foldchanges.view) = data_results$name
gene <- names(foldchanges.view)[abs(foldchanges.view) > 1.5]
kk <- enrichKEGG(gene         = gene,
                 organism     = 'ibu',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 minGSSize = 5)

choices =   data <- read.csv("data/kegg_pathway_list.csv", 
                             header = T, 
                             row.names = 1)
# List of choices for selectInput
mylist <- rownames(choices)
mylist[1]

mylist[1]

id <- gsub("ibu", "", input$Select)

pathview(gene.data  = foldchanges.view,
         pathway.id = "00051",
         species    = "ko",
         kegg.native = T,
         limit      = list(gene=max(abs(foldchanges.view)), cpd=1),
         out.suffix = "FructoseMannose_metabolism", 
         low = "#d73027", mid = "#fee090", high = "#4575b4")


pathview::pathview(gene.data = foldchanges.view, 
                   pathway.id = "00310", 
                   gene.idtype = "KEGG",
                   species = "ibu", 
                   #out.suffix="kegg",
                   kegg.native=T,
                   limit      = list(gene=max(abs(foldchanges.view)), cpd=1),
                   same.layer=T,
                   low = "#d73027", mid = "#fee090", high = "#4575b4",
                   kegg.dir = "kegg_maps")

file.rename(".", paste0("keggmaps/ibu", id.select, '.pathview.png', sep=''))



id.select <- gsub("ibu","", "ibu00310")



paste0("keggmaps/ibu", id.select, '.pathview.png', sep='')
file.move(paste0("keggmaps/ibu", id.select, '.pathview.png', sep=''))

#p.glu_vs_lysine <- plot_volcano_custom(dep, contrast = "Glucose.plus.Acetate_vs_Lysine", 
#
#                                label_size = 2, 
#                               add_names = TRUE) 
#return(p.glu_vs_lysine)


p.glu_vs_lysine <- EnhancedVolcano::EnhancedVolcano(data_results,
                                                    lab = data_results$name,  
                                                    x = 'log2FoldChange',
                                                    y = 'p.val',
                                                    FCcutoff = input$FCcutoff,
                                                    cutoffLineType = 'twodash',
                                                    cutoffLineWidth = 0.8,
                                                    transcriptPointSize = 1.5,
                                                    transcriptLabSize = 4.0,
                                                    colAlpha = 1,
                                                    #DrawConnectors= T,
                                                    legend=c('NS','Log (base 2) fold-change','adj.P value',
                                                             'adj.P value & Log (base 2) fold-change'),
                                                    legendPosition = 'right',
                                                    legendLabSize = 12)
                                                    








ps0 <- pseq()
data_filt <- filter_missval(ps0, thr = input$Filter_missval_Thres)
data_norm <- normalize_vsn(data_filt)
set.seed(2156)
# All possible imputation methods are printed in an error, if an invalid function name is given.
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")

dep <- add_rejections(data_diff_all_contrasts, alpha = input$alpha, lfc = input$Lfc)
data_results <- get_results(dep)
df_wide <- get_df_wide(dep)

foldchanges.view = data_results$Lysine_centered
names(foldchanges.view) = data_results$name
gene <- names(foldchanges.view)[abs(foldchanges.view) > 1.5]
kk <- enrichKEGG(gene         = gene,
                 organism     = 'ibu',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 minGSSize = 5)
kk

kk@result




overlap <- plot_frequency(data_prot)
overlap.fil <- overlap$data
overlap.fil.sub <- subset(overlap.fil, Var1!=0)

p <- ggplot(overlap.fil.sub, aes(x = Var1, y = Freq, fill = Var1)) + 
  geom_col() + scale_fill_grey(start = 0.8, end = 0.2) + 
  labs(title = "Protein identifications overlap", x = "Identified in number of samples", 
       y = "Number of proteins") + theme_DEP2() + theme(legend.position = "none")
plotly::ggplotly(p)

library(tidyverse)
df <- assay(data_prot) %>% data.frame() %>% rownames_to_column() %>% 
  gather(ID, bin, -rowname) %>% mutate(bin = ifelse(is.na(bin), 
                                                    0, 1))
df
stat <- df %>% group_by(rowname) %>% summarize(sum = sum(bin))
overlap$layers

DT::datatable(stat)
