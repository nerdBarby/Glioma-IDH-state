## ML Part -----------------------------------------------------------------
#calling Packages 
library(TCGAbiolinks)
library(sesame)
library(sesameData)
library(SummarizedExperiment)
library(dplyr)
library(data.table)

#methylation counts

methylation <- GDCquery(
  project = "TCGA-LGG",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  access = "open",
  data.type = "Methylation Beta Value",
)
methylationoutput = getResults(methylation)
GDCdownload(methylation)
gdcglio_meth= GDCprepare(methylation,summarizedExperiment = TRUE)
methdata= assay(gdcglio_meth)
################################
methylationdata=as.data.frame(methdata)
methylationdata <- na.omit(methylationdata)
#------------------------------------------------------------------------------
#METADATA
idh_status <- as.data.frame(gdcglio_meth@colData@listData[["paper_IDH.status"]])
lgm_clusters <- as.data.frame(gdcglio_meth@colData@listData[["paper_Pan.Glioma.DNA.Methylation.Cluster"]]) 
#into data frame
anno<-as.data.frame(methylationoutput$cases)
anno <- cbind(anno,idh_status)
anno <- cbind(anno,lgm_clusters)
colnames(anno) <- c("Sample_ids","IDH_status","LGM_stage")
#----------------------------------------------------------------------------
#Saving fileas 
write.csv(methylationdata,"counts_data.csv")
write.csv(anno,"coldata.csv")



###############################################################################################################################
# Read a large CSV file using fread
counts_data<-  fread("counts_data.csv")
# Step 2: Convert the first column to row names efficiently
setDF(counts_data)                      # Convert data.table to data.frame for rownames()
rownames(counts_data) <- counts_data[, 1]   # Set first column as row names
counts_data <- counts_data[, -1] 

#reading metadata
colData<- read.csv("coldata.csv")
colData<- colData[,-1]
#----------------------------------------------------------------------------------------------------
#Preprocessing 
mean_beta <- 0.3 # setting a threshold 

# Filter out the columns where the mean β value is >= 0.3
filtered_data <- counts_data %>%
  filter(if_all(everything(), ~ . < 0.3))
#saving 
write.csv(filtered_data,"filtered_data.csv")

#---------------------------------------------

#merging 
expression_data=fread("filtered_data.csv")
setDF(expression_data)
rownames(expression_data)=expression_data$V1
expression_data=expression_data[,-1]
#transposing 
texp=t(expression_data)
#---------------------------------------------
merged_df=cbind(texp,colData)
merged_df=merged_df[,-55704]
merged_df <- merged_df[ , !(names(merged_df) %in% "X")]
#removing NAs
cleaned_df <- merged_df[!is.na(merged_df$IDH_status), ]
head(cleaned_df)
#saving and heading to the model
write.csv(cleaned_df,"filtered_data.csv")




## BD Part -----------------------------------------------------------------
# Loading packages --------------------------------------------------------
library(readr)
library(org.Hs.eg.db)
library(DESeq2)
library(vsn)

# Reading Data ------------------------------------------------------------

meth_data=read.csv("C:/Users/marie/Downloads/HackBio_stage4/New folder/NEW/filtered_data.csv", row.names = 1)
pheno=read.csv("C:/Users/marie/Downloads/HackBio_stage4/New folder/NEW/coldata.csv", row.names = 1)

table(pheno$IDH_status)

#Removing nulls from pheno:
colnames(meth_data) = gsub("\\.", "-", colnames(meth_data))
pheno = na.omit(pheno)
data = meth_data[, colnames(meth_data) %in% pheno$Sample_ids]
meth = rownames(data)
data = apply(data, 2, function(x) as.numeric(gsub(",", "", x)))
rownames(data) = meth

save(data, meth_data, pheno, file = "C:/Users/marie/Downloads/HackBio_stage4/New folder/NEW/dataANDmetadata.RDATA")
load("C:/Users/marie/Downloads/HackBio_stage4/New folder/NEW/dataANDmetadata.RDATA")



# Differential Methylation Analysis (limma) ------------------------------------------------------------------
M_values <- log2(data / (1 - data))

library(limma)
# Define design matrix
type = as.character(pheno$IDH_status)
design = model.matrix(~0 + factor(type))
colnames(design) = levels(factor(type))

# creating the model and fitting it
contrast = makeContrasts(Mutant - WT, levels = design)
fit = lmFit(as.matrix(M_values), design)
fit2 = contrasts.fit(fit, contrast)
fit2 = eBayes(fit2)

# Get top differentially methylated CpG sites
DMPs=topTable(fit2, adjust.method='fdr', number=999999999,p.value=1,coef = 1)



# Volcano plot ------------------------------------------------------------
require("ggrepel")
require("tidyverse")


# Define the volcano plot function
my_vlocano_plot <- function(tt1, ct) {
  celltype = ct
  pax2wtGgplot1 = tt1 %>%
    rownames_to_column("meth") %>%
    mutate(threshold = ifelse(logFC <= -0.001 & adj.P.Val < 0.05, 'Down',
                              ifelse(logFC >= 0.001 & adj.P.Val < 0.05, 'Up',
                                     ifelse(adj.P.Val >= 0.05, 'Not sig', 'NA'))),
           name = ifelse(-log10(adj.P.Val) > 15, meth, 'NA'))
  
  # Create the ggplot object for the volcano plot
  p = ggplot(pax2wtGgplot1, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = threshold)) +
    scale_color_manual(values = c("blue", "gray", "red")) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text.x = element_text(color = "grey20", size = 15, angle = 45, hjust = .5, vjust = 0.5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_blank(),
          legend.text = element_text(size = 20, face = "plain"),
          strip.text = element_text(size = 20, face = "plain"),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(color = 0),
          axis.line = element_line(colour = "black")
    ) +
    
    ggtitle(celltype) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(
      data = subset(pax2wtGgplot1, name != 'NA'),
      aes(label = name),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )
  
  
  
  # Print the plot
  
  print(p)
}
my_vlocano_plot(DMPs,'Mutant vs. WT')



# Overrepresentation Analysis ---------------------------------------------
#BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", lib = "/content/drive/MyDrive/Colab Notebooks/R_lib")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19, lib = "/content/drive/MyDrive/Colab Notebooks/R_lib")

dmps = DMPs
ds = rownames(dmps)
dmps = cbind(ds, dmps)

# Get annotation data
annotation_data = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

mapped_genes = annotation_data[annotation_data$Name %in% cpg_ids, c("Name", "UCSC_RefGene_Name")]
merged_data = merge(dmps, mapped_genes, by.x = "ds", by.y = "Name", all.x = T)

merged_data$UCSC_RefGene_Name = sapply(strsplit( merged_data$UCSC_RefGene_Name, ";"), function(x) {
  paste(unique(x), collapse = ";")
})

rownames(merged_data) = merged_data$UCSC_RefGene_Name



#install.packages("enrichR")
#library(enrichR)

#options(repr.plot.width = 8, repr.plot.height = 6)

#dbs = c("Reactome_2022")

d1 = rownames(merged_data)[merged_data$adj.P.Val < 0.05]

enriched = enrichr(d1, dbs)

d = plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",
               title = paste0(length(d1), "_", "_", "_", dbs[1]))

print(d)



library(enrichR)

options(repr.plot.width=8,repr.plot.height=6)

dbs <- c("Reactome_2022","GO_Biological_Process_2021","GO_Cellular_Component_2021","GO_Molecular_Function_2021")

d1 = rownames(merged_data)[merged_data$adj.P.Val < 0.05 & merged_data$logFC > 0]
enriched <- enrichr(d1 , dbs)

d=plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",
             
             title=paste0(length(d1),'_','Up','_','_',dbs[1]))
print(d)

d=plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",
             
             title=paste0(length(d1),'_','Up','_','_',dbs[2]))

print(d)



d=plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",
             
             title=paste0(length(d1),'_','Up','_','_',dbs[3]))

print(d)

d=plotEnrich(enriched[[4]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",
             
             title=paste0(length(d1),'_','Up','_',dbs[4]))

print(d)

d2= rownames(merged_data)[merged_data$adj.P.Val < 0.05 & merged_data$logFC < 0]

if(length(d1)==0){skipped_cell=append(i,skipped_cell);next}

enriched <- enrichr(d2 , dbs)



d=plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",
             
             title=paste0(length(d2),'_','Down','_','_',dbs[1]))

print(d)

d=plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",
             
             title=paste0(length(d2),'_','DOWN','_','_',dbs[2]))

print(d)



d=plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",
             
             title=paste0(length(d2),'_','DOWN','_','_',dbs[3]))

print(d)

d=plotEnrich(enriched[[4]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",
             
             title=paste0(length(d2),'_','DOWN','_','_',dbs[4]))

print(d)