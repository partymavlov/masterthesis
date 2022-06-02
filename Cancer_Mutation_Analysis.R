options(repos = c(CRAN = "https://packagemanager.rstudio.com/cran/latest"))
library(plyr)
library(dostats)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(scales)
library(tidyverse)
library(data.table)
library("VennDiagram") 
library(ggVennDiagram)
library(hrbrthemes)
library(plotly)
library(htmlwidgets)
library(ggplotify)
library(patchwork)
library(VariantAnnotation)
library(vcfR)
library(biomaRt)
library(Gviz)
library(maftools)

### Load somatic mutation dataset
somatic_db = read.maf(maf= '/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/mc3.v0.2.8.CONTROLLED.maf.gz')
### Extract statistics and summary
getSampleSummary(somatic_db)
#Shows gene summary.
getGeneSummary(somatic_db)
#shows clinical data associated with samples
getClinicalData(somatic_db)
#Shows all fields in MAF
getFields(somatic_db)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = somatic_db, basename = '/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/mc3.v0.2.8.CONTROLLED.maf.summary_new_unfiltered')
pdf(file='/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/mc3.v0.2.8.CONTROLLED.maf.summary_new_unfiltered.pdf', width = 10,
    height = 6)
plotmafSummary(maf = unfilteredMAF, addStat = 'mean', dashboard = TRUE, fs = 0.7, titvRaw = FALSE)
dev.off()
pdf(file='/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/mc3.v0.2.8.CONTROLLED.maf.oncoplot_new_unfiltered.pdf', width = 10,
    height = 6)
oncoplot(maf = somatic_db, top = 10)
dev.off()
rm(somatic_db)
gc()

### Somatic mutation dataset filtering
filteredMAF <- subsetMaf(somatic_db, fields = c("Chromosome", "Start_Position", "End_Position", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", 
                                          "Tumor_Seq_Allele2",  "dbSNP_RS", "FILTER", "t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count", "n_alt_count", "Strand", "Variant_Classification",   "all_effects", "Allele", "Gene", "Feature", "SOMATIC", "ExAC_AF", "Entrez_Gene_Id"), query = "!FILTER == 'nonpreferredpair,wga'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'nonpreferredpair,oxog'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'nonpreferredpair'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'oxog'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'common_in_exac,nonpreferredpair'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'gapfiller,nonpreferredpair'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'native_wga_mix,nonpreferredpair'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'StrandBias,nonpreferredpair'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'StrandBias,oxog'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'oxog,wga'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'gapfiller,oxog'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'gapfiller,native_wga_mix,oxog'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'common_in_exac,oxog'")
filteredMAF <- subsetMaf(filteredMAF, query = "!FILTER == 'StrandBias,nonpreferredpair,oxog'")
### Extract statistics and summary of filtered dataset
getSampleSummary(filteredMAF)
#Shows gene summary.
getGeneSummary(filteredMAF)
#shows clinical data associated with samples
getClinicalData(filteredMAF)
#Shows all fields in MAF
getFields(filteredMAF)
#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = filteredMAF, basename = '/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/mc3.v0.2.8.CONTROLLED.maf.summary_new')
pdf(file='/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/mc3.v0.2.8.CONTROLLED.maf.summary_new.pdf', width = 10,
    height = 6)
plotmafSummary(maf = filteredMAF, addStat = 'mean', dashboard = TRUE, fs = 0.7, titvRaw = FALSE)
dev.off()

write.table(filteredMAF@data, "/Users/mcpftw/Documents/Master\ Bioinformatik/Master\ Thesis/Data/TCGA/filtered_controlled_MAF.txt", sep = "\t", quote = F, row.names = F)

unfilteredMAF<- read.maf(maf= '/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/mc3.v0.2.8.CONTROLLED.maf.summary_new_unfiltered_maftools.maf')
somatic_filtered <- read.csv("/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/filtered_controlled_MAF.txt", sep = "\t")

levels_var_class = c("Frameshift Deletion", "Frameshift Insertion", "Inframe Deletion",
                     "Inframe Insertion", "Missense Mutation", "Nonsense Mutation", "Nonstop Mutation",
                     "Splice Site", "Translation Start Site")
colnames_Var_class = c("Tumor_Sample_Barcode","Frameshift Deletion", "Frameshift Insertion", "Inframe Deletion",
                       "Inframe Insertion", "Missense Mutation", "Nonsense Mutation", "Nonstop Mutation",
                       "Splice Site", "Translation Start Site", "total")      
colnames_Gene_summary= c("Hugo_Symbol","Frameshift Deletion", "Frameshift Insertion", "Inframe Deletion",
                       "Inframe Insertion", "Missense Mutation", "Nonsense Mutation", "Nonstop Mutation",
                       "Splice Site", "Translation Start Site", "total", "MutatedSamples", "AlteredSamples")      

germline_header <- scanVcfHeader("/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/0fbc9ce6-aed8-4bc6-a24f-9cf8229654a1/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz")
germline_header <- info(germline_header)
germline_header <- rbind(rownames(germline_header), germline_header$Description)
germline_header <- t(germline_header)
mysvp = ScanVcfParam(info=c("END", "SVLEN", "SVTYPE", "AC", "AF", "AN"),geno = NA)
germline_vcf <- readVcf("/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/0fbc9ce6-aed8-4bc6-a24f-9cf8229654a1/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz", param = mysvp)
germline_filtered <- germline_vcf[!is.na(info(germline_vcf)$AF),]
write.table(germline_filtered, file = "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/TCGA_germline_filtered_AF_END.txt", sep = "\t", quote = F, row.names = F)
germline_filtered <- germline_filtered[,c(1,3,2)]
write.table(germline_filtered, file = "/Users/mcpftw/annovar/example/TCGA_germline_filtered_AF_END.avinput", sep = "\t", quote = F, row.names = F, col.names = F)

germline_filtered <- read.csv("/Users/mcpftw/annovar/example/TCGA_germline_END.avinput", sep = "\t")
germline_filtered <- germline_filtered[!is.na(germline_filtered$AF),]
germline_filtered$END[is.na(germline_filtered$END)] <- germline_filtered$POS[is.na(germline_filtered$END)]
germline_filtered$X <- NULL
write.table(germline_filtered, file = "/Users/mcpftw/annovar/example/TCGA_germline_filtered_AF_END.avinput", sep = "\t", quote = F, row.names = F, col.names = F)
germline_filtered <- read.csv("/Users/mcpftw/annovar/example/TCGA_germline_filtered_AF_END.avinput", sep = "\t", header = F)
germline_annotations <- read.csv("/Users/mcpftw/annovar/myannoTCGA.hg19_multianno.txt", sep = "\t")
germline_filtered$'Func.refGene' <- germline_annotations$Func.refGene
germline_filtered$'Gene.refGene' <- germline_annotations$Gene.refGene
germline_filtered$'GeneDetail.refGene' <- germline_annotations$GeneDetail.refGene
germline_filtered$'ExonicFunc.refGene' <- germline_annotations$ExonicFunc.refGene
germline_filtered$'AAChange.refGene' <- germline_annotations$AAChange.refGene
germline_filtered$'avsnp147' <- germline_annotations$avsnp147
write.table(germline_filtered, file = "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/TCGA_germline_annotated.txt", sep = "\t", quote = F, row.names = F)

germline_MAF <- annovarToMaf(annovar = "/Users/mcpftw/annovar/myannoTCGA.hg19_multianno.txt", refBuild = 'hg19', table = 'refGene', MAFobj = TRUE)
write.mafSummary(maf = germline_MAF, basename = '/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/TCGA_germline.summary_new')
pdf(file='/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/germline_maf.summary_new_v2.pdf', width = 10,
    height = 6)
plotmafSummary(maf = germline_MAF, addStat = 'mean', dashboard = TRUE, fs = 0.7, titvRaw = FALSE)
dev.off()
write.table(germline_MAF@data, "/Users/mcpftw/Documents/Master\ Bioinformatik/Master\ Thesis/Data/TCGA/germline_controlled_MAF.txt", sep = "\t", quote = F, row.names = F)
germline_MAF = read.maf(maf= '/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/TCGA_germline.summary_new_maftools.maf')


germline_filtered <- read.csv("/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/TCGA_germline_annotated.txt", sep = "\t")
top_genes <- as.data.frame(table(germline_filtered$Gene.refGene))
top_genes <- top_genes[order(-top_genes$Freq),]
top_types <- as.data.frame(table(germline_filtered$TYPE))
top_func <- as.data.frame(table(germline_filtered$Func.refGene))
top_exonic_func <- as.data.frame(table(germline_filtered$ExonicFunc.refGene))
top_AAchange <- as.data.frame(table(paste(germline_filtered$REF, germline_filtered$ALT)))
top_AAchange <- top_AAchange[order(-top_AAchange$Freq),]
top_AAchange <- top_AAchange[1:12,]
top_AAchange$Var1 <- gsub(" ", ">", top_AAchange$Var1)
aaplot <- ggplot(top_AAchange,aes(x=Var1, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(legend.position="none") +
  coord_flip() +
  labs(x=NULL,y=NULL) +
  scale_y_continuous(labels = comma) + 
  ggtitle("SNV class")
aaplot

top_func <- top_func[-1,]
top_func <- top_func[order(-top_func$Freq),]
funcplot <- ggplot(top_func, aes(x= reorder(Var1, -desc(Freq)), y=Freq, fill=Var1, label=Freq)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(legend.position="none") +
  coord_flip() +
  labs(x=NULL,y=NULL) +
  scale_y_continuous(labels = comma) +
  geom_text(size = 3, position = position_stack(vjust = 0.9)) +
  ggtitle("Variant Location")
funcplot

top_exonic_func <- top_exonic_func[-1,]
top_exonic_func <- top_exonic_func[order(-top_exonic_func$Freq),]
exonicfuncplot <- ggplot(top_exonic_func, aes(x= reorder(Var1, -desc(Freq)), y=Freq, fill=Var1, label=Freq)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(legend.position="none") +
  coord_flip() +
  labs(x=NULL,y=NULL) +
  scale_y_continuous(labels = comma) +
  geom_text(size = 3, position = position_stack(vjust = 0.9)) +
  ggtitle("Variant Classification")
exonicfuncplot

mutationtypeplot <- ggplot(top_types, aes(x= reorder(Var1,-desc(Freq)), y=Freq, fill=Var1, label=Freq)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(legend.position="none") +
  coord_flip() +
  labs(x=NULL,y=NULL) +
  scale_y_continuous(labels = comma) +
  geom_text(size = 3, position = position_stack(vjust = 0.9)) +
  ggtitle("Variant Type")
mutationtypeplot

top_genes <- top_genes[-1,]
top_genes <- top_genes[1:10,]
geneplot <-
  ggplot(top_genes, aes(
    x = reorder(Var1, -desc(Freq)),
    y = Freq,
    fill = Var1,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  ggtitle("Top 10 mutated genes")
geneplot

allele_counts <- as.data.frame(as.numeric(germline_filtered$AC))
#allele_test <- as.data.frame(allele_counts[1:100000, ])
allele_counts <- allele_counts[order(allele_counts$`as.numeric(germline_filtered$AC)`)]
alleleplot <- ggplot(allele_counts, aes(x= `as.numeric(germline_filtered$AC)`)) + 
                       geom_density() +
              xlim(0,50) + xlab("alleles") +ggtitle("Allele counts per mutation") + theme_minimal()
              #theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))
alleleplot

germline_plots <- ggarrange(funcplot, aaplot, geneplot, alleleplot,
                            labels = c("Mutation function distribution", "Amino acid substitution distribution", "Top 10 mutated genes", "Allele count distribution"),
                            ncol = 2, nrow = 2)



# 0-based to 1-based coordinates
germline_filtered$POS <- germline_filtered$POS + 1


germlineDbNsfp1 <- as.data.frame(cbind(germline_filtered$CHROM, germline_filtered$POS))
write.table(germlineDbNsfp1, file = "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/germline_forNSFP_pos_1based.txt", sep = "\t", quote = F, row.names = F, col.names = F)
germline_filtered <- germline_filtered[germline_filtered$TYPE == "SNP",]
germlineDbNsfp2 <- as.data.frame(cbind(germline_filtered$CHROM, germline_filtered$POS, germline_filtered$REF, germline_filtered$ALT))
write.table(germlineDbNsfp2, file = "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/germline_forNSFP_pos_SNP_only_1based.txt", sep = "\t", quote = F, row.names = F, col.names = F)

somatic_filtered <- read.csv("/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/filtered_controlled_MAF.txt", sep = "\t")
somaticDbNsfp1 <- as.data.frame(cbind(somatic_filtered$Chromosome, somatic_filtered$Start_Position))
write.table(somaticDbNsfp1, file = "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/somatic_forNSFP_pos.txt", sep = "\t", quote = F, row.names = F, col.names = F)
somatic_filtered <- somatic_filtered[somatic_filtered$Variant_Type == "SNP",]
somaticDbNsfp2 <- as.data.frame(cbind(somatic_filtered$Chromosome, somatic_filtered$Start_Position, somatic_filtered$Reference_Allele, somatic_filtered$Tumor_Seq_Allele2))
write.table(somaticDbNsfp2, file = "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/somatic_forNSFP_pos_SNP_only.txt", sep = "\t", quote = F, row.names = F, col.names = F)


somatic_NSFP <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/somatic_NSFP_pos.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "SIFT4G_converted_rankscore",
      "Polyphen2_HVAR_rankscore",
      "MetaSVM_rankscore",
      "MetaLR_rankscore",
      "MetaRNN_rankscore"
      
    )
  )


somatic_NSFP <- somatic_NSFP[somatic_NSFP$SIFT4G_converted_rankscore != "."]
somatic_NSFP$SIFT4G_converted_rankscore <- as.numeric(somatic_NSFP$SIFT4G_converted_rankscore)
somatic_NSFP <- somatic_NSFP[somatic_NSFP$Polyphen2_HVAR_rankscore != "."]
somatic_NSFP$Polyphen2_HVAR_rankscore <- as.numeric(somatic_NSFP$Polyphen2_HVAR_rankscore)
somatic_NSFP <- somatic_NSFP[somatic_NSFP$MetaSVM_rankscore != "."]
somatic_NSFP$MetaSVM_rankscore <- as.numeric(somatic_NSFP$MetaSVM_rankscore)
somatic_NSFP <- somatic_NSFP[somatic_NSFP$MetaLR_rankscore != "."]
somatic_NSFP$MetaLR_rankscore <- as.numeric(somatic_NSFP$MetaLR_rankscore)
somatic_NSFP <- somatic_NSFP[somatic_NSFP$MetaRNN_rankscore != "."]
somatic_NSFP$MetaRNN_rankscore <- as.numeric(somatic_NSFP$MetaRNN_rankscore)
somatic_NSFP = somatic_NSFP %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_rankscore`))
somatic_NSFP = somatic_NSFP %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaLR_rankscore`))

somatic_NSFP_raw <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/somatic_NSFP_pos.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "SIFT4G_score",
      "Polyphen2_HVAR_score",
      "MetaSVM_score",
      "MetaLR_score",
      "MetaRNN_score"
    )
  )

somatic_NSFP_raw$SIFT4G_score <- sapply(str_split(somatic_NSFP_raw$SIFT4G_score, ";"), tail, 1)
somatic_NSFP_raw$Polyphen2_HVAR_score <- sapply(str_split(somatic_NSFP_raw$Polyphen2_HVAR_score, ";"), tail, 1)
somatic_NSFP_raw$MetaSVM_score <- sapply(str_split(somatic_NSFP_raw$MetaSVM_score, ";"), tail, 1)
somatic_NSFP_raw$MetaLR_score <- sapply(str_split(somatic_NSFP_raw$MetaLR_score, ";"), tail, 1)
somatic_NSFP_raw$MetaRNN_score <- sapply(str_split(somatic_NSFP_raw$MetaRNN_score, ";"), tail, 1)
somatic_NSFP_raw <- somatic_NSFP_raw[somatic_NSFP_raw$SIFT4G_score != "."]
somatic_NSFP_raw <- somatic_NSFP_raw[somatic_NSFP_raw$Polyphen2_HVAR_score != "."]
somatic_NSFP_raw <- somatic_NSFP_raw[somatic_NSFP_raw$MetaSVM_score != "."]
somatic_NSFP_raw <- somatic_NSFP_raw[somatic_NSFP_raw$MetaLR_score != "."]
somatic_NSFP_raw <- somatic_NSFP_raw[somatic_NSFP_raw$MetaRNN_score != "."]
somatic_NSFP_raw$SIFT4G_score <- as.numeric(somatic_NSFP_raw$SIFT4G_score)
somatic_NSFP_raw$Polyphen2_HVAR_score <- as.numeric(somatic_NSFP_raw$Polyphen2_HVAR_score)
somatic_NSFP_raw$MetaSVM_score <- as.numeric(somatic_NSFP_raw$MetaSVM_score)
somatic_NSFP_raw$MetaLR_score <- as.numeric(somatic_NSFP_raw$MetaLR_score)
somatic_NSFP_raw$MetaRNN_score <- as.numeric(somatic_NSFP_raw$MetaRNN_score)
somatic_NSFP_raw = somatic_NSFP_raw %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_score`))


somatic_NSFP_pred <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/somatic_NSFP_pos.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "SIFT4G_pred",
      "Polyphen2_HVAR_pred",
      "MetaSVM_pred",
      "MetaLR_pred",
      "MetaRNN_pred"
    )
  )
somatic_NSFP_pred$SIFT4G_pred <- substr(somatic_NSFP_pred$SIFT4G_pred, nchar(somatic_NSFP_pred$SIFT4G_pred), nchar(somatic_NSFP_pred$SIFT4G_pred))
somatic_NSFP_pred$Polyphen2_HVAR_pred <- substr(somatic_NSFP_pred$Polyphen2_HVAR_pred, nchar(somatic_NSFP_pred$Polyphen2_HVAR_pred), nchar(somatic_NSFP_pred$Polyphen2_HVAR_pred))
somatic_NSFP_pred$MetaSVM_pred <- substr(somatic_NSFP_pred$MetaSVM_pred, nchar(somatic_NSFP_pred$MetaSVM_pred), nchar(somatic_NSFP_pred$MetaSVM_pred))
somatic_NSFP_pred$MetaLR_pred <- substr(somatic_NSFP_pred$MetaLR_pred, nchar(somatic_NSFP_pred$MetaLR_pred), nchar(somatic_NSFP_pred$MetaLR_pred))
somatic_NSFP_pred$MetaRNN_pred <- substr(somatic_NSFP_pred$MetaRNN_pred, nchar(somatic_NSFP_pred$MetaRNN_pred), nchar(somatic_NSFP_pred$MetaRNN_pred))
somatic_NSFP_pred <- somatic_NSFP_pred[somatic_NSFP_pred$SIFT4G_pred != "."]
somatic_NSFP_pred <- somatic_NSFP_pred[somatic_NSFP_pred$Polyphen2_HVAR_pred != "."]
somatic_NSFP_pred <- somatic_NSFP_pred[somatic_NSFP_pred$MetaSVM_pred != "."]
somatic_NSFP_pred <- somatic_NSFP_pred[somatic_NSFP_pred$MetaLR_pred != "."]
somatic_NSFP_pred <- somatic_NSFP_pred[somatic_NSFP_pred$MetaRNN_pred != "."]

somatic_NSFP_SNP_only <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/somatic_NSFP_pos_SNP_only.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "SIFT4G_converted_rankscore",
      "Polyphen2_HVAR_rankscore",
      "MetaSVM_rankscore",
      "MetaLR_rankscore",
      "MetaRNN_rankscore"
    )
  )
somatic_NSFP_SNP_only <- somatic_NSFP_SNP_only[somatic_NSFP_SNP_only$SIFT4G_converted_rankscore != "."]
somatic_NSFP_SNP_only <- somatic_NSFP_SNP_only[somatic_NSFP_SNP_only$Polyphen2_HVAR_rankscore != "."]
somatic_NSFP_SNP_only <- somatic_NSFP_SNP_only[somatic_NSFP_SNP_only$MetaSVM_rankscore != "."]
somatic_NSFP_SNP_only <- somatic_NSFP_SNP_only[somatic_NSFP_SNP_only$MetaLR_rankscore != "."]
somatic_NSFP_SNP_only <- somatic_NSFP_SNP_only[somatic_NSFP_SNP_only$MetaRNN_rankscore != "."]
somatic_NSFP_SNP_only$SIFT4G_converted_rankscore <- as.numeric(somatic_NSFP_SNP_only$SIFT4G_converted_rankscore)
somatic_NSFP_SNP_only$Polyphen2_HVAR_rankscore <- as.numeric(somatic_NSFP_SNP_only$Polyphen2_HVAR_rankscore)
somatic_NSFP_SNP_only$MetaSVM_rankscore <- as.numeric(somatic_NSFP_SNP_only$MetaSVM_rankscore)
somatic_NSFP_SNP_only$MetaLR_rankscore <- as.numeric(somatic_NSFP_SNP_only$MetaLR_rankscore)
somatic_NSFP_SNP_only$MetaRNN_rankscore <- as.numeric(somatic_NSFP_SNP_only$MetaRNN_rankscore)
somatic_NSFP_SNP_only = somatic_NSFP_SNP_only %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_rankscore`))

somatic_NSFP_raw_SNP <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/somatic_NSFP_pos_SNP_only.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "SIFT4G_score",
      "Polyphen2_HVAR_score",
      "MetaSVM_score",
      "MetaLR_score",
      "MetaRNN_score"
    )
  )
somatic_NSFP_raw_SNP$SIFT4G_score <- sapply(str_split(somatic_NSFP_raw_SNP$SIFT4G_score, ";"), tail, 1)
somatic_NSFP_raw_SNP$Polyphen2_HVAR_score <- sapply(str_split(somatic_NSFP_raw_SNP$Polyphen2_HVAR_score, ";"), tail, 1)
somatic_NSFP_raw_SNP$MetaSVM_score <- sapply(str_split(somatic_NSFP_raw_SNP$MetaSVM_score, ";"), tail, 1)
somatic_NSFP_raw_SNP$MetaLR_score <- sapply(str_split(somatic_NSFP_raw_SNP$MetaLR_score, ";"), tail, 1)
somatic_NSFP_raw_SNP$MetaRNN_score <- sapply(str_split(somatic_NSFP_raw_SNP$MetaRNN_score, ";"), tail, 1)
somatic_NSFP_raw_SNP <- somatic_NSFP_raw_SNP[somatic_NSFP_raw_SNP$SIFT4G_score != "."]
somatic_NSFP_raw_SNP <- somatic_NSFP_raw_SNP[somatic_NSFP_raw_SNP$Polyphen2_HVAR_score != "."]
somatic_NSFP_raw_SNP <- somatic_NSFP_raw_SNP[somatic_NSFP_raw_SNP$MetaSVM_score != "."]
somatic_NSFP_raw_SNP <- somatic_NSFP_raw_SNP[somatic_NSFP_raw_SNP$MetaLR_score != "."]
somatic_NSFP_raw_SNP <- somatic_NSFP_raw_SNP[somatic_NSFP_raw_SNP$MetaRNN_score != "."]
somatic_NSFP_raw_SNP$SIFT4G_score <- as.numeric(somatic_NSFP_raw_SNP$SIFT4G_score)
somatic_NSFP_raw_SNP$Polyphen2_HVAR_score <- as.numeric(somatic_NSFP_raw_SNP$Polyphen2_HVAR_score)
somatic_NSFP_raw_SNP$MetaSVM_score <- as.numeric(somatic_NSFP_raw_SNP$MetaSVM_score)
somatic_NSFP_raw_SNP$MetaLR_score <- as.numeric(somatic_NSFP_raw_SNP$MetaLR_score)
somatic_NSFP_raw_SNP$MetaRNN_score <- as.numeric(somatic_NSFP_raw_SNP$MetaRNN_score)
somatic_NSFP_raw_SNP = somatic_NSFP_raw_SNP %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_score`))

somatic_NSFP_pred_SNP <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/somatic_NSFP_pos_SNP_only.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "SIFT4G_pred",
      "Polyphen2_HVAR_pred",
      "MetaSVM_pred",
      "MetaLR_pred",
      "MetaRNN_pred"
    )
  )
somatic_NSFP_pred_SNP$SIFT4G_pred <- substr(somatic_NSFP_pred_SNP$SIFT4G_pred, nchar(somatic_NSFP_pred_SNP$SIFT4G_pred), nchar(somatic_NSFP_pred_SNP$SIFT4G_pred))
somatic_NSFP_pred_SNP$Polyphen2_HVAR_pred <- substr(somatic_NSFP_pred_SNP$Polyphen2_HVAR_pred, nchar(somatic_NSFP_pred_SNP$Polyphen2_HVAR_pred), nchar(somatic_NSFP_pred_SNP$Polyphen2_HVAR_pred))
somatic_NSFP_pred_SNP$MetaSVM_pred <- substr(somatic_NSFP_pred_SNP$MetaSVM_pred, nchar(somatic_NSFP_pred_SNP$MetaSVM_pred), nchar(somatic_NSFP_pred_SNP$MetaSVM_pred))
somatic_NSFP_pred_SNP$MetaLR_pred <- substr(somatic_NSFP_pred_SNP$MetaLR_pred, nchar(somatic_NSFP_pred_SNP$MetaLR_pred), nchar(somatic_NSFP_pred_SNP$MetaLR_pred))
somatic_NSFP_pred_SNP$MetaRNN_pred <- substr(somatic_NSFP_pred_SNP$MetaRNN_pred, nchar(somatic_NSFP_pred_SNP$MetaRNN_pred), nchar(somatic_NSFP_pred_SNP$MetaRNN_pred))
somatic_NSFP_pred_SNP <- somatic_NSFP_pred_SNP[somatic_NSFP_pred_SNP$SIFT4G_pred != "."]
somatic_NSFP_pred_SNP <- somatic_NSFP_pred_SNP[somatic_NSFP_pred_SNP$Polyphen2_HVAR_pred != "."]
somatic_NSFP_pred_SNP <- somatic_NSFP_pred_SNP[somatic_NSFP_pred_SNP$MetaSVM_pred != "."]
somatic_NSFP_pred_SNP <- somatic_NSFP_pred_SNP[somatic_NSFP_pred_SNP$MetaLR_pred != "."]
somatic_NSFP_pred_SNP <- somatic_NSFP_pred_SNP[somatic_NSFP_pred_SNP$MetaRNN_pred != "."]


germline_NSFP <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/germline_NSFP_pos.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "SIFT4G_converted_rankscore",
      "Polyphen2_HVAR_rankscore",
      "MetaSVM_rankscore",
      "MetaLR_rankscore",
      "MetaRNN_rankscore"
    )
  )
germline_NSFP <- germline_NSFP[germline_NSFP$SIFT4G_converted_rankscore != "."]
germline_NSFP <- germline_NSFP[germline_NSFP$Polyphen2_HVAR_rankscore != "."]
germline_NSFP <- germline_NSFP[germline_NSFP$MetaSVM_rankscore != "."]
germline_NSFP <- germline_NSFP[germline_NSFP$MetaLR_rankscore != "."]
germline_NSFP <- germline_NSFP[germline_NSFP$MetaRNN_rankscore != "."]
germline_NSFP$SIFT4G_converted_rankscore <- as.numeric(germline_NSFP$SIFT4G_converted_rankscore)
germline_NSFP$Polyphen2_HVAR_rankscore <- as.numeric(germline_NSFP$Polyphen2_HVAR_rankscore)
germline_NSFP$MetaSVM_rankscore <- as.numeric(germline_NSFP$MetaSVM_rankscore)
germline_NSFP$MetaLR_rankscore <- as.numeric(germline_NSFP$MetaLR_rankscore)
germline_NSFP$MetaRNN_rankscore <- as.numeric(germline_NSFP$MetaRNN_rankscore)
germline_NSFP = germline_NSFP %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_rankscore`))

germline_NSFP_raw <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/germline_NSFP_pos.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "SIFT4G_score",
      "Polyphen2_HVAR_score",
      "MetaSVM_score",
      "MetaLR_score",
      "MetaRNN_score"
    )
  )

germline_NSFP_raw$SIFT4G_score <- sapply(str_split(germline_NSFP_raw$SIFT4G_score, ";"), tail, 1)
germline_NSFP_raw$Polyphen2_HVAR_score <- sapply(str_split(germline_NSFP_raw$Polyphen2_HVAR_score, ";"), tail, 1)
germline_NSFP_raw$MetaSVM_score <- sapply(str_split(germline_NSFP_raw$MetaSVM_score, ";"), tail, 1)
germline_NSFP_raw$MetaLR_score <- sapply(str_split(germline_NSFP_raw$MetaLR_score, ";"), tail, 1)
germline_NSFP_raw$MetaRNN_score <- sapply(str_split(germline_NSFP_raw$MetaRNN_score, ";"), tail, 1)
germline_NSFP_raw <- germline_NSFP_raw[germline_NSFP_raw$SIFT4G_score != "."]
germline_NSFP_raw <- germline_NSFP_raw[germline_NSFP_raw$Polyphen2_HVAR_score != "."]
germline_NSFP_raw <- germline_NSFP_raw[germline_NSFP_raw$MetaSVM_score != "."]
germline_NSFP_raw <- germline_NSFP_raw[germline_NSFP_raw$MetaLR_score != "."]
germline_NSFP_raw <- germline_NSFP_raw[germline_NSFP_raw$MetaRNN_score != "."]
germline_NSFP_raw$SIFT4G_score <- as.numeric(germline_NSFP_raw$SIFT4G_score)
germline_NSFP_raw$Polyphen2_HVAR_score <- as.numeric(germline_NSFP_raw$Polyphen2_HVAR_score)
germline_NSFP_raw$MetaSVM_score <- as.numeric(germline_NSFP_raw$MetaSVM_score)
germline_NSFP_raw$MetaLR_score <- as.numeric(germline_NSFP_raw$MetaLR_score)
germline_NSFP_raw$MetaRNN_score <- as.numeric(germline_NSFP_raw$MetaRNN_score)
germline_NSFP_raw = germline_NSFP_raw %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_score`))

germline_NSFP_pred <- fread(
  "/Users/mcpftw/dbNSFP4.2a/germline_NSFP_pos.txt",
  select = c(
    "#chr",
    "pos(1-based)",
    "ref",
    "alt",
    "SIFT4G_pred",
    "Polyphen2_HVAR_pred",
    "MetaSVM_pred",
    "MetaLR_pred",
    "MetaRNN_pred"
  )
)
germline_NSFP_pred$SIFT4G_pred <- substr(germline_NSFP_pred$SIFT4G_pred, nchar(germline_NSFP_pred$SIFT4G_pred), nchar(germline_NSFP_pred$SIFT4G_pred))
germline_NSFP_pred$Polyphen2_HVAR_pred <- substr(germline_NSFP_pred$Polyphen2_HVAR_pred, nchar(germline_NSFP_pred$Polyphen2_HVAR_pred), nchar(germline_NSFP_pred$Polyphen2_HVAR_pred))
germline_NSFP_pred$MetaSVM_pred <- substr(germline_NSFP_pred$MetaSVM_pred, nchar(germline_NSFP_pred$MetaSVM_pred), nchar(germline_NSFP_pred$MetaSVM_pred))
germline_NSFP_pred$MetaLR_pred <- substr(germline_NSFP_pred$MetaLR_pred, nchar(germline_NSFP_pred$MetaLR_pred), nchar(germline_NSFP_pred$MetaLR_pred))
germline_NSFP_pred$MetaRNN_pred <- substr(germline_NSFP_pred$MetaRNN_pred, nchar(germline_NSFP_pred$MetaRNN_pred), nchar(germline_NSFP_pred$MetaRNN_pred))
germline_NSFP_pred <- germline_NSFP_pred[germline_NSFP_pred$SIFT4G_pred != "."]
germline_NSFP_pred <- germline_NSFP_pred[germline_NSFP_pred$Polyphen2_HVAR_pred != "."]
germline_NSFP_pred <- germline_NSFP_pred[germline_NSFP_pred$MetaSVM_pred != "."]
germline_NSFP_pred <- germline_NSFP_pred[germline_NSFP_pred$MetaLR_pred != "."]
germline_NSFP_pred <- germline_NSFP_pred[germline_NSFP_pred$MetaRNN_pred != "."]


germline_NSFP_SNP_only <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/germline_NSFP_pos_SNP_only.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "SIFT4G_converted_rankscore",
      "Polyphen2_HVAR_rankscore",
      "MetaSVM_rankscore",
      "MetaLR_rankscore",
      "MetaRNN_rankscore"
    )
  )
germline_NSFP_SNP_only <- germline_NSFP_SNP_only[germline_NSFP_SNP_only$SIFT4G_converted_rankscore != "."]
germline_NSFP_SNP_only <- germline_NSFP_SNP_only[germline_NSFP_SNP_only$Polyphen2_HVAR_rankscore != "."]
germline_NSFP_SNP_only <- germline_NSFP_SNP_only[germline_NSFP_SNP_only$MetaSVM_rankscore != "."]
germline_NSFP_SNP_only <- germline_NSFP_SNP_only[germline_NSFP_SNP_only$MetaLR_rankscore != "."]
germline_NSFP_SNP_only <- germline_NSFP_SNP_only[germline_NSFP_SNP_only$MetaRNN_rankscore != "."]
germline_NSFP_SNP_only$SIFT4G_converted_rankscore <- as.numeric(germline_NSFP_SNP_only$SIFT4G_converted_rankscore)
germline_NSFP_SNP_only$Polyphen2_HVAR_rankscore <- as.numeric(germline_NSFP_SNP_only$Polyphen2_HVAR_rankscore)
germline_NSFP_SNP_only$MetaSVM_rankscore <- as.numeric(germline_NSFP_SNP_only$MetaSVM_rankscore)
germline_NSFP_SNP_only$MetaLR_rankscore <- as.numeric(germline_NSFP_SNP_only$MetaLR_rankscore)
germline_NSFP_SNP_only$MetaRNN_rankscore <- as.numeric(germline_NSFP_SNP_only$MetaRNN_rankscore)
germline_NSFP_SNP_only = germline_NSFP_SNP_only %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_rankscore`))

germline_NSFP_raw_SNP <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/germline_NSFP_pos_SNP_only.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "SIFT4G_score",
      "Polyphen2_HVAR_score",
      "MetaSVM_score",
      "MetaLR_score",
      "MetaRNN_score"
    )
  )
germline_NSFP_raw_SNP$SIFT4G_score <- sapply(str_split(germline_NSFP_raw_SNP$SIFT4G_score, ";"), tail, 1)
germline_NSFP_raw_SNP$Polyphen2_HVAR_score <- sapply(str_split(germline_NSFP_raw_SNP$Polyphen2_HVAR_score, ";"), tail, 1)
germline_NSFP_raw_SNP$MetaSVM_score <- sapply(str_split(germline_NSFP_raw_SNP$MetaSVM_score, ";"), tail, 1)
germline_NSFP_raw_SNP$MetaLR_score <- sapply(str_split(germline_NSFP_raw_SNP$MetaLR_score, ";"), tail, 1)
germline_NSFP_raw_SNP$MetaRNN_score <- sapply(str_split(germline_NSFP_raw_SNP$MetaRNN_score, ";"), tail, 1)
germline_NSFP_raw_SNP <- germline_NSFP_raw_SNP[germline_NSFP_raw_SNP$SIFT4G_score != "."]
germline_NSFP_raw_SNP <- germline_NSFP_raw_SNP[germline_NSFP_raw_SNP$Polyphen2_HVAR_score != "."]
germline_NSFP_raw_SNP <- germline_NSFP_raw_SNP[germline_NSFP_raw_SNP$MetaSVM_score != "."]
germline_NSFP_raw_SNP <- germline_NSFP_raw_SNP[germline_NSFP_raw_SNP$MetaLR_score != "."]
germline_NSFP_raw_SNP <- germline_NSFP_raw_SNP[germline_NSFP_raw_SNP$MetaRNN_score != "."]
germline_NSFP_raw_SNP$SIFT4G_score <- as.numeric(germline_NSFP_raw_SNP$SIFT4G_score)
germline_NSFP_raw_SNP$Polyphen2_HVAR_score <- as.numeric(germline_NSFP_raw_SNP$Polyphen2_HVAR_score)
germline_NSFP_raw_SNP$MetaSVM_score <- as.numeric(germline_NSFP_raw_SNP$MetaSVM_score)
germline_NSFP_raw_SNP$MetaLR_score <- as.numeric(germline_NSFP_raw_SNP$MetaLR_score)
germline_NSFP_raw_SNP$MetaRNN_score <- as.numeric(germline_NSFP_raw_SNP$MetaRNN_score)
germline_NSFP_raw_SNP = germline_NSFP_raw_SNP %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_score`))

germline_NSFP_pred_SNP <- fread(
  "/Users/mcpftw/dbNSFP4.2a/germline_NSFP_pos.txt",
  select = c(
    "#chr",
    "pos(1-based)",
    "ref",
    "alt",
    "SIFT4G_pred",
    "Polyphen2_HVAR_pred",
    "MetaSVM_pred",
    "MetaLR_pred",
    "MetaRNN_pred"
  )
)
germline_NSFP_pred_SNP$SIFT4G_pred <- substr(germline_NSFP_pred_SNP$SIFT4G_pred, nchar(germline_NSFP_pred_SNP$SIFT4G_pred), nchar(germline_NSFP_pred_SNP$SIFT4G_pred))
germline_NSFP_pred_SNP$Polyphen2_HVAR_pred <- substr(germline_NSFP_pred_SNP$Polyphen2_HVAR_pred, nchar(germline_NSFP_pred_SNP$Polyphen2_HVAR_pred), nchar(germline_NSFP_pred_SNP$Polyphen2_HVAR_pred))
germline_NSFP_pred_SNP$MetaSVM_pred <- substr(germline_NSFP_pred_SNP$MetaSVM_pred, nchar(germline_NSFP_pred_SNP$MetaSVM_pred), nchar(germline_NSFP_pred_SNP$MetaSVM_pred))
germline_NSFP_pred_SNP$MetaLR_pred <- substr(germline_NSFP_pred_SNP$MetaLR_pred, nchar(germline_NSFP_pred_SNP$MetaLR_pred), nchar(germline_NSFP_pred_SNP$MetaLR_pred))
germline_NSFP_pred_SNP$MetaRNN_pred <- substr(germline_NSFP_pred_SNP$MetaRNN_pred, nchar(germline_NSFP_pred_SNP$MetaRNN_pred), nchar(germline_NSFP_pred_SNP$MetaRNN_pred))
germline_NSFP_pred_SNP <- germline_NSFP_pred_SNP[germline_NSFP_pred_SNP$SIFT4G_pred != "."]
germline_NSFP_pred_SNP <- germline_NSFP_pred_SNP[germline_NSFP_pred_SNP$Polyphen2_HVAR_pred != "."]
germline_NSFP_pred_SNP <- germline_NSFP_pred_SNP[germline_NSFP_pred_SNP$MetaSVM_pred != "."]
germline_NSFP_pred_SNP <- germline_NSFP_pred_SNP[germline_NSFP_pred_SNP$MetaLR_pred != "."]
germline_NSFP_pred_SNP <- germline_NSFP_pred_SNP[germline_NSFP_pred_SNP$MetaRNN_pred != "."]


colors <- c("somatic" = "red", "germline" = "blue")
SIFT_conv_rankscores_means <- as.data.frame(t(cbind(mean((somatic_NSFP$SIFT4G_converted_rankscore)), mean((germline_NSFP$SIFT4G_converted_rankscore)))))
SIFT_conv_rankscores_plot <- ggplot() + 
  geom_density(data = somatic_NSFP, aes(x = somatic_NSFP$SIFT4G_converted_rankscore, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP, aes(x = germline_NSFP$SIFT4G_converted_rankscore, y=..count../sum(..count..), color = "germline")) +
  geom_vline(data=SIFT_conv_rankscores_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.015)) + labs(x = "SIFT4G rankscore", y = "scaled density", title = "SIFT4G rankscore distribution", color = "Dataset") + scale_color_manual(values = colors)
SIFT_conv_rankscores_plot

SIFT_raw_means <- as.data.frame(t(cbind(mean((somatic_NSFP_raw$SIFT4G_score)), mean((germline_NSFP_raw$SIFT4G_score)))))
SIFT_raw_plot <- ggplot() + 
  geom_density(data = somatic_NSFP_raw, aes(x = somatic_NSFP_raw$SIFT4G_score, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP_raw, aes(x = germline_NSFP_raw$SIFT4G_score, y=..count../sum(..count..), color = "germline")) +
  geom_vline(data=SIFT_raw_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.035)) + labs(x = "SIFT4G raw score", y = "scaled density", title = "SIFT4G raw score distribution", color = "Dataset") + scale_color_manual(values = colors)
SIFT_raw_plot

SIFT_conv_rankscores_SNP_means <- as.data.frame(t(cbind(mean((somatic_NSFP_SNP_only$SIFT4G_converted_rankscore)), mean((germline_NSFP_SNP_only$SIFT4G_converted_rankscore)))))
SIFT_conv_rankscores_SNP_plot <- ggplot() + 
  geom_density(data = somatic_NSFP_SNP_only, aes(x = somatic_NSFP_SNP_only$SIFT4G_converted_rankscore, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP_SNP_only, aes(x = germline_NSFP_SNP_only$SIFT4G_converted_rankscore, y=..count../sum(..count..), color= "germline")) +
  geom_vline(data=SIFT_conv_rankscores_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.015)) + labs(x = "SIFT4G rankscore", y = "scaled density", title = "SIFT4G rankscore distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
SIFT_conv_rankscores_SNP_plot

SIFT_raw_SNP_means <- as.data.frame(t(cbind(mean((somatic_NSFP_raw_SNP$SIFT4G_score)), mean((germline_NSFP_raw$SIFT4G_score)))))
SIFT_raw_SNP_plot <- ggplot() + 
  geom_density(data = somatic_NSFP_raw_SNP, aes(x = somatic_NSFP_raw_SNP$SIFT4G_score, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP_raw_SNP, aes(x = germline_NSFP_raw_SNP$SIFT4G_score, y=..count../sum(..count..), color= "germline")) +
  geom_vline(data=SIFT_raw_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.035)) + labs(x = "SIFT4G raw score", y = "scaled density", title = "SIFT4G raw score distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
SIFT_raw_SNP_plot

Polyphen2_conv_rankscores_means <- as.data.frame(t(cbind(mean((somatic_NSFP$Polyphen2_HVAR_rankscore)), mean((germline_NSFP$Polyphen2_HVAR_rankscore)))))
Polyphen2_conv_rankscores_plot <- ggplot() + 
  geom_density(data = somatic_NSFP, aes(x = somatic_NSFP$Polyphen2_HVAR_rankscore, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP, aes(x = germline_NSFP$Polyphen2_HVAR_rankscore, y=..count../sum(..count..), color = "germline")) +
  geom_vline(data=Polyphen2_conv_rankscores_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.005)) + labs(x = "Polyphen2 rankscore", y = "scaled density", title = "Polyphen2 rankscore distribution", color = "Dataset") + scale_color_manual(values = colors)
Polyphen2_conv_rankscores_plot

Polyphen2_raw_means <- as.data.frame(t(cbind(mean((somatic_NSFP_raw$Polyphen2_HVAR_score)), mean((germline_NSFP_raw$Polyphen2_HVAR_score)))))
Polyphen2_raw_plot <- ggplot() + 
  geom_density(data = somatic_NSFP_raw, aes(x = somatic_NSFP_raw$Polyphen2_HVAR_score, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP_raw, aes(x = germline_NSFP_raw$Polyphen2_HVAR_score, y=..count../sum(..count..), color = "germline")) +
  geom_vline(data=Polyphen2_raw_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.015)) + labs(x = "Polyphen2 raw score", y = "scaled density", title = "Polyphen2 raw score distribution", color = "Dataset") + scale_color_manual(values = colors)
Polyphen2_raw_plot

Polyphen2_conv_rankscores_SNP_means <- as.data.frame(t(cbind(mean((somatic_NSFP_SNP_only$Polyphen2_HVAR_rankscore)), mean((germline_NSFP_SNP_only$Polyphen2_HVAR_rankscore)))))
Polyphen2_conv_rankscores_SNP_plot <- ggplot() + 
  geom_density(data = somatic_NSFP_SNP_only, aes(x = somatic_NSFP_SNP_only$Polyphen2_HVAR_rankscore, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP_SNP_only, aes(x = germline_NSFP_SNP_only$Polyphen2_HVAR_rankscore, y=..count../sum(..count..), color = "germline")) +
  geom_vline(data=Polyphen2_conv_rankscores_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.005)) + labs(x = "Polyphen2 rankscore", y = "scaled density", title = "Polyphen2 rankscore distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
Polyphen2_conv_rankscores_SNP_plot

Polyphen2_raw_SNP_means <- as.data.frame(t(cbind(mean((somatic_NSFP_raw_SNP$Polyphen2_HVAR_score)), mean((germline_NSFP_raw_SNP$Polyphen2_HVAR_score)))))
Polyphen2_raw_SNP_plot <- ggplot() + 
  geom_density(data = somatic_NSFP_raw_SNP, aes(x = somatic_NSFP_raw_SNP$Polyphen2_HVAR_score, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP_raw_SNP, aes(x = germline_NSFP_raw_SNP$Polyphen2_HVAR_score, y=..count../sum(..count..), color = "germline")) +
  geom_vline(data=Polyphen2_raw_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.015)) + labs(x = "Polyphen2 raw score", y = "scaled density", title = "Polyphen2 raw score distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
Polyphen2_raw_SNP_plot

MetaLR_rankscores_means <- as.data.frame(t(cbind(mean((somatic_NSFP$MetaLR_rankscore)), mean((germline_NSFP$MetaLR_rankscore)))))
MetaLR_rankscores_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP, aes(x = somatic_NSFP$MetaLR_rankscore, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP, aes(x = germline_NSFP$MetaLR_rankscore, y=..count../sum(..count..), color= "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaLR_rankscores_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.1)) + labs(x = "MetaLR rankscore", y = "normalized density", title = "MetaLR rankscore distribution", color = "Dataset") + scale_color_manual(values = colors)
MetaLR_rankscores_barplot

MetaLR_raw_means <- as.data.frame(t(cbind(mean((somatic_NSFP_raw$MetaLR_score)), mean((germline_NSFP_raw$MetaLR_score)))))
MetaLR_raw_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP_raw, aes(x = somatic_NSFP_raw$MetaLR_score, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP_raw, aes(x = germline_NSFP_raw$MetaLR_score, y=..count../sum(..count..), color= "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaLR_raw_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.25)) + labs(x = "MetaLR raw score", y = "normalized density", title = "MetaLR raw score distribution", color = "Dataset") + scale_color_manual(values = colors)
MetaLR_raw_barplot

MetaLR_rankscores_SNP_means <- as.data.frame(t(cbind(mean((somatic_NSFP_SNP_only$MetaLR_rankscore)), mean((germline_NSFP_SNP_only$MetaLR_rankscore)))))
MetaLR_rankscores_SNP_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP_SNP_only, aes(x = somatic_NSFP_SNP_only$MetaLR_rankscore, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP_SNP_only, aes(x = germline_NSFP_SNP_only$MetaLR_rankscore, y=..count../sum(..count..), color = "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaLR_rankscores_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.1)) + labs(x = "MetaLR rankscore", y = "normalized density", title = "MetaLR rankscore distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
MetaLR_rankscores_SNP_barplot

MetaLR_raw_SNP_means <- as.data.frame(t(cbind(mean((somatic_NSFP_raw_SNP$MetaLR_score)), mean((germline_NSFP_raw_SNP$MetaLR_score)))))
MetaLR_raw_SNP_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP_raw_SNP, aes(x = somatic_NSFP_raw_SNP$MetaLR_score, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP_raw_SNP, aes(x = germline_NSFP_raw_SNP$MetaLR_score, y=..count../sum(..count..), color = "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaLR_raw_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.25)) + labs(x = "MetaLR raw score", y = "normalized density", title = "MetaLR raw score distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
MetaLR_raw_SNP_barplot

MetaSVM_rankscores_means <- as.data.frame(t(cbind(mean((somatic_NSFP$MetaSVM_rankscore)), mean((germline_NSFP$MetaSVM_rankscore)))))
MetaSVM_rankscores_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP, aes(x = somatic_NSFP$MetaSVM_rankscore, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP, aes(x = germline_NSFP$MetaSVM_rankscore, y=..count../sum(..count..), color = "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaSVM_rankscores_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.1)) + labs(x = "MetaSVM rankscore", y = "normalized density", title = "MetaSVM rankscore distribution", color = "Dataset") + scale_color_manual(values = colors)
MetaSVM_rankscores_barplot

MetaSVM_raw_means <- as.data.frame(t(cbind(mean((somatic_NSFP_raw$MetaSVM_score)), mean((germline_NSFP_raw$MetaSVM_score)))))
MetaSVM_raw_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP_raw, aes(x = somatic_NSFP_raw$MetaSVM_score, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP_raw, aes(x = germline_NSFP_raw$MetaSVM_score, y=..count../sum(..count..), color = "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaSVM_raw_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(-1.5,1.5)) + ylim(c(0, 0.1)) + labs(x = "MetaSVM raw score", y = "normalized density", title = "MetaSVM raw score distribution", color = "Dataset") + scale_color_manual(values = colors)
MetaSVM_raw_barplot

MetaSVM_rankscores_SNP_means <- as.data.frame(t(cbind(mean((somatic_NSFP_SNP_only$MetaSVM_rankscore)), mean((germline_NSFP_SNP_only$MetaSVM_rankscore)))))
MetaSVM_rankscores_SNP_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP_SNP_only, aes(x = somatic_NSFP_SNP_only$MetaSVM_rankscore, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP_SNP_only, aes(x = germline_NSFP_SNP_only$MetaSVM_rankscore, y=..count../sum(..count..), color = "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaSVM_rankscores_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.1)) + labs(x = "MetaSVM rankscore", y = "normalized density", title = "MetaSVM rankscore distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
MetaSVM_rankscores_SNP_barplot

MetaSVM_raw_SNP_means <- as.data.frame(t(cbind(mean((somatic_NSFP_raw_SNP$MetaSVM_score)), mean((germline_NSFP_raw_SNP$MetaSVM_score)))))
MetaSVM_raw_SNP_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP_raw_SNP, aes(x = somatic_NSFP_raw_SNP$MetaSVM_score, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP_raw_SNP, aes(x = germline_NSFP_raw_SNP$MetaSVM_score, y=..count../sum(..count..), color = "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaSVM_raw_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(-1.5,1.5)) + ylim(c(0, 0.1)) + labs(x = "MetaSVM raw score", y = "normalized density", title = "MetaSVM raw score distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
MetaSVM_raw_SNP_barplot

MetaRNN_rankscores_means <- as.data.frame(t(cbind(mean((somatic_NSFP$MetaRNN_rankscore)), mean((germline_NSFP$MetaRNN_rankscore)))))
MetaRNN_rankscores_plot <- ggplot() + 
  geom_density(data = somatic_NSFP, aes(x = somatic_NSFP$MetaRNN_rankscore, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP, aes(x = germline_NSFP$MetaRNN_rankscore, y=..count../sum(..count..), color = "germline")) +
  geom_vline(data=MetaRNN_rankscores_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1))+ ylim(c(0, 0.012)) +  labs(x = "MetaRNN rankscore", y = "scaled density", title = "MetaRNN rankscore distribution", color = "Dataset") + scale_color_manual(values = colors)
MetaRNN_rankscores_plot

MetaRNN_raw_means <- as.data.frame(t(cbind(mean((somatic_NSFP_raw$MetaRNN_score)), mean((germline_NSFP_raw$MetaRNN_score)))))
MetaRNN_raw_plot <- ggplot() + 
  geom_density(data = somatic_NSFP_raw, aes(x = somatic_NSFP_raw$MetaRNN_score, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP_raw, aes(x = germline_NSFP_raw$MetaRNN_score, y=..count../sum(..count..), color = "germline")) +
  geom_vline(data=MetaRNN_raw_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1))+ ylim(c(0, 0.01)) +  labs(x = "MetaRNN raw score", y = "scaled density", title = "MetaRNN raw score distribution", color = "Dataset") + scale_color_manual(values = colors)
MetaRNN_raw_plot

MetaRNN_rankscores_SNP_means <- as.data.frame(t(cbind(mean((somatic_NSFP_SNP_only$MetaRNN_rankscore)), mean((germline_NSFP_SNP_only$MetaRNN_rankscore)))))
MetaRNN_rankscores_SNP_plot <- ggplot() + 
  geom_density(data = somatic_NSFP_SNP_only, aes(x = somatic_NSFP_SNP_only$MetaRNN_rankscore, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP_SNP_only, aes(x = germline_NSFP_SNP_only$MetaRNN_rankscore, y=..count../sum(..count..), color = "germline")) +
  geom_vline(data=MetaRNN_rankscores_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.012)) + labs(x = "MetaRNN rankscore", y = "scaled density", title = "MetaRNN rankscore distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
MetaRNN_rankscores_SNP_plot

MetaRNN_raw_SNP_means <- as.data.frame(t(cbind(mean((somatic_NSFP_raw_SNP$MetaRNN_score)), mean((germline_NSFP_raw_SNP$MetaRNN_score)))))
MetaRNN_raw_SNP_plot <- ggplot() + 
  geom_density(data = somatic_NSFP_raw_SNP, aes(x = somatic_NSFP_raw_SNP$MetaRNN_score, y=..count../sum(..count..), color= "somatic")) +
  geom_density(data = germline_NSFP_raw_SNP, aes(x = germline_NSFP_raw_SNP$MetaRNN_score, y=..count../sum(..count..), color = "germline")) +
  geom_vline(data=MetaRNN_raw_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.01)) + labs(x = "MetaRNN raw score", y = "scaled density", title = "MetaRNN raw score distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
MetaRNN_raw_SNP_plot

MetaRNN_rankscores_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP, aes(x = somatic_NSFP$MetaRNN_rankscore, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP, aes(x = germline_NSFP$MetaRNN_rankscore, y=..count../sum(..count..), color = "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaRNN_rankscores_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.1)) + labs(x = "MetaRNN rankscore", y = "normalized density", title = "MetaRNN rankscore distribution", color = "Dataset") + scale_color_manual(values = colors)
MetaRNN_rankscores_barplot

MetaRNN_raw_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP_raw, aes(x = somatic_NSFP_raw$MetaRNN_score, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP_raw, aes(x = germline_NSFP_raw$MetaRNN_score, y=..count../sum(..count..), color = "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaRNN_raw_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.2)) + labs(x = "MetaRNN raw score", y = "normalized density", title = "MetaRNN raw score distribution", color = "Dataset") + scale_color_manual(values = colors)
MetaRNN_raw_barplot

MetaRNN_rankscores_SNP_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP_SNP_only, aes(x = somatic_NSFP_SNP_only$MetaRNN_rankscore, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP_SNP_only, aes(x = germline_NSFP_SNP_only$MetaRNN_rankscore, y=..count../sum(..count..), color = "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaRNN_rankscores_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.1)) + labs(x = "MetaRNN rankscore", y = "normalized density", title = "MetaRNN rankscore distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
MetaRNN_rankscores_SNP_barplot

MetaRNN_raw_SNP_barplot <- ggplot() + 
  geom_histogram(data = somatic_NSFP_raw_SNP, aes(x = somatic_NSFP_raw_SNP$MetaRNN_score, y=..count../sum(..count..), color= "somatic"), binwidth=.05,  alpha=.5, position="identity") +
  geom_histogram(data = germline_NSFP_raw_SNP, aes(x = germline_NSFP_raw_SNP$MetaRNN_score, y=..count../sum(..count..), color = "germline"), binwidth=.05,  alpha=.5, position="identity") +
  geom_vline(data=MetaRNN_raw_SNP_means, aes(xintercept=V1, color = c("somatic", "germline")),
             linetype="dashed", size=0.75, show.legend = F) + theme_bw() + xlim(c(0,1)) + ylim(c(0, 0.2)) + labs(x = "MetaRNN raw score", y = "normalized density", title = "MetaRNN raw score distribution (SNPs only)", color = "Dataset") + scale_color_manual(values = colors)
MetaRNN_raw_SNP_barplot

algorithm <- c(rep("SIFT4G", 3), rep("Polyphen2", 3), rep("MetaSVM", 3), rep("MetaLR", 3), rep("MetaRNN", 3))
predictions <- rep(c("Damaging", "Possibly damaging", "Tolerated"), 5)
vals_somatic <- c(sum(somatic_NSFP_pred$SIFT4G_pred == "D"), sum(somatic_NSFP_pred$SIFT4G_pred == "P"), sum(somatic_NSFP_pred$SIFT4G_pred == "T"), 
            sum(somatic_NSFP_pred$Polyphen2_HVAR_pred == "D"), sum(somatic_NSFP_pred$Polyphen2_HVAR_pred == "P"), sum(somatic_NSFP_pred$Polyphen2_HVAR_pred == "B"),
            sum(somatic_NSFP_pred$MetaSVM_pred == "D"), sum(somatic_NSFP_pred$MetaSVM_pred == "P"), sum(somatic_NSFP_pred$MetaSVM_pred == "T"),
            sum(somatic_NSFP_pred$MetaLR_pred == "D"), sum(somatic_NSFP_pred$MetaLR_pred == "P"), sum(somatic_NSFP_pred$MetaLR_pred == "T"),
            sum(somatic_NSFP_pred$MetaRNN_pred == "D"), sum(somatic_NSFP_pred$MetaRNN_pred == "P"), sum(somatic_NSFP_pred$MetaRNN_pred == "T"))
data_somatic <- data.frame(algorithm, predictions, vals_somatic)
somatic_prediction_plot <- ggplot(data_somatic, aes(fill=predictions, y=vals_somatic, x=algorithm)) + 
  geom_bar(position="fill", stat="identity") + ylab("percentage") + theme_bw() + ggtitle("Score-based predictions in somatic dataset")
somatic_prediction_plot

vals_somatic_SNP <- c(sum(somatic_NSFP_pred_SNP$SIFT4G_pred == "D"), sum(somatic_NSFP_pred_SNP$SIFT4G_pred == "P"), sum(somatic_NSFP_pred_SNP$SIFT4G_pred == "T"), 
                  sum(somatic_NSFP_pred_SNP$Polyphen2_HVAR_pred == "D"), sum(somatic_NSFP_pred_SNP$Polyphen2_HVAR_pred == "P"), sum(somatic_NSFP_pred_SNP$Polyphen2_HVAR_pred == "B"),
                  sum(somatic_NSFP_pred_SNP$MetaSVM_pred == "D"), sum(somatic_NSFP_pred_SNP$MetaSVM_pred == "P"), sum(somatic_NSFP_pred_SNP$MetaSVM_pred == "T"),
                  sum(somatic_NSFP_pred_SNP$MetaLR_pred == "D"), sum(somatic_NSFP_pred_SNP$MetaLR_pred == "P"), sum(somatic_NSFP_pred_SNP$MetaLR_pred == "T"),
                  sum(somatic_NSFP_pred_SNP$MetaRNN_pred == "D"), sum(somatic_NSFP_pred_SNP$MetaRNN_pred == "P"), sum(somatic_NSFP_pred_SNP$MetaRNN_pred == "T"))
data_somatic_SNP <- data.frame(algorithm, predictions, vals_somatic_SNP)
somatic_prediction_SNP_plot <- ggplot(data_somatic_SNP, aes(fill=predictions, y=vals_somatic_SNP, x=algorithm)) + 
  geom_bar(position="fill", stat="identity") + ylab("percentage") + theme_bw() + ggtitle("Score-based predictions in somatic dataset (SNPs only)")
somatic_prediction_SNP_plot

vals_germline <- c(sum(germline_NSFP_pred$SIFT4G_pred == "D"), sum(germline_NSFP_pred$SIFT4G_pred == "P"), sum(germline_NSFP_pred$SIFT4G_pred == "T"), 
                  sum(germline_NSFP_pred$Polyphen2_HVAR_pred == "D"), sum(germline_NSFP_pred$Polyphen2_HVAR_pred == "P"), sum(germline_NSFP_pred$Polyphen2_HVAR_pred == "B"),
                  sum(germline_NSFP_pred$MetaSVM_pred == "D"), sum(germline_NSFP_pred$MetaSVM_pred == "P"), sum(germline_NSFP_pred$MetaSVM_pred == "T"),
                  sum(germline_NSFP_pred$MetaLR_pred == "D"), sum(germline_NSFP_pred$MetaLR_pred == "P"), sum(germline_NSFP_pred$MetaLR_pred == "T"),
                  sum(germline_NSFP_pred$MetaRNN_pred == "D"), sum(germline_NSFP_pred$MetaRNN_pred == "P"), sum(germline_NSFP_pred$MetaRNN_pred == "T"))
data_germline <- data.frame(algorithm, predictions, vals_germline)
germline_prediction_plot <- ggplot(data_germline, aes(fill=predictions, y=vals_germline, x=algorithm)) + 
  geom_bar(position="fill", stat="identity") + ylab("percentage") + theme_bw() + ggtitle("Score-based predictions in germline dataset")
germline_prediction_plot

vals_germline_SNP <- c(sum(germline_NSFP_pred_SNP$SIFT4G_pred == "D"), sum(germline_NSFP_pred_SNP$SIFT4G_pred == "P"), sum(germline_NSFP_pred_SNP$SIFT4G_pred == "T"), 
                  sum(germline_NSFP_pred_SNP$Polyphen2_HVAR_pred == "D"), sum(germline_NSFP_pred_SNP$Polyphen2_HVAR_pred == "P"), sum(germline_NSFP_pred_SNP$Polyphen2_HVAR_pred == "B"),
                  sum(germline_NSFP_pred_SNP$MetaSVM_pred == "D"), sum(germline_NSFP_pred_SNP$MetaSVM_pred == "P"), sum(germline_NSFP_pred_SNP$MetaSVM_pred == "T"),
                  sum(germline_NSFP_pred_SNP$MetaLR_pred == "D"), sum(germline_NSFP_pred_SNP$MetaLR_pred == "P"), sum(germline_NSFP_pred_SNP$MetaLR_pred == "T"),
                  sum(germline_NSFP_pred_SNP$MetaRNN_pred == "D"), sum(germline_NSFP_pred_SNP$MetaRNN_pred == "P"), sum(germline_NSFP_pred_SNP$MetaRNN_pred == "T"))
data_germline_SNP <- data.frame(algorithm, predictions, vals_germline_SNP)
germline_prediction_SNP_plot <- ggplot(data_germline_SNP, aes(fill=predictions, y=vals_germline_SNP, x=algorithm)) + 
  geom_bar(position="fill", stat="identity") + ylab("percentage") + theme_bw() + ggtitle("Score-based predictions in germline dataset (SNPs only)")
germline_prediction_SNP_plot


dataset <- c(rep("somatic", 2), rep("somatic (SNPs only)", 2), rep("germline", 2), rep("germline (SNPs only)", 2))
predictions <- rep(c("Damaging", "Tolerated"), 4)
vals_sift <- c(sum(somatic_NSFP_pred$SIFT4G_pred == "D"), sum(somatic_NSFP_pred$SIFT4G_pred == "T"), 
               sum(somatic_NSFP_pred_SNP$SIFT4G_pred == "D"), sum(somatic_NSFP_pred_SNP$SIFT4G_pred == "T"),
               sum(germline_NSFP_pred$SIFT4G_pred == "D"), sum(germline_NSFP_pred$SIFT4G_pred == "T"),
               sum(germline_NSFP_pred_SNP$SIFT4G_pred == "D"), sum(germline_NSFP_pred_SNP$SIFT4G_pred == "T"))
data_sift <- data.frame(dataset, predictions, vals_sift)
sift_prediction_plot <- ggplot(data_sift, aes(fill=predictions, y=vals_sift, x=dataset)) + 
  geom_bar(position="fill", stat="identity") + ylab("percentage") + theme_bw() + ggtitle("Score-based predictions of SIFT4G") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
sift_prediction_plot

dataset <- c(rep("somatic", 3), rep("somatic (SNPs only)", 3), rep("germline", 3), rep("germline (SNPs only)", 3))
predictions <- rep(c("Damaging", "Possibly damaging", "Benign"), 4)
vals_polyphen <- c(sum(somatic_NSFP_pred$Polyphen2_HVAR_pred == "D"), sum(somatic_NSFP_pred$Polyphen2_HVAR_pred == "P"), sum(somatic_NSFP_pred$Polyphen2_HVAR_pred == "B"), 
               sum(somatic_NSFP_pred_SNP$Polyphen2_HVAR_pred == "D"), sum(somatic_NSFP_pred_SNP$Polyphen2_HVAR_pred == "P"), sum(somatic_NSFP_pred_SNP$Polyphen2_HVAR_pred == "B"),
               sum(germline_NSFP_pred$Polyphen2_HVAR_pred == "D"), sum(germline_NSFP_pred$Polyphen2_HVAR_pred == "P"), sum(germline_NSFP_pred$Polyphen2_HVAR_pred == "B"),
               sum(germline_NSFP_pred_SNP$Polyphen2_HVAR_pred == "D"), sum(germline_NSFP_pred_SNP$Polyphen2_HVAR_pred == "P"), sum(germline_NSFP_pred_SNP$Polyphen2_HVAR_pred == "B"))
data_polyphen <- data.frame(dataset, predictions, vals_polyphen)
polyphen_prediction_plot <- ggplot(data_polyphen, aes(fill=predictions, y=vals_polyphen, x=dataset)) + 
  geom_bar(position="fill", stat="identity") + ylab("percentage") + theme_bw() + ggtitle("Score-based predictions of Polyphen2") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
polyphen_prediction_plot

dataset <- c(rep("somatic", 2), rep("somatic (SNPs only)", 2), rep("germline", 2), rep("germline (SNPs only)", 2))
predictions <- rep(c("Damaging", "Tolerated"), 4)
vals_svm <- c(sum(somatic_NSFP_pred$MetaSVM_pred == "D"), sum(somatic_NSFP_pred$MetaSVM_pred == "T"), 
               sum(somatic_NSFP_pred_SNP$MetaSVM_pred == "D"), sum(somatic_NSFP_pred_SNP$MetaSVM_pred == "T"),
               sum(germline_NSFP_pred$MetaSVM_pred == "D"), sum(germline_NSFP_pred$MetaSVM_pred == "T"),
               sum(germline_NSFP_pred_SNP$MetaSVM_pred == "D"), sum(germline_NSFP_pred_SNP$MetaSVM_pred == "T"))
data_svm <- data.frame(dataset, predictions, vals_svm)
svm_prediction_plot <- ggplot(data_svm, aes(fill=predictions, y=vals_svm, x=dataset)) + 
  geom_bar(position="fill", stat="identity") + ylab("percentage") + theme_bw() + ggtitle("Score-based predictions of MetaSVM") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
svm_prediction_plot

dataset <- c(rep("somatic", 2), rep("somatic (SNPs only)", 2), rep("germline", 2), rep("germline (SNPs only)", 2))
predictions <- rep(c("Damaging", "Tolerated"), 4)
vals_rnn <- c(sum(somatic_NSFP_pred$SIFT4G_pred == "D"), sum(somatic_NSFP_pred$SIFT4G_pred == "T"), 
               sum(somatic_NSFP_pred_SNP$SIFT4G_pred == "D"), sum(somatic_NSFP_pred_SNP$SIFT4G_pred == "T"),
               sum(germline_NSFP_pred$SIFT4G_pred == "D"), sum(germline_NSFP_pred$SIFT4G_pred == "T"),
               sum(germline_NSFP_pred_SNP$SIFT4G_pred == "D"), sum(germline_NSFP_pred_SNP$SIFT4G_pred == "T"))
data_rnn <- data.frame(dataset, predictions, vals_rnn)
rnn_prediction_plot <- ggplot(data_rnn, aes(fill=predictions, y=vals_rnn, x=dataset)) + 
  geom_bar(position="fill", stat="identity") + ylab("percentage") + theme_bw() + ggtitle("Score-based predictions of MetaRNN") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
rnn_prediction_plot

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/NSFP_multiplot_new2.pdf', 
    width = 8,
    height = 10)
counter_somatic <- nrow(somatic_NSFP)
counter_germline <- nrow(germline_NSFP)
data_multi <- data.frame(algorithm = c(rep("SIFT4G", counter_somatic), rep("Polyphen2", counter_somatic), rep("MetaSVM", counter_somatic), rep("MetaLR", counter_somatic), rep("MetaRNN", counter_somatic), rep("SIFT4G", counter_germline), rep("Polyphen2", counter_germline), rep("MetaSVM", counter_germline), rep("MetaLR", counter_germline), rep("MetaRNN", counter_germline)),
                                 dataset= c(rep("somatic", (5*counter_somatic)), rep("germline", (5*counter_germline))), scores = c(somatic_NSFP$SIFT4G_converted_rankscore, 
  somatic_NSFP$Polyphen2_HVAR_rankscore,somatic_NSFP$MetaSVM_rankscore,somatic_NSFP$MetaLR_rankscore,somatic_NSFP$MetaRNN_rankscore, germline_NSFP$SIFT4G_converted_rankscore, germline_NSFP$Polyphen2_HVAR_rankscore, germline_NSFP$MetaSVM_rankscore, germline_NSFP$MetaLR_rankscore, germline_NSFP$MetaRNN_rankscore))
multiplot <- ggplot(data_multi, aes(x = scores, y = algorithm, color = dataset, fill = dataset)) +
  geom_density_ridges() +
  theme_ridges() + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), name = "rank score") +
  scale_fill_manual(values = c("#D55E0050", "#0072B250")) +
  scale_color_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(
    override.aes = list(
      fill = c("#D55E00A0", "#0072B2A0"),
      color = NA, point_color = NA)
  )
  ) +
  theme(text = element_text(size=10))+
  ggtitle("Rank score distributions in somatic and germline datasets") +
  theme_ridges(center = TRUE)
multiplot

dev.off()

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/NSFP_multiplot_SNP_new.pdf', 
    width = 8,
    height = 10)
counter_somatic_SNP <- nrow(somatic_NSFP_SNP_only)
counter_germline_SNP <- nrow(germline_NSFP_SNP_only)
data_multi_SNP <- data.frame(algorithm = c(rep("SIFT4G", counter_somatic_SNP), rep("Polyphen2", counter_somatic_SNP), rep("MetaSVM", counter_somatic_SNP), rep("MetaLR", counter_somatic_SNP), rep("MetaRNN", counter_somatic_SNP), rep("SIFT4G", counter_germline_SNP), rep("Polyphen2", counter_germline_SNP), rep("MetaSVM", counter_germline_SNP), rep("MetaLR", counter_germline_SNP), rep("MetaRNN", counter_germline_SNP)),
  dataset= c(rep("somatic", (5*counter_somatic_SNP)), rep("germline", (5*counter_germline_SNP))), 
  scores = c(somatic_NSFP_SNP_only$SIFT4G_converted_rankscore, somatic_NSFP_SNP_only$Polyphen2_HVAR_rankscore,somatic_NSFP_SNP_only$MetaSVM_rankscore,somatic_NSFP_SNP_only$MetaLR_rankscore,somatic_NSFP_SNP_only$MetaRNN_rankscore, germline_NSFP_SNP_only$SIFT4G_converted_rankscore, germline_NSFP_SNP_only$Polyphen2_HVAR_rankscore, germline_NSFP_SNP_only$MetaSVM_rankscore, germline_NSFP_SNP_only$MetaLR_rankscore, germline_NSFP_SNP_only$MetaRNN_rankscore))
multiplot_SNP <- ggplot(data_multi_SNP, aes(x = scores, y = algorithm, color = dataset, fill = dataset)) +
  geom_density_ridges() +
  theme_ridges() + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), name = "rank score") +
  scale_fill_manual(values = c("#D55E0050", "#0072B250")) +
  scale_color_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(
    override.aes = list(
      fill = c("#D55E00A0", "#0072B2A0"),
      color = NA, point_color = NA)
  )
  ) +
  ggtitle("Rank score distributions in somatic and germline datasets (SNPs only)") +
  theme_ridges(center = TRUE)
multiplot_SNP
dev.off()

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/NSFP_multiplot_raw_new.pdf', 
    width = 8,
    height = 10)
counter_somatic_raw <- nrow(somatic_NSFP_raw)
counter_germline_raw <- nrow(germline_NSFP_raw)
data_multi_raw <- data.frame(algorithm = c(rep("SIFT4G", counter_somatic_raw), rep("Polyphen2", counter_somatic_raw), rep("MetaSVM", counter_somatic_raw), rep("MetaLR", counter_somatic_raw), rep("MetaRNN", counter_somatic_raw), rep("SIFT4G", counter_germline_raw), rep("Polyphen2", counter_germline_raw), rep("MetaSVM", counter_germline_raw), rep("MetaLR", counter_germline_raw), rep("MetaRNN", counter_germline_raw)),
                         dataset= c(rep("somatic", (5*counter_somatic_raw)), rep("germline", (5*counter_germline_raw))), scores = c(somatic_NSFP_raw$SIFT4G_score, 
                                                                                                                                    somatic_NSFP_raw$Polyphen2_HVAR_score,somatic_NSFP_raw$MetaSVM_score,somatic_NSFP_raw$MetaLR_score,somatic_NSFP_raw$MetaRNN_score, germline_NSFP_raw$SIFT4G_score, germline_NSFP_raw$Polyphen2_HVAR_score, germline_NSFP_raw$MetaSVM_score, germline_NSFP_raw$MetaLR_score, germline_NSFP_raw$MetaRNN_score))
multiplot_raw <- ggplot(data_multi_raw, aes(x = scores, y = algorithm, color = dataset, fill = dataset)) +
  geom_density_ridges() +
  theme_ridges() + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), name = "raw score") +
  scale_fill_manual(values = c("#D55E0050", "#0072B250")) +
  scale_color_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(
    override.aes = list(
      fill = c("#D55E00A0", "#0072B2A0"),
      color = NA, point_color = NA)
  )
  ) +
  ggtitle("Raw score distributions in \n somatic and germline datasets") +
  theme_ridges(center = TRUE)
multiplot_raw
dev.off()

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/NSFP_multiplot_SNP_raw_new.pdf', 
    width = 8,
    height = 10)
counter_somatic_SNP_raw <- nrow(somatic_NSFP_raw_SNP)
counter_germline_SNP_raw <- nrow(germline_NSFP_raw_SNP)
data_multi_SNP_raw <- data.frame(algorithm = c(rep("SIFT4G", counter_somatic_SNP_raw), rep("Polyphen2", counter_somatic_SNP_raw), rep("MetaSVM", counter_somatic_SNP_raw), rep("MetaLR", counter_somatic_SNP_raw), rep("MetaRNN", counter_somatic_SNP_raw), rep("SIFT4G", counter_germline_SNP_raw), rep("Polyphen2", counter_germline_SNP_raw), rep("MetaSVM", counter_germline_SNP_raw), rep("MetaLR", counter_germline_SNP_raw), rep("MetaRNN", counter_germline_SNP_raw)),
                             dataset= c(rep("somatic", (5*counter_somatic_SNP_raw)), rep("germline", (5*counter_germline_SNP_raw))), 
                             scores = c(somatic_NSFP_raw_SNP$SIFT4G_score, somatic_NSFP_raw_SNP$Polyphen2_HVAR_score,somatic_NSFP_raw_SNP$MetaSVM_score,somatic_NSFP_raw_SNP$MetaLR_score,somatic_NSFP_raw_SNP$MetaRNN_score, germline_NSFP_raw_SNP$SIFT4G_score, germline_NSFP_raw_SNP$Polyphen2_HVAR_score, germline_NSFP_raw_SNP$MetaSVM_score, germline_NSFP_raw_SNP$MetaLR_score, germline_NSFP_raw_SNP$MetaRNN_score))
multiplot_SNP_raw <- ggplot(data_multi_SNP_raw, aes(x = scores, y = algorithm, color = dataset, fill = dataset)) +
  geom_density_ridges() +
  theme_ridges() + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), name = "raw score") +
  scale_fill_manual(values = c("#D55E0050", "#0072B250")) +
  scale_color_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(
    override.aes = list(
      fill = c("#D55E00A0", "#0072B2A0"),
      color = NA, point_color = NA)
  )
  ) +
  ggtitle("Raw score distributions in somatic \n and germline datasets (SNPs only)") +
  theme_ridges(center = TRUE)
multiplot_SNP_raw
dev.off()

###
###Testing for statistical significance between datasets###
###

pvals_rankscores = t.test(somatic_NSFP, germline_NSFP)
pvals_rankscores = t.test(somatic_NSFP$MetaRNN_rankscore, germline_NSFP$MetaRNN_rankscore)
pvals_rankscores_bonf = p.adjust(pvals_rankscores)

###
### Finding best fit of data distribution ###
###
scores = c(somatic_NSFP_ext$MetaRNN_rankscore, germline_NSFP_ext$MetaRNN_rankscore)
#dataset = c(rep("somatic", 5000), rep("germline", 5000))
dataset = c(rep("somatic", nrow(somatic_NSFP_ext)), rep("germline", nrow(germline_NSFP_ext)))
X = data.frame(scores, dataset)
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/Statistical Significance Analysis/violin_boxplot_MetaRNN_rank.pdf', 
    width = 6,
    height = 6)
ggplot(data = X, aes(x=dataset, y=scores, fill=dataset)) +
  geom_violin(width=0.4) +
  geom_boxplot(width=0.1, color="black", alpha=0.5) +
#  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  theme(legend.position="none")+
  ylab("MetaRNN rank score")

dev.off()

scores_raw = c(somatic_NSFP_ext$MetaRNN_score, germline_NSFP_ext$MetaRNN_score)
#dataset = c(rep("somatic", 5000), rep("germline", 5000))
dataset = c(rep("somatic", nrow(somatic_NSFP_ext)), rep("germline", nrow(germline_NSFP_ext)))
X_raw = data.frame(scores_raw, dataset)
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/Statistical Significance Analysis/violin_boxplot_MetaRNN_raw.pdf', 
    width = 6,
    height = 6)
ggplot(data = X_raw, aes(x=dataset, y=scores_raw, fill=dataset)) +
  geom_violin(width=0.4) +
  geom_boxplot(width=0.1, color="black", alpha=0.5) +
  #  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  theme(legend.position="none")+
  ylab("MetaRNN raw score")
dev.off()


with(X_sub, shapiro.test(scores[dataset == "somatic"]))
with(X_sub, shapiro.test(scores[dataset == "germline"]))
ks_ranscores <- ks.test(X$scores[X$dataset == "somatic"], X$scores[X$dataset == "germline"], alternative = "two.sided")

library(fitdistrplus)
library(logspline)
descdist(somatic_NSFP$MetaRNN_rankscore, discrete = FALSE)
descdist(germline_NSFP$MetaRNN_rankscore, discrete = FALSE)
descdist(X$scores, discrete = F)
fit.norm_somatic <- fitdist(somatic_NSFP$MetaRNN_rankscore, "norm")
fit.norm_germline <- fitdist(germline_NSFP$MetaRNN_rankscore, "norm")
t.test()
fit.unif <- fitdist(X$scores, "unif")
fit.unif_sub <- fitdist(X_sub$scores, "unif")
fit.norm <- fitdist(X$scores, "norm")
plot(fit.norm)
plot(fit.unif_sub)

res.ftest <- var.test(scores ~ dataset, data = X)
res.ftest

t_test = t.test(scores~dataset, data = X, var.equal = T, alternative = "two.sided")



scores_raw = c(somatic_NSFP_raw$MetaRNN_score, germline_NSFP_raw$MetaRNN_score)
#dataset = c(rep("somatic", 5000), rep("germline", 5000))
dataset_raw = c(rep("somatic", counter_somatic_raw), rep("germline", counter_germline_raw))
X_raw = data.frame(scores_raw, dataset_raw)
X_raw = X_raw[sample(nrow(X_raw)),]
X_raw_sub <- X_raw[1:5000,]

with(X_raw_sub, shapiro.test(scores_raw[dataset_raw == "somatic"]))
with(X_raw_sub, shapiro.test(scores_raw[dataset_raw == "germline"]))



ggboxplot(X_raw, x = "dataset_raw", y = "scores_raw", 
          color = "dataset_raw", palette = c("#00AFBB", "#E7B800"),
          ylab = "MetaRNN raw score", xlab = "dataset")

ks_raw <- ks.test(X_raw$scores_raw[X_raw$dataset_raw == "somatic"], X_raw$scores_raw[X_raw$dataset_raw == "germline"], alternative = "two.sided")
descdist(X_raw$scores_raw, discrete = F)
fit.unif_raw <- fitdist(X_raw$scores_raw, "unif")
#fit.norm_raw <- fitdist(X_raw$scores_raw, "norm")
#plot(fit.norm_raw)
plot(fit.unif_raw)
res.ftest_raw <- var.test(scores_raw ~ dataset_raw, data = X_raw)
res.ftest_raw
wilcox_test_raw = wilcox.test(scores_raw ~ dataset_raw, data = X_raw, alternative = "two.sided")
wilcox_test_raw


###
### Actual test ###
###


### Shared variants included
counter_somatic_ext <- nrow(somatic_NSFP_ext)
counter_germline_ext <- nrow(germline_NSFP_ext)
scores = c(somatic_NSFP_ext$MetaRNN_rankscore, germline_NSFP_ext$MetaRNN_rankscore)
#dataset = c(rep("somatic", 5000), rep("germline", 5000))
dataset = c(rep("somatic", counter_somatic_ext), rep("germline", counter_germline_ext))

### Shared variants excluded
counter_somatic_ext <- nrow(somatic_NSFP_extended)
counter_germline_ext <- nrow(germline_NSFP_extended)
scores = c(somatic_NSFP_extended$MetaRNN_rankscore, germline_NSFP_extended$MetaRNN_rankscore)
#dataset = c(rep("somatic", 5000), rep("germline", 5000))
dataset = c(rep("somatic", counter_somatic_ext), rep("germline", counter_germline_ext))

X= data.frame(scores, dataset)
#X = X[sample(nrow(X)),]
#X_sub <- X[1:5000,]

t_test_ori = t.test(scores~dataset, data = X, var.equal = T, alternative = "two.sided")
ks_test_ori <- ks.test(X$scores[X$dataset == "somatic"], X$scores[X$dataset == "germline"], alternative = "two.sided")
wilcox_test_ori = wilcox.test(scores ~ dataset, data = X, alternative = "two.sided")


### Bootstrapping with 1000 tests ###


t_tests_pvals = c()
t_tests_nullvals = c()
ks_tests_pvals = c()
ks_tests_nullvals = c()
wilcox_tests_pvals = c()
wilcox_tests_nullvals = c()
wilcox_tests_stats = c()
for (i in 1:1000) {
  X$dataset <- sample(X$dataset)
  print(i)
  t_test_temp = t.test(scores~dataset, data = X, var.equal = T, alternative = "two.sided")
  t_tests_pvals = c(t_tests_pvals, t_test_temp$p.value)
  temp_val1 = t_test_temp$estimate[1]
  temp_val2 = t_test_temp$estimate[2]  
  t_tests_nullvals = c(t_tests_nullvals, abs(temp_val1 - temp_val2))
  ks_test_temp = ks.test(X$scores[X$dataset == "somatic"], X$scores[X$dataset == "germline"], alternative = "two.sided")
  ks_tests_pvals = c(ks_tests_pvals, ks_test_temp$p.value)
  ks_tests_nullvals = c(ks_tests_nullvals, ks_test_temp$statistic)
  wilcox_temp = wilcox.test(scores ~ dataset, data = X, alternative = "two.sided")
  wilcox_tests_pvals = c(wilcox_tests_pvals, wilcox_temp$p.value)
  wilcox_tests_nullvals = c(wilcox_tests_nullvals, wilcox_temp$null.value)
  wilcox_tests_stats = c(wilcox_tests_stats, wilcox_temp$statistic)
}

###
### Plots from bootstrapping ###
###

### T-test ###
t_tests_pvals = data.frame(t_tests_pvals)
t_tests_pvals_plot <- ggplot() + 
  geom_histogram(data = t_tests_pvals, aes(x= t_tests_pvals)) +
  geom_vline(xintercept= t_test_ori$p.value, size=0.75, show.legend = F, color = "red") + 
  geom_text(aes(x=t_test_ori$p.value, label="\np-value of t-test with true labels", y=15), colour="red", angle=90, text=element_text(size=11)) +
  theme_bw()  +  labs(x = "p-values", y = "count", title = "T-test p-values from random dataset resampling (1000 runs)" )
t_tests_pvals_plot

t_tests_nullvals = data.frame(t_tests_nullvals)
t_test_ori_diff = abs(t_test_ori$estimate[1] - t_test_ori$estimate[2])
t_tests_nullvals_plot <- ggplot() + 
  geom_histogram(data = t_tests_nullvals, aes(x= t_tests_nullvals), bins = 10) +
  geom_vline(xintercept= t_test_ori_diff, size=0.75, show.legend = F, color = "red") + 
  geom_text(aes(x=t_test_ori_diff, label="\nt-test estimate with true labels", y=300), colour="red", angle=90, text=element_text(size=11)) +
  theme_bw()  +  labs(x = "differences in means", y = "count", title = "T-test means differences from random dataset resampling (1000 runs)" )
t_tests_nullvals_plot

### KS-test ###
ks_tests_pvals = data.frame(ks_tests_pvals)
ks_tests_pvals_plot <- ggplot() + 
  geom_histogram(data = ks_tests_pvals, aes(x= ks_tests_pvals)) +
  geom_vline(xintercept= ks_test_ori$p.value, size=0.75, show.legend = F, color = "red") + 
  geom_text(aes(x=ks_test_ori$p.value, label="\np-value of KS-test with true labels", y=15), colour="red", angle=90, text=element_text(size=11)) +
  theme_bw()  +  labs(x = "p-values", y = "count", title = "KS-test p-values from random dataset resampling (1000 runs)" )
ks_tests_pvals_plot

ks_tests_nullvals = data.frame(ks_tests_nullvals)
ks_tests_nullvals_plot <- ggplot() + 
  geom_histogram(data = ks_tests_nullvals, aes(x= ks_tests_nullvals)) +
  geom_vline(xintercept= ks_test_ori$statistic, size=0.75, show.legend = F, color = "red") + 
  geom_text(aes(x=ks_test_ori$statistic, label="\nKS-test statistic with true labels", y=300), colour="red", angle=90, text=element_text(size=11)) +
  theme_bw()  +  labs(x = "KS-statistic", y = "count", title = "KS-test statistic from dataset random resampling (1000 runs)" )
ks_tests_nullvals_plot

### Wilcox test ###

wilcox_tests_pvals = data.frame(wilcox_tests_pvals)
wilcox_tests_pvals_plot <- ggplot() + 
  geom_histogram(data = wilcox_tests_pvals, aes(x= wilcox_tests_pvals)) +
  geom_vline(xintercept= wilcox_test_ori$p.value, size=0.75, show.legend = F, color = "red") + 
  geom_text(aes(x=wilcox_test_ori$p.value, label="\np-value of Wilcox-test with true labels", y=15), colour="red", angle=90, text=element_text(size=11)) +
  theme_bw()  +  labs(x = "p-values", y = "count", title = "Wilcox-test p-values from random dataset resampling (1000 runs)" )
wilcox_tests_pvals_plot

wilcox_tests_nullvals = data.frame(wilcox_tests_nullvals)
wilcox_tests_nullvals_plot <- ggplot() + 
  geom_histogram(data = wilcox_tests_nullvals, aes(x= wilcox_tests_nullvals)) +
  geom_vline(xintercept= wilcox_test_ori$null.value, size=0.75, show.legend = F, color = "red") + 
  geom_text(aes(x=wilcox_test_ori$null.value, label="\nWilcox-test mu-location with true labels", y=300), colour="red", angle=90, text=element_text(size=11)) +
  theme_bw()  +  labs(x = "KS-null-value", y = "count", title = "KS-test null-value from dataset random resampling (1000 runs)" )
wilcox_tests_nullvals_plot

wilcox_tests_stats = data.frame(wilcox_tests_stats)
wilcox_tests_stat_plot <- ggplot() + 
  geom_histogram(data = wilcox_tests_stats, aes(x= wilcox_tests_stats)) +
  geom_vline(xintercept= wilcox_test_ori$statistic, size=0.75, show.legend = F, color = "red") + 
  geom_text(aes(x=wilcox_test_ori$statistic, label="\nWilcox-test statistic with true labels", y=300), colour="red", angle=90, text=element_text(size=11)) +
  theme_bw()  +  labs(x = "Wilcox-statistic", y = "count", title = "Wilcox-test statistic from dataset random resampling (1000 runs)" )
wilcox_tests_stat_plot




###
### Dataset overlap analysis
###


germline_NSFP_pos <- fread(
  "/Users/mcpftw/dbNSFP4.2a/germline_NSFP_pos.txt",
  select = c(
    "#chr",
    "pos(1-based)",
    "ref",
    "alt"
  )
)

somatic_NSFP_pos <- fread(
  "/Users/mcpftw/dbNSFP4.2a/somatic_NSFP_pos.txt",
  select = c(
    "#chr",
    "pos(1-based)",
    "ref",
    "alt"
  )
)

somatic_NSFP_ext <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/somatic_NSFP_pos.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "MetaRNN_score",
      "MetaRNN_rankscore",
      "MutPred_score",
      "MutPred_rankscore",
      "MutPred_Top5features",
      "GO_biological_process",
      "GO_cellular_component",
      "GO_molecular_function"
      
    )
  )


somatic_NSFP_ext <- somatic_NSFP_ext[somatic_NSFP_ext$MetaRNN_rankscore != "."]
somatic_NSFP_ext$MetaRNN_rankscore <- as.numeric(somatic_NSFP_ext$MetaRNN_rankscore)
somatic_NSFP_ext$MetaRNN_score <- sapply(str_split(somatic_NSFP_ext$MetaRNN_score, ";"), tail, 1)
somatic_NSFP_ext <- somatic_NSFP_ext[somatic_NSFP_ext$MetaRNN_score != "."]
somatic_NSFP_ext$MetaRNN_score <- as.numeric(somatic_NSFP_ext$MetaRNN_score)
somatic_NSFP_ext <- somatic_NSFP_ext[somatic_NSFP_ext$MutPred_score != "."]
somatic_NSFP_ext$MutPred_score <- as.numeric(somatic_NSFP_ext$MutPred_score)
somatic_NSFP_ext <- somatic_NSFP_ext[somatic_NSFP_ext$MutPred_rankscore != "."]
somatic_NSFP_ext$MutPred_rankscore <- as.numeric(somatic_NSFP_ext$MutPred_rankscore)
somatic_NSFP_ext = somatic_NSFP_ext %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_score`))


germline_NSFP_ext <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/germline_NSFP_pos.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "MetaRNN_score",
      "MetaRNN_rankscore",
      "MutPred_score",
      "MutPred_rankscore",
      "MutPred_Top5features",
      "GO_biological_process",
      "GO_cellular_component",
      "GO_molecular_function"
    )
  )

germline_NSFP_ext <- germline_NSFP_ext[germline_NSFP_ext$MetaRNN_rankscore != "."]
germline_NSFP_ext$MetaRNN_rankscore <- as.numeric(germline_NSFP_ext$MetaRNN_rankscore)
germline_NSFP_ext$MetaRNN_score <- sapply(str_split(germline_NSFP_ext$MetaRNN_score, ";"), tail, 1)
germline_NSFP_ext <- germline_NSFP_ext[germline_NSFP_ext$MetaRNN_score != "."]
germline_NSFP_ext$MetaRNN_score <- as.numeric(germline_NSFP_ext$MetaRNN_score)
germline_NSFP_ext <- germline_NSFP_ext[germline_NSFP_ext$MutPred_score != "."]
germline_NSFP_ext$MutPred_score <- as.numeric(germline_NSFP_ext$MutPred_score)
germline_NSFP_ext <- germline_NSFP_ext[germline_NSFP_ext$MutPred_rankscore != "."]
germline_NSFP_ext$MutPred_rankscore <- as.numeric(germline_NSFP_ext$MutPred_rankscore)
germline_NSFP_ext = germline_NSFP_ext %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_score`))

### Extract shared mutations and remove them from somatic and germline sets
germline_NSFP_extended <- anti_join(germline_NSFP_ext, somatic_NSFP_ext, by = c("#chr", "pos(1-based)", "ref", "alt"))
somatic_NSFP_extended <- anti_join(somatic_NSFP_ext, germline_NSFP_ext, by = c("#chr", "pos(1-based)", "ref", "alt"))
shared_extended <- inner_join(somatic_NSFP_ext, germline_NSFP_ext, by = c("#chr", "pos(1-based)", "ref", "alt"))
names(shared_extended)[names(shared_extended) == 'MutPred_score.x'] <- 'MutPred_score'
names(shared_extended)[names(shared_extended) == 'MutPred_Top5features.x'] <- 'MutPred_Top5features'
rm(germline_NSFP_ext, somatic_NSFP_ext)

### Shared mutations Venn Diagram 
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/venn_shared_mutations_new.pdf', 
    width = 8,
    height = 4)
grid.newpage()
draw.pairwise.venn(area1 = nrow(somatic_NSFP_ext),                        
                   area2 = nrow(germline_NSFP_ext),
                   cross.area = nrow(shared_extended),
                   fill = c("red", "blue"),
                   category = c("somatic", "germline"))
dev.off()

### Shared mutations additional statistics

### Chromosome distribution
shared_chromosome_counts <- as.data.frame(table(shared_extended$`#chr`))
library(gtools)
library(jam)
chrOrder<-c(1:22, "X","Y","M")
#chrOrder<-c("M", "Y", "X", 22:1)
shared_chromosome_counts$Var1<-factor(shared_chromosome_counts$Var1, levels=chrOrder)
shared_chromosome_counts$Var1
shared_chromosome_counts <- shared_chromosome_counts[order(shared_chromosome_counts$Var1),]

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/shared_mutations_chr.pdf', 
    width = 8,
    height = 4)
shared_chromosome_plot <-
  ggplot(shared_chromosome_counts, aes(
    x = Var1,
    y = Freq,
    fill = Var1,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
#  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("variants") +
  xlab("chromosome") +
  geom_text(size = 3, angle = 90, position = position_stack(vjust = 0.8)) 
#top_func_somatic_plot
plot(shared_chromosome_plot)
dev.off()

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/shared_mutations_MetaRNN_distrib.pdf', 
    width = 8,
    height = 4)
shared_MetaRNN_rankscores_mean <- as.data.frame(mean((shared_extended$MutPred_rankscore.x)))
shared_MetaRNN_rankscores_plot <- ggplot() + 
  geom_density(data = shared_extended, aes(x = shared_extended$MetaRNN_rankscore.x, y=..count../sum(..count..), color = "red")) +
  geom_vline(data=shared_MetaRNN_rankscores_mean, aes(xintercept=`mean((shared_extended$MutPred_rankscore.x))`, color = "blue"),
             linetype="dashed", size=0.75, show.legend = F) +
  theme_minimal() + xlim(c(0,1))+ ylim(c(0, 0.005)) + 
  theme(legend.position="none")+
  labs(x = "MetaRNN rank score", y = "scaled density") 
#+ scale_color_manual(values = colors)
plot(shared_MetaRNN_rankscores_plot)
dev.off()

### Save extended datasets
write.table(somatic_NSFP_extended, "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/somatic_NSFP_extended.txt", sep = "\t", quote = F, row.names = F)
write.table(germline_NSFP_extended, "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/germline_NSFP_extended.txt", sep = "\t", quote = F, row.names = F)
write.table(shared_extended, "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/shared_NSFP_extended.txt", sep = "\t", quote = F, row.names = F)




### Load extended datasets
somatic_NSFP_extended <- read.csv("/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/somatic_NSFP_extended.txt", sep = ",")
germline_NSFP_extended <- read.csv("/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/germline_NSFP_extended_mut.txt", sep = ",")
shared_extended <- read.csv("/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/shared_NSFP_extended.txt", sep = "\t")

somatic_NSFP_extended$MetaRNN_score <- as.numeric(somatic_NSFP_extended$MetaRNN_score)
somatic_NSFP_extended$MetaRNN_rankscore <- as.numeric(somatic_NSFP_extended$MetaRNN_rankscore)
somatic_NSFP_extended$MutPred_score <- as.numeric(somatic_NSFP_extended$MutPred_score)
somatic_NSFP_extended$MutPred_rankscore <- as.numeric(somatic_NSFP_extended$MutPred_rankscore)

germline_NSFP_extended$MetaRNN_score <- as.numeric(germline_NSFP_extended$MetaRNN_score)
germline_NSFP_extended$MetaRNN_rankscore <- as.numeric(germline_NSFP_extended$MetaRNN_rankscore)
germline_NSFP_extended$MutPred_score <- as.numeric(germline_NSFP_extended$MutPred_score)
germline_NSFP_extended$MutPred_rankscore <- as.numeric(germline_NSFP_extended$MutPred_rankscore)

shared_extended$MetaRNN_score.x <- as.numeric(shared_extended$MetaRNN_score.x)
shared_extended$MetaRNN_rankscore.x <- as.numeric(shared_extended$MetaRNN_rankscore.x)
shared_extended$MutPred_score.x <- as.numeric(shared_extended$MutPred_score.x)
names(shared_extended)[names(shared_extended) == 'MutPred_score.x'] <- 'MutPred_score'
names(shared_extended)[names(shared_extended) == 'MutPred_Top5features.x'] <- 'MutPred_Top5features'
shared_extended$MutPred_rankscore.x <- as.numeric(shared_extended$MutPred_rankscore.x)


###
### GO Enrichment Analysis
###


### Somatic

### Additional filtering step for significant scores

somatic_NSFP_extended <- somatic_NSFP_extended[somatic_NSFP_extended$MetaRNN_rankscore > mean(somatic_NSFP_extended$MetaRNN_rankscore),]
germline_NSFP_extended <- germline_NSFP_extended[germline_NSFP_extended$MetaRNN_rankscore > mean(germline_NSFP_extended$MetaRNN_rankscore),]
shared_extended <- shared_extended[shared_extended$MetaRNN_rankscore.x > mean(shared_extended$MetaRNN_rankscore.x) | shared_extended$MetaRNN_rankscore.y > mean(shared_extended$MetaRNN_rankscore.y),]

### or

somatic_NSFP_extended <- somatic_NSFP_extended[somatic_NSFP_extended$MetaRNN_rankscore > 0.75,]
germline_NSFP_extended <- germline_NSFP_extended[germline_NSFP_extended$MetaRNN_rankscore > 0.75,]
shared_extended <- shared_extended[shared_extended$MetaRNN_rankscore.x > 0.75 | shared_extended$MetaRNN_rankscore.y > 0.75,]
###

somatic_GOprocess = as.data.frame(somatic_NSFP_extended$GO_biological_process)
somatic_GOproc = unlist(str_split(somatic_GOprocess$`somatic_NSFP_extended$GO_biological_process`, ";"))
somatic_GOprocess = as.data.frame((table(somatic_GOproc)))
somatic_GOprocess = somatic_GOprocess[-1,]
somatic_GOprocess = somatic_GOprocess[order(-somatic_GOprocess$Freq), ]
somatic_GOcompartment = as.data.frame(((somatic_NSFP_extended$GO_cellular_component)))
somatic_GOcomp = unlist(str_split(somatic_GOcompartment$`((somatic_NSFP_extended$GO_cellular_component))`, ";"))
somatic_GOcompartment = as.data.frame(table(somatic_GOcomp))
somatic_GOcompartment = somatic_GOcompartment[-1,]
somatic_GOcompartment = somatic_GOcompartment[order(-somatic_GOcompartment$Freq), ]
somatic_GOfunction = as.data.frame(((somatic_NSFP_extended$GO_molecular_function)))
somatic_GOfunc = unlist(str_split(somatic_GOfunction$`((somatic_NSFP_extended$GO_molecular_function))`, ";"))
somatic_GOfunction = as.data.frame(table(somatic_GOfunc))
somatic_GOfunction = somatic_GOfunction[-1,]
somatic_GOfunction = somatic_GOfunction[order(-somatic_GOfunction$Freq), ]


### GO top 10 plots

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/GO Enrichment/somatic_cut_0.75.pdf', 
    width = 8,
    height = 4)
top_func_somatic <- somatic_GOfunction[1:10,]
#top_func_somatic$Var1[top_func_somatic$Freq == 34081] <- "test"
top_func_somatic_plot <-
  ggplot(top_func_somatic, aes(
    x = reorder(somatic_GOfunc, -desc(Freq)),
    y = Freq,
    fill = somatic_GOfunc,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.8)) +
  ggtitle("Top 10 GO molecular functions in somatic dataset")
#top_func_somatic_plot
plot(top_func_somatic_plot)

top_func_somatic_pie = ggplot(top_func_somatic, aes(x="", y=Freq, fill=somatic_GOfunc)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
#top_func_somatic_pie
plot(top_func_somatic_pie)

###

top_comp_somatic <- somatic_GOcompartment[1:10,]
#top_func_somatic$Var1[top_func_somatic$Freq == 34081] <- "test"
top_comp_somatic_plot <-
  ggplot(top_comp_somatic, aes(
    x = reorder(somatic_GOcomp, -desc(Freq)),
    y = Freq,
    fill = somatic_GOcomp,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.8)) +
  ggtitle("Top 10 GO cellular compartments in somatic dataset")
plot(top_comp_somatic_plot)

top_comp_somatic_pie = ggplot(top_comp_somatic, aes(x="", y=Freq, fill=somatic_GOcomp)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
plot(top_comp_somatic_pie)

###

top_proc_somatic <- somatic_GOprocess[1:10,]
#top_func_somatic$Var1[top_func_somatic$Freq == 34081] <- "test"
top_proc_somatic_plot <-
  ggplot(top_proc_somatic, aes(
    x = reorder(somatic_GOproc, -desc(Freq)),
    y = Freq,
    fill = somatic_GOproc,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.8)) +
  ggtitle("Top 10 GO biological processes in somatic dataset")
plot(top_proc_somatic_plot)

top_proc_somatic_pie = ggplot(top_proc_somatic, aes(x="", y=Freq, fill=somatic_GOproc)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
plot(top_proc_somatic_pie)

dev.off()

### Germline

germline_GOprocess = as.data.frame(germline_NSFP_extended$GO_biological_process)
germline_GOproc = unlist(str_split(germline_GOprocess$`germline_NSFP_extended$GO_biological_process`, ";"))
germline_GOprocess = as.data.frame((table(germline_GOproc)))
germline_GOprocess = germline_GOprocess[-1,]
germline_GOprocess = germline_GOprocess[order(-germline_GOprocess$Freq), ]
germline_GOcompartment = as.data.frame(((germline_NSFP_extended$GO_cellular_component)))
germline_GOcomp = unlist(str_split(germline_GOcompartment$`((germline_NSFP_extended$GO_cellular_component))`, ";"))
germline_GOcompartment = as.data.frame(table(germline_GOcomp))
germline_GOcompartment = germline_GOcompartment[-1,]
germline_GOcompartment = germline_GOcompartment[order(-germline_GOcompartment$Freq), ]
germline_GOfunction = as.data.frame(((germline_NSFP_extended$GO_molecular_function)))
germline_GOfunc = unlist(str_split(germline_GOfunction$`((germline_NSFP_extended$GO_molecular_function))`, ";"))
germline_GOfunction = as.data.frame(table(germline_GOfunc))
germline_GOfunction = germline_GOfunction[-1,]
germline_GOfunction = germline_GOfunction[order(-germline_GOfunction$Freq), ]


### GO top 10 plots

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/GO Enrichment/germline_cut_0.75.pdf', 
    width = 8,
    height = 4)

top_func_germline <- germline_GOfunction[1:10,]
top_func_germline_plot <-
  ggplot(top_func_germline, aes(
    x = reorder(germline_GOfunc, -desc(Freq)),
    y = Freq,
    fill = germline_GOfunc,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.8)) +
  ggtitle("Top 10 GO molecular functions in germline dataset")
plot(top_func_germline_plot)

top_func_germline_pie = ggplot(top_func_germline, aes(x="", y=Freq, fill=germline_GOfunc)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
plot(top_func_germline_pie)

###

top_comp_germline <- germline_GOcompartment[1:10,]
top_comp_germline_plot <-
  ggplot(top_comp_germline, aes(
    x = reorder(germline_GOcomp, -desc(Freq)),
    y = Freq,
    fill = germline_GOcomp,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.8)) +
  ggtitle("Top 10 GO cellular compartments in germline dataset")
plot(top_comp_germline_plot)

top_comp_germline_pie = ggplot(top_comp_germline, aes(x="", y=Freq, fill=germline_GOcomp)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
plot(top_comp_germline_pie)

###

top_proc_germline <- germline_GOprocess[1:10,]
top_proc_germline_plot <-
  ggplot(top_proc_germline, aes(
    x = reorder(germline_GOproc, -desc(Freq)),
    y = Freq,
    fill = germline_GOproc,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.8)) +
  ggtitle("Top 10 GO biological processes in germline dataset")
plot(top_proc_germline_plot)

top_proc_germline_pie = ggplot(top_proc_germline, aes(x="", y=Freq, fill=germline_GOproc)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
plot(top_proc_germline_pie)

dev.off()

### Shared mutations

shared_GOprocess = as.data.frame(shared_extended$GO_biological_process.x)
shared_GOproc = unlist(str_split(shared_GOprocess$`shared_extended$GO_biological_process.x`, ";"))
shared_GOprocess = as.data.frame((table(shared_GOproc)))
shared_GOprocess = shared_GOprocess[-1,]
shared_GOprocess = shared_GOprocess[order(-shared_GOprocess$Freq), ]
shared_GOcompartment = as.data.frame(((shared_extended$GO_cellular_component.x)))
shared_GOcomp = unlist(str_split(shared_GOcompartment$`((shared_extended$GO_cellular_component.x))`, ";"))
shared_GOcompartment = as.data.frame(table(shared_GOcomp))
shared_GOcompartment = shared_GOcompartment[-1,]
shared_GOcompartment = shared_GOcompartment[order(-shared_GOcompartment$Freq), ]
shared_GOfunction = as.data.frame(((shared_extended$GO_molecular_function.x)))
shared_GOfunc = unlist(str_split(shared_GOfunction$`((shared_extended$GO_molecular_function.x))`, ";"))
shared_GOfunction = as.data.frame(table(shared_GOfunc))
shared_GOfunction = shared_GOfunction[-1,]
shared_GOfunction = shared_GOfunction[order(-shared_GOfunction$Freq), ]


### GO top 10 plots

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/GO Enrichment/shared_cut.pdf', 
    width = 8,
    height = 4)

top_func_shared <- shared_GOfunction[1:10,]
top_func_shared_plot <-
  ggplot(top_func_shared, aes(
    x = reorder(shared_GOfunc, -desc(Freq)),
    y = Freq,
    fill = shared_GOfunc,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.8)) +
  ggtitle("Top 10 GO molecular functions in shared dataset")
plot(top_func_shared_plot)

top_func_shared_pie = ggplot(top_func_shared, aes(x="", y=Freq, fill=shared_GOfunc)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
plot(top_func_shared_pie)

###

top_comp_shared <- shared_GOcompartment[1:10,]
top_comp_shared_plot <-
  ggplot(top_comp_shared, aes(
    x = reorder(shared_GOcomp, -desc(Freq)),
    y = Freq,
    fill = shared_GOcomp,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.8)) +
  ggtitle("Top 10 GO cellular compartments in shared dataset")
plot(top_comp_shared_plot)

top_comp_shared_pie = ggplot(top_comp_shared, aes(x="", y=Freq, fill=shared_GOcomp)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
plot(top_comp_shared_pie)

###

top_proc_shared <- shared_GOprocess[1:10,]
top_proc_shared_plot <-
  ggplot(top_proc_shared, aes(
    x = reorder(shared_GOproc, -desc(Freq)),
    y = Freq,
    fill = shared_GOproc,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.8)) +
  ggtitle("Top 10 GO biological processes in shared dataset")
plot(top_proc_shared_plot)

top_proc_shared_pie = ggplot(top_proc_shared, aes(x="", y=Freq, fill=shared_GOproc)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
plot(top_proc_shared_pie)


dev.off()

### Combine GO summaries
# biological processes 
colnames(somatic_GOprocess)[1] = "GO_proc"
somatic_GOprocess_top = somatic_GOprocess[1:10,]

somatic_GOprocess_top = rbind(somatic_GOprocess_top, somatic_GOprocess[somatic_GOprocess$GO_proc == "neutrophil degranulation",])

somatic_GOprocess_top$'dataset' = "somatic"
total_somatic_GOproc = sum(somatic_GOprocess$Freq)
colnames(germline_GOprocess)[1] = "GO_proc"
germline_GOprocess_top = germline_GOprocess[1:10,]
germline_GOprocess_top$'dataset' = "germline"
total_germline_GOproc = sum(germline_GOprocess$Freq)
colnames(shared_GOprocess)[1] = "GO_proc"
shared_GOprocess_top = shared_GOprocess[1:10,]

shared_GOprocess_top = rbind(shared_GOprocess_top, shared_GOprocess[shared_GOprocess$GO_proc == "regulation of transcription by RNA polymerase II",])
shared_GOprocess_top$'dataset' = "shared"
total_shared_GOproc = sum(shared_GOprocess$Freq)
dataset_GOproc = rbind(somatic_GOprocess_top, germline_GOprocess_top, shared_GOprocess_top)
dataset_GOproc$'percent' = 0.0
dataset_GOproc$percent[dataset_GOproc$dataset == "somatic"] = dataset_GOproc$Freq[dataset_GOproc$dataset == "somatic"]/total_somatic_GOproc
dataset_GOproc$percent[dataset_GOproc$dataset == "germline"] = dataset_GOproc$Freq[dataset_GOproc$dataset == "germline"]/total_germline_GOproc
dataset_GOproc$percent[dataset_GOproc$dataset == "shared"] = dataset_GOproc$Freq[dataset_GOproc$dataset == "shared"]/total_shared_GOproc

dataset_GOproc = dataset_GOproc[!dataset_GOproc$GO_proc == "extracellular matrix organization",]
dataset_GOproc = dataset_GOproc[!dataset_GOproc$GO_proc == "positive regulation of transcription, DNA-templated",]

# cellular compartments
colnames(somatic_GOcompartment)[1] = "GO_comp"
somatic_GOcompartment_top = somatic_GOcompartment[1:10,]
somatic_GOcompartment_top$'dataset' = "somatic"
total_somatic_GOcomp = sum(somatic_GOcompartment$Freq)
colnames(germline_GOcompartment)[1] = "GO_comp"
germline_GOcompartment_top = germline_GOcompartment[1:10,]
germline_GOcompartment_top$'dataset' = "germline"
total_germline_GOcomp = sum(germline_GOcompartment$Freq)
colnames(shared_GOcompartment)[1] = "GO_comp"
shared_GOcompartment_top = shared_GOcompartment[1:10,]
shared_GOcompartment_top$'dataset' = "shared"
total_shared_GOcomp = sum(shared_GOcompartment$Freq)
dataset_GOcomp = rbind(somatic_GOcompartment_top, germline_GOcompartment_top, shared_GOcompartment_top)
dataset_GOcomp$'percent' = 0.0
dataset_GOcomp$percent[dataset_GOcomp$dataset == "somatic"] = dataset_GOcomp$Freq[dataset_GOcomp$dataset == "somatic"]/total_somatic_GOcomp
dataset_GOcomp$percent[dataset_GOcomp$dataset == "germline"] = dataset_GOcomp$Freq[dataset_GOcomp$dataset == "germline"]/total_germline_GOcomp
dataset_GOcomp$percent[dataset_GOcomp$dataset == "shared"] = dataset_GOcomp$Freq[dataset_GOcomp$dataset == "shared"]/total_shared_GOcomp

# molecular functions
colnames(somatic_GOfunction)[1] = "GO_func"
somatic_GOfunction_top = somatic_GOfunction[1:10,]
somatic_GOfunction_top$'dataset' = "somatic"
total_somatic_GOfunc = sum(somatic_GOfunction$Freq)
colnames(germline_GOfunction)[1] = "GO_func"
germline_GOfunction_top = germline_GOfunction[1:10,]
germline_GOfunction_top$'dataset' = "germline"
total_germline_GOfunc = sum(germline_GOfunction$Freq)
colnames(shared_GOfunction)[1] = "GO_func"
shared_GOfunction_top = shared_GOfunction[1:10,]
shared_GOfunction_top$'dataset' = "shared"
total_shared_GOfunc = sum(shared_GOfunction$Freq)
dataset_GOfunc = rbind(somatic_GOfunction_top, germline_GOfunction_top, shared_GOfunction_top)
dataset_GOfunc$'percent' = 0.0
dataset_GOfunc$percent[dataset_GOfunc$dataset == "somatic"] = dataset_GOfunc$Freq[dataset_GOfunc$dataset == "somatic"]/total_somatic_GOfunc
dataset_GOfunc$percent[dataset_GOfunc$dataset == "germline"] = dataset_GOfunc$Freq[dataset_GOfunc$dataset == "germline"]/total_germline_GOfunc
dataset_GOfunc$percent[dataset_GOfunc$dataset == "shared"] = dataset_GOfunc$Freq[dataset_GOfunc$dataset == "shared"]/total_shared_GOfunc

### Plot summarized GO results (3 in 1)
colors <- c("somatic" = "red", "germline" = "blue", "shared" = "green")
Go_plot_features_combo  = function (dF, level, filtStage){
  feature = dF[,1]
  freq = dF[,2]
  logFreq = log10(freq)
  dataset = dF[,3]
  percentage = round((dF[,4]*100), digits = 2)
  pdf(paste0('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/GO Enrichment/', level, '_', filtStage, '_combined_bar.pdf'), 
      width = 10,
      height = 8)
  MutPred_plot <-
    ggplot(dF, aes(
      #  x = reorder(feature, -dplyr::desc(freq)),
      x = reorder(feature, dplyr::desc(feature)),
      y = logFreq,
      fill = dataset,
      label = paste0(freq, " (", percentage, "%)")
    )) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    coord_flip() +
    labs(x = NULL, y = NULL, color= "dataset") +
    scale_y_continuous(labels = comma) +
    ylab(paste("count and fraction of all", level, "terms in the dataset \n(axis log-transformed)")) +
    xlab(level) +
    geom_text(size = 2, position = position_stack(vjust = 0.70)) +
    theme(axis.text.y = element_text(angle = 40, vjust = 1, hjust=1, size = 8)) +
    ggtitle(paste("Top 10", level, "GO terms \n in all three datasets", paste0("(",filtStage, " filtering stage)")))
  
  plot(MutPred_plot)
  
  dev.off()
  pdf(paste0('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/GO Enrichment/', level, '_combined_pie_small.pdf'), 
      width = 8,
      height = 4)
  
  MutPred_pie = ggplot(dF, aes(x=dataset, y=logFreq, fill=feature, label=freq)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) 
  #    theme_void()
  plot(MutPred_pie)
  
  dev.off()
}
### Execute modified function

# first filtering level
Go_plot_features_combo(dataset_GOproc, "biological process", "first")
Go_plot_features_combo(dataset_GOcomp, "cellular compartment", "first")
Go_plot_features_combo(dataset_GOfunc, "molecular function", "first")

# second filtering level
Go_plot_features_combo(dataset_GOproc, "biological process", "second")
Go_plot_features_combo(dataset_GOcomp, "cellular compartment", "second")
Go_plot_features_combo(dataset_GOfunc, "molecular function", "second")

# third filtering level
Go_plot_features_combo(dataset_GOproc, "biological process", "third")
Go_plot_features_combo(dataset_GOcomp, "cellular compartment", "third")
Go_plot_features_combo(dataset_GOfunc, "molecular function", "third")

###
### Most mutated genes analysis
###

### Data loading


somatic_NSFP_ext <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/somatic_NSFP_pos.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "genename",
      "MetaRNN_score",
      "MetaRNN_rankscore",
      "MutPred_score",
      "MutPred_rankscore",
      "MutPred_Top5features",
      "GO_biological_process",
      "GO_cellular_component",
      "GO_molecular_function"
      
    )
  )


somatic_NSFP_ext <- somatic_NSFP_ext[somatic_NSFP_ext$MetaRNN_rankscore != "."]
somatic_NSFP_ext <- somatic_NSFP_ext[somatic_NSFP_ext$genename != "."]
somatic_NSFP_ext$genename <- sapply(str_split(somatic_NSFP_ext$genename, ";"), tail, 1)
somatic_NSFP_ext$MetaRNN_rankscore <- as.numeric(somatic_NSFP_ext$MetaRNN_rankscore)
somatic_NSFP_ext$MetaRNN_score <- sapply(str_split(somatic_NSFP_ext$MetaRNN_score, ";"), tail, 1)
somatic_NSFP_ext <- somatic_NSFP_ext[somatic_NSFP_ext$MetaRNN_score != "."]
somatic_NSFP_ext$MetaRNN_score <- as.numeric(somatic_NSFP_ext$MetaRNN_score)
somatic_NSFP_ext <- somatic_NSFP_ext[somatic_NSFP_ext$MutPred_score != "."]
somatic_NSFP_ext$MutPred_score <- as.numeric(somatic_NSFP_ext$MutPred_score)
somatic_NSFP_ext <- somatic_NSFP_ext[somatic_NSFP_ext$MutPred_rankscore != "."]
somatic_NSFP_ext$MutPred_rankscore <- as.numeric(somatic_NSFP_ext$MutPred_rankscore)
somatic_NSFP_ext = somatic_NSFP_ext %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_score`))


germline_NSFP_ext <-
  fread(
    "/Users/mcpftw/dbNSFP4.2a/germline_NSFP_pos.txt",
    select = c(
      "#chr",
      "pos(1-based)",
      "ref",
      "alt",
      "genename",
      "MetaRNN_score",
      "MetaRNN_rankscore",
      "MutPred_score",
      "MutPred_rankscore",
      "MutPred_Top5features",
      "GO_biological_process",
      "GO_cellular_component",
      "GO_molecular_function"
    )
  )

germline_NSFP_ext <- germline_NSFP_ext[germline_NSFP_ext$MetaRNN_rankscore != "."]
germline_NSFP_ext <- germline_NSFP_ext[germline_NSFP_ext$genename != "."]
germline_NSFP_ext$genename <- sapply(str_split(germline_NSFP_ext$genename, ";"), tail, 1)
germline_NSFP_ext$MetaRNN_rankscore <- as.numeric(germline_NSFP_ext$MetaRNN_rankscore)
germline_NSFP_ext$MetaRNN_score <- sapply(str_split(germline_NSFP_ext$MetaRNN_score, ";"), tail, 1)
germline_NSFP_ext <- germline_NSFP_ext[germline_NSFP_ext$MetaRNN_score != "."]
germline_NSFP_ext$MetaRNN_score <- as.numeric(germline_NSFP_ext$MetaRNN_score)
germline_NSFP_ext <- germline_NSFP_ext[germline_NSFP_ext$MutPred_score != "."]
germline_NSFP_ext$MutPred_score <- as.numeric(germline_NSFP_ext$MutPred_score)
germline_NSFP_ext <- germline_NSFP_ext[germline_NSFP_ext$MutPred_rankscore != "."]
germline_NSFP_ext$MutPred_rankscore <- as.numeric(germline_NSFP_ext$MutPred_rankscore)
germline_NSFP_ext = germline_NSFP_ext %>% group_by(`#chr`, `pos(1-based)`, ref, alt) %>% dplyr::slice(which.max(`MetaRNN_score`))

### Extract shared mutations and remove them from somatic and germline sets
germline_NSFP_extended <- anti_join(germline_NSFP_ext, somatic_NSFP_ext, by = c("#chr", "pos(1-based)", "ref", "alt"))
somatic_NSFP_extended <- anti_join(somatic_NSFP_ext, germline_NSFP_ext, by = c("#chr", "pos(1-based)", "ref", "alt"))
shared_extended <- inner_join(somatic_NSFP_ext, germline_NSFP_ext, by = c("#chr", "pos(1-based)", "ref", "alt"))
names(shared_extended)[names(shared_extended) == 'MutPred_score.x'] <- 'MutPred_score'
names(shared_extended)[names(shared_extended) == 'MutPred_Top5features.x'] <- 'MutPred_Top5features'
names(shared_extended)[names(shared_extended) == 'MetaRNN_rankscore.x'] <- 'MetaRNN_rankscore'
names(shared_extended)[names(shared_extended) == 'genename.x'] <- 'genename'
rm(germline_NSFP_ext, somatic_NSFP_ext)

### Extract top mutated genes with frequency and average score

gene_count_somatic = as.data.frame(table(somatic_NSFP_extended$genename))
colnames(gene_count_somatic) = c("genename", "Freq")
gene_count_germline = as.data.frame(table(germline_NSFP_extended$genename))
colnames(gene_count_germline) = c("genename", "Freq")
gene_count_shared = as.data.frame(table(shared_extended$genename))
colnames(gene_count_shared) = c("genename", "Freq")
gene_score_somatic = aggregate( MetaRNN_rankscore ~ genename, somatic_NSFP_extended, mean)
gene_score_germline = aggregate( MetaRNN_rankscore ~ genename, germline_NSFP_extended, mean)
gene_score_shared = aggregate( MetaRNN_rankscore ~ genename, shared_extended, mean)
genes_somatic = join(gene_count_somatic, gene_score_somatic, by="genename")
genes_germline = join(gene_count_germline, gene_score_germline, by="genename")
genes_shared = join(gene_count_shared, gene_score_shared, by="genename")

genes_somatic = genes_somatic[order(-genes_somatic$Freq),]
genes_germline = genes_germline[order(-genes_germline$Freq),]
genes_shared = genes_shared[order(-genes_shared$Freq),]

### Plot results

# Somatic 
top_genes_somatic = genes_somatic[1:10,]
top_genes_somatic_scores = data.frame(genename = somatic_NSFP_extended$genename, MetaRNN_rankscore = somatic_NSFP_extended$MetaRNN_rankscore)
top_genes_somatic_scores = top_genes_somatic_scores[top_genes_somatic_scores$genename %in% top_genes_somatic$genename,]
top_genes_somatic_scores$genename <- factor(top_genes_somatic_scores$genename, levels=c("USH2A", "LRP1B", "NEB", "CSMD3", "SYNE1", "RYR2", "OBSCN", "FLG", "MUC16", "TTN"))
colors <- rainbow(length(top_genes_somatic$genename))
names(colors)  <- top_genes_somatic$genename
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/Mutated Genes Analysis/most_mutated_genes_somatic.pdf', 
    width = 8,
    height = 4)
most_mutated_somatic_plot <-
  ggplot(top_genes_somatic, aes(
    x = reorder(genename, -dplyr::desc(Freq)),
    y = Freq,
    fill = genename,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("acummulated mutations") +
  xlab("gene") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors)
most_mutated_somatic_boxplot <- 
  ggplot(top_genes_somatic_scores, aes(
    x = genename,
    y = MetaRNN_rankscore,
    fill = genename,
  )) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("MetaRNN rank scores") +
  xlab("gene")+
  scale_fill_manual(values = colors)
cowplot::plot_grid(most_mutated_somatic_plot + ggtitle(""),
                   most_mutated_somatic_boxplot + ggtitle("") + theme(axis.text.y = element_blank(),
                                                                      axis.ticks.y = element_blank(),
                                                                      axis.title.y = element_blank()) ,
                   nrow=1, labels = c("A", "B"), align = "h", axis = "lr", rel_heights = c(1,1))
dev.off()

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/Mutated Genes Analysis/genes_scores_freq_density_somatic.pdf', 
    width = 8,
    height = 4)
ggplot(data=genes_somatic, aes(x=MetaRNN_rankscore, y=Freq)) +
  geom_point(alpha=1/5)+
  coord_trans(y = "sqrt") +
  theme_minimal() + 
  labs(x = NULL, y = NULL) +
  ylab("Mutation frequency") +
  xlab("MetaRNN rank scores") 
#+  ggtitle("Scatter plot of mutation frequency versus MetaRNN rank scores of genes in somatic dataset")
dev.off()

# Germline

top_genes_germline = genes_germline[1:10,]
top_genes_germline_scores = data.frame(genename = germline_NSFP_extended$genename, MetaRNN_rankscore = germline_NSFP_extended$MetaRNN_rankscore)
top_genes_germline_scores = top_genes_germline_scores[top_genes_germline_scores$genename %in% top_genes_germline$genename,]
top_genes_germline_scores$genename <- factor(top_genes_germline_scores$genename, levels=c("PKD1", "SYNE2", "PLEC", "SYNE1", "NEB", "MUC4", "OBSCN", "AHNAK2", "MUC16", "TTN"))
colors <- rainbow(length(top_genes_germline$genename))
names(colors)  <- top_genes_germline$genename
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/Mutated Genes Analysis/most_mutated_genes_germline.pdf', 
    width = 8,
    height = 4)
most_mutated_germline_plot <-
  ggplot(top_genes_germline, aes(
    x = reorder(genename, -dplyr::desc(Freq)),
    y = Freq,
    fill = genename,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("acummulated mutations") +
  xlab("gene") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors)
most_mutated_germline_boxplot <- 
  ggplot(top_genes_germline_scores, aes(
    x = genename,
    y = MetaRNN_rankscore,
    fill = genename,
  )) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("MetaRNN rank scores") +
  xlab("gene")+
  scale_fill_manual(values = colors)
cowplot::plot_grid(most_mutated_germline_plot + ggtitle(""),
                   most_mutated_germline_boxplot + ggtitle("") + theme(axis.text.y = element_blank(),
                                                                      axis.ticks.y = element_blank(),
                                                                      axis.title.y = element_blank()) ,
                   nrow=1, labels = c("A", "B"), align = "h", axis = "lr", rel_heights = c(1,1))
dev.off()

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/Mutated Genes Analysis/genes_scores_freq_density_germline.pdf', 
    width = 8,
    height = 4)
ggplot(data=genes_germline, aes(x=MetaRNN_rankscore, y=Freq)) +
  geom_point(alpha=1/5)+
  coord_trans(y = "sqrt") +
  theme_minimal() + 
  labs(x = NULL, y = NULL) +
  ylab("Mutation frequency") +
  xlab("MetaRNN rank scores") 
#+  ggtitle("Scatter plot of mutation frequency versus MetaRNN rank scores of genes in somatic dataset")
dev.off()

# Shared

top_genes_shared = genes_shared[1:10,]
top_genes_shared_scores = data.frame(genename = shared_extended$genename, MetaRNN_rankscore = shared_extended$MetaRNN_rankscore)
top_genes_shared_scores = top_genes_shared_scores[top_genes_shared_scores$genename %in% top_genes_shared$genename,]
top_genes_shared_scores$genename <- factor(top_genes_shared_scores$genename, levels=c("NEB", "SYNE1", "MUC17", "HRNR", "AHNAK2", "OBSCN", "MUC4", "MUC16", "FLG", "TTN"))
colors <- rainbow(length(top_genes_shared$genename))
names(colors)  <- top_genes_shared$genename
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/Mutated Genes Analysis/most_mutated_genes_shared.pdf', 
    width = 8,
    height = 4)
most_mutated_shared_plot <-
  ggplot(top_genes_shared, aes(
    x = reorder(genename, -dplyr::desc(Freq)),
    y = Freq,
    fill = genename,
    label = Freq
  )) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("acummulated mutations") +
  xlab("gene") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = colors)
most_mutated_shared_boxplot <- 
  ggplot(top_genes_shared_scores, aes(
    x = genename,
    y = MetaRNN_rankscore,
    fill = genename,
  )) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(labels = comma) +
  ylab("MetaRNN rank scores") +
  xlab("gene")+
  scale_fill_manual(values = colors)
cowplot::plot_grid(most_mutated_shared_plot + ggtitle(""),
                   most_mutated_shared_boxplot + ggtitle("") + theme(axis.text.y = element_blank(),
                                                                      axis.ticks.y = element_blank(),
                                                                      axis.title.y = element_blank()) ,
                   nrow=1, labels = c("A", "B"), align = "h", axis = "lr", rel_heights = c(1,1))
dev.off()

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/Mutated Genes Analysis/genes_scores_freq_density_shared.pdf', 
    width = 8,
    height = 4)
ggplot(data=genes_shared, aes(x=MetaRNN_rankscore, y=Freq)) +
  geom_point(alpha=1/5)+
  coord_trans(y = "sqrt") +
  theme_minimal() + 
  labs(x = NULL, y = NULL) +
  ylab("Mutation frequency") +
  xlab("MetaRNN rank scores") 
#+  ggtitle("Scatter plot of mutation frequency versus MetaRNN rank scores of genes in shared dataset")
dev.off()

###
### MutPred predictions ###
###

regexp <- "([:space:])([[:digit:]]+)(\\.*)([[:digit:]]*)"

### Somatic MutPred

somatic_MutPred_pvals <-  vector(mode="list", nrow(somatic_NSFP_extended))
somatic_MutPred_actionable <-  vector(mode="list", nrow(somatic_NSFP_extended))
somatic_MutPred_confident <-  vector(mode="list", nrow(somatic_NSFP_extended))
somatic_MutPred_v_confident <-  vector(mode="list", nrow(somatic_NSFP_extended))

somatic_NSFP_extended$MutPred_Top5features <- as.list(sapply(somatic_NSFP_extended$MutPred_Top5features, str_split, ";"))
somatic_MutPred_pvals <- lapply(somatic_NSFP_extended$MutPred_Top5features, lapply, str_extract, regexp)
somatic_NSFP_extended$MutPred_pvals <- somatic_MutPred_pvals
somatic_NSFP_extended$MutPred_pvals <- lapply(somatic_NSFP_extended$MutPred_pvals, as.numeric)


for (i in 1:nrow(somatic_NSFP_extended)) {
  if (i %% 1000 == 0){
    print(i)
  }
  if (somatic_NSFP_extended$MutPred_score[[i]] > 0.5) {
    for (j in 1:length(somatic_NSFP_extended$MutPred_Top5features[[i]])){
      if (somatic_NSFP_extended$MutPred_score[[i]] > 0.5 & somatic_NSFP_extended$MutPred_pvals[[i]][j] < 0.05) {
        somatic_MutPred_actionable[[i]][j] = somatic_NSFP_extended$MutPred_Top5features[[i]][j]
      }
      if (somatic_NSFP_extended$MutPred_score[[i]] > 0.75 & somatic_NSFP_extended$MutPred_pvals[[i]][j] < 0.05) {
        somatic_MutPred_confident[[i]][j] <- somatic_NSFP_extended$MutPred_Top5features[[i]][j]
      }
      if (somatic_NSFP_extended$MutPred_score[[i]] > 0.75 & somatic_NSFP_extended$MutPred_pvals[[i]][j] < 0.01) {
        somatic_MutPred_v_confident[[i]][j] <- somatic_NSFP_extended$MutPred_Top5features[[i]][j]
      }
      
    }
  }
}

somatic_NSFP_extended$MutPred_actionable = somatic_MutPred_actionable
somatic_NSFP_extended$MutPred_confident = somatic_MutPred_confident
somatic_NSFP_extended$MutPred_v_confident = somatic_MutPred_v_confident
rm(somatic_MutPred_pvals, somatic_MutPred_actionable, somatic_MutPred_confident, somatic_MutPred_v_confident)
gc()

### Germline MutPred

germline_MutPred_pvals <-  vector(mode="list", nrow(germline_NSFP_extended))
germline_MutPred_actionable <-  vector(mode="list", nrow(germline_NSFP_extended))
germline_MutPred_confident <-  vector(mode="list", nrow(germline_NSFP_extended))
germline_MutPred_v_confident <-  vector(mode="list", nrow(germline_NSFP_extended))

germline_NSFP_extended$MutPred_Top5features <- as.list(sapply(germline_NSFP_extended$MutPred_Top5features, str_split, ";"))
germline_MutPred_pvals <- lapply(germline_NSFP_extended$MutPred_Top5features, lapply, str_extract, regexp)
germline_NSFP_extended$MutPred_pvals <- germline_MutPred_pvals
germline_NSFP_extended$MutPred_pvals <- lapply(germline_NSFP_extended$MutPred_pvals, as.numeric)

for (i in 1:nrow(germline_NSFP_extended)) {
  if (i %% 1000 == 0){
    print(i)
  }
  if (germline_NSFP_extended$MutPred_score[[i]] > 0.5) {
    for (j in 1:length(germline_NSFP_extended$MutPred_Top5features[[i]])){
      #      print(i)
      #      print(j)
      #      print(germline_NSFP_extended$MutPred_pvals[[i]][j])
      #      print(germline_NSFP_extended$MutPred_Top5features[[i]][j])
      if (germline_NSFP_extended$MutPred_score[[i]] > 0.5 & germline_NSFP_extended$MutPred_pvals[[i]][j] < 0.05) {
        germline_MutPred_actionable[[i]][j] = germline_NSFP_extended$MutPred_Top5features[[i]][j]
      }
      if (germline_NSFP_extended$MutPred_score[[i]] > 0.75 & germline_NSFP_extended$MutPred_pvals[[i]][j] < 0.05) {
        germline_MutPred_confident[[i]][j] <- germline_NSFP_extended$MutPred_Top5features[[i]][j]
      }
      if (germline_NSFP_extended$MutPred_score[[i]] > 0.75 & germline_NSFP_extended$MutPred_pvals[[i]][j] < 0.01) {
        germline_MutPred_v_confident[[i]][j] <- germline_NSFP_extended$MutPred_Top5features[[i]][j]
      }
      
    }
  }
}
germline_NSFP_extended$MutPred_actionable = germline_MutPred_actionable
germline_NSFP_extended$MutPred_confident = germline_MutPred_confident
germline_NSFP_extended$MutPred_v_confident = germline_MutPred_v_confident

rm(germline_MutPred_actionable, germline_MutPred_confident, germline_MutPred_v_confident, germline_MutPred_pvals)
gc()

### Shared MutPred

shared_MutPred_pvals <-  vector(mode="list", nrow(shared_extended))
shared_MutPred_actionable <-  vector(mode="list", nrow(shared_extended))
shared_MutPred_confident <-  vector(mode="list", nrow(shared_extended))
shared_MutPred_v_confident <-  vector(mode="list", nrow(shared_extended))

shared_extended$MutPred_Top5features.x <- as.list(sapply(shared_extended$MutPred_Top5features.x, str_split, ";"))
shared_MutPred_pvals <- lapply(shared_extended$MutPred_Top5features.x, lapply, str_extract, regexp)
shared_extended$MutPred_pvals <- shared_MutPred_pvals 
shared_extended$MutPred_pvals <- lapply(shared_extended$MutPred_pvals, as.numeric) 

for (i in 1:nrow(shared_extended)) {
  if (i %% 1000 == 0){
    print(i)
  }
  if (shared_extended$MutPred_score.x[[i]] > 0.5) {
    for (j in 1:length(shared_extended$MutPred_Top5features.x[[i]])){
      if (shared_extended$MutPred_score.x[[i]] > 0.5 & shared_extended$MutPred_pvals[[i]][j] < 0.05) {
        shared_MutPred_actionable[[i]][j] = shared_extended$MutPred_Top5features.x[[i]][j]
      }
      if (shared_extended$MutPred_score.x[[i]] > 0.75 & shared_extended$MutPred_pvals[[i]][j] < 0.05) {
        shared_MutPred_confident[[i]][j] <- shared_extended$MutPred_Top5features.x[[i]][j]
      }
      if (shared_extended$MutPred_score.x[[i]] > 0.75 & shared_extended$MutPred_pvals[[i]][j] < 0.01) {
        shared_MutPred_v_confident[[i]][j] <- shared_extended$MutPred_Top5features.x[[i]][j]
      }
      
    }
  }
}

shared_extended$MutPred_actionable = shared_MutPred_actionable
shared_extended$MutPred_confident = shared_MutPred_confident
shared_extended$MutPred_v_confident = shared_MutPred_v_confident

rm(shared_MutPred_actionable, shared_MutPred_confident, shared_MutPred_pvals, shared_MutPred_v_confident)
gc()

### Write results
somatic_NSFP_extended <- apply(somatic_NSFP_extended,2,as.character)
write.csv(somatic_NSFP_extended, "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/somatic_NSFP_extended_mut.txt", sep = "\t", quote = F, row.names = F)
germline_NSFP_extended <- apply(germline_NSFP_extended,2,as.character)
write.csv(germline_NSFP_extended, "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/germline_NSFP_extended_mut.txt", sep = "\t", quote = F, row.names = F)
shared_extended <- apply(shared_extended,2,as.character)
write.csv(shared_extended, "/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/TCGA/shared_NSFP_extended_mut.txt", sep = "\t", quote = F, row.names = F)


### MutPred Prediction analysis

### Func to extract ordered table of features

MutPred_features_toTable <- function (MutPred_list) {
  MutPred_list = unlist(MutPred_list)
  MutPred_list = sapply(str_split(MutPred_list, "\\("), head, 1)
  MutPred_list = sapply(str_split(MutPred_list, " at "), head, 1)
  MutPred_list = trimws(MutPred_list)
  MutPred_list = as.data.frame((table(MutPred_list)))
  MutPred_list =  MutPred_list[order(-MutPred_list$Freq),]
  return(MutPred_list)
}

### Execute func

somatic_Mut_Act = MutPred_features_toTable(somatic_NSFP_extended$MutPred_actionable)
somatic_Mut_Conf = MutPred_features_toTable(somatic_NSFP_extended$MutPred_confident)
somatic_Mut_V_conf = MutPred_features_toTable(somatic_NSFP_extended$MutPred_v_confident)

germline_Mut_Act = MutPred_features_toTable(germline_NSFP_extended$MutPred_actionable)
germline_Mut_Conf = MutPred_features_toTable(germline_NSFP_extended$MutPred_confident)
germline_Mut_V_conf = MutPred_features_toTable(germline_NSFP_extended$MutPred_v_confident)

shared_Mut_Act = MutPred_features_toTable(shared_extended$MutPred_actionable)
shared_Mut_Conf = MutPred_features_toTable(shared_extended$MutPred_confident)
shared_Mut_V_conf = MutPred_features_toTable(shared_extended$MutPred_v_confident)

### Plot results func

MutPred_plot_features  = function (dataset, dt_name, level){
  feature = dataset[,1]
  freq = dataset[,2]
  pdf(paste('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/', dt_name, '_', level, '_bar.pdf'), 
      width = 8,
      height = 4)
  
  MutPred_plot <-
    ggplot(dataset, aes(
      x = reorder(feature, -desc(freq)),
      y = freq,
      fill = feature,
      label = freq
    )) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(legend.position = "none") +
    coord_flip() +
    labs(x = NULL, y = NULL) +
    scale_y_continuous(labels = comma) +
    ylab("") +
    geom_text(size = 3, position = position_stack(vjust = 0.8)) +
    ggtitle(paste("All", level, "MutPred structural features in", dt_name,  "dataset"))
  
  plot(MutPred_plot)
  
  dev.off()
  
  pdf(paste('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/', dt_name, '_', level, '_pie.pdf'), 
      width = 8,
      height = 4)
  
  MutPred_pie = ggplot(dataset, aes(x="", y=Freq, fill=feature)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void()
  plot(MutPred_pie)
  
  dev.off()
}

### Execute func

MutPred_plot_features(somatic_Mut_Act, "somatic", "actionable")
MutPred_plot_features(somatic_Mut_Conf, "somatic", "confident")
MutPred_plot_features(somatic_Mut_V_conf, "somatic", "very confident")

MutPred_plot_features(germline_Mut_Act, "germline", "actionable")
MutPred_plot_features(germline_Mut_Conf, "germline", "confident")
MutPred_plot_features(germline_Mut_V_conf, "germline", "very confident")

MutPred_plot_features(shared_Mut_Act, "shared", "actionable")
MutPred_plot_features(shared_Mut_Conf, "shared", "confident")
MutPred_plot_features(shared_Mut_V_conf, "shared", "very confident")

### Plot summarized results (3 in 1)
colors <- c("somatic" = "red", "germline" = "blue", "shared" = "green")
MutPred_plot_features_combo  = function (dF, level){
  feature = dF[,1]
  freq = dF[,2]
  logFreq = log10(freq)
  dataset = dF[,3]
  percentage = round((dF[,4]*100), digits = 2)
  pdf(paste0('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/', level, '_combined_bar.pdf'), 
      width = 8,
      height = 6)
  MutPred_plot <-
    ggplot(dF, aes(
      #  x = reorder(feature, -dplyr::desc(freq)),
      x = reorder(feature, dplyr::desc(feature)),
      y = logFreq,
      fill = dataset,
      label = paste0(freq, " (", percentage, "%)")
    )) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    coord_flip() +
    labs(x = NULL, y = NULL, color= "dataset") +
    scale_y_continuous(labels = comma) +
    ylab("count and fraction of the dataset \n(axis log-transformed)") +
    xlab("feature change") +
    geom_text(size = 2, position = position_stack(vjust = 0.70)) +
    #    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 10)) +
    ggtitle(paste("Total", level, "MutPred structural and functional feature changes"))
  
  plot(MutPred_plot)
  
  dev.off()
  pdf(paste0('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/', level, '_combined_pie_small.pdf'), 
      width = 8,
      height = 4)
  
  MutPred_pie = ggplot(dF, aes(x=dataset, y=logFreq, fill=feature, label=freq)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) 
  #    theme_void()
  plot(MutPred_pie)
  
  dev.off()
}
### Execute modified function

MutPred_plot_features_combo(dataset_Mut_V_conf, "very confident")
MutPred_plot_features_combo(dataset_Mut_Conf, "confident")
MutPred_plot_features_combo(dataset_Mut_Act, "actionable")

### All three confidence levels in one plot
library(gridExtra)
library(cowplot)
MutPred_plot_features_combined  = function (dF, level){
  feature = dF[,1]
  freq = dF[,2]
  logFreq = log10(freq)
  dataset = dF[,3]
  percentage = round((dF[,4]*100), digits = 2)
  
  MutPred_plot <-
    ggplot(dF, aes(
      #  x = reorder(feature, -dplyr::desc(freq)),
      x = feature,
      y = logFreq,
      fill = dataset,
      label = paste0(freq, " (", percentage, "%)")
    )) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    #    coord_flip() +
    labs(x = NULL, y = NULL, color= "dataset") +
    scale_y_continuous(labels = comma) +
    ylab("count and fraction of the dataset \n(axis log-transformed)") +
    xlab("feature change") +
    geom_text(size = 2, position = position_stack(vjust = 0.50), angle = 90) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 10)) +
    ggtitle(paste("Total", level, "MutPred structural and functional feature changes"))
  
  return(MutPred_plot)
}
pdf(paste0('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/multi_combined_bar_small.pdf'), 
    width = 8,
    height = 18)
pvconf = MutPred_plot_features_combined(dataset_Mut_V_conf, "very confident")
pconf = MutPred_plot_features_combined(dataset_Mut_Conf, "confident")
pact = MutPred_plot_features_combined(dataset_Mut_Act, "actionable")
#grid.arrange(pvconf, pconf, pact, nrow=3)
cowplot::plot_grid(pact + theme(axis.text.x = element_blank(),
                                axis.ticks.x = element_blank(),
                                axis.title.x = element_blank()) 
                   + ggtitle(""),
                   pconf + theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.title.x = element_blank())
                   + ggtitle(""),
                   pvconf + ggtitle(""),
                   ncol=1, labels = c("actionable", "confident", "very confident"), align = "v", axis = "lr", rel_heights = c(1,1,1.45))
dev.off()

### Additional MutPred statistics

MutPred_feature_counter = function(feature_list){
  filtered_list = feature_list[-which(sapply(feature_list, is.null))]
  return(sum(lengths(filtered_list)))
}

MutPred_feature_counter(somatic_NSFP_extended$MutPred_actionable)
MutPred_feature_counter(somatic_NSFP_extended$MutPred_confident)
MutPred_feature_counter(somatic_NSFP_extended$MutPred_v_confident)

MutPred_feature_counter(germline_NSFP_extended$MutPred_actionable)
MutPred_feature_counter(germline_NSFP_extended$MutPred_confident)
MutPred_feature_counter(germline_NSFP_extended$MutPred_v_confident)

MutPred_feature_counter(shared_extended$MutPred_actionable)
MutPred_feature_counter(shared_extended$MutPred_confident)
MutPred_feature_counter(shared_extended$MutPred_v_confident)

### Additional MutPred histogram of lengths
colors <- c("somatic" = "red", "germline" = "blue", "shared" = "green")

#Actionable
pdf(paste('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Quantitative Analysis/actionable_bar.pdf'), 
    width = 8,
    height = 4)
Mut_som_act = as.data.frame(lengths(somatic_NSFP_extended$MutPred_actionable))
Mut_germ_act = as.data.frame(lengths(germline_NSFP_extended$MutPred_actionable))
Mut_shared_act = as.data.frame(lengths(shared_extended$MutPred_actionable))
MutPred_actionable_histogram <- ggplot() + 
  geom_histogram(data = Mut_som_act, aes(x = `lengths(somatic_NSFP_extended$MutPred_actionable)`, color= "somatic", fill = "somatic"), binwidth=.1, alpha=.5, position="identity") +
  geom_histogram(data = Mut_germ_act, aes(x = `lengths(germline_NSFP_extended$MutPred_actionable)`, color = "germline", fill = "germline"), binwidth=.1,  alpha=.5, position="identity") +
  geom_histogram(data = Mut_shared_act, aes(x = `lengths(shared_extended$MutPred_actionable)`, color = "shared", fill = "shared"), binwidth=.1,  alpha=.5, position="identity") +
  theme_bw() + xlim(c(0.5,5.5)) + ylim(c(0, 875000)) + labs(x = "actionable hypotheses per mutation", y = "mutations", title = "Distribution of actionable hypotheses across datasets", color = "Dataset", fill = "Dataset")
plot(MutPred_actionable_histogram)
dev.off()

#Confident
pdf(paste('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Quantitative Analysis/confident_bar.pdf'), 
    width = 8,
    height = 4)
Mut_som_conf = as.data.frame(lengths(somatic_NSFP_extended$MutPred_confident))
Mut_germ_conf = as.data.frame(lengths(germline_NSFP_extended$MutPred_confident))
Mut_shared_conf = as.data.frame(lengths(shared_extended$MutPred_confident))
MutPred_confident_histogram <- ggplot() + 
  geom_histogram(data = Mut_som_conf, aes(x = `lengths(somatic_NSFP_extended$MutPred_confident)`, color= "somatic", fill = "somatic"), binwidth=.1, alpha=.5, position="identity") +
  geom_histogram(data = Mut_germ_conf, aes(x = `lengths(germline_NSFP_extended$MutPred_confident)`, color = "germline", fill = "germline"), binwidth=.1,  alpha=.5, position="identity") +
  geom_histogram(data = Mut_shared_conf, aes(x = `lengths(shared_extended$MutPred_confident)`, color = "shared", fill = "shared"), binwidth=.1,  alpha=.5, position="identity") +
  theme_bw() + xlim(c(0.5,5.5)) + ylim(c(0, 250000)) + labs(x = "confident hypotheses per mutation", y = "mutations", title = "Distribution of confident hypotheses across datasets", color = "Dataset", fill = "Dataset")
plot(MutPred_confident_histogram)
dev.off()

#Very confident
pdf(paste('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Quantitative Analysis/v_confident_bar.pdf'), 
    width = 8,
    height = 4)
Mut_som_v_conf = as.data.frame(lengths(somatic_NSFP_extended$MutPred_v_confident))
Mut_germ_v_conf = as.data.frame(lengths(germline_NSFP_extended$MutPred_v_confident))
Mut_shared_v_conf = as.data.frame(lengths(shared_extended$MutPred_v_confident))
MutPred_v_conf_histogram <- ggplot() + 
  geom_histogram(data = Mut_som_v_conf, aes(x = `lengths(somatic_NSFP_extended$MutPred_v_confident)`, color= "somatic", fill = "somatic"), binwidth=.1, alpha=.5, position="identity") +
  geom_histogram(data = Mut_germ_v_conf, aes(x = `lengths(germline_NSFP_extended$MutPred_v_confident)`, color = "germline", fill = "germline"), binwidth=.1,  alpha=.5, position="identity") +
  geom_histogram(data = Mut_shared_v_conf, aes(x = `lengths(shared_extended$MutPred_v_confident)`, color = "shared", fill = "shared"), binwidth=.1,  alpha=.5, position="identity") +
  theme_bw() + xlim(c(0.5,5.5)) + ylim(c(0, 100000)) + labs(x = "very confident hypotheses per mutation", y = "mutations", title = "Distribution of very confident hypotheses across datasets", color = "Dataset", fill = "Dataset")
plot(MutPred_v_conf_histogram)
dev.off()


#Actionable features

somatic_Mut_Act$'dataset' = "somatic"
germline_Mut_Act$'dataset' = "germline"
shared_Mut_Act$'dataset' = "shared"
dataset_Mut_Act = rbind(somatic_Mut_Act, germline_Mut_Act, shared_Mut_Act)
dataset_Mut_Act = subset(dataset_Mut_Act, !is.na(MutPred_list))
dataset_Mut_Act = as.data.frame(group_by(dataset_Mut_Act, dataset) %>% mutate(percent = Freq/sum(Freq)))


#Confident features

somatic_Mut_Conf$'dataset' = "somatic"
germline_Mut_Conf$'dataset' = "germline"
shared_Mut_Conf$'dataset' = "shared"
dataset_Mut_Conf = rbind(somatic_Mut_Conf, germline_Mut_Conf, shared_Mut_Conf)
dataset_Mut_Conf = subset(dataset_Mut_Conf, !is.na(MutPred_list))
dataset_Mut_Conf = as.data.frame(group_by(dataset_Mut_Conf, dataset) %>% mutate(percent = Freq/sum(Freq)))

#Very confident features

somatic_Mut_V_conf$'dataset' = "somatic"
germline_Mut_V_conf$'dataset' = "germline"
shared_Mut_V_conf$'dataset' = "shared"
dataset_Mut_V_conf = rbind(somatic_Mut_V_conf, germline_Mut_V_conf, shared_Mut_V_conf)
dataset_Mut_V_conf = subset(dataset_Mut_V_conf, !is.na(MutPred_list))
dataset_Mut_V_conf = as.data.frame(group_by(dataset_Mut_V_conf, dataset) %>% mutate(percent = Freq/sum(Freq)))
#set1 = c("Gain of catalytic residue", "Gain of disorder", "Gain of glycosylation", "Gain of methylation", "Gain of MoRF binding", "Gain of solvent accessibility",  "Loss of catalytic residue", "Loss of disorder", "Loss of helix", "Loss of methylation", "Loss of MoRF binding", "Loss of sheet", "Loss of stability")
set1_1 = c("Gain of catalytic residue", "Gain of disorder", "Gain of methylation", "Gain of MoRF binding",  "Loss of catalytic residue", "Loss of MoRF binding", "Loss of stability")
set2_1 = c("Gain of glycosylation", "Gain of solvent accessibility", "Loss of disorder",  "Loss of helix", "Loss of methylation", "Loss of sheet", "Gain of loop", "Gain of relative solvent accessibility", "Gain of sheet")
set3_1 = c("Gain of helix", "Gain of phosphorylation", "Gain of stability", "Gain of ubiquitination", "Loss of glycosylation",  "Loss of loop", "Loss of phosphorylation", "Loss of relative solvent accessibility", "Loss of solvent accessibility", "Loss of ubiquitination")
dataset_Mut_V_conf1 = dataset_Mut_V_conf[dataset_Mut_V_conf$MutPred_list %in% set1_1,]
dataset_Mut_V_conf2 = dataset_Mut_V_conf[dataset_Mut_V_conf$MutPred_list %in% set2_1,]
dataset_Mut_V_conf3 = dataset_Mut_V_conf[dataset_Mut_V_conf$MutPred_list %in% set3_1,]

pdf(paste('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Quantitative Analysis/v_confident_features_bar_high.pdf'), 
    width = 6,
    height = 12)

MutPred_v_conf_feat_histogram1 <-
  ggplot(dataset_Mut_V_conf1, aes(x = MutPred_list, y = Freq, color= dataset, fill = dataset, label = Freq)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  #  coord_flip() +
  labs(x = "very confident hypotheses", y = "mutations", title = "Distribution of very confident hypotheses across datasets (1/3)", color = "Dataset", fill = "Dataset") +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = 'black')

plot(MutPred_v_conf_feat_histogram1)
dev.off()

pdf(paste('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Quantitative Analysis/v_confident_features_bar_mid.pdf'), 
    width = 6,
    height = 12)

MutPred_v_conf_feat_histogram2 <-
  ggplot(dataset_Mut_V_conf2, aes(x = MutPred_list, y = Freq, color= dataset, fill = dataset, label = Freq)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  #  coord_flip() +
  labs(x = "very confident hypotheses", y = "mutations", title = "Distribution of very confident hypotheses across datasets (2/3)", color = "Dataset", fill = "Dataset") +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = 'black')

plot(MutPred_v_conf_feat_histogram2)
dev.off()

pdf(paste('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Quantitative Analysis/v_confident_features_bar_low.pdf'), 
    width = 6,
    height = 12)

MutPred_v_conf_feat_histogram3 <-
  ggplot(dataset_Mut_V_conf3, aes(x = MutPred_list, y = Freq, color= dataset, fill = dataset, label = Freq)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  #  coord_flip() +
  labs(x = "very confident hypotheses", y = "mutations", title = "Distribution of very confident hypotheses across datasets (3/3)", color = "Dataset", fill = "Dataset") +
  scale_y_continuous(labels = comma) +
  ylab("") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = 'black')

plot(MutPred_v_conf_feat_histogram3)
dev.off()

###MutPred feature heatmap


x <- c("Catalytic residue",
       "Disorder",
       "Glycosylation",
       "Helix",
       "Loop",
       "Methylation",
       "MoRF binding",
       "Phosphorylation",
       "Relative solvent accessibility",
       "Sheet",
       "Solvent accessibility",
       "Stability",
       "Ubiquitination")
y <- c("Gain in somatic (%)", "Loss in somatic (%)", "Gain in germline (%)", "Loss in germline (%)", "Gain in shared (%)", "Loss in shared (%)")
data <- expand.grid(X=x, Y=y)
data$Z <- c(c(14.17, 23.89, 2.34, 0.73, 1.55, 6.87, 8.63, 0.63, 1.70, 1.58, 2.56, 0.08, 0.65),
            c(6.59, 3.59, 0.52, 2.16, 0.32, 2.43, 5.48, 0.05, 0.36, 2.41, 0.16, 10.48, 0.09),
            c(14.15, 24.15, 2.55, 0.73, 1.76, 9.43, 7.66, 0.42, 1.35, 1.48, 2.10, 0.14, 0.81),
            c(6.40, 2.98, 0.53, 2.20, 0.26, 2.06, 5.98, 0.06, 0.20, 2.25, 0.14, 10.11, 0.09),
            c(13.76, 11.87, 2.72, 0.84, 1.78, 7.32, 9.23, 0.53, 1.88, 1.61, 3.28, 0.07, 0.73),
            c(5.68, 5.04, 0.59, 2.61, 0.28, 2.30, 19.52, 0.04, 0.22, 2.16, 0.35, 5.49, 0.07))

data <- data %>%
  mutate(text = paste("feature: ", X, "\n", "change: ", Y, "\n", "Percentage of feature changes in dataset: ",round(Z,2), "% \n"))


MutHeat <- ggplot(data, aes(X, Y, fill= Z, text=text)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  geom_text(aes(label = round(Z, 2))) +
  labs(x = "structural feature", y = "change", title = "Distribution of very confident structual changes", fill = "Percentage") +
  scale_fill_gradient(low="white", high="red") 
MutHeat
MutHeat_ia <- ggplotly(MutHeat, tooltip="text")
MutHeat_ia

saveWidget(MutHeat_ia, file="/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred Enrichment/heatmap1.html")


###MutPred feature heatmap with fold-changes


x <- c("Catalytic residue",
       "Disorder",
       "Glycosylation",
       "Helix",
       "Loop",
       "Methylation",
       "MoRF binding",
       "Phosphorylation",
       "Relative solvent accessibility",
       "Sheet",
       "Solvent accessibility",
       "Stability",
       "Ubiquitination")
y <- c("Gain in somatic (%)", "Loss in somatic (%)", "Fold-change in somatic", "Gain in germline (%)", "Loss in germline (%)", "Fold-change in germline", "Gain in shared (%)", "Loss in shared (%)", "Fold-change in shared")
data <- expand.grid(X=x, Y=y)
data$Z <- c(c(14.17, 23.89, 2.34, 0.73, 1.55, 6.87, 8.63, 0.63, 1.70, 1.58, 2.56, 0.08, 0.65),
            c(6.59, 3.59, 0.52, 2.16, 0.32, 2.43, 5.48, 0.05, 0.36, 2.41, 0.16, 10.48, 0.09),
            c(2.15, 6.65, 4.47, -2.97, 4.88, 2.83, 1.58, 11.97, 4.70, -1.53, 16.17, -135.55, 7.29),
            c(14.15, 24.15, 2.55, 0.73, 1.76, 9.43, 7.66, 0.42, 1.35, 1.48, 2.10, 0.14, 0.81),
            c(6.40, 2.98, 0.53, 2.20, 0.26, 2.06, 5.98, 0.06, 0.20, 2.25, 0.14, 10.11, 0.09),
            c(2.21, 8.12, 4.84, -3.00, 6.71, 4.57, 1.28, 6.80, 6.73, -1.52, 15.30, -70.39, 9.25),
            c(13.76, 11.87, 2.72, 0.84, 1.78, 7.32, 9.23, 0.53, 1.88, 1.61, 3.28, 0.07, 0.73),
            c(5.68, 5.04, 0.59, 2.61, 0.28, 2.30, 19.52, 0.04, 0.22, 2.16, 0.35, 5.49, 0.07),
            c(2.42, 2.36, 4.62, -3.10, 6.35, 3.18, -2.11, 12.67, 8.38, -1.34, 9.36, -78.20, 10.40))

data <- data %>%
  mutate(text = paste0("feature: ", X, "\n", "change: ", Y, "\n", "value: ",round(Z,2), "\n"))

rng = range(data$Z)


pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/heatmap_static_percentage_fc.pdf', 
    width = 10,
    height = 8)
MutHeat2 <- ggplot(data, aes(X, Y, fill= Z, text=text)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  geom_text(aes(label = round(Z, 2))) +
  labs(x = "structural feature", y = "change", title = "Distribution of very confident structual feature changes", fill = "value") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(floor(rng[1]), ceiling(rng[2])))
#  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=2, limits=c(-130, 20))
#scale_fill_gradient2(low="white", high="red", midpoint = 2)
MutHeat2
dev.off()

MutHeat_ia2 <- ggplotly(MutHeat2, tooltip="text")
MutHeat_ia2

saveWidget(MutHeat_ia2, file="/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred Enrichment/heatmap_perc_fc.html")

###MutPred statistical significance analysis

regexp <- "([:space:])([[:digit:]]+)(\\.*)([[:digit:]]*)"
MutPred_features_pvals <- function (MutPred_list) {
  MutPred_list = unlist(MutPred_list)
  MutPred_terms = sapply(str_split(MutPred_list, "\\("), head, 1)
  MutPred_terms = sapply(str_split(MutPred_terms, " at "), head, 1)
  MutPred_terms = trimws(MutPred_terms)
  MutPred_pvals = lapply(MutPred_list, lapply, str_extract, regexp)
  MutPred_pvals = unlist(MutPred_pvals)
  MutPred_pvals = trimws(MutPred_pvals)
  MutPred_pvals = as.numeric(MutPred_pvals)
  MutPred_processed = data.frame(MutPred_terms, MutPred_pvals)
  MutPred_processed = subset(MutPred_processed, !is.na(MutPred_terms))
  return(MutPred_processed)
}

somatic_MutPred_terms_pvals_v_conf = MutPred_features_pvals(somatic_NSFP_extended$MutPred_v_confident)
germline_MutPred_terms_pvals_v_conf = MutPred_features_pvals(somatic_NSFP_extended$MutPred_v_confident)
shared_MutPred_terms_pvals_v_conf = MutPred_features_pvals(shared_extended$MutPred_v_confident)

MutPred_features_scores_pvals <- function (dataset) {
  output = data.frame(MutPred_term = character(), MutPred_score = double(), MutPred_pval = double())
  for (i in 1:nrow(dataset)) {
    if(length(dataset$MutPred_v_confident[[i]]) != 0){
      MutPred_score = dataset$MutPred_score[i]
      MutPred_terms = dataset$MutPred_v_confident[[i]]
      for (j in 1:length(MutPred_terms)) {
        MutPred_term = sapply(str_split(MutPred_terms[j], "\\("), head, 1)
        MutPred_term = sapply(str_split(MutPred_term, " at "), head, 1)
        MutPred_term = trimws(MutPred_term)
        MutPred_pval = dataset$MutPred_pvals[[i]][j]
        output = rbind(output, c(MutPred_term, MutPred_score, MutPred_pval))
      }
    }
  }
  output = subset(output, !is.na(output[,1]))
  output$MutPred_score = as.numeric(output$MutPred_score)
  output$MutPred_pval = as.numeric(output$MutPred_pval)
  colnames(output) = c("MutPred_term", "MutPred_score", "MutPred_pval")
  return(output)
}

somatic_MutPred_terms_scores_pvals_v_conf = MutPred_features_scores_pvals(somatic_NSFP_extended)
germline_MutPred_terms_scores_pvals_v_conf = MutPred_features_scores_pvals(germline_NSFP_extended)
shared_MutPred_terms_scores_pvals_v_conf = MutPred_features_scores_pvals(shared_extended)


###Find best fit for data

#Describe distributions
library(fitdistrplus)
library(logspline)
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/cullen_frey_somatic_scores.pdf', 
    width = 8,
    height = 6)
descdist(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_score, discrete = FALSE)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/cullen_frey_somatic_pvals.pdf', 
    width = 8,
    height = 6)
descdist(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_pval, discrete = FALSE)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/cullen_frey_germline_scores.pdf', 
    width = 8,
    height = 6)
descdist(germline_MutPred_terms_scores_pvals_v_conf$MutPred_score, discrete = FALSE)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/cullen_frey_germline_pvals.pdf', 
    width = 8,
    height = 6)
descdist(germline_MutPred_terms_scores_pvals_v_conf$MutPred_pval, discrete = FALSE)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/cullen_frey_shared_scores.pdf', 
    width = 8,
    height = 6)
descdist(shared_MutPred_terms_scores_pvals_v_conf$MutPred_score, discrete = FALSE)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/cullen_frey_shared_pvals.pdf', 
    width = 8,
    height = 6)
descdist(shared_MutPred_terms_scores_pvals_v_conf$MutPred_pval, discrete = FALSE)
dev.off()
#Fit to normal distribution
fit.norm_somatic_score <- fitdist(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_score, "norm")
fit.norm_somatic_pval <- fitdist(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_pval, "norm")
fit.norm_germline_score <- fitdist(germline_MutPred_terms_scores_pvals_v_conf$MutPred_score, "norm")
fit.norm_germline_pval <- fitdist(germline_MutPred_terms_scores_pvals_v_conf$MutPred_pval, "norm")
fit.norm_shared_score <- fitdist(shared_MutPred_terms_scores_pvals_v_conf$MutPred_score, "norm")
fit.norm_shared_pval <- fitdist(shared_MutPred_terms_scores_pvals_v_conf$MutPred_pval, "norm")

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_norm_somatic_scores.pdf', 
    width = 8,
    height = 6)
plot(fit.norm_somatic_score)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_norm_somatic_pvals.pdf', 
    width = 8,
    height = 6)
plot(fit.norm_somatic_pval)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_norm_germline_scores.pdf', 
    width = 8,
    height = 6)
plot(fit.norm_germline_score)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_norm_germline_pvals.pdf', 
    width = 8,
    height = 6)
plot(fit.norm_germline_pval)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_norm_shared_scores.pdf', 
    width = 8,
    height = 6)
plot(fit.norm_shared_score)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_norm_shared_pvals.pdf', 
    width = 8,
    height = 6)
plot(fit.norm_shared_pval)
dev.off()
#Fit to uniform distribution
fit.unif_somatic_score <- fitdist(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_score, "unif")
fit.unif_somatic_pval <- fitdist(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_pval, "unif")
fit.unif_germline_score <- fitdist(germline_MutPred_terms_scores_pvals_v_conf$MutPred_score, "unif")
fit.unif_germline_pval <- fitdist(germline_MutPred_terms_scores_pvals_v_conf$MutPred_pval, "unif")
fit.unif_shared_score <- fitdist(shared_MutPred_terms_scores_pvals_v_conf$MutPred_score, "unif")
fit.unif_shared_pval <- fitdist(shared_MutPred_terms_scores_pvals_v_conf$MutPred_pval, "unif")

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_unif_somatic_scores.pdf', 
    width = 8,
    height = 6)
plot(fit.unif_somatic_score)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_unif_somatic_pvals.pdf', 
    width = 8,
    height = 6)
plot(fit.unif_somatic_pval)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_unif_germline_scores.pdf', 
    width = 8,
    height = 6)
plot(fit.unif_germline_score)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_unif_germline_pvals.pdf', 
    width = 8,
    height = 6)
plot(fit.unif_germline_pval)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_unif_shared_scores.pdf', 
    width = 8,
    height = 6)
plot(fit.unif_shared_score)
dev.off()
pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/fit_unif_shared_pvals.pdf', 
    width = 8,
    height = 6)
plot(fit.unif_shared_pval)
dev.off()

###Perform ori tests
feature_set = c(set1_1, set2_1, set3_1)
feature_set = feature_set[order(feature_set)]
t_tests = data.frame(feature = character(), pval_score_1 = double(), pval_score_2 = double(), pval_score_3 = double(),
                     pval_pv_1 = double(), pval_pv_2 = double(), pval_pv_3 = double())
ks_tests = data.frame(feature = character(), pval_score_1 = double(), pval_score_2 = double(), pval_score_3 = double(),
                      pval_pv_1 = double(), pval_pv_2 = double(), pval_pv_3 = double())
wilcox_tests = data.frame(feature = character(), pval_score_1 = double(), pval_score_2 = double(), pval_score_3 = double(),
                          pval_pv_1 = double(), pval_pv_2 = double(), pval_pv_3 = double())
for (feature in feature_set) {
  #t-test
  t_test_score_1 = t.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_score[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                          germline_MutPred_terms_scores_pvals_v_conf$MutPred_score[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                          var.equal = T, alternative = "two.sided")
  t_test_score_2 = t.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_score[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                          shared_MutPred_terms_scores_pvals_v_conf$MutPred_score[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                          var.equal = T, alternative = "two.sided")
  t_test_score_3 = t.test(shared_MutPred_terms_scores_pvals_v_conf$MutPred_score[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                          germline_MutPred_terms_scores_pvals_v_conf$MutPred_score[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                          var.equal = T, alternative = "two.sided")
  t_test_pval_1 = t.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_pval[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                         germline_MutPred_terms_scores_pvals_v_conf$MutPred_pval[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                         var.equal = T, alternative = "two.sided")
  t_test_pval_2 = t.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_pval[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                         shared_MutPred_terms_scores_pvals_v_conf$MutPred_pval[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                         var.equal = T, alternative = "two.sided")
  t_test_pval_3 = t.test(shared_MutPred_terms_scores_pvals_v_conf$MutPred_pval[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                         germline_MutPred_terms_scores_pvals_v_conf$MutPred_pval[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                         var.equal = T, alternative = "two.sided")
  t_tests = rbind(t_tests, c(feature, t_test_score_1$p.value, t_test_score_2$p.value, t_test_score_3$p.value, t_test_pval_1$p.value, t_test_pval_2$p.value, t_test_pval_3$p.value))
  colnames(t_tests) = c("feature", "somatic vs. germline, based on score",  "somatic vs. shared, based on score",  "shared vs. germline, based on score", "somatic vs. germline, based on p-value",  "somatic vs. shared, based on p-value",  "shared vs. germline, based on p-value")
  #ks-test
  ks_test_score_1 = ks.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_score[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                            germline_MutPred_terms_scores_pvals_v_conf$MutPred_score[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                            alternative = "two.sided")
  ks_test_score_2 = ks.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_score[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                            shared_MutPred_terms_scores_pvals_v_conf$MutPred_score[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                            alternative = "two.sided")
  ks_test_score_3 = ks.test(shared_MutPred_terms_scores_pvals_v_conf$MutPred_score[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                            germline_MutPred_terms_scores_pvals_v_conf$MutPred_score[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                            alternative = "two.sided")
  ks_test_pval_1 = ks.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_pval[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                           germline_MutPred_terms_scores_pvals_v_conf$MutPred_pval[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                           alternative = "two.sided")
  ks_test_pval_2 = ks.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_pval[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                           shared_MutPred_terms_scores_pvals_v_conf$MutPred_pval[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                           alternative = "two.sided")
  ks_test_pval_3 = ks.test(shared_MutPred_terms_scores_pvals_v_conf$MutPred_pval[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                           germline_MutPred_terms_scores_pvals_v_conf$MutPred_pval[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                           alternative = "two.sided")
  ks_tests = rbind(ks_tests, c(feature, ks_test_score_1$p.value, ks_test_score_2$p.value, ks_test_score_3$p.value, ks_test_pval_1$p.value, ks_test_pval_2$p.value, ks_test_pval_3$p.value))
  colnames(ks_tests) = c("feature", "somatic vs. germline, based on score",  "somatic vs. shared, based on score",  "shared vs. germline, based on score", "somatic vs. germline, based on p-value",  "somatic vs. shared, based on p-value",  "shared vs. germline, based on p-value")
  #wilcox-test
  w_test_score_1 = wilcox.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_score[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                               germline_MutPred_terms_scores_pvals_v_conf$MutPred_score[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                               alternative = "two.sided")
  w_test_score_2 = wilcox.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_score[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                               shared_MutPred_terms_scores_pvals_v_conf$MutPred_score[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                               alternative = "two.sided")
  w_test_score_3 = wilcox.test(shared_MutPred_terms_scores_pvals_v_conf$MutPred_score[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                               germline_MutPred_terms_scores_pvals_v_conf$MutPred_score[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                               alternative = "two.sided")
  w_test_pval_1 = wilcox.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_pval[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                              germline_MutPred_terms_scores_pvals_v_conf$MutPred_pval[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                              alternative = "two.sided")
  w_test_pval_2 = wilcox.test(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_pval[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                              shared_MutPred_terms_scores_pvals_v_conf$MutPred_pval[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                              alternative = "two.sided")
  w_test_pval_3 = wilcox.test(shared_MutPred_terms_scores_pvals_v_conf$MutPred_pval[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                              germline_MutPred_terms_scores_pvals_v_conf$MutPred_pval[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature], 
                              alternative = "two.sided")
  wilcox_tests = rbind(wilcox_tests, c(feature, w_test_score_1$p.value, w_test_score_2$p.value, w_test_score_3$p.value, w_test_pval_1$p.value, w_test_pval_2$p.value, w_test_pval_3$p.value))
  colnames(wilcox_tests) = c("feature", "somatic vs. germline, based on score",  "somatic vs. shared, based on score",  "shared vs. germline, based on score", "somatic vs. germline, based on p-value",  "somatic vs. shared, based on p-value",  "shared vs. germline, based on p-value")
  
}

###Generate Heatmaps

#T-tests
x <- t_tests$feature
y <- colnames(t_tests)[2:7]
data <- expand.grid(X=x, Y=y)
data$Z <- as.numeric(unlist(t_tests[,-1]))
data <- data %>%
  mutate(text = paste("feature: ", X, "\n", "test setting: ", Y, "\n", "p-value: ", Z, "\n"))
rng = c(min(data$Z), max(data$Z))

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/heatmap_t_tests.pdf', 
    width = 16,
    height = 8)
MutHeat_t_tests <- ggplot(data, aes(X, Y, fill= Z, text=text)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 12),axis.text.y = element_text(angle = 60, hjust=1, size = 12)) +
  geom_text(aes(label = round(Z, 6), angle = 60)) +
  labs(x = "feature change", y = "test setting", title = "Distribution of t-test p-values on feature changes", fill = "p-value") +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0.1,    
                       limits=c(0, 1))
MutHeat_t_tests
dev.off()

MutHeat_t_tests_ia <- ggplotly(MutHeat_t_tests, tooltip="text")
MutHeat_t_tests_ia

saveWidget(MutHeat_t_tests_ia, file="/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred Enrichment/heatmap_pvals_t_tests.html")

#KS-tests
x <- ks_tests$feature
y <- colnames(ks_tests)[2:7]
data <- expand.grid(X=x, Y=y)
data$Z <- as.numeric(unlist(ks_tests[,-1]))
data <- data %>%
  mutate(text = paste("feature: ", X, "\n", "test setting: ", Y, "\n", "p-value: ", Z, "\n"))
rng = c(min(data$Z), max(data$Z))

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/heatmap_ks_tests.pdf', 
    width = 16,
    height = 8)
MutHeat_ks_tests <- ggplot(data, aes(X, Y, fill= Z, text=text)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 12),axis.text.y = element_text(angle = 60, hjust=1, size = 12)) +
  geom_text(aes(label = round(Z, 6), angle = 60)) +
  labs(x = "feature change", y = "test setting", title = "Distribution of KS-test p-values on feature changes", fill = "p-value") +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0.1,    
                       limits=c(0, 1))
MutHeat_ks_tests
dev.off()

MutHeat_ks_tests_ia <- ggplotly(MutHeat_ks_tests, tooltip="text")
MutHeat_ks_tests_ia

saveWidget(MutHeat_ks_tests_ia, file="/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred Enrichment/heatmap_pvals_ks_tests.html")

#Wilcoxon-tests
x <- wilcox_tests$feature
y <- colnames(wilcox_tests)[2:7]
data <- expand.grid(X=x, Y=y)
data$Z <- as.numeric(unlist(wilcox_tests[,-1]))
data <- data %>%
  mutate(text = paste("feature: ", X, "\n", "test setting: ", Y, "\n", "p-value: ", Z, "\n"))
rng = c(min(data$Z), max(data$Z))

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/heatmap_wilcox_tests.pdf', 
    width = 16,
    height = 8)
MutHeat_wilcox_tests <- ggplot(data, aes(X, Y, fill= Z, text=text)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 12),axis.text.y = element_text(angle = 60, hjust=1, size = 12)) +
  geom_text(aes(label = round(Z, 6), angle = 60)) +
  labs(x = "feature change", y = "test setting", title = "Distribution of Wilcoxon test p-values on feature changes", fill = "p-value") +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0.1,    
                       limits=c(0, 1))
MutHeat_wilcox_tests
dev.off()

MutHeat_wilcox_tests_ia <- ggplotly(MutHeat_wilcox_tests, tooltip="text")
MutHeat_wilcox_tests_ia

saveWidget(MutHeat_wilcox_tests_ia, file="/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred Enrichment/heatmap_pvals_wilcox_tests.html")

### Multiple testing with 1000 reps ###
t_tests_list = list()
ks_tests_list = list()
wilcox_tests_list = list()

for (i in 1:1000) {
  print(i)
  t_tests = data.frame(feature = character(), pval_score_1 = double(), pval_score_2 = double(), pval_score_3 = double(),
                       pval_pv_1 = double(), pval_pv_2 = double(), pval_pv_3 = double())
  ks_tests = data.frame(feature = character(), pval_score_1 = double(), pval_score_2 = double(), pval_score_3 = double(),
                        pval_pv_1 = double(), pval_pv_2 = double(), pval_pv_3 = double())
  wilcox_tests = data.frame(feature = character(), pval_score_1 = double(), pval_score_2 = double(), pval_score_3 = double(),
                            pval_pv_1 = double(), pval_pv_2 = double(), pval_pv_3 = double())
  for (feature in feature_set) {
    somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term <- sample(somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term)
    germline_MutPred_terms_scores_pvals_v_conf$MutPred_term <- sample(germline_MutPred_terms_scores_pvals_v_conf$MutPred_term)
    shared_MutPred_terms_scores_pvals_v_conf$MutPred_term <- sample(shared_MutPred_terms_scores_pvals_v_conf$MutPred_term)
    som_scores = somatic_MutPred_terms_scores_pvals_v_conf$MutPred_score[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature] 
    germ_scores = germline_MutPred_terms_scores_pvals_v_conf$MutPred_score[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature]
    sha_scores = shared_MutPred_terms_scores_pvals_v_conf$MutPred_score[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature]
    som_pvals = somatic_MutPred_terms_scores_pvals_v_conf$MutPred_pval[somatic_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature]
    germ_pvals = germline_MutPred_terms_scores_pvals_v_conf$MutPred_pval[germline_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature]
    sha_pvals = shared_MutPred_terms_scores_pvals_v_conf$MutPred_pval[shared_MutPred_terms_scores_pvals_v_conf$MutPred_term == feature]
    #t-test
    t_test_score_1 = t.test(som_scores,
                            germ_scores,
                            var.equal = T, alternative = "two.sided")
    t_test_score_2 = t.test(som_scores,
                            sha_scores,
                            var.equal = T, alternative = "two.sided")
    t_test_score_3 = t.test(sha_scores,
                            germ_scores,
                            var.equal = T, alternative = "two.sided")
    t_test_pval_1 = t.test(som_pvals,
                           germ_pvals,
                           var.equal = T, alternative = "two.sided")
    t_test_pval_2 = t.test(som_pvals,
                           sha_pvals,
                           var.equal = T, alternative = "two.sided")
    t_test_pval_3 = t.test(sha_pvals,
                           germ_pvals,
                           var.equal = T, alternative = "two.sided")
    t_tests = rbind(t_tests, c(feature, t_test_score_1$p.value, t_test_score_2$p.value, t_test_score_3$p.value, t_test_pval_1$p.value, t_test_pval_2$p.value, t_test_pval_3$p.value))
    colnames(t_tests) = c("feature", "somatic vs. germline, based on score",  "somatic vs. shared, based on score",  "shared vs. germline, based on score", "somatic vs. germline, based on p-value",  "somatic vs. shared, based on p-value",  "shared vs. germline, based on p-value")
    #ks-test
    ks_test_score_1 = ks.test(som_scores,
                              germ_scores,
                              alternative = "two.sided")
    ks_test_score_2 = ks.test(som_scores,
                              sha_scores,
                              alternative = "two.sided")
    ks_test_score_3 = ks.test(sha_scores,
                              germ_scores,
                              alternative = "two.sided")
    ks_test_pval_1 = ks.test(som_pvals,
                             germ_pvals,
                             alternative = "two.sided")
    ks_test_pval_2 = ks.test(som_pvals,
                             sha_pvals,
                             alternative = "two.sided")
    ks_test_pval_3 = ks.test(sha_pvals,
                             germ_pvals,
                             alternative = "two.sided")
    ks_tests = rbind(ks_tests, c(feature, ks_test_score_1$p.value, ks_test_score_2$p.value, ks_test_score_3$p.value, ks_test_pval_1$p.value, ks_test_pval_2$p.value, ks_test_pval_3$p.value))
    colnames(ks_tests) = c("feature", "somatic vs. germline, based on score",  "somatic vs. shared, based on score",  "shared vs. germline, based on score", "somatic vs. germline, based on p-value",  "somatic vs. shared, based on p-value",  "shared vs. germline, based on p-value")
    #wilcox-test
    w_test_score_1 = wilcox.test(som_scores,
                                 germ_scores,
                                 alternative = "two.sided")
    w_test_score_2 = wilcox.test(som_scores,
                                 sha_scores,
                                 alternative = "two.sided")
    w_test_score_3 = wilcox.test(sha_scores,
                                 germ_scores,
                                 alternative = "two.sided")
    w_test_pval_1 = wilcox.test(som_pvals,
                                germ_pvals,
                                alternative = "two.sided")
    w_test_pval_2 = wilcox.test(som_pvals,
                                sha_pvals,
                                alternative = "two.sided")
    w_test_pval_3 = wilcox.test(sha_pvals,
                                germ_pvals,
                                alternative = "two.sided")
    wilcox_tests = rbind(wilcox_tests, c(feature, w_test_score_1$p.value, w_test_score_2$p.value, w_test_score_3$p.value, w_test_pval_1$p.value, w_test_pval_2$p.value, w_test_pval_3$p.value))
    colnames(wilcox_tests) = c("feature", "somatic vs. germline, based on score",  "somatic vs. shared, based on score",  "shared vs. germline, based on score", "somatic vs. germline, based on p-value",  "somatic vs. shared, based on p-value",  "shared vs. germline, based on p-value")
  }
  t_tests_list[[i]] <- t_tests
  ks_tests_list[[i]] <- ks_tests
  wilcox_tests_list[[i]] <- wilcox_tests
}  

### Merging
t_tests_merged = wilcox_tests
ks_tests_merged = wilcox_tests
wilcox_tests_merged = t_tests

for (i in 2:ncol(t_tests)) {
  for (j in 1:nrow(t_tests)) {
    t_tests_merged[[i]][j] <- mean(sapply(t_tests_list, function(x){
      return(as.numeric(x[[i]][j]))
    }))
    
  }
}

for (i in 2:ncol(ks_tests)) {
  for (j in 1:nrow(ks_tests)) {
    ks_tests_merged[[i]][j] <- mean(sapply(ks_tests_list, function(x){
      return(as.numeric(x[[i]][j]))
    }))
    
  }
}

for (i in 2:ncol(wilcox_tests)) {
  for (j in 1:nrow(wilcox_tests)) {
    wilcox_tests_merged[[i]][j] <- mean(sapply(wilcox_tests_list, function(x){
      return(as.numeric(x[[i]][j]))
    }))
  }
}

###Generate merged Heatmaps

#T-tests
x <- t_tests_merged$feature
y <- colnames(t_tests_merged)[2:7]
data <- expand.grid(X=x, Y=y)
data$Z <- as.numeric(unlist(t_tests_merged[,-1]))
data <- data %>%
  mutate(text = paste("feature: ", X, "\n", "test setting: ", Y, "\n", "p-value: ", Z, "\n"))
rng = c(min(data$Z), max(data$Z))

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/heatmap_t_tests_merged.pdf', 
    width = 16,
    height = 8)
MutHeat_t_tests_merged <- ggplot(data, aes(X, Y, fill= Z, text=text)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 12),axis.text.y = element_text(angle = 60, hjust=1, size = 12)) +
  geom_text(aes(label = round(Z, 6), angle = 60)) +
  labs(x = "feature change", y = "test setting", title = "Distribution of t-test p-values on feature changes with random resampling (1000 runs, averaged)", fill = "p-value") +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0.1,    
                       limits=c(0, 1))
MutHeat_t_tests_merged
dev.off()

MutHeat_t_tests_ia_merged <- ggplotly(MutHeat_t_tests_merged, tooltip="text")
MutHeat_t_tests_ia_merged

saveWidget(MutHeat_t_tests_ia_merged, file="/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred Enrichment/heatmap_pvals_t_tests_merged.html")

#KS-tests
x <- ks_tests_merged$feature
y <- colnames(ks_tests_merged)[2:7]
data <- expand.grid(X=x, Y=y)
data$Z <- as.numeric(unlist(ks_tests_merged[,-1]))
data <- data %>%
  mutate(text = paste("feature: ", X, "\n", "test setting: ", Y, "\n", "p-value: ", Z, "\n"))
rng = c(min(data$Z), max(data$Z))

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/heatmap_ks_tests_merged.pdf', 
    width = 16,
    height = 8)
MutHeat_ks_tests_merged <- ggplot(data, aes(X, Y, fill= Z, text=text)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 12),axis.text.y = element_text(angle = 60, hjust=1, size = 12)) +
  geom_text(aes(label = round(Z, 6), angle = 60)) +
  labs(x = "feature change", y = "test setting", title = "Distribution of KS-test p-values on feature changes with random resampling (1000 runs, averaged)", fill = "p-value") +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0.1,    
                       limits=c(0, 1))
MutHeat_ks_tests_merged
dev.off()

MutHeat_ks_tests_ia_merged <- ggplotly(MutHeat_ks_tests_merged, tooltip="text")
MutHeat_ks_tests_ia_merged

saveWidget(MutHeat_ks_tests_ia_merged, file="/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred Enrichment/heatmap_pvals_ks_tests_merged.html")

#Wilcoxon-tests
x <- wilcox_tests_merged$feature
y <- colnames(wilcox_tests_merged)[2:7]
data <- expand.grid(X=x, Y=y)
data$Z <- as.numeric(unlist(wilcox_tests_merged[,-1]))
data <- data %>%
  mutate(text = paste("feature: ", X, "\n", "test setting: ", Y, "\n", "p-value: ", Z, "\n"))
rng = c(min(data$Z), max(data$Z))

pdf('/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred enrichment/Statistical Significance/heatmap_wilcox_tests_merged.pdf', 
    width = 16,
    height = 8)
MutHeat_wilcox_tests_merged <- ggplot(data, aes(X, Y, fill= Z, text=text)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 12),axis.text.y = element_text(angle = 60, hjust=1, size = 12)) +
  geom_text(aes(label = round(Z, 6), angle = 60)) +
  labs(x = "feature change", y = "test setting", title = "Distribution of Wilcoxon test p-values on feature changes with random resampling (1000 runs, averaged)", fill = "p-value") +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0.1,    
                       limits=c(0, 1))
MutHeat_wilcox_tests_merged
dev.off()

MutHeat_wilcox_tests_ia_merged <- ggplotly(MutHeat_wilcox_tests_merged, tooltip="text")
MutHeat_wilcox_tests_ia_merged

saveWidget(MutHeat_wilcox_tests_ia_merged, file="/Users/mcpftw/Documents/Master Bioinformatik/Master Thesis/Data/Plots/MutPred Enrichment/heatmap_pvals_wilcox_tests_merged.html")

