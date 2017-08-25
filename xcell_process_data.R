my_filepath <- "/Users/raimofranke/Downloads/xcell2017-2/data/"
#my_filepath <- "E:/My Dropbox/HZI/xcell2017/data/"

#xcelldata_1
xcell_raw <- read_xcell("1", my_filepath)
xcell_raw_edited <- edit_df("1", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("1", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)


xcell_median_1 <- calculate_median_curves("1", xcell_norm)
xcell_median_norm_1 <- normalize_dmso("1", xcell_median_1)
xcell_sp <- smooth_splines(xcell_median_norm_1)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_1"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_2
xcell_raw <- read_xcell("2", my_filepath)
xcell_raw_edited <- edit_df("2", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("2", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)

#removal of outliers
xcell_norm[,"2_OkadaicAcid.2"] <- NULL


xcell_median_2 <- calculate_median_curves("2", xcell_norm)
xcell_median_norm_2 <- normalize_dmso("2", xcell_median_2)
xcell_sp <- smooth_splines(xcell_median_norm_2)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_3
xcell_raw <- read_xcell("3", my_filepath)
xcell_raw_edited <- edit_df("3", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("3", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)


xcell_median_3 <- calculate_median_curves("3", xcell_norm)
xcell_median_norm_3 <- normalize_dmso("3", xcell_median_3)
xcell_sp <- smooth_splines(xcell_median_norm_3)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_4
xcell_raw <- read_xcell("4", my_filepath)
xcell_raw_edited <- edit_df("4", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("4", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)


xcell_median_4 <- calculate_median_curves("4", xcell_norm)
xcell_median_norm_4 <- normalize_dmso("4", xcell_median_4)
xcell_sp <- smooth_splines(xcell_median_norm_4)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_4"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_5
xcell_raw <- read_xcell("5", my_filepath)
xcell_raw_edited <- edit_df("5", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("5", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)


xcell_median_5 <- calculate_median_curves("5", xcell_norm)
xcell_median_norm_5 <- normalize_dmso("5", xcell_median_5)
xcell_sp <- smooth_splines(xcell_median_norm_5)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_5"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_6
xcell_raw <- read_xcell("6", my_filepath)
xcell_raw_edited <- edit_df("6", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("6", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)

#removal of outliers
xcell_norm[,"6_ArchazolidB.4"] <- NULL

xcell_median_6 <- calculate_median_curves("6", xcell_norm)
xcell_median_norm_6 <- normalize_dmso("6", xcell_median_6)
xcell_sp <- smooth_splines(xcell_median_norm_6)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_6"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_7
xcell_raw <- read_xcell("7", my_filepath)
xcell_raw_edited <- edit_df("7", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("7", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)

#removal of outliers
xcell_norm[,"7_H89.3"] <- NULL
#Mevastatin has to be removed completely, only two replicates and one had defective electrode
xcell_norm[,"7_Mevastatin.1"] <- NULL
xcell_norm[,"7_Mevastatin.2"] <- NULL

xcell_median_7 <- calculate_median_curves("7", xcell_norm)
xcell_median_norm_7 <- normalize_dmso("7", xcell_median_7)
xcell_sp <- smooth_splines(xcell_median_norm_7)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_7"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_8
xcell_raw <- read_xcell("8", my_filepath)
xcell_raw_edited <- edit_df("8", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("8", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)

#removal of outliers

#H89 has to be removed completely, only two replicates and one had defective electrode
xcell_norm[,"8_H89.1"] <- NULL
xcell_norm[,"8_H89.2"] <- NULL

xcell_median_8 <- calculate_median_curves("8", xcell_norm)
xcell_median_norm_8 <- normalize_dmso("8", xcell_median_8)
xcell_sp <- smooth_splines(xcell_median_norm_8)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_8"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_9
xcell_raw <- read_xcell("9", my_filepath)
xcell_raw_edited <- edit_df("9", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("9", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)


xcell_median_9 <- calculate_median_curves("9", xcell_norm)
xcell_median_norm_9 <- normalize_dmso("9", xcell_median_9)
xcell_sp <- smooth_splines(xcell_median_norm_9)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_9"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_10
xcell_raw <- read_xcell("10", my_filepath)
xcell_raw_edited <- edit_df("10", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("10", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)


xcell_median_10 <- calculate_median_curves("10", xcell_norm)
xcell_median_norm_10 <- normalize_dmso("10", xcell_median_10)
xcell_sp <- smooth_splines(xcell_median_norm_10)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_10"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_11
xcell_raw <- read_xcell("11", my_filepath)
xcell_raw_edited <- edit_df("11", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("11", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)

#removal of outliers
xcell_norm[,"11_Alsterpaullone.2"] <- NULL
xcell_norm[,"11_Podophyllotoxin.4"] <- NULL
xcell_norm[,"11_Simvastatin.4"] <- NULL
xcell_norm[,"11_Velcade.3"] <- NULL

xcell_median_11 <- calculate_median_curves("11", xcell_norm)
xcell_median_norm_11 <- normalize_dmso("11", xcell_median_11)
xcell_sp <- smooth_splines(xcell_median_norm_11)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_11"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_12
xcell_raw <- read_xcell("12", my_filepath)
xcell_raw_edited <- edit_df("12", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("12", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)


xcell_median_12 <- calculate_median_curves("12", xcell_norm)
xcell_median_norm_12 <- normalize_dmso("12", xcell_median_12)
xcell_sp <- smooth_splines(xcell_median_norm_12)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_12"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_13
xcell_raw <- read_xcell("13", my_filepath)
xcell_raw_edited <- edit_df("13", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("13", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)


xcell_median_13 <- calculate_median_curves("13", xcell_norm)
xcell_median_norm_13 <- normalize_dmso("13", xcell_median_13)
xcell_sp <- smooth_splines(xcell_median_norm_13)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_13"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_14
xcell_raw <- read_xcell("14", my_filepath)
xcell_raw_edited <- edit_df("14", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("14", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)


xcell_median_14 <- calculate_median_curves("14", xcell_norm)
xcell_median_norm_14 <- normalize_dmso("14", xcell_median_14)
xcell_sp <- smooth_splines(xcell_median_norm_14)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_14"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_15
xcell_raw <- read_xcell("15", my_filepath)
xcell_raw_edited <- edit_df("15", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("15", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)

#removal of outliers
xcell_norm[,"15_SB203580.4"] <- NULL

xcell_median_15 <- calculate_median_curves("15", xcell_norm)
xcell_median_norm_15 <- normalize_dmso("15", xcell_median_15)
xcell_sp <- smooth_splines(xcell_median_norm_15)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_15"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

##########
#xcelldata_16
xcell_raw <- read_xcell("16", my_filepath)
xcell_raw_edited <- edit_df("16", xcell_raw)
xcell_norm <- normalize_xcell(xcell_raw_edited)
results.xcell <- do_median_polish("16", xcell_norm)

#For a selection of samples with RANGEresiduals >0.5
selection.xcell<- subset(results.xcell, results.xcell[,"RNGvector"] > 0.5)
print(selection.xcell)

xcell_median_16 <- calculate_median_curves("16", xcell_norm)
xcell_median_norm_16 <- normalize_dmso("16", xcell_median_16)
xcell_sp <- smooth_splines(xcell_median_norm_16)

pheatmap(xcell_sp, scale = "column", cluster_cols = FALSE, clustering_distance_rows = "correlation")

library(gplots)
rowv <- as.dendrogram(hclust(dist(xcell_sp),method="complete"))
heatmap.2(xcell_sp, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75),
          main = paste0("xcelldata_16"), scale="column", labCol=F, key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.8)

