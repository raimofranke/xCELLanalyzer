#working directory
setwd("./DMSO_test/")

#compounds
DMSO_test.raw<-read.csv2(file="DMSO_test_raw.csv", header=T)


###############################################
#Normalization
x<-as.matrix(DMSO_test.raw[,2:97])
norm<-x[1,]
norm<-as.vector(norm)
DMSO_test.norm<-x/rep(norm, each = nrow(x))
DMSO_test.norm<-DMSO_test.norm[2:163,] #remove first row (last measurement before compound addition)

    
          
#set measurement
DMSO_test.norm<-DMSO_test.norm[1:150,]
DMSO_test.norm <- as.data.frame(DMSO_test.norm)

#boxplot
postscript("figure_1.eps", width = 860, height = 600)
par(mar=c(5,3,2,2)+0.1)
boxplot(DMSO_test.norm,  ylab = "NCI", xlab = "well position", cex.axis=0.4, las=2, col = "lightgray")
dev.off()


# plot the TCRPs
my_timepoints <- DMSO_test.raw[2:151,]$t - DMSO_test.raw[2,]$t
rownames(DMSO_test.norm) <- my_timepoints

pdf(file = "DMSO_controls.pdf", paper = "a4r")
par(mfrow = c(2,3))
for (i in 1: ncol(DMSO_test.norm)){
plot(DMSO_test.norm[,i], type="l", col="blue",
     main=colnames(DMSO_test.norm)[i],cex.main=0.8, ylim=c(0.8,2), ylab="NCI", xlab="t [h]")
}
dev.off()

#calculate some medians
my_matrix <- as.matrix(DMSO_test.norm)
col_medians <- apply(my_matrix, 2, median)

#wilcox test
col_medians["G7"]<- NA
col_medians["H2"] <- NA
col_medians["H5"] <- NA
col_medians["H10"] <- NA

row_A <- col_medians[1:12]
row_B <- col_medians[13:24]
row_C <- col_medians[25:36]
row_D <- col_medians[37:48]
row_E <- col_medians[49:60]
row_F <- col_medians[61:72]
row_G <- col_medians[73:84]
row_G["G6"] <- NA
row_G["G7"] <- NA
row_H <- col_medians[85:96]
row_H["H5"] <- NA
row_H["H10"] <- NA

outer_rows <- c(row_A, row_H)
inner_rows <- c(row_B, row_C, row_D, row_E, row_F, row_G)
median(outer_rows, na.rm = T)
median(inner_rows, na.rm = T)

wilcox.test(outer_rows, inner_rows, na.rm = T, alternative = "two.sided")
