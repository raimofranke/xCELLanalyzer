#libraries
library(class)
library(gplots)
library(MASS)


#working directory
setwd("E:/xCelligence/data/bettina/11082014")

#compounds
bettina11082014.raw<-read.csv2(file="bettina11082014_raw.csv", header=T)
bettina11082014.anno<-read.csv2(file="bettina11082014_anno_blinded.csv", header=T)
bettina11082014.match<-read.csv2(file="bettina11082014_match_blinded.csv", header=F) 
bettina11082014.match<-bettina11082014.match$V1

###############################################
#Normalization
x<-as.matrix(bettina11082014.raw[,2:97])
norm<-x[1,]
norm<-as.vector(norm)
bettina11082014.norm<-x/rep(norm, each = nrow(x))
bettina11082014.norm<-bettina11082014.norm[2:895,] #remove first row (last measurement before compound addition)


#rownamens, colnames
row.names(bettina11082014.norm)<-bettina11082014.raw$t[2:895]-bettina11082014.raw[2,1]
ma<-match(colnames(bettina11082014.norm), bettina11082014.anno$Position)
colnames(bettina11082014.norm)=bettina11082014.anno$Compound[ma]



#set measurement
bettina11082014.norm<-bettina11082014.norm[1:800,]


#Median Polishing to identify and remove outliers

RNGvector<-vector(mode="numeric")
result.coleff<-c()

#setwd("C:/Users/rmf/PowerFolders/xCell_data/reference_data/results/")
pdf("bettina11082014.pdf", paper="a4")

i<-0
repeat{
  i<-i+1
  
  
  ma<-grep(bettina11082014.match[i], colnames(bettina11082014.norm)) #matching using reg expr #colnames NOT in succesive order 1,2,3
  z<-bettina11082014.norm[,ma]
  z<-as.data.frame(z)
  z<-z[,order(names(z))]
  mp<-medpolish(z)
  boxplot(z)
  boxplot(mp$residuals, main="median polish residuals")
  mp.residuals<-data.frame(mp$residuals)
  
  rng<-diff(apply(mp.residuals, 2, range)) #calculate range for each column of mp.residuals
  RNGvector<-c(RNGvector,rng)
  result.coleff<-c(result.coleff,mp$col)
  
  
  plot(x=rownames(z), y=z[,1], type="l", col="blue", main=colnames(z), sub=("blue, green,red,orange"), ylim=c(0.2,4))
  lines(x=rownames(z), y=z[,2], type="l", col="green")
  #lines(x=rownames(z), y=z[,3], type="l", col="red")
  #lines(x=rownames(z), y=z[,4], type="l", col="orange")
  
  
  if (i==21) break
}

dev.off()

RNG<-as.data.frame(RNGvector)
result.coleff<-as.data.frame(result.coleff)
results.bettina11082014<-cbind(RNG, result.coleff)
write.csv2(results.bettina11082014, "medianpolish_bettina11082014.csv")

#For a selection of samples with IQRresiduals >0.4
selection.bettina11082014<- subset(results.bettina11082014, results.bettina11082014[,"RNGvector"] > 0.3)
print(selection.bettina11082014)


#removal of outliers
bettina11082014.norm<-as.data.frame(bettina11082014.norm)

bettina11082014.norm[,"12_DMSO.3"]<-NULL
bettina11082014.norm[,"12_GAL1.1"]<-NULL
bettina11082014.norm[,"12_GAL2.3"]<-NULL
bettina11082014.norm[,"12_GAL4.2"]<-NULL
bettina11082014.norm[,"12_TOFA.3"]<-NULL






#calculate median curves
#matching using reg expr
bettina11082014.median<-matrix(ncol=24, nrow=800)
colnames(bettina11082014.median)<-bettina11082014.match
rownames(bettina11082014.median)<-row.names(bettina11082014.norm)
i<-0
repeat{
  i<-i+1
  
  name<-bettina11082014.match[i]
  ma<-grep(bettina11082014.match[i], colnames(bettina11082014.norm))
  z<-bettina11082014.norm[,ma]
  bettina11082014.median[,i]<-apply(z, 1, median) #rowMeans(z)   #median of rows of a matrix    #rowMeans(x)
  
  if (i==24) break
}



#Normalization with DMSO median TCRP
bettina11082014.median <- t(bettina11082014.median)
bettina11082014.median <- scale(bettina11082014.median, center=bettina11082014.median["12_DMSO",], scale=F)
bettina11082014.median <- t(bettina11082014.median)
temp <- as.data.frame(bettina11082014.median)
temp[,"12_DMSO"] <- NULL
bettina11082014.median <- as.matrix(temp)

#smoothing splines für bettina11082014
xbettina11082014<-matrix(ncol=22, nrow=23)
row.names(xbettina11082014)<-colnames(bettina11082014.median)
t<-rownames(bettina11082014.median)
t<-as.numeric(t)


i<-0
repeat{
  i<-i+1
  temp<-smooth.spline(x=t, y= bettina11082014.median[,i], nknots=20)
  xbettina11082014[i,]<-temp$fit$coef
  if (i==23) break
}


###########################################################################################################
#heatmap and hierarchical clustering
#alone
rowv <- as.dendrogram(hclust(dist(xbettina11082014),method="complete"))
heatmap.2(xbettina11082014, dendrogram=c("row"), Colv=F, Rowv=rowv, hclustfun=d, col=redgreen(75), main="bettina11082014", scale="column", labCol=F,
          key=T, symkey=FALSE, density.info="none", trace="none", cexRow=0.45)
