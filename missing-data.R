if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
BiocManager::install("ArrayExpress")
BiocManager::install("GEOquery")

#UCI 1 random forest
library(GEOquery)
url="https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer/breast-cancer.data"
okunan = read.table(url, sep = ",", na.strings = "?")
head(okunan)
veri = data.frame("Class"=factor(okunan[,1]),
                  "age"=as.factor(okunan[,2]),
                  "menopause"=as.factor(okunan[,3]),
                  "tumor-size"=as.factor(okunan[,4]),
                  "inv-nodes"=as.factor(okunan[,5]),
                  "node-caps"=as.factor(okunan[,6]),
                  "deg-malig"=as.integer(okunan[,7]),
                  "breast"=as.factor(okunan[,8]),
                  "breast-quad"=as.factor(okunan[,9]),
                  "irradiat"=as.factor(okunan[,10]))
head(veri)
library(VIM)
aggr(veri, prop=F, numbers=T)
library(randomForest)
set.seed(222)
veri.imputed=rfImpute(irradiat~.,veri)
head(veri.imputed)


#UCI 2 knn
library(GEOquery)
url="https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer/breast-cancer.data"
okunan = read.table(url, sep = ",", na.strings = "?")
head(okunan)
veri = data.frame("Class"=factor(okunan[,1]),
                  "age"=as.factor(okunan[,2]),
                  "menopause"=as.factor(okunan[,3]),
                  "tumor-size"=as.factor(okunan[,4]),
                  "inv-nodes"=as.factor(okunan[,5]),
                  "node-caps"=as.factor(okunan[,6]),
                  "deg-malig"=as.integer(okunan[,7]),
                  "breast"=as.factor(okunan[,8]),
                  "breast-quad"=as.factor(okunan[,9]),
                  "irradiat"=as.factor(okunan[,10]))
head(veri)
library(VIM)
aggr(veri, prop=F, numbers=T)
library(impute)
set.seed(2222)
veri.imputed = impute.knn(veri, k=10)
head(veri.imputed)


#GEO 1 KNN
library(GEOquery)
gds=getGEO("GDS2635")
eset=GDS2eSet(gds,do.log2 = TRUE)
veri1=exprs(eset)
kayip=nrow(which(is.na(exprs(eset)),arr.ind = TRUE))
show(kayip)
dim(eset)
veri1=exprs(eset)
veri1[1:5,1:5]
veri=t(exprs(eset))
veri[1:5,1:4]
which(is.na(veri),arr.ind = TRUE)
veri[7,53299]
library(impute)
yeni_veri <- impute.knn(as.matrix(veri))
veri = yeni_veri$data
veri[7,53299]

modelpca = prcomp(veri)
summary(modelpca)
plot(modelpca,type="l",main="Temel Bileşenler")

#GEO 2 random forest
library(GEOquery)
gds=getGEO("GDS2635")
eset=GDS2eSet(gds,do.log2 = TRUE)
veri1=exprs(eset)
kayip=nrow(which(is.na(exprs(eset)),arr.ind = TRUE))
show(kayip)
dim(eset)
veri1=exprs(eset)
veri1[1:5,1:5]
veri=t(exprs(eset))
veri[1:5,1:4]
which(is.na(veri),arr.ind = TRUE)
veri[7,53299]
library(randomForest)
set.seed(222)
veri.imputed=rfImpute(1:5~.,veri)
head(veri.imputed)
veri[7,53299]

modelpca = prcomp(veri)
summary(modelpca)
plot(modelpca,type="l",main="Temel Bileşenler")

