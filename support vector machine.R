library(GEOquery)
gds=getGEO("GDS3837")

eset=GDS2eSet(gds, do.log2 = TRUE)
eset


colnames(pData(eset))


durum=pData(eset)$disease.state

table(durum)

dim(eset)

annotation(eset)<-"hgu133plus2.db"
library(genefilter)
filtrelenmis=nsFilter(eset)$eset


dim(filtrelenmis)

sonveri=data.frame(t(exprs(filtrelenmis)))


set.seed(1234)
sınır=floor(.80*length(durum))
print(sınır)


ind=sample.int(n=length(durum),size=sınır,replace = F)

veriegitim=sonveri[ind,]
sınıfegitim=durum[ind]
veritest=sonveri[-ind,]
sınıftest=durum[-ind]


library(caret)
set.seed(1234)

modelNB=train(y=sınıfegitim, x=veriegitim[,1:100], method="nb", importance=TRUE)


print(modelNB)

confusionMatrix(modelNB)



tahmin=predict(modelNB,veritest)
con=table(tahmin,sınıftest)
dogruluk=sum(diag(con))/sum(con)
print(dogruluk)


library(pROC)
rocveri=roc(as.numeric(sınıftest),as.numeric(tahmin))
plot(rocveri,col="blue")

auc(rocveri)

onemligenler=varImp(modelNB,scale=FALSE)
show(onemligenler)


plot(onemligenler, top=20)

detach("package:hgu133plus2.db",unload = TRUE)
detach("package:org.Hs.eg.db",unload = TRUE)

library(hgu133plus2.db)
select(hgu133plus2.db,c("203256_at","238915_at"),c("SYMBOL","GENENAME"))

#######################################
#destek vektör makinesi

library(GEOquery)
gds=getGEO("GDS1209")
eset=GDS2eSet(gds,do.log2 = TRUE)

boxplot(eset)

library(affyPLM)
eset=normalize(eset)
boxplot(eset)

dim(eset)

colnames(pData(eset))

durum=pData(eset)$disease.state
print(durum)

annotation(eset)<-"hgu133plus2.db"
library(genefilter)
filtrelenmis=nsFilter(eset)$eset


dim(filtrelenmis)

sonveri=data.frame(t(exprs(filtrelenmis)))


set.seed(1234)
sınır=floor(.80*length(durum))
print(sınır)


ind=sample.int(n=length(durum),size=sınır,replace = F)

veriegitim=sonveri[ind,]
sınıfegitim=durum[ind]
veritest=sonveri[-ind,]
sınıftest=durum[-ind]


library(caret)
set.seed(1234)

modelSVM=train(y=sınıfegitim, x=veriegitim[,1:100], method="svmLinear", importance=TRUE)


print(modelSVM)

confusionMatrix(modelSVM)



tahmin=predict(modelSVM,veritest)
con=table(tahmin,sınıftest)
dogruluk=sum(diag(con))/sum(con)
print(dogruluk)


library(pROC)
rocveri=roc(as.numeric(sınıftest),as.numeric(tahmin))
plot(rocveri,col="blue")

auc(rocveri)

onemligenler=varImp(modelSVM,scale=FALSE)
show(onemligenler)


plot(onemligenler, top=20)


detach("package:hgu133plus2.db",unload = TRUE)
detach("package:org.Hs.eg.db",unload = TRUE)

library(hgu133plus2.db)
select(hgu133plus2.db,c("203362_s_at","201930_at","214512_s_at"),c("SYMBOL","GENENAME"))