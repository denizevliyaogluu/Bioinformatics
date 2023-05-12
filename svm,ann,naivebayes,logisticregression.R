#SVM, ANN, naive bayes, Lojistik Regresyon
#Rastgele orman algoritması ile sınıflandırma
#veri kümesinden  boostrap yöntemi ile örnekler seçilir ve bu örneklere dayalı
#olarak sınıflandırma ağaçları oluşturulur.
#Küçük hücreli akciğer kanseri ile ilgili çalışma

#çalışmda 23 akciğer kanseri, 42 normal dokuya sahip veri kümesi üzerinde çalışılacak

#Veri GEO platformanunda GDS4794 ile indirilerek kullanılacak.

library(GEOquery)
gds=getGEO("GDS4794")

#GDS doprudan analiz edilemediğinden veri kümesinin expressionSet nesnesine 
#önüştürülmektedir

eset=GDS2eSet(gds,do.log2 = TRUE)
print(eset)

#Sınıf özelliklerini öğrenmek üzere pData(eset) kullanılır
colnames(pData(eset))

durum=pData(eset)$disease.state

#sınıf değeri arasında yer alan "small cell lung cancer" uzun açıklamadır
#bunu SCLC biçiminde kısaltmak için aşağıdaki işlem uygulanır

levels(durum)[levels(durum)=="small cell lung cancer"]="SCLC"

#kayıp verinin bulunması analizlerde sorunlara neden olabilir. 
#vei kümesindeki eksik verilerin tamamlanması gerekmektedir.
#kayığ verinni olup olmadığını araştırmak için aşağıdaki yol izlenir

kayıp=nrow(which(is.na(exprs(eset)),arr.ind = TRUE))

#veri kümesinde görüldüğü gibi eksik veriler bulunmaktadır 
#kayıp verilerin yerine yeni değerler tahmin ederek yerleştirmek için inpute paketinin 
#inpute.knn fonksiyonu kullanılır.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")

veri=t(exprs(eset))
library(impute)
veri.imputed<-impute.knn(as.matrix(veri))
kayıpyok=veri.imputed$data


dim(eset)

#filtreleme işlemini yerine getirdikten sonra veri kümesi boyutlarının küçüldüğü görülür

library(genefilter)
varFiltered=varFilter(t(kayıpyok),var.cutoff = 0.9)
dim(varFiltered)

#veri kümesinin transpozesi alınarak analiz için hazır hale getirilir.

sonveri=data.frame(t(varFiltered))

#verinin bölünmesi ½70 eğitim % 30 test

set.seed(123)
sınır=floor(.7*length(durum))
ind=sample.int(n=length(durum),size=sınır,replace = F)

veriegitim=sonveri[ind,]
sınıfegitim=durum[ind]
sınıftest=durum[-ind]
veritest=sonveri[-ind,]

#ayırma işlemi sonucunda 45 gözlemden oluşan eğitim, 20 gözlemden oluşan test 
#veri kümeleri oluşturulmuştur

dim(veriegitim)

dim(veritest)

#svm
library(e1071)
model_svm <- svm(veriegitim,method='svmRadial',metric='ROC')
show(model_svm)

öngörü=predict(model_svm,veritest)
con=table(öngörü,sınıftest)
dogruluk=sum(diag(con))/sum(con)
print(dogruluk)

library(pROC)
rocveri=roc(as.numeric(sınıftest),as.numeric(öngörü))
plot(rocveri,col="blue")

auc(rocveri)

#önemli genleri ortaya çıkartmamız lazım
varImpPlot(model_svm)


#naive bayes
library(e1071)
set.seed(33)
model_nb <- naiveBayes(x=veriegitim,y=sınıfegitim)
show(model_nb)

öngörü=predict(model_nb,veritest)
con=table(öngörü,sınıftest)
dogruluk=sum(diag(con))/sum(con)
print(dogruluk)

library(pROC)
rocveri=roc(as.numeric(sınıftest),as.numeric(öngörü))
plot(rocveri,col="blue")

auc(rocveri)

#önemli genleri ortaya çıkartmamız lazım
varImpPlot(model_nb)

#logistic regression
#Logistic Regression
log_model <- glm(sınıfegitim ~ ., data = veriegitim, family = binomial())
summary(log_model)
show(log_model)

öngörü=predict(log_model,veritest)
con=table(öngörü,sınıftest)
dogruluk=sum(diag(con))/sum(con)
print(dogruluk)

library(pROC)
rocveri=roc(as.numeric(sınıftest),as.numeric(öngörü))
plot(rocveri,col="blue")

auc(rocveri)

#ann 
library(nnet)
ann_model <- nnet(veriegitim, sınıfegitim, size=c(0.1,0.5),maxit=1)
ann_tahmin <- predict(ann_model, veritest,type="class")
ann_tahmin_sınıf <- ifelse(ann_tahmin >= 0.5, "SCLC", "Normal")
mean(ann_tahmin_sınıf == sınıftest)


