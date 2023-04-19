#veri yukleme
df <- data(iris)
head(iris)

#veri setinin %90'ı 
ran <- sample(1:nrow(iris), 0.9 * nrow(iris))

#normalizasyon islemi
nor <- function(x) {(x-min(x))/max(x)-min(x)}

#veri setinin ilk 4 sutunununda normalizasyonu calistirma
iris_norm <- as.data.frame(lapply(iris[,c(1,2,3,4)],nor))
summary(iris_norm)

#egitim seti
iris_train <- iris_norm[ran,]

#test seti
iris_test <- iris_norm[-ran,]

#egitim veri setinin 5.sutununu cikarma
iris_target_category <- iris[ran,5]

iris_test_category <- iris[-ran,5]

#class paketi yukleme
library(class)

#knn islevi calistirma
pr <- knn(iris_train, iris_test, cl=iris_target_category, k=1)

#confusion matrix olusturma
tab <- table(pr, iris_test_category)

#dogruluk hesaplama
accuracy <- function(x){sum(diag(x)/(sum(rowSums(x))))*100}
accuracy(tab)

