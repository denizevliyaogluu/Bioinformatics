#ornek1 - lof
library(pec)
data("cost")
veri=cost
head(veri)
veri=veri[,c(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-15)]

library(DMwR2)
aykirideger=lofactor(veri,k=10)
plot(density(aykirideger))

aykırı=order(aykirideger,decreasing = T)[1:4]
show(aykırı) 
veri[aykırı,]

n=nrow(veri)
nokta=rep(".",n)
nokta[aykırı]="*"
col=rep("black",n)
col[aykırı]="red"
  pairs(veri[1:4],pch=nokta, col=col)
  veri=veri[-aykırı,]
  
veri[aykırı,]


#ornek2
library(GEOquery)
gds=getGEO("GDS3470")
eset=GDS2eSet(gds,do.log2 = TRUE)
veri=exprs(eset)
head(veri)
veri[1:5,1:4]
library(mvoutlier)
veri=t(data.frame(veri))
sonuc=pcout(veri,makeplot = TRUE)
order(sonuc$x.dist1,decreasing = TRUE)



#ornek3 sınıflandırma algoritması ile
library(GEOquery)
gds=getGEO("GDS3470")

eset=GDS2eSet(gds,do.log2 = TRUE)
veri=t(exprs(eset))
veri[1:5,1:4]
durum=pData(eset)$diseae.state
library(randomForest)
set.seed(22)
aykırı=randomForest(veri,durum, proximity = TRUE)
order(outlier(aykırı),decreasing = TRUE)[1:5]



#ornek4 copa aykırı değer analizi
BiocManager::install("copa")
library(GEOquery)
library(copa)
gds=getGEO("GDS3470")
eset=GDS2eSet(gds,do.log2 = TRUE)
pData(eset)[,1:2]
c1=abs(3-as.numeric(pData(eset)[,2]))
c1
sonuc=copa(eset,c1,pct=0.99)
tableCopa(sonuc)
summaryCopa(sonuc,3)
