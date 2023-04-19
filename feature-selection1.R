#ANATASYON PAKETİ KULLANARAK ÖZELLİK SEÇİMİ
#CMA PAKETİ KULLANARAK ÖZNİTELİK SEÇİMİ
library(GEOquery)
okunan=getGEO("GSE5057")
eset=okunan[[1]]
veri=as.matrix(t(exprs(eset)))
colname(pData(eset))
durum=pData(eset)$characteristics_ch1.2
install.packages("CMA")
BiocManager::install("CMA")
library(CMA)
set.seed(111)
ogrenme=GenerateLearningsets(y=durum,method="CV",fold=5,strat=TRUE)
data(golub)
golubY <- golub[,1]
golubX <- as.matrix(golub[,-1])
set.seed(111)
five <- GenerateLearningsets(y=golubY,method="CV",fold=5,strat=TRUE)
selttest<-GeneSelection(golubX, golubY, learningsets=five, method="t.test")
show(selttest)
toplist(selttest,k=10,iter=1)
plot(selttest,iter=1)
library(GEOquery)
okunan=getGEO("GSE5057")
save.image("path")
okunan=getGEO("GSE5057")


#GENEFILTER KULLANARAK ÖZNİTELİK SEÇİMİ
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("genefilter")
library(GEOquery)
gds=getGEO("GDS5913")
eset=GDS2eSet(gds, do.log2 = TRUE)
library(genefilter)
dim(eset)
filtrele=varFilter(eset,var.cutoff=0.9)
dim(filtrele)


#ENTROPİYE DAYALI
temp=tempfile()
download.file("https://www.cancerdata.org/system/files/publications/OSCC_Maastro_AZM_nomogram_manuscript_CSV.zip",temp)
veri=read.table(unz(temp,"OSCC_Maastro_AZM_nomogram_manuscript_CSV.csv"), sep=";",header=TRUE)
unlink(temp)
str(veri, give.head = TRUE)
veri=veri[,-22]
durum=as.factor(veri$Death_status)
library(FSelector)
onem=information.gain(durum~.,veri)
onem[order(-onem),1,drop=FALSE]
subset=cutoff.k(onem,5)
oznitelik=as.simple.formula(subset,"Durum")
print(oznitelik)


