library(GEOquery)
veri=getGEO("GDS4900")
show(veri)
dim(veri)
ozellik?veri[[1]]
temp=tempfile()
download.file("https://cancerdata.org/system/files/publications/OSCC_Maastro_AZM_nomogram_manuscript_CSV.zip",temp)
veri_oku=read.table(unz(temp,"OSCC_Maastro_AZM_nomogram_manuscript_CSV.csv"),sep=";",header=TRUE)
#hangi ozniteliklerin oldugunu gormek istiyorum
str(veri_oku, give.head=FALSE)
veriler=veri_oku[,-22]
veriler=veri_oku[,-13]
sinif=as.factor(veriler$Death_status)
head(sinif)
#r dilinde ozellik secim islemesi fSelector dedigimiz paket icerisinden gerceklesmektedir.
#install.packages('FSelector')
install.packages('FSelector')
library(FSelector)
#bilgi kazancı algoritmasını kullanarak gerçekleştirmek istiyorum
#information.gain fonksiyonu çağırılarak özellik seçiminde ilgili özelliklere ulaşmış olacağız
secilen=information.gain(sinif~.,veriler)
#secilen[order(-secilen),1,drop=FALSE)
library(FSelector)
secilen=information.gain(sinif~.,veriler)

#Özellik Secimi
#Entropiye dayalı öznitelik secimi
#Entropi tabanlı öznitelik secimi islemlerinde FSelector isimli, R paketi yuklenir.
#Bu paketin icerisinde, bilgi kazancı kazanç oranı yöntemleri ile özellik seçim işlemi gerçekleştirilmektedir. libra
#Bilgi kazancı algoritmasını kullanabilmek için information.gain(), kazanç oranı gain.ratio() 



#AĞIZ VE YUTAK KANSERİ VERŞ KÜMESİ ÜZERİNDE ÖZELLİK SEÇİMİ

#Bu veri kümesi üzeerinde bilgi kazancı algoritması kullanılarak özellik seçimi gerçekleştirilecek
#ilgili data CancerData.org adresinde yüklü olan geçici dizinden indirilmesi
temp=tempfile()
download.file("https://cancerdata.org/system/files/publications/OSCC_Maastro_AZM_nomogram_manuscript_CSV.zip",temp)
#geçici dizine indirilen sıkıştırılmış dosyanın okutulabilmesi için açılması gerekmektedir.
#Açma işlemi aşağıdaki gösterildiği biçimde unz() fonksiyonu ile yapılır.
#Dosya okunduktan sonra oluşturulan geçici dizin yok etmek için unlink() fonksiyonu kullanılır.

veri=read.table(unz(temp,"OSCC_Maastro_AZM_nomogram_manuscript_CSV.csv"), sep=";",header=TRUE)
unlink(temp)
#Veri kümesinin hangi özniteliklere sahip olduğunu göstermek için 
str(veri, give.head = TRUE)
veri=veri[,-22]
durum=as.factor(veri$Death_status) #sınıf değişkeni olarak seçiyoruz
library(FSelector)
onem=information.gain(durum~.,veri)
onem[order(-onem),1,drop=FALSE]
subset=cutoff.k(onem,5)
oznitelik=as.simple.formula(subset,"Durum")
print(oznitelik)


#MELONOMİ VERİ KÜMESİ ÜZERİNDE ÖZELLİKLE SEÇİMİ
#genefilter filtre kullanılarak öznitelik seçimi yapalım
#expressionSet bunun özelliği olan genefilter öznitelik seçimini görelim
install.packages("geneFilter")
BiocManager::install("geneFilter")
library(GEOquery)
gds=getGEO("GDS1375")
eset=GDS2eSet(gds, do.log2 = TRUE)
library(genefilter)
BiocManager::install("geneFilter")
library(genefilter)
dim(eset)
filtreke=varFilter(eset,var.cutoff=0.9)
dim(filtrele)
#anatasyon paketi kullanılarak özellik seçimi yapılmak isteniyor
#anatasyon
#CMA paketi kullanılarak öznitelik seçme işlemi
library(GEOquery)
okunan=getGEO("GSE21510")
eset=okunan[[1]]
veri=as.matrix(t(exprs(eset)))
colname(pData(eset))
durum=pData(eset)$characteristics_ch1.2
install.packages("CMA")
BiocManager::install("CMA")
library(CMA)
set.seed(111)
#capraz dogrulamaislemi icin kullanilir;
ogrenme=GenerateLearningsets(y=durum,method="CV",fold=5,strat=TRUE)
data(golub)
#extract class labels
golubY <- golub[,1]
#extract gene expression from first 10 genes
golubX <- as.matrix(golub[,-1])
#generate five different learningsets
set.seed(111)
five <- GenerateLearningsets(y=golubY,method="CV",fold=5,strat=TRUE)
#simple t-test
selttest<-GeneSelection(golubX, golubY, learningsets=five, method="t.test")
show(selttest)
toplist(selttest,k=10,iter=1)
plot(selttest,iter=1)
library(GEOquery)
okunan=getGEO("GSE21510")
save.image("path")
okunan=getGEO("GSE21510")