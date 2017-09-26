#####加载所需的程序包

library(ggfortify) #pca时使用
library(limma)#差异分析时使用
library(ggplot2)#画图时使用
library(pheatmap)


#####数据读取
getwd() #查看当前路径
setwd("/Users/normacy/Desktop/Rdemo/Rtest")#设置工作路径

expData <- read.csv("exp.csv",row.names=1)
sampleInfo <- read.csv("infor.csv",row.names=1)
geneInfo <- read.csv("00geneinf.csv",row.names=1)

dim(expData)
dim(sampleInfo)
dim(geneInfo)

head(expData)
head(sampleInfo)

summary(expData)
summary(sampleInfo)

jpeg(file="histAge.jpeg")
hist(sampleInfo$Age)
dev.off()

jpeg(file="tenGeneBoxplot.jpeg")
expSample <- expData[1:10][1:10]
boxplot(expSample)
dev.off()

#从sampleInfo里筛选
sampleInfo <-merge(sampleInfo,expData,by="row.names")
sampleInfo <- sampleInfo[1:9]
rownames(sampleInfo)<-sampleInfo$Row.nemes
sampleInfo <- sampleInfo[,-1]

######质控
#detect p值控制 (无相关信息)
#PCA分析
pca.exp <- prcomp(expData)
summary(pca.exp)#查看每个主成分所占的比例
screeplot(pca.exp,type="lines")#画碎石图

jpeg(file="samplePca.jpeg")
autoplot(prcomp(expData), data = expData, shape = FALSE, label.size = 3) #一行代码画图，看离群值
dev.off()

jpeg(file="genePca.jpeg")
autoplot(prcomp(t(expData)), data = expData, shape = FALSE, label.size = 3) #一行代码画图，看离群值
dev.off()


#线性回归校正协变量
mmage=model.matrix(~-1+sampleInfo$race+sampleInfo$educ+sampleInfo$Cogdx+ sampleInfo$PMI+ sampleInfo$Neu.+ sampleInfo$Neu..1)
lm.exp=apply(t(expData),1,function(x){residuals(lm(x~mmage))})#这里老是报错，不知道为啥
expData <- lm.exp

#批次矫正（无相关信息）
#Quantile normalization

#####数据分析
#按照性别t.test(limma) 这里的参考文献是
group_list=sampleInfo$sex
design <- model.matrix(~factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(t(expData))
v <- voom(t(expData),design,normalize="quantile")
#fit <- lmFit(t(expData),design)
fit <- lmFit(v,design)
fit2 <- eBayes(fit)
tempOutput = topTable(fit2, coef=2, n=Inf)
DEG_voom = na.omit(tempOutput)
write.table(DEG_voom,"DEG_voom.txt",row.names=T, sep="\t",quote=F)

#火山图
png(file="volcano.png")
threshold <- as.factor((DEG_voom$logFC>2 | DEG_voom$logFC< -2)&DEG_voom$adj.P.Val<0.05)
r03 = ggplot(DEG_voom,aes(logFC,-1*log10(adj.P.Val ),colour=threshold))
r03 <- r03 + xlab("adjust.p.value") + ylab("log2FC") + ggtitle("Volcano plot")
r03 + geom_point()
dev.off()

#热图
png(file="heatmap.png")
pheatmap(expData,fontsize=9, fontsize_row=6, labRow=F)
dev.off()

#按照age correlation
c <- apply(expData,2,function(x){cor.test(x,sampleInfo$Age)})#这个还没想好怎么提取结果

#前十个基因两两做correlation
expSample <- expData[1:10][1:10]
jpeg(file="tenGeneCorrelation.jpeg")
panel.hist <- function(x, ...)
 {
     usr <- par("usr"); on.exit(par(usr))
     par(usr = c(usr[1:2], 0, 1.5) )
     h <- hist(x, plot = FALSE)
     breaks <- h$breaks; nB <- length(breaks)
     y <- h$counts; y <- y/max(y)
     rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
 }
 
 panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
 {
     usr <- par("usr"); on.exit(par(usr))
     par(usr = c(0, 1, 0, 1))
     r <- abs(cor(x, y))
     txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
     if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
     text(0.5, 0.5, txt, cex = cex.cor * r)
 }
 pairs(expSample, diag.panel = panel.hist,
 upper.panel = panel.smooth, lower.panel = panel.cor)
dev.off()

