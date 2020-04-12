# 基因芯片研究突变/未突变的套细胞淋巴瘤患者：B细胞淋巴细胞

# 本文字符格式基于UTF-8编码

# 【废弃的前数据-不合适】：1.下载数据：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47811
# 1.下载数据：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36000
# 2.读取下载的数据
# 2.1 先下载并加载必要的package，已安装并加载的，可以不用重复这些步骤
# Bioconductor用于生物芯片数据分析和基因组数据分析的软件包
source("https://bioconductor.org/biocLite.R")
# 2.1.1 Bioconductor version 3.7版本以前，使用bioCLite下载
biocLite("affy")
biocLite("simpleaffy")
biocLite("pheatmap")
biocLite("RColorBrewer")
biocLite("affyPLM")
biocLite("GOstats")

# 2.1.2 Bioconductor version 3.8 其中biocLite已在3.8废弃，改用BiocManager::install
BiocManager::install("affy")
BiocManager::install("simpleaffy")
BiocManager::install("pheatmap")
BiocManager::install("RColorBrewer")
BiocManager::install("affyPLM")
BiocManager::install("GOstats")
BiocManager::install("survival")


# 2.2 设置工作目录——放置源文件
setwd("C:/libtobrain/R/R_Package/GSE36000_RAW")

# 2.3 将GSE36000_RAW.tar解压到GSE36000_RAW目录下
# 2.4 读取基因数据
dir.files <- list.files(path = "C:/libtobrain/R/R_Package/GSE36000_RAW/", pattern = "*.CEL.gz")
dir.files
library(parallel)
library(BiocGenerics)
library(Biobase)
library(affy)
affy.data <- ReadAffy(filenames = dir.files)


# 2.4.1 查看数据类型
data.class(affy.data)
class(affy.data)
# 2.4.2 查看芯片的基本信息
show(affy.data)
# 2.4.3 查看芯片的样本名
sampleNames(affy.data)

# 2.4.4 重命名样品名
# GSM879118	MCL_IGHV_UNMUT（R1394）       mcl.ighv.unmut.1
# GSM879119	MCL_IGHV_UNMUT（R1328）       mcl.ighv.unmut.2
# GSM879120	MCL_IGHV_UNMUT（R1329）       mcl.ighv.unmut.3
# GSM879121	MCL_IGHV_UNMUT（R1306）       mcl.ighv.unmut.4
# GSM879122	MCL_IGHV_UNMUT（R1332）       mcl.ighv.unmut.5
# GSM879123	MCL_IGHV_MUT（R1399_M）       mcl.ighv.mut.1
# GSM879124	MCL_IGHV_MUT（R1400）         mcl.ighv.mut.2
# GSM879125	MCL_IGHV_UNMUT（R1587）       mcl.ighv.unmut.6
# GSM879126	MCL_IGHV_UNMUT（R1680）       mcl.ighv.unmut.7
# GSM879127	MCL_IGHV_UNMUT（268-01-5TR）  mcl.ighv.unmut.8
# GSM879128	MCL_IGHV_MUT（269-01-5TR）    mcl.ighv.mut.3
# GSM879129	MCL_IGHV_UNMUT（043-01-4TR）  mcl.ighv.unmut.9
# GSM879130	MCL_IGHV_MUT（R1302）         mcl.ighv.mut.4
# GSM879131	MCL_IGHV_MUT（R1304）         mcl.ighv.mut.5
# GSM879132	MCL_IGHV_MUT（R1305）         mcl.ighv.mut.6
# GSM879133	MCL_IGHV_MUT（R1628）         mcl.ighv.mut.7
# GSM879134	MCL_IGHV_MUT（R1629）         mcl.ighv.mut.8
# GSM879135	MCL_IGHV_MUT（R1341）         mcl.ighv.mut.9
# GSM879136	MCL_IGHV_MUT（R1333）         mcl.ighv.mut.10
# GSM879137	MCL_IGHV_MUT（R1585）         mcl.ighv.mut.11
# GSM879138	MCL_IGHV_MUT（R1589）         mcl.ighv.mut.12
# GSM879139	MCL_IGHV_MUT（R1591）         mcl.ighv.mut.13
# GSM879140	MCL_IGHV_UNMUT（R1595）       mcl.ighv.unmut.10
# GSM879141	MCL_IGHV_MUT（R1640）         mcl.ighv.mut.14
# GSM879142	MCL_IGHV_MUT（R1626）         mcl.ighv.mut.15
# GSM879143	MCL_IGHV_MUT（R1610）         mcl.ighv.mut.16
# GSM879144	MCL_IGHV_MUT（02_201）        mcl.ighv.mut.17
# GSM879145	MCL_IGHV_MUT（R1339）         mcl.ighv.mut.18
# GSM879146	MCL_IGHV_MUT（R1301）         mcl.ighv.mut.19
# GSM879147	MCL_IGHV_UNMUT（R1338）       mcl.ighv.unmut.11
# GSM879148	MCL_IGHV_UNMUT（R1762）       mcl.ighv.unmut.12
# GSM879149	MCL_IGHV_MUT（R1938）         mcl.ighv.mut.20
# GSM879150	MCL_IGHV_MUT（R1940）         mcl.ighv.mut.21
# GSM879151	MCL_IGHV_MUT（R1942）         mcl.ighv.mut.22
# GSM879152	MCL_IGHV_MUT（R1941）         mcl.ighv.mut.23
# GSM879153	MCL_IGHV_MUT（R1970）         mcl.ighv.mut.24
# GSM879154	MCL_IGHV_UNMUT（R1897）       mcl.ighv.unmut.13
# GSM879155	MCL_IGHV_UNMUT（R1388）       mcl.ighv.unmut.14

c.all <- c("mcl.ighv.unmut.1","mcl.ighv.unmut.2","mcl.ighv.unmut.3","mcl.ighv.unmut.4","mcl.ighv.unmut.5",
           "mcl.ighv.mut.1","mcl.ighv.mut.2",
           "mcl.ighv.unmut.6","mcl.ighv.unmut.7","mcl.ighv.unmut.8",
           "mcl.ighv.mut.3",
           "mcl.ighv.unmut.9",
           "mcl.ighv.mut.4","mcl.ighv.mut.5", "mcl.ighv.mut.6","mcl.ighv.mut.7", "mcl.ighv.mut.8","mcl.ighv.mut.9", "mcl.ighv.mut.10","mcl.ighv.mut.11", "mcl.ighv.mut.12","mcl.ighv.mut.13",
           "mcl.ighv.unmut.10",
           "mcl.ighv.mut.14","mcl.ighv.mut.15", "mcl.ighv.mut.16","mcl.ighv.mut.17","mcl.ighv.mut.18", "mcl.ighv.mut.19",
           "mcl.ighv.unmut.11","mcl.ighv.unmut.12",
           "mcl.ighv.mut.20","mcl.ighv.mut.21", "mcl.ighv.mut.22","mcl.ighv.mut.23","mcl.ighv.mut.24",
           "mcl.ighv.unmut.13","mcl.ighv.unmut.14")

c.all.status <- c("unmut","unmut","unmut","unmut","unmut",
                  "mut","mut",
                  "unmut","unmut","unmut",
                  "mut",
                  "unmut",
                  "mut","mut", "mut","mut", "mut","mut", "mut","mut", "mut","mut",
                  "unmut",
                  "mut","mut", "mut","mut","mut", "mut",
                  "unmut","unmut",
                  "mut","mut", "mut","mut","mut",
                  "unmut","unmut")
sampleNames(affy.data) <- c.all
sampleNames(affy.data)

# 2.5读取芯片灰度图
# 2.5.1 设置工作目录——放置结果
setwd("C:/libtobrain/R/R_Package/workSpaces/Version5/result2/")
# 2.5.2 将基因芯片进行解析读取，存放为result目录下的pdf文件，即芯片灰度图
pdf(file = "芯片灰度图.pdf")
image(affy.data[,1])
dev.off()

# 2.6 质量分析与控制
# 2.6.1 质量控制，化权重图、残差图等
library(preprocessCore)
library(gcrma)
library(affyPLM)

# 对探针数据集做线性拟合
Pset <- fitPLM(affy.data)
# 根据计算结果，画权重图
pdf(file = "权重图.pdf")
image(Pset, type="weights", which=1, main="Weights")
dev.off()
# 根据计算结果，画残差图
pdf(file = "残差图.pdf")
image(Pset,  type="resids", which=1, main="Residuals")
dev.off()
# 根据计算结果，画残差符号图
pdf(file = "残差符号图.pdf")
image(Pset, type="sign.resids", which=1, main="Residuals.sign")
dev.off()

# 2.6.2 将基因芯片进行解析读取，获取质量分析报告
library(genefilter)
library(gcrma)
library(simpleaffy)
affy.qc <- qc(affy.data)
pdf(file = "质量分析报告.pdf", width = 10, height =20)
# 图形化显示分析报告
plot(affy.qc)
# 剔除这些质量不合格的数据：mut3/unmut3,unmut8
dev.off()

# 2.6.3质量控制，绘制RLE和NUSE箱线图
# 载入一组颜色
library(RColorBrewer)
library(graph)
colors <- brewer.pal(12,"Set3")
# 绘制RLE箱线图，如果数据有部分数据远离0值，那么就需要剔除这些质量不合格的数据，这里无须剔除
Mbox(Pset,ylim=c(-1,1),col=colors,main="RLE",las=2)
# 绘制NUSE箱线图，如果数据有部分数据远离1值，那么就需要剔除这些质量不合格的数据，这里无须剔除
boxplot(Pset,ylim=c(0.95,1.2),col=colors,main="NUSE",las=2)

# 2.6.4 RNA降解，质量分析，如果图形中5'趋向3'斜率较低，那么就需要把这些不够明显区别的剔除,这里看出都明显有斜率，这里无须剔除
affy.deg <- AffyRNAdeg(affy.data)
plotAffyRNAdeg(affy.deg, col = colors)
legend("topleft", rownames(pData(affy.data)), col = colors, lwd = 1, inset = 0.05, cex = 0.5)

# 2.6.5 将以上的筛选重新组成实验组
# GSM879118	MCL_IGHV_UNMUT（R1394）       mcl.ighv.unmut.1 
# GSM879119	MCL_IGHV_UNMUT（R1328）       mcl.ighv.unmut.2
# GSM879120	MCL_IGHV_UNMUT（R1329）       mcl.ighv.unmut.3    x
# GSM879121	MCL_IGHV_UNMUT（R1306）       mcl.ighv.unmut.4
# GSM879122	MCL_IGHV_UNMUT（R1332）       mcl.ighv.unmut.5
# GSM879123	MCL_IGHV_MUT（R1399_M）       mcl.ighv.mut.1
# GSM879124	MCL_IGHV_MUT（R1400）         mcl.ighv.mut.2
# GSM879125	MCL_IGHV_UNMUT（R1587）       mcl.ighv.unmut.6
# GSM879126	MCL_IGHV_UNMUT（R1680）       mcl.ighv.unmut.7 
# GSM879127	MCL_IGHV_UNMUT（268-01-5TR）  mcl.ighv.unmut.8    x
# GSM879128	MCL_IGHV_MUT（269-01-5TR）    mcl.ighv.mut.3      x
# GSM879129	MCL_IGHV_UNMUT（043-01-4TR）  mcl.ighv.unmut.9
# GSM879130	MCL_IGHV_MUT（R1302）         mcl.ighv.mut.4 
# GSM879131	MCL_IGHV_MUT（R1304）         mcl.ighv.mut.5
# GSM879132	MCL_IGHV_MUT（R1305）         mcl.ighv.mut.6
# GSM879133	MCL_IGHV_MUT（R1628）         mcl.ighv.mut.7
# GSM879134	MCL_IGHV_MUT（R1629）         mcl.ighv.mut.8
# GSM879135	MCL_IGHV_MUT（R1341）         mcl.ighv.mut.9
# GSM879136	MCL_IGHV_MUT（R1333）         mcl.ighv.mut.10
# GSM879137	MCL_IGHV_MUT（R1585）         mcl.ighv.mut.11
# GSM879138	MCL_IGHV_MUT（R1589）         mcl.ighv.mut.12
# GSM879139	MCL_IGHV_MUT（R1591）         mcl.ighv.mut.13
# GSM879140	MCL_IGHV_UNMUT（R1595）       mcl.ighv.unmut.10
# GSM879141	MCL_IGHV_MUT（R1640）         mcl.ighv.mut.14
# GSM879142	MCL_IGHV_MUT（R1626）         mcl.ighv.mut.15
# GSM879143	MCL_IGHV_MUT（R1610）         mcl.ighv.mut.16
# GSM879144	MCL_IGHV_MUT（02_201）        mcl.ighv.mut.17
# GSM879145	MCL_IGHV_MUT（R1339）         mcl.ighv.mut.18
# GSM879146	MCL_IGHV_MUT（R1301）         mcl.ighv.mut.19
# GSM879147	MCL_IGHV_UNMUT（R1338）       mcl.ighv.unmut.11
# GSM879148	MCL_IGHV_UNMUT（R1762）       mcl.ighv.unmut.12
# GSM879149	MCL_IGHV_MUT（R1938）         mcl.ighv.mut.20
# GSM879150	MCL_IGHV_MUT（R1940）         mcl.ighv.mut.21
# GSM879151	MCL_IGHV_MUT（R1942）         mcl.ighv.mut.22
# GSM879152	MCL_IGHV_MUT（R1941）         mcl.ighv.mut.23
# GSM879153	MCL_IGHV_MUT（R1970）         mcl.ighv.mut.24
# GSM879154	MCL_IGHV_UNMUT（R1897）       mcl.ighv.unmut.13
# GSM879155	MCL_IGHV_UNMUT（R1388）       mcl.ighv.unmut.14


c.experiment  <- c("mcl.ighv.unmut.1","mcl.ighv.unmut.2","mcl.ighv.unmut.4","mcl.ighv.unmut.5","mcl.ighv.unmut.6","mcl.ighv.unmut.7","mcl.ighv.unmut.9","mcl.ighv.unmut.10",
                   "mcl.ighv.unmut.11","mcl.ighv.unmut.12", "mcl.ighv.unmut.13","mcl.ighv.unmut.14",
                   "mcl.ighv.mut.1","mcl.ighv.mut.2","mcl.ighv.mut.4","mcl.ighv.mut.5", "mcl.ighv.mut.6","mcl.ighv.mut.7", "mcl.ighv.mut.8","mcl.ighv.mut.9", "mcl.ighv.mut.10",
                   "mcl.ighv.mut.11", "mcl.ighv.mut.12","mcl.ighv.mut.13","mcl.ighv.mut.14","mcl.ighv.mut.15", "mcl.ighv.mut.16","mcl.ighv.mut.17","mcl.ighv.mut.18", "mcl.ighv.mut.19",
                   "mcl.ighv.mut.20","mcl.ighv.mut.21", "mcl.ighv.mut.22","mcl.ighv.mut.23","mcl.ighv.mut.24"
                   )

c.experiment.status <- c("unmut","unmut","unmut","unmut","unmut","unmut",
                  "unmut","unmut","unmut","unmut","unmut","unmut",
                  "mut","mut", "mut","mut", "mut","mut", "mut","mut", "mut","mut",
                  "mut","mut", "mut","mut", "mut","mut", "mut","mut", "mut","mut",
                  "mut","mut", "mut"
                  )
length(c.experiment)
length(c.experiment.status)
affy.experiment <- affy.data[,c.experiment]
colnames(affy.experiment)

# 2.7 将原始数据组成数据框，对数据进行一个分类的整理（为差异分析做准备）
affy.experiment.frame <- data.frame(SampleID = c.experiment, Disease = c.experiment.status)

sampleNames(affy.experiment)

c.unmut <- c("mcl.ighv.unmut.1","mcl.ighv.unmut.2","mcl.ighv.unmut.4","mcl.ighv.unmut.5","mcl.ighv.unmut.6","mcl.ighv.unmut.7","mcl.ighv.unmut.9","mcl.ighv.unmut.10",
             "mcl.ighv.unmut.11","mcl.ighv.unmut.12", "mcl.ighv.unmut.13","mcl.ighv.unmut.14")
c.mut <- c("mcl.ighv.mut.1","mcl.ighv.mut.2","mcl.ighv.mut.4","mcl.ighv.mut.5", "mcl.ighv.mut.6","mcl.ighv.mut.7", "mcl.ighv.mut.8","mcl.ighv.mut.9", "mcl.ighv.mut.10",
           "mcl.ighv.mut.11", "mcl.ighv.mut.12","mcl.ighv.mut.13","mcl.ighv.mut.14","mcl.ighv.mut.15", "mcl.ighv.mut.16","mcl.ighv.mut.17","mcl.ighv.mut.18", "mcl.ighv.mut.19",
           "mcl.ighv.mut.20","mcl.ighv.mut.21", "mcl.ighv.mut.22","mcl.ighv.mut.23","mcl.ighv.mut.24")
length(c.unmut)
length(c.mut)


affy.experiment.unmut <- affy.experiment[,c.unmut]
affy.experiment.mut <- affy.experiment[,c.mut]

# 2.8 预处理，（预处理）背景校正、标准化、汇总
# 预处理的方式有ms、expresso、gcrma
bgcorrect.methods()
normalize.methods(affy.experiment)
express.summary.stat.methods()

# 2.8.1 通过聚类分析，查看数据质量，预处理
# 使用gcrma算法来预处理数据，也可以通过mas5、rma算法
affy.experiment.gcrma <- gcrma(affy.experiment)
affy.experiment.gcrma.exprs <- exprs(affy.experiment.gcrma)

# 计算样品两两之间的Pearson相关系数
pearson_cor <- cor(affy.experiment.gcrma.exprs)
# 得到Pearson距离的下三角矩阵
dist.lower <- as.dist(1 - pearson_cor)
# 聚类分析、画图
hc <- hclust(dist.lower,"ave")

pdf(file = "聚类分析.pdf", width = 10, height =10)
plot(hc)
dev.off()


# PCA
# samplenames <- sub(pattern = "\\.CEL", replacement = "", colnames(affy.data.gcrma.eset))
groups <- factor(affy.experiment.frame[,2])
# BiocManager::install("affycoretools")
# BiocManager::install("affycoretools")
library("affycoretools")
affycoretools::plotPCA(affy.experiment.gcrma.exprs, addtext=c.experiment, groups=groups, groupnames=levels(c.experiment.status))


# 2.8.2 使用expresso，进行背景校正、标准化处理和汇总
affy.experiment.expresso.eset <- expresso(affy.experiment, bgcorrect.method = "mas", normalize.method = "constant", pmcorrect.method = "mas",
                                          summary.method = "mas")
affy.experiment.expresso.exprs = exprs(affy.experiment.expresso.eset) 

# 2.8.3 直接采用MAS5算法进行数据预处理
affy.experiment.mas5.eset = mas5(affy.experiment)  
# 用exprs()从eset中获取数字化的表达谱矩阵
affy.experiment.mas5.exprs = exprs(affy.experiment.mas5.eset) 


# 基因芯片差异分析
# 3.1 差异分析方法：limma
# 3.1.1 选取差异表达基因，使用Bioconductor中的limma包
# 将（未突变、突变）状态数据框数据读取为因子
library(limma)
status <- factor(affy.experiment.frame[, "Disease"])
# 构建实验设计矩阵（即未突变、突变两种差异对比）
design <- model.matrix(~-1+status)
head(design)
# 构建对比模型，比较两个实验条件下表达数据,其中contrasts可以参考design的列名
contrast.matrix <- makeContrasts(contrasts = "statusmut - statusunmut", levels=design)

# 3.1.2 线性模型拟合
fit <- lmFit(affy.experiment.gcrma.exprs , design)
colnames(affy.experiment.gcrma.exprs)
# 根据对比模型进行差值计算 
fit1 <- contrasts.fit(fit, contrast.matrix)
# 3.1.3 贝叶斯检验
fit2 <- eBayes(fit1)
View(fit2)

# 3.1.4 生成所有基因的检验结果报告
result <- topTable(fit2, coef="statusmut - statusunmut", n=nrow(fit2), lfc=log2(2))
nrow(result)  # 1446
# 用P.Value进行筛选，得到全部差异表达基因
result <- result[result[, "P.Value"]<0.01,]
result2 <- result[result[, "P.Value"]<0.001,]
# 显示一部分报告结果
head(result)
head(result2)
# 分析不同p.value值下，解析出来匹配的数据行数
nrow(result)  # 909
nrow(result2) # 701
write.table(result, file="statusmut_statusunmut+0.01.txt", quote=F, sep="\t")
write.table(result2, file="statusmut_statusunmut+0.001.txt", quote=F, sep="\t")

# 3.2 寻找在不同条件下存在差异表达的基因
# 3.2.1 对表达谱矩阵进行对数化处理
affy.experiment.mas5.exprs.log <- log(affy.experiment.mas5.exprs, 2)
probeset.id <- row.names(affy.experiment.mas5.exprs.log)
# 保存数据
write.table(affy.experiment.mas5.exprs.log, file="affy.experiment.mas5.exprs.log.txt", quote=F, sep="\t")

# 3.2.2 计算每个基因在样本中的平均表达值
unmut.mean <- apply(affy.experiment.mas5.exprs.log[, c.unmut], 1, mean)
mut.mean <- apply(affy.experiment.mas5.exprs.log[, c.mut], 1, mean)
# 计算差异基因的相对表达量 (log(fold change))
unmut.to.mut.mean <- unmut.mean - mut.mean

# 使用cbind，将数据表拼接在一起
unmut.to.mut.mean.data <- cbind(affy.experiment.mas5.exprs.log, unmut.mean, mut.mean, unmut.to.mut.mean)
head(unmut.to.mut.mean.data)
# 保存数据
write.table(unmut.to.mut.mean.data, file="unmut.to.mut.mean.data.txt", quote=F, sep="\t")

# 3.2.3 倍数法，寻找在不同条件下表达量存在2倍以上差异的基因(探针组) (log(fold change)>1 or <-1)
unmut.to.mut.mean.fc.probesets <- names(unmut.to.mut.mean[unmut.to.mut.mean >1 | unmut.to.mut.mean<(-1)])

# 3.2.4 用t test检验某个基因在突变与未突变基因中的表达是否存在差异
dataset.unmut <- affy.experiment.mas5.exprs.log[1, c.unmut]
dataset.mut <- affy.experiment.mas5.exprs.log[1, c.mut]
test.gene <- t.test(dataset.unmut, dataset.mut, "two.sided")
test.gene
test.gene$p.value
length(dataset.mut)
length(dataset.unmut)
?t.test

# 3.2.5 使用apply函数，对t test检验每个基因在未突变与突变基因中的表达是否存在差异,计算p值
dim(affy.experiment.mas5.exprs.log)
p.value.experiment.genes <- apply(affy.experiment.mas5.exprs.log, 1, function(x) { t.test(x[1:12], x[13:35]) $p.value } )

# 3.2.6 用mas5calls函数来检测每张芯片中每个基因是否表达,P/M/A（有/临界值/没有）
affy.experiment.mas5calls <- mas5calls(affy.experiment)
affy.experiment.mas5calls.exprs = exprs(affy.experiment.mas5calls)
# 保存数据
write.table(affy.experiment.mas5calls.exprs, file="affy.experiment.mas5calls.exprs.txt", quote=F, sep="\t")
affy.experiment.PMA <- apply(affy.experiment.mas5calls.exprs, 1, paste, collapse="")
head(affy.experiment.PMA)

# 3.2.7 把至少在一个芯片中有表达的基因选出来
genes.present <- names(affy.experiment.PMA[affy.experiment.PMA != "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"])

# 54675
length(affy.experiment.PMA)
# 40673
length(genes.present)
affy.experiment.mas5.exprs.log.present <- affy.experiment.mas5.exprs.log[genes.present,]

# 3.2.8 设置FDR（False Discovery Rate）方法，来校正P值
p.adjust.methods
raw.pvals.present <- p.value.experiment.genes[genes.present]
fdr.pvals.present <- p.adjust(raw.pvals.present, method="fdr")

# 3.2.9 对FDR p值按从小到大排序
fdr.pvals.present.sorted <- fdr.pvals.present[order(fdr.pvals.present)]
length(fdr.pvals.present.sorted)
# 输出FDR p值最小的10个基因
fdr.pvals.present.sorted[1:10]
fdr.pvals.present.sorted[40663:40673]
# 将结果整理成表格
expression.plus.pvals <- cbind(affy.experiment.mas5.exprs.log.present, raw.pvals.present, fdr.pvals.present)
write.table(expression.plus.pvals, "mas5_expression.plus.pvals.txt", sep="\t", quote=F)

# 3.2.10 找出那些至少在一个芯片中有表达，且在不同组织中表达有差异（FDR p值<0.23)的基因(探针组)
DE.fdr.probesets <- names(fdr.pvals.present[fdr.pvals.present < 0.05])
DE.fc.probesets <- names(unmut.to.mut.mean[unmut.to.mut.mean >1 | unmut.to.mut.mean<(-1)])
DE.probesets <- intersect(DE.fdr.probesets,DE.fc.probesets)

# 获取那些至少在一个芯片中有表达且在不同组织中表达有差异的基因的表达量
unmut.to.mut.mean.mas5_DE <- unmut.to.mut.mean.data[DE.probesets, c.experiment]

# 3.2.11 构造热图
library(pheatmap)
pdf(file = "差异表达基因热图.pdf",width = 14,height = 35)
# pheatmap(unmut.to.mut.mean.mas5_DE,color=colorRampPalette(c("green","black","red"))(100),clustering_distance_rows = "correlation",clustering_distance_cols = "euclidean",clustering_method="complete")
pheatmap(unmut.to.mut.mean.mas5_DE,color=colorRampPalette(c("green","black","red"))(100),
         clustering_distance_rows = "correlation",clustering_distance_cols = "euclidean",clustering_method="complete")
dev.off()

rnames.mas5_DE <- as.matrix(rownames(unmut.to.mut.mean.mas5_DE))
colnames(rnames.mas5_DE) <- "probeset.id"

affy.experiment.mas5.DE.mat <- cbind(rnames.mas5_DE,unmut.to.mut.mean.mas5_DE)
#保存为txt文件
write.table(affy.experiment.mas5.DE.mat, "unmut.to.mut.mean.mas5_DE.txt", sep="\t", quote=F,row.names=FALSE)


# 4.基因功能分析
# 4.1 注释工具包
# 4.1.1 加载注释工具包
# BiocManager::install(affydb)
library(XML)
library(annotate)
library(org.Hs.eg.db)

# 4.1.2 获得基因芯片注释包名称
affy.experiment@annotation
# 4.1.3 加载注释包"hgu133plus2.db"
affy.experiment.affydb <- annPkgName(affy.experiment@annotation, type="db")
library(affy.experiment.affydb, character.only=TRUE)
# 根据每个探针组的ID,获取那些至少在一个芯片中有表达且在不同组织中表达有差异的基因名、Entrez ID
unmut.to.mut.mean.symbols <- getSYMBOL(rownames(unmut.to.mut.mean.mas5_DE),affy.experiment.affydb)
unmut.to.mut.mean.EntrezID <- getEG(rownames(unmut.to.mut.mean.mas5_DE),affy.experiment.affydb)
unmut.to.mut.mean.mas5.DE.anno=cbind(unmut.to.mut.mean.EntrezID,unmut.to.mut.mean.symbols,unmut.to.mut.mean.mas5_DE)
# 保存为txt文件
write.table(unmut.to.mut.mean.mas5.DE.anno, "unmut.to.mut.mean.mas5.DE.anno.txt", sep="\t", quote=F,row.names = FALSE)

##4.2	对芯片数据进行聚类分析
#打开Cluster 3.0, 首先点File/Open data file 选取"brain_mas5_DE.txt"，
#接着点Adjust Data -> Center genes -> Mean(即将每一行的所有值减去这一行的平均值，使这一行的平均值为0。)
#然后对基因和样本都进行层次聚类分析，聚类结果保存在*.cdt文件里。
#Hierarchical -> Genes -> Cluster -> Correlation(Centred)
#Hierarchical -> Arrays -> Cluster -> Euclidean distance
#Clustering method -> Complete linkage

#打开Java TreeView，然后点File/open 选取*.cdt文件进行可视化分析。

# 5 GO富集分析
# 加载所需R包
library(Matrix)
library(Category)
library(graph)
library(GOstats)
library("GO.db")
# 5.1 获取基因芯片所有探针组与差异表达基因的EntrezID
# 提取芯片中affy.data所有探针组对应的EntrezID，注意保证uniq
entrezUniverse <- unique(getEG(rownames(affy.experiment),affy.experiment.affydb));
# 提取所有差异表达基因及其对应的EntrezID，去除na值，注意保证uniq
entrezSelected <- unique(unmut.to.mut.mean.EntrezID[!is.na(unmut.to.mut.mean.EntrezID)]);
# 设置GO富集分析的所有参数
params <- new("GOHyperGParams", geneIds = entrezSelected, universeGeneIds = entrezUniverse, 
              annotation = affy.experiment.affydb, ontology = "BP", pvalueCutoff = 0.001, conditional = FALSE, testDirection = "over");
# 5.2 对所有的GO term根据params参数做超几何检验，并进行汇总
hgOver <- hyperGTest(params);
bp <- summary(hgOver) ;

# 5.3 同时生成所有GO term的检验结果文件，每个GOterm都有指向官方网站的链接，可以获得其详细信息
htmlReport(hgOver,  file='unmut.to.mut.mean.mas5_DE_go.html') ;
head(bp)
# 保存数据
write.table(bp, "unmut.to.mut.mean.mas5_DE_go.txt", sep="\t", quote=F,row.names = FALSE)

# 根据每个探针组的ID获取对应基因Gene Symbol，并作为新的一列
result$symbols <- getSYMBOL(rownames(result), affy.experiment.affydb)
# 根据探针ID获取对应基因Entrez ID
result$EntrezID <- getEG(rownames(result), affy.experiment.affydb)
# 显示前几行
head(result)
nrow(result)
