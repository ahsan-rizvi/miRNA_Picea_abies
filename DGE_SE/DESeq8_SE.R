#!/pdc/software/23.12/eb/software/R/4.4.1-cpeGNU-23.12/bin/Rscript
library(DESeq2)
library(ggplot2)

args                  = commandArgs(TRUE)
temp                  = args[1] # address of file with gene counts of experiment
a6                    = read.table(temp, header = TRUE)
a7                    = read.table(args[2], header = TRUE, sep = ',') #Design
a6$CV2                <-NULL
a6$Mean               <-NULL
a6$Median             <-NULL
a6$Variance           <-NULL
a81                   = (a6[,-1])#Remove gene name to form matrics
data_matrix           = data.matrix(a81)
rownames(data_matrix) = (a6[,1])#Extract gene/miRNA name
colnames(data_matrix) = a7$SubmittedID #Put coloumn name as stages of experiments
rownames(a7)          = a7$SubmittedID
head(data_matrix)
a7
method1='BH'
alpha1=0.01

#------------------------------------------
df=data_matrix
#df2=df[apply(df[,-1], 1, function(x) !all(x>0)),]
df2=df[rowSums(df[, -1]) > 0, ]
head(df2)
print(length(rownames(df)))
print(length(rownames(df2)))

#print(rownames(df2))

file2="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/Express_SE.txt"
write.table(rownames(df2), file=file2, quote=F,sep="\t", col.names = F, row.names = F)

#quit()


#------------------------------------------
d          <- DESeqDataSetFromMatrix(countData=data_matrix,colData=a7,design = ~stage)
d2=d
filter     = rowMedians(counts(d),normalized=TRUE)
d          <- DESeq(d)
png(file = '/cfs/klemming/home/a/ahsan2/script_pdc/DGE/PCA_SE.png',width = 500,height = 500)
data1=getVarianceStabilizedData(d)
DESeq2::plotPCA(varianceStabilizingTransformation(d,blind=T),intgroup="stage")+theme_bw() + # change theme
geom_point(size=5) + # add point size
theme( plot.background = element_blank() , panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,axis.text=element_text(size=22),axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22) ,legend.text=element_text(size=22)  ,panel.border = element_rect(size = 3), legend.title=element_blank())
dev.off()

th=0.01
print('S2-S1')
res1       <- results(d,contrast=c("stage","S2","S1"),pAdjustMethod = method1,alpha=alpha1,lfcThreshold=0.5,filter = filter)
summary(res1)
S2_S1=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$padj<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = Name_S2_S1
print('S3-S2')
res1       <- results(d,contrast=c("stage","S3","S2"),pAdjustMethod = method1,alpha=alpha1,lfcThreshold=0.5,filter = filter)
summary(res1)
S3_S2=res1
res2       = res1[complete.cases(res1), ]
Name_S3_S2 = rownames(res2)[res2$padj<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S3_S2)
#------------------------------------------
print('S4-S3')
res1       <- results(d,contrast=c("stage","S4","S3"),pAdjustMethod = method1,alpha=alpha1,lfcThreshold=0.5,filter = filter)
summary(res1)
S4_S3=res1
res2       = res1[complete.cases(res1), ]
Name_S4_S3 = rownames(res2)[res2$padj<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S4_S3)
#------------------------------------------
print('S5-S4')
res1       <- results(d,contrast=c("stage","S5","S4"),pAdjustMethod = method1,alpha=alpha1,lfcThreshold=0.5,filter = filter)
summary(res1)
S5_S4=res1
res2       = res1[complete.cases(res1), ]
Name_S5_S4 = rownames(res2)[res2$padj<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S5_S4)
#------------------------------------------
print('S6-S5')
res1       <- results(d,contrast=c("stage","S6","S5"),pAdjustMethod = method1,alpha=alpha1,lfcThreshold=0.5,filter = filter)
summary(res1)
S6_S5=res1
res2       = res1[complete.cases(res1), ]
Name_S6_S5 = rownames(res2)[res2$padj<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S6_S5)
#------------------------------------------
print('S7-S6')
res1       <- results(d,contrast=c("stage","S7","S6"),pAdjustMethod = method1,alpha=alpha1,lfcThreshold=0.5,filter = filter)
summary(res1)
S7_S6=res1
res2       = res1[complete.cases(res1), ]
Name_S7_S6 = rownames(res2)[res2$padj<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S7_S6)
#------------------------------------------
print('S8-S7')
res1       <- results(d,contrast=c("stage","S8","S7"),pAdjustMethod = method1,alpha=alpha1,lfcThreshold=0.5,filter = filter)
summary(res1)
S8_S7=res1
res2       = res1[complete.cases(res1), ]
Name_S8_S7 = rownames(res2)[res2$padj<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S8_S7)
##############################################




##############################################
d=cbind(S2_S1$log2FoldChange,S3_S2$log2FoldChange,S4_S3$log2FoldChange,S5_S4$log2FoldChange,S6_S5$log2FoldChange,S7_S6$log2FoldChange,S8_S7$log2FoldChange)
d[is.na(d)] <- 0
rownames(d)=rownames(res1)
colnames(d)=c('S2_S1','S3_S2','S4_S3','S5_S4','S6_S5','S7_S6','S8_S7')
file3="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/SE_Log2FC_miRNA.csv"
write.table(d, file=file3, quote=F,sep="\t", col.names = T, row.names = T)
##############################################


##############################################
Uniq       = unique(temp)
file2="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/DGE_SE.txt"
write.table(Uniq, file=file2, quote=F,sep="\t", col.names = F, row.names = F)
colnames(data1)=a7[,2]
file3="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/DGE_SE_VST.txt"
write.table(data1, file=file3, quote=F,sep="\t", col.names = T, row.names = T)
##############################################



