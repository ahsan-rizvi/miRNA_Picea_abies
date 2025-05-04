#!/pdc/software/23.12/eb/software/R/4.4.1-cpeGNU-23.12/bin/Rscript
#plotPCA(rld5Family,intgroup = c('Sample','Type'),returnData = FALSE)
library(DESeq2)
library(ggplot2)
#-----------------------------------------
args                  = commandArgs(TRUE)
temp                  = args[1] # address of file with gene counts of experiment
a6                    = read.table(temp, header = TRUE)
a7                    = read.table(args[2], header = TRUE, sep = ';') #Design
a6$CV2                <-NULL
a6$Mean               <-NULL
a6$Median             <-NULL
a6$Variance           <-NULL
a81                   = (a6[,-1])#Remove gene name to form matrics
data_matrix           = data.matrix(a81)
rownames(data_matrix) = (a6[,1])#Extract gene/miRNA name
#-----------------------------------------
a7$SamplesWithPath_UPPMAX <- NULL
colnames(data_matrix) = a7$NGI_ID #Put coloumn name as stages of experiments
rownames(a7)          = a7$NGI_ID
a7$NGI_ID             <-NULL
print('----------------MATRIX------------------------------------------------------')
data_matrix <- data_matrix[, -41]
data_matrix <- data_matrix[, -52]
head(data_matrix)
print('----------------DIAMENSION MATRIX-------------------------------------------')
dim(data_matrix)
print('----------------DESIGN------------------------------------------------------')
row.names.remove <- c("P11601_145", "P11601_203")
a7<-a7[!(row.names(a7) %in% row.names.remove), ]
a7
print('----------------DESeq2------------------------------------------------------')
d      <- DESeqDataSetFromMatrix(countData=data_matrix,colData=a7,design = ~Sample3)
filter = rowMedians(counts(d),normalized=TRUE)
d      <- DESeq(d)
data1=getVarianceStabilizedData(d)

rdds   <- rlog(d)
print('----------------PCA------------- -------------------------------------------')
#pdf(file = '/cfs/klemming/home/a/ahsan2/script_pdc/DGE/PCA_ZE.pdf',width = 16,height = 8)
file22=paste(args[3], "PCA_ZE.png", sep="")

#png(file = '/cfs/klemming/home/a/ahsan2/script_pdc/DGE/PCA_ZE.png',width = 500,height = 500)
png(file = file22,width = 500,height = 500)
pcaData <- DESeq2::plotPCA(varianceStabilizingTransformation(d,blind=T),intgroup=c('Sample','Type'),returnData=TRUE)  #+ geom_text(aes(label=name),vjust=2)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2,color=factor(rdds$Sample),shape=factor(rdds$Type))) +
    geom_point(size=5) + # geom_text(aes(label=name),vjust=2) +
    theme_bw() +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() + theme(legend.title=element_blank()) +
    theme( plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=22),axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22),legend.text=element_text(size=22),panel.border = element_rect(size = 3))
   
#theme( plot.background = element_blank() , panel.grid.major = element_blank(), panel.grid.minor = element_blank() ,axis.text=element_text(size=22),axis.title.x = element_text(size = 22),axis.title.y = element_text(size = 22) ,legend.text=element_text(size=22)  ,panel.border = element_rect(size = 3), legend.title=element_blank())



dev.off()
print('----------------------------------------------------------------------------')

#------------------------------------------
method1='none'
#------------------------------------------
print('Z2.SD-Z1.SD')
res1       <- results(d,contrast=c("Sample3","Z2-SD","Z1-SD"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z2_Z1_SD=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = Name_S2_S1
#------------------------------------------

#------------------------------------------
print('Z3.SD-Z2.SD')
res1       <- results(d,contrast=c("Sample3","Z3-SD","Z2-SD"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z3_Z2_SD=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1) 
#------------------------------------------


#------------------------------------------
print('Z4.FMG-Z3.SD')
res1       <- results(d,contrast=c("Sample3","Z4-FMG","Z3-SD"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z4_FMG_Z3_SD=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)


print('Z4.ZE-Z3.SD')
res1       <- results(d,contrast=c("Sample3","Z4-ZE","Z3-SD"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z4_ZE_Z3_SD=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
#------------------------------------------


#------------------------------------------
print('Z5.FMG-Z4.FMG')
res1       <- results(d,contrast=c("Sample3","Z5-FMG","Z4-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z5_Z4_FMG=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)

print('Z5.ZE-Z4.ZE')
res1       <- results(d,contrast=c("Sample3","Z5-ZE","Z4-ZE"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z5_Z4_ZE=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
#------------------------------------------


#------------------------------------------
print('Z6.FMG-Z5.FMG')
res1       <- results(d,contrast=c("Sample3","Z6-FMG","Z5-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z6_Z5_FMG=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)

print('Z6.ZE-Z5.ZE')
res1       <- results(d,contrast=c("Sample3","Z6-ZE","Z5-ZE"),pAdjustMethod = method1,,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z6_Z5_ZE=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
#------------------------------------------


#------------------------------------------
print('Z7.FMG-Z6.FMG')
res1       <- results(d,contrast=c("Sample3","Z7-FMG","Z6-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z7_Z6_FMG=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)

print('Z7.ZE-Z6.ZE')
res1       <- results(d,contrast=c("Sample3","Z7-ZE","Z6-ZE"),pAdjustMethod = method1,,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z7_Z6_ZE=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
#------------------------------------------


#------------------------------------------
print('Z8.FMG-Z7.FMG')
res1       <- results(d,contrast=c("Sample3","Z8-FMG","Z7-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z8_Z7_FMG=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)

print('Z8.ZE-Z7.ZE')
res1       <- results(d,contrast=c("Sample3","Z8-ZE","Z7-ZE"),pAdjustMethod = method1,,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z8_Z7_ZE=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
#------------------------------------------



#------------------------------------------
print('Z9.FMG-Z8.FMG')
res1       <- results(d,contrast=c("Sample3","Z9-FMG","Z8-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z9_Z8_FMG=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)

print('Z9.ZE-Z8.ZE')
res1       <- results(d,contrast=c("Sample3","Z9-ZE","Z8-ZE"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z9_Z8_ZE=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
#------------------------------------------


#------------------------------------------
print('Z9.FMG-Z8.FMG')
res1       <- results(d,contrast=c("Sample3","Z9-FMG","Z8-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z9_Z8_FMG=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)


print('Z9.ZE-Z8.ZE')
res1       <- results(d,contrast=c("Sample3","Z9-ZE","Z8-ZE"),pAdjustMethod = method1,,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z9_Z8_ZE=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
#------------------------------------------



#------------------------------------------
print('Z10.FMG-Z9.FMG')
res1       <- results(d,contrast=c("Sample3","Z10-FMG","Z9-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z10_Z9_FMG=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.05 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)

print('Z10.ZE-Z9.ZE')
res1       <- results(d,contrast=c("Sample3","Z10-ZE","Z9-ZE"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z10_Z9_ZE=res1 #

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
#------------------------------------------


##############################################
d=cbind(Z2_Z1_SD$log2FoldChange   ,Z3_Z2_SD$log2FoldChange, Z4_FMG_Z3_SD$log2FoldChange,Z4_ZE_Z3_SD$log2FoldChange ,Z5_Z4_FMG$log2FoldChange, Z5_Z4_ZE$log2FoldChange  ,Z6_Z5_FMG$log2FoldChange, Z6_Z5_ZE$log2FoldChange ,Z7_Z6_FMG$log2FoldChange, Z7_Z6_ZE$log2FoldChange  ,Z8_Z7_FMG$log2FoldChange, Z8_Z7_ZE$log2FoldChange  ,Z9_Z8_FMG$log2FoldChange, Z9_Z8_ZE$log2FoldChange, Z10_Z9_FMG$log2FoldChange, Z10_Z9_ZE$log2FoldChange)

d[is.na(d)] <- 0

rownames(d)=rownames(res1)

colnames(d)=c('Z2_Z1_SD'   ,'Z3_Z2_SD', 'Z4_FMG_Z3_SD','Z4_ZE_Z3_SD' ,'Z5_Z4_FMG', 'Z5_Z4_ZE'  ,'Z6_Z5_FMG', 'Z6_Z5_ZE' ,'Z7_Z6_FMG', 'Z7_Z6_ZE'  ,'Z8_Z7_FMG', 'Z8_Z7_ZE'  ,'Z9_Z8_FMG', 'Z9_Z8_ZE', 'Z10_Z9_FMG', 'Z10_Z9_ZE')

#head(d)
file33=paste(args[3], "ZE_Log2FC_miRNA.csv", sep="")

#file3="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/ZE_Log2FC_miRNA.csv"
write.table(d, file=file33, quote=F,sep="\t", col.names = T, row.names = T)
##############################################


Uniq       = unique(temp)
#file2="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/DGE_ZE.txt"
#print(file2)
file22=paste(args[3], "DGE_ZE.txt", sep="")
write.table(Uniq, file=file22, quote=F,sep="\t", col.names = F, row.names = F)



colnames(data1)=a7[,3]
#head(data1)
print('----------------------------------------------------------------------------')
#colnames(data1)
#print('----------------------------------------------------------------------------')



#head(data1)
file22=paste(args[3], "DGE_ZE_VST.tsv", sep="")
#file3="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/DGE_ZE_VST.txt"
write.table(data1, file=file22, quote=F,sep="\t", col.names = T, row.names = T)





