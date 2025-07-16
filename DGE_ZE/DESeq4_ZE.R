#!/pdc/software/24.11/eb/software/R/4.4.2-cpeGNU-24.11/bin/Rscript
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
#----------------------------------------------

df=data_matrix
#df2=df[apply(df[,-1], 1, function(x) !all(x>0)),]
df2=df[rowSums(df[, -1]) > 0, ]
print(length(rownames(df2)))

df_SD=df2[,c(1,2,3,4,5,6,7,8,50)]
df_SD2=df_SD[rowSums(df_SD[, -1]) > 0, ]
print(length(rownames(df_SD2)))
file2="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/Express_SD.txt"
write.table(rownames(df_SD2), file=file2, quote=F,sep="\t", col.names = F, row.names = F)


df_FG=df2[,c(9,11,13,15,17,19,22,24,26,28,31,33,35,37,39,41,44,47,51,52,53,54)]
#head(df_FG)
df_FG2=df_FG[rowSums(df_FG[, -1]) > 0, ]
print(length(rownames(df_FG2)))
file2="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/Express_FG.txt"
write.table(rownames(df_FG2), file=file2, quote=F,sep="\t", col.names = F, row.names = F)


df_ZE=df2[,c(10,12,14,16,18,20,21,23,25,27,29,30,32,34,36,38,40,42,43,45,46,48)]
#head(df_ZE)
df_ZE2=df_ZE[rowSums(df_ZE[, -1]) > 0, ]
print(length(rownames(df_ZE2)))
file2="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/Express_ZE.txt"
write.table(rownames(df_ZE2), file=file2, quote=F,sep="\t", col.names = F, row.names = F)




#quit()
#----------------------------------------------

print('----------------DESeq2------------------------------------------------------')
d      <- DESeqDataSetFromMatrix(countData=data_matrix,colData=a7,design = ~Sample3)
filter = rowMedians(counts(d),normalized=TRUE)
d      <- DESeq(d)
data1=getVarianceStabilizedData(d)

rdds   <- rlog(d)


print('----------------PCA------------- -------------------------------------------')
png(file = '/cfs/klemming/home/a/ahsan2/script_pdc/DGE/PCA_ZE.png',width = 500,height = 500)
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

#quit()

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
tempSD     = Name_S2_S1
#------------------------------------------

#------------------------------------------
print('Z3.SD-Z2.SD')
res1       <- results(d,contrast=c("Sample3","Z3-SD","Z2-SD"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z3_Z2_SD=res1



res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]

temp       = union(temp,Name_S2_S1) 
tempSD       = union(tempSD,Name_S2_S1)
#------------------------------------------


#------------------------------------------
print('Z4.FMG-Z3.SD')
res1       <- results(d,contrast=c("Sample3","Z4-FMG","Z3-SD"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z4_FMG_Z3_SD=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
#tempFG     = Name_S2_S1

print('Z4.ZE-Z3.SD')
res1       <- results(d,contrast=c("Sample3","Z4-ZE","Z3-SD"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z4_ZE_Z3_SD=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempSD       = union(tempSD,Name_S2_S1)
#------------------------------------------


#------------------------------------------
print('Z5.FMG-Z4.FMG')
res1       <- results(d,contrast=c("Sample3","Z5-FMG","Z4-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z5_Z4_FMG=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempFG=Name_S2_S1

print('Z5.ZE-Z4.ZE')
res1       <- results(d,contrast=c("Sample3","Z5-ZE","Z4-ZE"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z5_Z4_ZE=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempZE=Name_S2_S1
#------------------------------------------


#------------------------------------------
print('Z6.FMG-Z5.FMG')
res1       <- results(d,contrast=c("Sample3","Z6-FMG","Z5-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z6_Z5_FMG=res1


res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempFG       = union(tempFG,Name_S2_S1)


print('Z6.ZE-Z5.ZE')
res1       <- results(d,contrast=c("Sample3","Z6-ZE","Z5-ZE"),pAdjustMethod = method1,,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z6_Z5_ZE=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempZE       = union(tempZE,Name_S2_S1)
#------------------------------------------


#------------------------------------------
print('Z7.FMG-Z6.FMG')
res1       <- results(d,contrast=c("Sample3","Z7-FMG","Z6-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z7_Z6_FMG=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempFG       = union(tempFG,Name_S2_S1)

print('Z7.ZE-Z6.ZE')
res1       <- results(d,contrast=c("Sample3","Z7-ZE","Z6-ZE"),pAdjustMethod = method1,,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z7_Z6_ZE=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempZE       = union(tempZE,Name_S2_S1)
#------------------------------------------


#------------------------------------------
print('Z8.FMG-Z7.FMG')
res1       <- results(d,contrast=c("Sample3","Z8-FMG","Z7-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z8_Z7_FMG=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempFG       = union(tempFG,Name_S2_S1)

print('Z8.ZE-Z7.ZE')
res1       <- results(d,contrast=c("Sample3","Z8-ZE","Z7-ZE"),pAdjustMethod = method1,,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z8_Z7_ZE=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempZE       = union(tempZE,Name_S2_S1)
#------------------------------------------



#------------------------------------------
print('Z9.FMG-Z8.FMG')
res1       <- results(d,contrast=c("Sample3","Z9-FMG","Z8-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z9_Z8_FMG=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempFG       = union(tempFG,Name_S2_S1)

print('Z9.ZE-Z8.ZE')
res1       <- results(d,contrast=c("Sample3","Z9-ZE","Z8-ZE"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z9_Z8_ZE=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempZE       = union(tempZE,Name_S2_S1)
#------------------------------------------


#------------------------------------------
print('Z9.FMG-Z8.FMG')
res1       <- results(d,contrast=c("Sample3","Z9-FMG","Z8-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z9_Z8_FMG=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempFG       = union(tempFG,Name_S2_S1)


print('Z9.ZE-Z8.ZE')
res1       <- results(d,contrast=c("Sample3","Z9-ZE","Z8-ZE"),pAdjustMethod = method1,,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z9_Z8_ZE=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempZE       = union(tempZE,Name_S2_S1)
#------------------------------------------



#------------------------------------------
print('Z10.FMG-Z9.FMG')
res1       <- results(d,contrast=c("Sample3","Z10-FMG","Z9-FMG"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z10_Z9_FMG=res1

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.05 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempFG       = union(tempFG,Name_S2_S1)

print('Z10.ZE-Z9.ZE')
res1       <- results(d,contrast=c("Sample3","Z10-ZE","Z9-ZE"),pAdjustMethod = method1,alpha=0.01,lfcThreshold=0.5,filter = filter)
summary(res1)
Z10_Z9_ZE=res1 #

res2       = res1[complete.cases(res1), ]
Name_S2_S1 = rownames(res2)[res2$pvalue<0.01 & abs(res2$log2FoldChange)>0.5]
temp       = union(temp,Name_S2_S1)
tempZE       = union(tempZE,Name_S2_S1)
#------------------------------------------


##############################################
d=cbind(Z2_Z1_SD$log2FoldChange   ,Z3_Z2_SD$log2FoldChange, Z4_FMG_Z3_SD$log2FoldChange,Z4_ZE_Z3_SD$log2FoldChange ,Z5_Z4_FMG$log2FoldChange, Z5_Z4_ZE$log2FoldChange  ,Z6_Z5_FMG$log2FoldChange, Z6_Z5_ZE$log2FoldChange ,Z7_Z6_FMG$log2FoldChange, Z7_Z6_ZE$log2FoldChange  ,Z8_Z7_FMG$log2FoldChange, Z8_Z7_ZE$log2FoldChange  ,Z9_Z8_FMG$log2FoldChange, Z9_Z8_ZE$log2FoldChange, Z10_Z9_FMG$log2FoldChange, Z10_Z9_ZE$log2FoldChange)

d[is.na(d)] <- 0

rownames(d)=rownames(res1)

colnames(d)=c('Z2_Z1_SD'   ,'Z3_Z2_SD', 'Z4_FMG_Z3_SD','Z4_ZE_Z3_SD' ,'Z5_Z4_FMG', 'Z5_Z4_ZE'  ,'Z6_Z5_FMG', 'Z6_Z5_ZE' ,'Z7_Z6_FMG', 'Z7_Z6_ZE'  ,'Z8_Z7_FMG', 'Z8_Z7_ZE'  ,'Z9_Z8_FMG', 'Z9_Z8_ZE', 'Z10_Z9_FMG', 'Z10_Z9_ZE')

#head(d)

file3="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/ZE_Log2FC_miRNA.csv"
write.table(d, file=file3, quote=F,sep="\t", col.names = T, row.names = T)
##############################################


Uniq       = unique(temp)
UniqSD       = unique(tempSD)
UniqFG       = unique(tempFG)
UniqZE       = unique(tempZE)
file2="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/DGE_ZE.txt"
file4="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/DGE_SD.txt"
file5="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/DGE_FG.txt"
#print(tempSD)
#print(tempZE)
#print(tempFG)

#print(file2)
write.table(UniqZE, file=file2, quote=F,sep="\t", col.names = F, row.names = F)
write.table(UniqSD, file=file4, quote=F,sep="\t", col.names = F, row.names = F)
write.table(UniqFG, file=file5, quote=F,sep="\t", col.names = F, row.names = F)

colnames(data1)=a7[,3]
#head(data1)
print('----------------------------------------------------------------------------')
#colnames(data1)
#print('----------------------------------------------------------------------------')



#head(data1)
file3="/cfs/klemming/home/a/ahsan2/script_pdc/DGE/DGE_ZE_VST.txt"
write.table(data1, file=file3, quote=F,sep="\t", col.names = T, row.names = T)




quit()





#-------------------
name='Picea_mature_1078'
print(name)
data2=data1[name,]
header=c('Z1-SD','Z2-SD','Z3-SD','Z4-FMG','Z4-ZE','Z5-FMG','Z5-ZE','Z6-FMG','Z6-ZE','Z7-FMG','Z7-ZE', 'Z8-FMG','Z8-ZE', 'Z9-FMG','Z9-ZE', 'Z10-FMG','Z10-ZE')
data3=c(median(data2[1:3]),median(data2[4:6]),median(data2[7],data2[8],data2[50]), median(data2[9],data2[11],data2[13]), median(data2[10],data2[12],data2[14]), median( data2[15],data2[17],data2[19] ), median( data2[16],data2[18],data2[20],data2[21] ) , median( data2[22],data2[24],data2[54] ) , median( data2[23],data2[25] ) , median( data2[26],data2[28],data2[51] ), median( data2[27],data2[29],data2[30] ),  median( data2[31],data2[33],data2[35] ),  median( data2[32],data2[34],data2[36] ),  median( data2[37],data2[39] ), median( data2[38],data2[40] ), median( data2[41],data2[44], data2[47], data2[52], data2[53] ), median( data2[42],data2[43], data2[45], data2[46], data2[48], data2[49] ) )
header
data3
#-------------------


#-------------------
name='Picea_mature_2515'
print(name)
data2=data1[name,]
header=c('Z1-SD','Z2-SD','Z3-SD','Z4-FMG','Z4-ZE','Z5-FMG','Z5-ZE','Z6-FMG','Z6-ZE','Z7-FMG','Z7-ZE', 'Z8-FMG','Z8-ZE', 'Z9-FMG','Z9-ZE', 'Z10-FMG','Z10-ZE')
data3=c(median(data2[1:3]),median(data2[4:6]),median(data2[7],data2[8],data2[50]), median(data2[9],data2[11],data2[13]), median(data2[10],data2[12],data2[14]), median( data2[15],data2[17],data2[19] ), median( data2[16],data2[18],data2[20],data2[21] ) , median( data2[22],data2[24],data2[54] ) , median( data2[23],data2[25] ) , median( data2[26],data2[28],data2[51] ), median( data2[27],data2[29],data2[30] ),  median( data2[31],data2[33],data2[35] ),  median( data2[32],data2[34],data2[36] ),  median( data2[37],data2[39] ), median( data2[38],data2[40] ), median( data2[41],data2[44], data2[47], data2[52], data2[53] ), median( data2[42],data2[43], data2[45], data2[46], data2[48], data2[49] ) )
header
data3
#-------------------

#-------------------
name='Picea_mature_13'
print(name)
data2=data1[name,]
header=c('Z1-SD','Z2-SD','Z3-SD','Z4-FMG','Z4-ZE','Z5-FMG','Z5-ZE','Z6-FMG','Z6-ZE','Z7-FMG','Z7-ZE', 'Z8-FMG','Z8-ZE', 'Z9-FMG','Z9-ZE', 'Z10-FMG','Z10-ZE')
data3=c(median(data2[1:3]),median(data2[4:6]),median(data2[7],data2[8],data2[50]), median(data2[9],data2[11],data2[13]), median(data2[10],data2[12],data2[14]), median( data2[15],data2[17],data2[19] ), median( data2[16],data2[18],data2[20],data2[21] ) , median( data2[22],data2[24],data2[54] ) , median( data2[23],data2[25] ) , median( data2[26],data2[28],data2[51] ), median( data2[27],data2[29],data2[30] ),  median( data2[31],data2[33],data2[35] ),  median( data2[32],data2[34],data2[36] ),  median( data2[37],data2[39] ), median( data2[38],data2[40] ), median( data2[41],data2[44], data2[47], data2[52], data2[53] ), median( data2[42],data2[43], data2[45], data2[46], data2[48], data2[49] ) )
header
data3
#-------------------

#-------------------
name='Picea_mature_674'
print(name)
data2=data1[name,]
header=c('Z1-SD','Z2-SD','Z3-SD','Z4-FMG','Z4-ZE','Z5-FMG','Z5-ZE','Z6-FMG','Z6-ZE','Z7-FMG','Z7-ZE', 'Z8-FMG','Z8-ZE', 'Z9-FMG','Z9-ZE', 'Z10-FMG','Z10-ZE')
data3=c(median(data2[1:3]),median(data2[4:6]),median(data2[7],data2[8],data2[50]), median(data2[9],data2[11],data2[13]), median(data2[10],data2[12],data2[14]), median( data2[15],data2[17],data2[19] ), median( data2[16],data2[18],data2[20],data2[21] ) , median( data2[22],data2[24],data2[54] ) , median( data2[23],data2[25] ) , median( data2[26],data2[28],data2[51] ), median( data2[27],data2[29],data2[30] ),  median( data2[31],data2[33],data2[35] ),  median( data2[32],data2[34],data2[36] ),  median( data2[37],data2[39] ), median( data2[38],data2[40] ), median( data2[41],data2[44], data2[47], data2[52], data2[53] ), median( data2[42],data2[43], data2[45], data2[46], data2[48], data2[49] ) )
header
data3
#-------------------


#-------------------
name='Picea_mature_1450'
print(name)
data2=data1[name,]
header=c('Z1-SD','Z2-SD','Z3-SD','Z4-FMG','Z4-ZE','Z5-FMG','Z5-ZE','Z6-FMG','Z6-ZE','Z7-FMG','Z7-ZE', 'Z8-FMG','Z8-ZE', 'Z9-FMG','Z9-ZE', 'Z10-FMG','Z10-ZE')
data3=c(median(data2[1:3]),median(data2[4:6]),median(data2[7],data2[8],data2[50]), median(data2[9],data2[11],data2[13]), median(data2[10],data2[12],data2[14]), median( data2[15],data2[17],data2[19] ), median( data2[16],data2[18],data2[20],data2[21] ) , median( data2[22],data2[24],data2[54] ) , median( data2[23],data2[25] ) , median( data2[26],data2[28],data2[51] ), median( data2[27],data2[29],data2[30] ),  median( data2[31],data2[33],data2[35] ),  median( data2[32],data2[34],data2[36] ),  median( data2[37],data2[39] ), median( data2[38],data2[40] ), median( data2[41],data2[44], data2[47], data2[52], data2[53] ), median( data2[42],data2[43], data2[45], data2[46], data2[48], data2[49] ) )
header
data3
#-------------------



#-------------------
name='Picea_mature_2480'
print(name)
data2=data1[name,]
header=c('Z1-SD','Z2-SD','Z3-SD','Z4-FMG','Z4-ZE','Z5-FMG','Z5-ZE','Z6-FMG','Z6-ZE','Z7-FMG','Z7-ZE', 'Z8-FMG','Z8-ZE', 'Z9-FMG','Z9-ZE', 'Z10-FMG','Z10-ZE')
data3=c(median(data2[1:3]),median(data2[4:6]),median(data2[7],data2[8],data2[50]), median(data2[9],data2[11],data2[13]), median(data2[10],data2[12],data2[14]), median( data2[15],data2[17],data2[19] ), median( data2[16],data2[18],data2[20],data2[21] ) , median( data2[22],data2[24],data2[54] ) , median( data2[23],data2[25] ) , median( data2[26],data2[28],data2[51] ), median( data2[27],data2[29],data2[30] ),  median( data2[31],data2[33],data2[35] ),  median( data2[32],data2[34],data2[36] ),  median( data2[37],data2[39] ), median( data2[38],data2[40] ), median( data2[41],data2[44], data2[47], data2[52], data2[53] ), median( data2[42],data2[43], data2[45], data2[46], data2[48], data2[49] ) )
header
data3
#-------------------


#-------------------
name='Picea_mature_3502'
print(name)
data2=data1[name,]
header=c('Z1-SD','Z2-SD','Z3-SD','Z4-FMG','Z4-ZE','Z5-FMG','Z5-ZE','Z6-FMG','Z6-ZE','Z7-FMG','Z7-ZE', 'Z8-FMG','Z8-ZE', 'Z9-FMG','Z9-ZE', 'Z10-FMG','Z10-ZE')
data3=c(median(data2[1:3]),median(data2[4:6]),median(data2[7],data2[8],data2[50]), median(data2[9],data2[11],data2[13]), median(data2[10],data2[12],data2[14]), median( data2[15],data2[17],data2[19] ), median( data2[16],data2[18],data2[20],data2[21] ) , median( data2[22],data2[24],data2[54] ) , median( data2[23],data2[25] ) , median( data2[26],data2[28],data2[51] ), median( data2[27],data2[29],data2[30] ),  median( data2[31],data2[33],data2[35] ),  median( data2[32],data2[34],data2[36] ),  median( data2[37],data2[39] ), median( data2[38],data2[40] ), median( data2[41],data2[44], data2[47], data2[52], data2[53] ), median( data2[42],data2[43], data2[45], data2[46], data2[48], data2[49] ) )
header
data3
#-------------------


#-------------------
name='Picea_mature_1078'
print(name)
data2=data1[name,]
header=c('Z1-SD','Z2-SD','Z3-SD','Z4-FMG','Z4-ZE','Z5-FMG','Z5-ZE','Z6-FMG','Z6-ZE','Z7-FMG','Z7-ZE', 'Z8-FMG','Z8-ZE', 'Z9-FMG','Z9-ZE', 'Z10-FMG','Z10-ZE')
data3=c(median(data2[1:3]),median(data2[4:6]),median(data2[7],data2[8],data2[50]), median(data2[9],data2[11],data2[13]), median(data2[10],data2[12],data2[14]), median( data2[15],data2[17],data2[19] ), median( data2[16],data2[18],data2[20],data2[21] ) , median( data2[22],data2[24],data2[54] ) , median( data2[23],data2[25] ) , median( data2[26],data2[28],data2[51] ), median( data2[27],data2[29],data2[30] ),  median( data2[31],data2[33],data2[35] ),  median( data2[32],data2[34],data2[36] ),  median( data2[37],data2[39] ), median( data2[38],data2[40] ), median( data2[41],data2[44], data2[47], data2[52], data2[53] ), median( data2[42],data2[43], data2[45], data2[46], data2[48], data2[49] ) )
header
data3
#-------------------
