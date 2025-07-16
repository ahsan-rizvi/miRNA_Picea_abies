#!/pdc/software/23.12/eb/software/R/4.4.1-cpeGNU-23.12/bin/Rscript
library(ggplot2)
library(topGO)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
library(treemap)
#-----------------------------------------
args                  = commandArgs(TRUE)
temp                  = args[1] # address of file with gene counts of experiment
go_spruce_genes <- read.csv(temp, header = TRUE)
geneID2GO_spruce <- unstack(go_spruce_genes[,c(1,2)])
a81                   = (go_spruce_genes[,-2])#Remove gene name to form matrics
names(a81)            = (go_spruce_genes[,-1])
temp                  = args[2]
a7                    = read.table(temp, header = F)
goi                   = a7$V1 # genes of interest
geneList <- factor(as.integer(names(a81) %in% (goi) ))
names(geneList) <- names(a81)
#---------------------------------------

print('Input to TopGO .......................................................')
GOdata <- new("topGOdata", 
              ontology = "BP",  #BP, MF or CC
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO_spruce)
allGO <- usedGO(GOdata)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#print('Fisher test results ..................................................')
#resultFisher
#print('......................................................................')
allRes <- GenTable(GOdata, 
                   weightFisher = resultFisher, 
                   orderBy = "padj"
                   ) %>% mutate_at(c('weightFisher'), as.numeric) %>% mutate(padj = round(p.adjust(weightFisher, method = "BH"), digits = 4)) %>%   arrange(padj)
allResSig <- allRes %>% filter(padj < 0.01) %>% mutate(scores = setNames(-log10(weightFisher), GO.ID))

if(!nrow(allResSig) > 0){
  message("No significant GO hits")
} else{
  message("Some GO are enriched! Saving significant results!")
}
#png(file = '/cfs/klemming/home/a/ahsan2/script_pdc/GO_Enrich.png',width = 2600,height = 1300)
pdf(file = '/cfs/klemming/home/a/ahsan2/script_pdc/GO_Enrich.pdf',width = 26,height = 13)
treemap(allResSig, 
        index = "Term", 
        vSize = "padj", 
        type = "value", 
        vColor = "padj", 
        palette = "RdYlBu", #Greens Purples YlOrBr RdYlBu YlOrBr
        fontsize.labels=c(85,85), 
        title = "",
        inflate.labels = FALSE, 
        lowerbound.cex.labels = 0, 
        bg.labels = "NA", 
        position.legend = "none", file = FALSE)

allRes2 <- GenTable(GOdata,
                   weightFisher = resultFisher,
                   orderBy = "padj"
                   )
write.csv(allResSig, "/cfs/klemming/home/a/ahsan2/script_pdc/GO_Enrich.csv")
print(allResSig)

#print(resultFisher)
#print(allResSig$Term)
