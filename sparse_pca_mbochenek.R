suppressPackageStartupMessages({
  library(devtools)
  library(Biobase)
  library(limma)
  library(edge)
  library(genefilter)
  library(qvalue)
  library(tidyverse)
  library(corpcor)
  library(data.table)
  library(jackstraw)
})
  
  library(devtools)
  library(Biobase)
  library(limma)
  library(edge)
  library(genefilter)
  library(qvalue)
  library(tidyverse)
  library(data.table)
  library(corpcor)

con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
summary(bottomly.eset)
class(bottomly.eset)
save(bottomly.eset, file="bottomly.Rdata")

load(file="bottomly.Rdata")

edata <- as.matrix(exprs(bottomly.eset))
edata <- log2(as.matrix(edata) + 1)
edata <- edata[rowMeans(edata) > 10, ]

#Homework Problem 1.
#Make one heatmap of the aforementioned Bottomly data with the following options: 
#a) both rows and columns are clustered, 
#b) show a dendrogram only on the columns., and 
#c) scale in the column direction. Send only one heatmap. 

library(RColorBrewer)
library(gplots)
my_palette <- colorRampPalette(c("blue", "white", "orange"))(n = 299)

pdf("problem1_mbochenek.pdf")
heatmap.2(edata,
          main = "Bottomly et al.", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          dendrogram="column",  # show a dendrogram only on the columns
          scale = "column",     # scale in the column direction
          Rowv=TRUE,
          Colv=TRUE)          # both rows and columns are clustered
dev.off()

# Singular value decomposition (SVD)
edata <- t(scale(t(edata), scale=FALSE, center=TRUE))
svd.out <- svd(edata)

PC = data.table(svd.out$v,pData(bottomly.eset))

#Homework Problem 2.
#Explore different combinations of PCs in scatter plots while coloring the data points by the genetic strains. 
#Find a combination of PCs that separate the strains well. Send only one scatterplot

#As importance of components decreases with their number, we want to search only in top PC. Thus, we analyze first 6.
ggplot(PC) + geom_point(aes(x=V1, y=V2, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V1, y=V3, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V1, y=V4, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V1, y=V5, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V1, y=V6, col=as.factor(strain)))

ggplot(PC) + geom_point(aes(x=V2, y=V3, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V2, y=V4, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V2, y=V5, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V2, y=V6, col=as.factor(strain)))

ggplot(PC) + geom_point(aes(x=V3, y=V4, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V3, y=V5, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V3, y=V6, col=as.factor(strain)))

ggplot(PC) + geom_point(aes(x=V4, y=V5, col=as.factor(strain)))
ggplot(PC) + geom_point(aes(x=V4, y=V6, col=as.factor(strain)))

ggplot(PC) + geom_point(aes(x=V5, y=V6, col=as.factor(strain)))

#the answer is 3rd and 4th component
pdf("problem2_mbochenek.pdf", height=7, width=10)
ggplot(PC) + geom_point(aes(x=V3, y=V4, col=as.factor(strain))) + ggtitle("PCs separating strains")
dev.off()

#However if we want to examine all PC combinations:
# PC_df <- as.data.frame(PC)
# for (i in 1:21){
#   for (j in 1:21) {
#     if (j > i) {
#       print(ggplot(PC_df) + geom_point(aes(x=PC_df[, i], y=PC_df[, j], col=as.factor(strain))) + xlab(paste("V", i)) + ylab(paste("V", j))  )
#     }
#   }
# }

#Homework Problem 3.
#Make a scatter plot of the top 2 left singular vectors.

pdf("problem3_mbochenek.pdf", height=7, width=10)
ggplot(LSV) + geom_point(aes(x=V1, y=V2)) + ggtitle("The top 2 left singular vectors")
dev.off()

#Homework Problem 4.
#Make one figure that contains violin plots of the top 5 left singular vectors (loadings). 
#Turn the top 5 left singular vectors into a data.table (or a data.frame) and ggplot2 to plot them altogether. Do not send 5 figures!

library(dplyr)
library(tidyr)

LSV_5 <- data.table(svd.out$u[, 1:5])
LSV_5 <- gather(LSV_5, key="LSV_number", value="Values")

pdf("problem4_mbochenek.pdf", height = 7, width = 10)
ggplot(LSV_5, aes(x=0, y=Values)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_grid(~LSV_number) + xlab("") + ylab("Vector values") + 
  ggtitle("The top 5 left singular vectors")
dev.off()






