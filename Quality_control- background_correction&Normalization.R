

# Install Bioconductor and packages for Microarray analysis

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')  
BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("sva")
BiocManager::install("AnnotationDbi")
BiocManager::install("simpleaffy")
aBiocManager::install("hgu133plus2.db")
BiocManager::install("gcrma")
BiocManager::install("affyQCReport")
BiocManager::install("arrayQualityMetrics")

# Data from GEO website
 library(affy)
rawdata <- ReadAffy(celfile.path ="/Users/vaidehipc/Desktop/Microarray_dataanlysis/GSE32323_RAW" )
data <- as.data.frame(exprs(rawdata))


# Quality control 
library("arrayQualityMetrics")
arrayQualityMetrics(expressionset = rawdata, outdir="report_rawdataset2", force = TRUE,do.logtransform = TRUE)

# Create a PCA plot to check the outliers
library(ggplot2)
matrix <- t(exprs(rma(rawdata)))
pcs <- prcomp(matrix, center = F, scale = F)
pcs_out <- as.data.frame(pcs$x)
group <- factor(c(rep("control",17), rep("cancer", 17)))

percent <- round(pcs$sdev / sum(pcs$sdev) * 100, 2)
percent <- paste( colnames(pcs$x), "(", paste( as.character(percent), "%", ")", sep="") )

pca_plot <- ggplot(pcs_out, aes(x=PC1, y = PC2, color = group, label = as.character(c(1:34))))
pca_plot <- pca_plot + geom_point() + stat_ellipse() + xlab(percent[1]) + ylab(percent[2]) + geom_text(aes(label=as.character(c(1:34))), hjust=0, vjust=0, size = 3) + theme(text = element_text(size=10))
pca_plot


# Remove the Outliers
position <- c(1,11,34)
clean <- rawdata[-position,]


# Background Correction and Normalization

First Method is Mas5 for background correction

library(affy)
Back <- mas5(clean, normalize = T, sc = 500)

# Write the expression in csv file
write.csv(exprs(Back), "mas5_cleandata.csv", row.names = F)

mas5_clean <- exprs(Back)

# Normal mas5 data plot
var <- c(1:31)
colnames(mas5_clean) <- var
boxplot(mas5_clean , xlab = " ", col=rainbow(31), main = "mas5 background correction of the clean dataset", cex.main=0.5, cex.axis=0.5)

#log2 Normalized plot

log2_mas5_cleandata <- log2(mas5_clean)
colnames(log2_mas5_cleandata) <- var
boxplot(log2_mas5_cleandata, xlab = "", col=rainbow(31), main = "log2 Normalization of mas5 background correction", cex.main=0.5, cex.axis=0.5)

# Second method is RMA()
rma_rawdata <- rma(rawdata, normalize = T, verbose = FALSE)
rma_cleandata <- rma(clean, normalize = T, verbose = FALSE)

RMA_RAW plot;
var_raw <- c(1:34)
express_raw <- exprs(rma_rawdata)
colnames(express_raw) <- var_raw
boxplot(express_raw , xlab = " ", col="lightblue",main = "rma background correction of the raw dataset", cex.main=1, cex.axis=1)

RMA_cleandata plot;
express_clean <- exprs(rma_cleandata)
colnames(express_clean) <- var
boxplot(express_clean , xlab = " ", col="lightblue",main = "rma background correction of the clean dataset", cex.main=1, cex.axis=1)

# Gcrma package
library(gcrma)
set <- gcrma(clean, verbose = FALSE)
gcrma_cleandata <- exprs(set)
colnames(gcrma_cleandata) <- var
boxplot(gcrma_cleandata,xlab= "", col = cm.colors(31), main= "gcrma of the clean dataset", cex.main=1, cex.axis=1)


