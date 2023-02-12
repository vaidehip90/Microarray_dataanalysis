# Microarray_data_analysis

Microarray technology has enables the researcher to monitor the expression levels of the thousands of the gene simulaneously. A microarray is a glass slide on which DNA molecules are fixed on an ordered manner at localation called spots or probes. The DNA in a spot is short strech of oligo-nucleotides that corresponds to the gene.

For this project, I used the GEO datasets; GSE32323. The study states " Screening for Epigenetically Masked Genes in Colorectal Cancer using 5-aza-2’-deoxycytidine treatment, Microarray and Gene Expression Profile".Gene expression profiles for 17 pairs of cancer and non-cancerous tissues from colorectal cancer patients were measured by Affymetrix HG-U133 Plus 2.0 arrays. Five cell lines were also used to investigate genes upregulated after 5-aza-2'-deoxycytidine treatment. Normalization was performed by robust multi-array average (RMA) method by Gene Expression Console (affymetrix). The normalization procedure was separately performed for each data set of clinical samples and the cellline samples. The normalized gene expression levels were presented as log2-transformed values by RMA.

Reference: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32323

 # Microarray data analysis steps includes;

1) Data Pre-processing and Quality Control

2) Data preprocessing : Background correction and normalization.

3) Visualization : Perform Principal Component Analysis (PCA)

4) Annotation

5) Gene filtering

6) Analysis with Limma

