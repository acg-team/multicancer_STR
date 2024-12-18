
#devtools::install_github("heatherjzhou/PCAForQTL")

setwd("C:\\Users\\Public\\Files\\STR_projects\\TCGA_STR\\gene_exp")

crc_tpm <- read.table("crc_tumor.tsv", sep= "\t", header = TRUE, row.names = "gene_id")
stad_tpm <- read.table("stad_tumor.tsv", sep= "\t", header = TRUE, row.names = "gene_id")
ucec_tpm <- read.table("ucec_tumor.tsv", sep= "\t", header = TRUE, row.names = "gene_id")

crc_prcompResult<-prcomp(t(crc_tpm),center=TRUE,scale.=TRUE)
stad_prcompResult<-prcomp(t(stad_tpm),center=TRUE,scale.=TRUE)
ucec_prcompResult<-prcomp(t(ucec_tpm),center=TRUE,scale.=TRUE)

crc_PCs <- crc_prcompResult$x 
stad_PCs <- stad_prcompResult$x 
ucec_PCs <- ucec_prcompResult$x 

crc_resultRunElbow<-PCAForQTL::runElbow(prcompResult = crc_prcompResult)
stad_resultRunElbow<-PCAForQTL::runElbow(prcompResult = stad_prcompResult)
ucec_resultRunElbow<-PCAForQTL::runElbow(prcompResult = ucec_prcompResult)
print(crc_resultRunElbow)
# CRC:25; STAD: 20; UCEC:17

K_crc<-crc_resultRunElbow
K_stad <- stad_resultRunElbow
K_ucec <- ucec_resultRunElbow
PCAForQTL::makeScreePlot(crc_prcompResult,labels=c("Elbow"),values=c(K_crc),
                         titleText="Gene expression CRC")
PCAForQTL::makeScreePlot(stad_prcompResult,labels=c("Elbow"),values=c(K_stad),
                         titleText="Gene expression STAD")
PCAForQTL::makeScreePlot(ucec_prcompResult,labels=c("Elbow"),values=c(K_ucec),
                         titleText="Gene expression UCEC")

crc_PCsTop<-crc_PCs[,1:K_crc]
crc_PCsTop<-scale(crc_PCsTop)
write.table(crc_PCsTop, file = "crc_cov.tsv", quote = FALSE, sep = "\t")

stad_PCsTop<-stad_PCs[,1:K_stad]
stad_PCsTop<-scale(stad_PCsTop)
write.table(stad_PCsTop, file = "stad_cov.tsv", quote = FALSE, sep = "\t")

ucec_PCsTop<-ucec_PCs[,1:K_ucec]
ucec_PCsTop<-scale(ucec_PCsTop)
write.table(ucec_PCsTop, file = "ucec_cov.tsv", quote = FALSE, sep = "\t")



