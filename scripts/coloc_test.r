#install.packages("coloc")

library(coloc)

setwd("C:\\Users\\Public\\Files\\STR_projects\\TCGA_STR\\scripts\\coloc")
gene_coloc <- read.csv("ucec_gene_coloc.csv")
me <- read.csv("ucec_me_coloc.csv")

comb <- merge(gene_coloc, me, by = "str_id", all = TRUE)
comb <- na.omit(comb)
length(unique(comb$str_id))

results_list <- list()
for (i in 1:nrow(comb)){
  row <- comb[i,]
  
  estr_data <- list(
    beta = row$coef.x,
    varbeta = row$stderror.x**2,
    MAF = row$maf.x,
    N = row$size.x,
    type = "quant",
    snp = row$str_id
  )
  
  mstr_data <- list(
    beta = row$coef.y,
    varbeta = row$stderror.y**2,
    MAF = row$maf.y,
    N = row$size.y,
    type = "quant",
    snp = row$str_id
  )
  
  coloc_result <- coloc.abf(estr_data, mstr_data)
  coloc_summary <- coloc_result$summary
  
  results_list[[i]] <- data.frame(
    str_id = row$str_id,
    PP.H0 = coloc_summary["PP.H0.abf"],
    PP.H1 = coloc_summary["PP.H1.abf"],
    PP.H2 = coloc_summary["PP.H2.abf"],
    PP.H3 = coloc_summary["PP.H3.abf"],
    PP.H4 = coloc_summary["PP.H4.abf"]
  )
}

coloc_results_df <- do.call(rbind, results_list)

length(unique(coloc_results_df[coloc_results_df$PP.H4 > 0.8, "str_id"]))

coloc_results_df$gene_eg <- comb$gene.x
coloc_results_df$gene_em <- comb$gene.y
coloc_results_df$cg <- comb$cg

write.csv(coloc_results_df, "ucec_coloc_res.csv", row.names = FALSE)

