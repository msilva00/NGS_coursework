df = read.csv("Galaxy159-Annotate_DESeq2.tabular", sep = "\t", header=F)


colnames(df) = c("GeneID", "BaseMean", "log2_FC", "StdErr", "WaldStats",
                 "p_value", "adj_pval", "chrom", "start", "end", "strand",
                 "flag", "gene_name")
head(df)

dim(df)
filtered_resdata = df[df$p_value<0.05,]
dim(filtered_resdata)

length(filtered_resdata$adj_pval)
filtered_resdata <- filtered_resdata[(filtered_resdata$log2_FC>1|filtered_resdata$log2_FC < -1),]
length(filtered_resdata$log2_FC)
length(filtered_resdata$log2_FC)/length(df$log2_FC)

write.csv(filtered_resdata, "quad_data_results.csv", row.names = FALSE)


test = DESeq(filtered_resdata)

