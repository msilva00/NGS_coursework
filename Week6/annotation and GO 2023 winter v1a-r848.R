library(biomaRt)
library(AnnotationDbi)
library(ComplexHeatmap)
library(gplots)
library(org.Hs.eg.db) ; #                      Entrez based annotation package 
library(EnsDb.Hsapiens.v86); #                 Ensembl based annotation package
library(TxDb.Hsapiens.UCSC.hg38.knownGene) ; # UCSC generated annotation package

library(tidyr)
library(purrr)
library(stringr)
library(topGO)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(circlize)

sess_dir="E:/NGS_UCSC/Week6/"
#getting the RNA-seq PMA datase stripped of annotation
# outdir<-"G:/working folder/output/with siRNA/"
outdir = "E:/NGS_UCSC/Week5/"
get1<-c("il6tnfvscontrol","r848gvscontrol")
toget<-paste0(outdir,"limma_transcripts_r848vscontrol.csv")
print(toget)

touse_ = read.csv(toget, header = TRUE)
touse_ = touse_[c('gene_id', 'logFC', 'AveExpr', 'P.Value', 'adj.P.val')]
head(touse_)
write.csv(touse_,paste(sess_dir,"limma_transcripts_r848vscontrol_raw.csv",sep=''))

touse <- touse_[,-1]
rownames(touse) <- touse_[,1]
head(touse)
workfile = cbind(rownames(touse),touse[,c("logFC","adj.P.val")])
colnames(workfile)<-c("gene_id", "logFC", "adj.P.Val")
rownames(workfile)<-NULL
# touse<-read.csv(toget,header = TRUE, row.names = 1)
# head(touse)
# head(rownames(touse))

# workfile<-cbind(rownames(touse),touse[,c("logFC","adj.P.val")])
# colnames(workfile)<-c("gene_id", "logFC", "adj.P.Val")
# rownames(workfile)<-NULL

head(workfile)
#---------------------------------------------------
#annotation using annotation.dbi

columns(org.Hs.eg.db)
columns(TxDb.Hsapiens.UCSC.hg38.knownGene)
columns(EnsDb.Hsapiens.v86)


txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene
txkeys<-head(keys(txdb,keytype="GENEID"))
txkeys

ensdb<-EnsDb.Hsapiens.v86
ensid1<-head(keys(ensdb,keytype="GENEID"))
ensid2<-head(keys(ensdb,keytype="ENTREZID"))
ensid1
ensid2


entrezadd<-function(to_use){
  head(to_use)
  toclip<-to_use$gene_id
  clipped <-strtrim(toclip, 15)
  
  to_use$ENSEMBL<-clipped
  clipped
  
  to_use$ENTREZID <- mapIds(ensdb,
                            keys=clipped,
                            column="ENTREZID",
                            keytype="GENEID",
                            multiVals="first")
  to_use$start <- mapIds(ensdb,
                         keys=clipped,
                         column="TXSEQSTART",
                         keytype="GENEID",
                         multiVals="first")
  
  to_use$end <- mapIds(ensdb,
                         keys=clipped,
                         column="TXSEQEND",
                         keytype="GENEID",
                         multiVals="first")
  
  to_use$gene_type <- mapIds(ensdb,
                             keys=clipped,
                             column= "GENEBIOTYPE",
                             keytype="GENEID",
                             multiVals="first")
  
  to_use$gene_name <- mapIds(ensdb,
                          keys=clipped,
                          column="SYMBOL",
                          keytype="GENEID",
                          multiVals="first")
  to_use$full_name <- mapIds(org.Hs.eg.db,
                             keys=clipped,
                             column="GENENAME",
                             keytype="ENSEMBL",
                             multiVals="first")
  
  to_use$ONTOLOGY<-mapIds(org.Hs.eg.db,
                          keys=clipped,
                          column="ONTOLOGY",
                          keytype="ENSEMBL",
                          multiVals="first")
  
  to_use$GO<-mapIds(org.Hs.eg.db,
                    keys=clipped,
                    column="GO",
                    keytype="ENSEMBL",
                    multiVals="first")
  
  return(to_use)
}

workfile<-entrezadd(workfile)
head(workfile)

notnacount<-sum(!is.na(workfile$ENTREZID))
notnacount/length(workfile$gene_id)*100

#---------------------------------------------------
#BiomaRt

listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
length(datasets$dataset)
datasets[71:81,]


ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
filters[1:10,]
attributes = listAttributes(ensembl)
length(attributes$name)
attributes[1:16,]

searchDatasets(mart = ensembl, pattern = "sapiens")
searchFilters(mart = ensembl, pattern = "id_version")
searchAttributes(mart = ensembl, pattern = "entrez")

filterval<-"ensembl_gene_id_version"

attval<-c("chromosome_name","start_position","end_position","strand","ensembl_gene_id_version","ensembl_gene_id","entrezgene_id","transcript_biotype","external_gene_name","description")
affyids=row.names(touse)

biomval<- getBM(attributes= attval, 
                filters = filterval, 
                values = affyids, 
                mart = ensembl)

head(biomval);length(biomval$chromosome_name)
tester<-!duplicated(biomval$ensembl_gene_id_version)
biomvalb<-biomval[tester,]
sum(is.na(biomval$entrezgene_id))/length(biomval$entrezgene_id)*100
biomvalb<-biomvalb[!is.na(biomvalb$entrezgene_id),]
head(biomvalb);length(biomvalb$chromosome_name)/length(touse$logFC)*100

#------------------------------------------------------
#setting up for methods comparison
methcomp<-data.frame(GO = numeric(3), KEGG = numeric (3), GSEA = numeric(3))
rownames(methcomp)<-c("MF","CC","BP")

#------------------------------------------------------
#gene ontology (GO) analysis
toget<-paste0(outdir,"limma_transcripts_",get1[3],".csv")
hma<-data.frame(read.csv(toget,header=TRUE))
hma[1,]

toclip<-hma$gene_id
clipped <-strtrim(toclip, 15)

hma$ENSEMBL<-clipped
clipped[1:10]
hma$ENTREZID <- mapIds(org.Hs.eg.db,
                          keys=clipped,
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")
length(hma$gene_id)
hma<-hma[hma$gene_type=="protein_coding",]
hma<-unique(hma)
length(hma$gene_id)

gene.df<-hma[,c("gene_id","ENSEMBL","ENTREZID","gene_name","logFC","adj.P.val")]
gene.df[1:10,]

ego2MF <- enrichGO(gene         = gene.df$ENSEMBL,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

ego2CC <- enrichGO(gene         = gene.df$ENSEMBL,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

ego2BP <- enrichGO(gene         = gene.df$ENSEMBL,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)


methcomp[1,1]<-length(ego2MF$Description)
methcomp[2,1]<-length(ego2CC$Description)
methcomp[3,1]<-length(ego2BP$Description)

methcomp


ego2MF$Description[1:10]
ego2CC$Description[1:10]
ego2BP$Description[1:10]

ego2BP[1,]

symbolget<-function(usefile,kegget){
  gene_symbol<-vector()
  tofind<-strsplit(as.character(kegget$geneID),split = "/",fixed = FALSE)
  tosearch<-usefile$ENSEMBL
  for(j in 1:length(kegget$Count)){
    hits<-match(tofind[[j]],tosearch)
    gene_symbol<-append(gene_symbol,paste(as.vector(usefile$gene_name[hits]),collapse = " "))
    print(j)
  }
  gene_symbol<-as.data.frame(gene_symbol)
  kegget<-cbind(kegget,gene_symbol)
  return(kegget)
}

ego2MFs<-symbolget(gene.df,ego2MF)
ego2CCs<-symbolget(gene.df,ego2CC)
ego2BPs<-symbolget(gene.df,ego2BP)
ego2BPs[1,]
ego2BPs$Description

barplot(ego2BP, showCategory=10)
dotplot(ego2BP, showCategory=10)
plotGOgraph(ego2BP)

#-----------------------------------------
#retrieving the gene information from all of the DEG
genid<-c()

for(a in 1:length(get1)){
  tmp<-data.frame(gene_id = character(), gene_name = character(),full_name = character())
  toget<-paste0(outdir,"limma_transcripts_",get1[a],".csv")
  print(toget)
  hma<-data.frame(read.csv(toget,header=TRUE))
  head(hma)
  length(hma$gene_id)
  hma<-hma[hma$gene_type=="protein_coding",]
  
  tmp<-hma[,c("gene_id", "gene_name")]
  genid<-rbind(genid,tmp)
  genid<-unique(genid)
}
length(genid$gene_id)
#Adding the log2FC data from the treatment groups

for(a in 1:length(get1)){
  toget<-paste0(outdir,"limma_transcripts_",get1[a],".csv",sep = "")
  print(toget)
  hma<-data.frame(read.csv(toget,header=TRUE))
  matchera<-match(genid$gene_id,hma$gene_id)
  genid<-cbind(genid,hma[matchera,"logFC"])
}
genid[1,]
#adding the group names for the log2FC
colnames(genid)<-append(colnames(genid[1:2]),get1)
colnames(genid)

toclip<-genid$gene_id
clipped <-strtrim(toclip, 15)
clipped[1:20]
genid$ENSEMBL<-clipped

genid$ENTREZID <- mapIds(org.Hs.eg.db,
                    keys=genid$ENSEMBL,
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")

head(genid)


#creating a matrix with only the log2FC in it
heatmat<-as.matrix(genid[,get1])
heatmat[is.na(heatmat)]<- 0
head(heatmat)
#---------------------------------------------------------------
#controlled graph with heatmap.2
to_lookat<-ego2BPs[c(24,12,30),]
to_lookat$Description

gene_lookat<-function(usefile,to_lookat){
  
  for(j in 1:length(to_lookat$geneID)){
    print(paste(j,to_lookat[j,"Description"],to_lookat[j,"Count"]))
    splitter<-as.vector(to_lookat$geneID[j])
    tofind<-unlist(strsplit(splitter,split = "/",fixed = FALSE))
    tosearch<-usefile$ENSEMBL
    hits<-match(tofind,tosearch)
    genes_get<-genid[hits,]
    rownames(genes_get)<-NULL
    genes_get$description<-to_lookat[j,"Description"]
    rownames(genes_get[1,])<-to_lookat[j,"Description"]
    if(j==1){
      genes_out<-genes_get
    }else{
      genes_out<-rbind(genes_out, genes_get)
    }
    
  }
  return(genes_out) 
}

resultsfile<-gene_lookat(genid,to_lookat)
head(resultsfile)

write.csv(resultsfile,paste0(outdir,"GO_heatmap_genes.csv"), row.names = FALSE)


head(resultsfile)
tolabel<-resultsfile$description
tolabel[duplicated(tolabel)]<-""
tolabel[1]<-""
tolabel
tosep<-which(tolabel!="")
tosep<-tosep-1
tosep
heatres<-as.matrix(resultsfile[,c("il6tnfvscontrol","r848gvscontrol")])
colnames(heatres)<-c("il6tnf","r848")
head(heatres)

#numbered heatmap
heatmap.2(heatres,
          trace="none",
          key=TRUE,
          keysize = 0.75,
          cexCol = 2,
          cexRow = 0.4,
          margins = c(6,6),
          rowsep = tosep,
          sepcolor = "black",
          density.info = "none",
          Rowv =  FALSE,
          Colv=FALSE,
          dendrogram = "none",
          col =  colorRampPalette(c("darkblue","white","darkred"))(256),
          main = to_lookat$Description)

#heatmap with genes on rows
heatmap.2(heatres,
          trace="none",
          key=TRUE,
          keysize = 0.75,
          margins = c(6,6),
          cexCol = 2,
          cexRow = 0.4,
          labRow = resultsfile$gene_name,
          rowsep = tosep,
          sepcolor = "black",
          density.info = "none",
          Rowv =  FALSE,
          Colv=FALSE,
          dendrogram = "none",
          col =  colorRampPalette(c("darkblue","white","darkred"))(256),
          main = to_lookat$Description)

#--------------------------------------------------------------------------------------
#Kyoto Encyclopedia of Genes and Genomes (KEGG) method of looking for gene enrichment

keggfind<- function(search){
  tofind<- search$ENSEMBL
  ENTREZ_ID <- mapIds(org.Hs.eg.db,
                      keys=tofind,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")
  
  rawmatrix<-data.frame(ENTREZ_ID = ENTREZ_ID, LogFC = search$logFC, FDR = search$adj.P.val)
  nrow(rawmatrix)
  sigGenes <- na.exclude(rawmatrix)
  nrow(sigGenes)
  kk <- enrichKEGG(gene = sigGenes$ENTREZ_ID, 
                   organism = 'hsa', 
                   pvalueCutoff = 0.05)
  
  return(kk)
}

symbolget2<-function(workfile,kegget){
  gene_symbol<-vector()
  tofind<-strsplit(as.character(kegget$geneID),split = "/",fixed = FALSE)
  tosearch<-workfile$ENTREZID
  for(j in 1:length(kegget$Count)){
    hits<-match(tofind[[j]],tosearch)
    gene_symbol<-append(gene_symbol,paste(as.vector(workfile$gene_name[hits]),collapse = " "))
  }
  gene_symbol<-as.data.frame(gene_symbol)
  kegget<-cbind(kegget,gene_symbol)
  return(kegget)
}

kk_trial<-keggfind(gene.df)

methcomp[3,2]<-length(kk_trial$Description)
methcomp

kk_trial$Description

kk_trial[5,]
kk_trial<-symbolget2(gene.df,kk_trial)
kk_trial[5,]

write.csv(kk_trial,paste0(outdir,"PMA KEGG analysis.csv"), row.names = FALSE)

browseKEGG(kk_trial, kk_trial[5,1])

pathview(gene.data = gene.df$ENTREZID, 
         pathway.id = kk_trial[5,1], 
         species = "Homo sapiens", 
         map.symbol = T)



#selected categories [3] Human T-cell leukemia virus 1 infection [5] Acute myeloid leukemia 
kk_trial$Description
to_lookat<-kk_trial[c(3,5),]
to_lookat$Description

gene_lookat_ENTREZ<-function(usefile,to_lookat){
  
  for(j in 1:length(to_lookat$geneID)){
    print(paste(j,to_lookat[j,"Description"],to_lookat[j,"Count"]))
    splitter<-as.vector(to_lookat$geneID[j])
    tofind<-unlist(strsplit(splitter,split = "/",fixed = FALSE))
    tosearch<-usefile$ENTREZID
    hits<-match(tofind,tosearch)
    genes_get<-genid[hits,]
    rownames(genes_get)<-NULL
    genes_get$description<-to_lookat[j,"Description"]
    rownames(genes_get[1,])<-to_lookat[j,"Description"]
    
    if(j==1){
      genes_out<-genes_get
    }else{
      genes_out<-rbind(genes_out, genes_get)
    }
    
  }
  return(genes_out) 
}

resultskk<-gene_lookat_ENTREZ(genid,to_lookat)
length(resultskk$description)

head(resultskk)
tolabel<-resultskk$description
tolabel[duplicated(tolabel)]<-""
tolabel[1]<-""
tolabel
tosep<-which(tolabel!="")
tosep<-tosep-1
tosep
heatres<-as.matrix(resultsfile[,c("il6tnfvscontrol","r848gvscontrol")])
colnames(heatres)<-c("il6tnf","r848")
head(heatres)
#numbered heatmap
heatmap.2(heatres,
          trace="none",
          key=TRUE,
          keysize = 0.75,
          cexCol = 2,
          cexRow = 0.75,
          margins = c(6,6),
          rowsep = tosep,
          sepcolor = "black",
          density.info = "none",
          Rowv =  FALSE,
          Colv=FALSE,
          dendrogram = "none",
          col =  colorRampPalette(c("darkblue","white","darkred"))(256))

#heatmap with genes on rows
heatmap.2(heatres,
          trace="none",
          key=TRUE,
          keysize = 0.75,
          margins = c(6,6),
          cexCol = 2,
          cexRow = 0.75,
          labRow = resultsfile$gene_name,
          rowsep = tosep,
          sepcolor = "black",
          density.info = "none",
          Rowv =  FALSE,
          Colv=FALSE,
          dendrogram = "none",
          col =  colorRampPalette(c("darkblue","white","darkred"))(256),
          main = to_lookat$Description)

#----------------------------------------------------------------
#Gene Set Enrichment Analysis (GSEA) method of ontology evaluation

#setting up matrix needed
d<-gene.df[,c("ENTREZID", "logFC")]
d<-d[!is.na(d[,1]),]
geneList = d[,2]
geneList
names(geneList) = as.character(d[,1])
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

#Running GSEA analysis
egoGSEA_MF <- gseGO(geneList  = geneList,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "MF",
                    minGSSize    = 120,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)

egoGSEA_CC <- gseGO(geneList  = geneList,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "CC",
                    minGSSize    = 120,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)

egoGSEA_BP <- gseGO(geneList  = geneList,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    minGSSize    = 120,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)

methcomp[1,3]<-length(egoGSEA_MF$Description)
methcomp[2,3]<-length(egoGSEA_CC$Description)
methcomp[3,3]<-length(egoGSEA_BP$Description)

methcomp

#egoGSEA_MF$Description
egoGSEA_CC$Description
egoGSEA_BP$Description[1:20]
head(egoGSEA_BP)


symbolget3<-function(usefile,kegget){
  #getting gene names for core_enrichment
  core_enrichment_names<-vector()
  tofind<-strsplit(as.character(kegget$core_enrichment),split = "/",fixed = FALSE)
  length(kegget$core_enrichment)
  tosearch<-usefile$ENTREZID
  j=1
  head(kegget)
  for(j in 1:length(kegget$core_enrichment)){
    hits<-match(tofind[[j]],tosearch)
    core_enrichment_names<-append(core_enrichment_names,paste(as.vector(usefile$gene_name[hits]),collapse = " "))
    print(j)
  }
  
  gene_symbol<-as.data.frame(core_enrichment_names)
  kegget<-cbind(kegget,gene_symbol)
  return(kegget)
}

#egoGSEA_MFs<-symbolget3(gene.df,egoGSEA_MF)
egoGSEA_CCs<-symbolget3(gene.df,egoGSEA_CC)
egoGSEA_BPs<-symbolget3(gene.df,egoGSEA_BP)


egoGSEA_BPs[1:10,c("ID","Description")]
egoGSEA_BPs[1,]

gseaplot(egoGSEA_BP, geneSetID = "GO:0002376")


#---------------------------------------------------------------------------
#GSEA heatmap
egoGSEA_BP$Description
#looking for [1] "immune system process" 
egoGSEA_BPs$Description
to_lookat<-egoGSEA_BPs[c(1),]
to_lookat$Description
to_lookat[1,]

gene_lookat_GSEA<-function(usefile,to_lookat){
  
  for(j in 1:length(to_lookat$core_enrichment)){
    print(paste(j,to_lookat[j,"Description"],to_lookat[j,"Count"]))
    #print(to_lookat[j,"Description"])
    splitter<-as.vector(to_lookat$core_enrichment[j])
    tofind<-unlist(strsplit(splitter,split = "/",fixed = FALSE))
    tosearch<-usefile$ENTREZID
    hits<-match(tofind,tosearch)
    genes_get<-genid[hits,]
    rownames(genes_get)<-NULL
    genes_get$description<-to_lookat[j,"Description"]
    rownames(genes_get[1,])<-to_lookat[j,"Description"]
    if(j==1){
      genes_out<-genes_get
    }else{
      genes_out<-rbind(genes_out, genes_get)
    }
    
  }
  return(genes_out) 
}

resultsfile<-gene_lookat_GSEA(genid,to_lookat)
tolabel<-resultsfile$description
tolabel[duplicated(tolabel)]<-""
tolabel[1]<-""
tolabel
tosep<-which(tolabel!="")
tosep<-tosep-1
tosep
heatres<-as.matrix(resultsfile[,c("ATRAvsDMSO","PMAvsDMSO","vitDvsDMSO")])
colnames(heatres)<-c("ATRA", "PMA","vitD")
heatres

#numbered heatmap
heatmap.2(heatres,
          trace="none",
          key=TRUE,
          keysize = 0.75,
          cexCol = 2,
          cexRow = 0.75,
          margins = c(6,6),
          rowsep = tosep,
          sepcolor = "black",
          density.info = "none",
          Rowv =  FALSE,
          Colv=FALSE,
          dendrogram = "none",
          col =  colorRampPalette(c("darkblue","white","darkred"))(256),
          main = to_lookat$Description)

#heatmap with genes on rows
heatmap.2(heatres,
          trace="none",
          key=TRUE,
          keysize = 0.75,
          margins = c(6,6),
          cexCol = 2,
          cexRow = 0.75,
          labRow = resultsfile$gene_name,
          rowsep = tosep,
          sepcolor = "black",
          density.info = "none",
          Rowv =  FALSE,
          Colv=FALSE,
          dendrogram = "none",
          col =  colorRampPalette(c("darkblue","white","darkred"))(256),
          main = to_lookat$Description)

#------------------------------------------------------------------------
#sanity check on GSEA ouput data:
#aquiring GO numbers from annotation.dbi

genid$GO<- mapIds(org.Hs.eg.db,
                  keys=clipped,
                  column="GO",
                  keytype="ENSEMBL",
                  multiVals = "list")
length(genid$GO)
genid[1,]

#Filtering genid data frame for chosen groups
go_look<-c(egoGSEA_BPs$ID[1])
go_look
found<-grep(go_look, genid$GO)
found; length(found)
go_found<-unique(genid[found,])
go_found$description<-go_look
go_found<-go_found[,c("gene_id","gene_name","ATRAvsDMSO","PMAvsDMSO","vitDvsDMSO","ENSEMBL","ENTREZID" ,"description")]

colnames(resultsfile)
colnames(go_found)

resultsfile<-rbind(resultsfile,go_found)

#return to gsea heatmap setup

#-----------------------------------------------------------------------
#This is for annotating the initial peaks file in assignment 7
# once run, you can use it for running the Gene Ontology algorithm
#on the primary output of the Diffbind program.
#------------------------------------------------------------------------

outdir<-"G:/working folder/tpa/output/"
get1<-c("HL60_H3K27ac_edgeR","HL60_H3K4me1_edgeR")

toget<-paste0(outdir,get1[1],".txt")
touse<-makeGRangesFromDataFrame(read.table(toget, header = TRUE),
                                keep.extra.columns=TRUE,
                                ignore.strand=FALSE,
                                seqnames.field="seqnames",
                                start.field="start",
                                end.field="end",
                                strand.field="strand")
head(touse)

gtf <- rtracklayer::import("D:/genomes/human/gencode_v34/gencode.v34.annotation.gff3")

match<-as.data.frame(findOverlaps(touse,gtf, ignore.strand = TRUE))
match<-match[!duplicated(match[,1]),]; #removing duplicate peaks
match
hma<-gtf[match$subjectHits,c("gene_id","gene_type","gene_name")]
hma$log2fc<-touse$log2fc[match$queryHits]
hma$FDR<-touse$FDR[match$queryHits]
hma$ENSEMBL<-strtrim(hma$gene_id, 15)
hma$ENTREZID <- mapIds(org.Hs.eg.db,
                       keys=hma$ENSEMBL,
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
gene.df<-hma[,c("gene_id","ENSEMBL","ENTREZID","gene_name","log2fc","FDR")]
head(gene.df)


