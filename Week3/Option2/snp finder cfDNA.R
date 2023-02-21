dir<-"E:/NGS_UCSC/Week3/Option2/"

samp_set<-c("ERR2838558")
tgene<- c(" brca2"," kras"," egfr"," pik3ca"," pten" )

ref_set<-c("ref brca2.csv",
           "ref kras.csv",
           "ref egfr.csv",
           "ref pik3ca.csv",
           "ref pten.csv")

sample<-read.csv(paste(dir,samp_set[1],".csv",sep = ""), header = TRUE)
# names(sample) = c("chrom","pos","ref","alt")
head(sample)
length(sample$chrom)

reference<-read.csv(paste(dir,ref_set[1],sep = ""), header = TRUE)
reference
length(reference$name)

hits<-match(sample$pos-1,reference$chromStart)
hits

sample_tokeep<-sample[!is.na(hits),]
sample_tokeep
ref_tokeep<-reference[hits[!is.na(hits)],]
ref_tokeep
hit_set<-cbind(sample_tokeep,ref_tokeep)
hit_set

#--------------------------------------
#setting up to access clinical significance annotation - see lecture 6
#this takes about a minute to run
library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
#---------------------------------------

#collecting SNP from multiple genes in a loop

for(i in 1:length(samp_set)){
  print(samp_set[i])
  for(j in 1:length(ref_set)){
    print(tgene[j])
    sample<-read.csv(paste(dir,samp_set[i],".csv",sep = ""), header = TRUE)
    print(length(sample$chrom))
    reference<-read.csv(paste(dir,ref_set[j],sep = ""), header = TRUE)
    reference$gene<-tgene[j]
    hits<-match(sample$pos-1,reference$chromStart)
    sample_tokeep<-sample[!is.na(hits),]
    ref_tokeep<-reference[hits[!is.na(hits)],]
    hit_set<-cbind(sample_tokeep,ref_tokeep)
    print(hit_set)
    if(nrow(hit_set)==0){
      print("no SNP")
      next
    }
    if(j==1){
      output<-hit_set
    }
    output<-rbind(output,hit_set)
  }
  #-----------------------------------------------
  #going to web to get clinical significance information and adding to file
  #this takes a few seconds to run
  clinsig<-getBM(attributes=c( "refsnp_id","clinical_significance"),
                 filters="snp_filter", values=output$name,
                 mart=ensembl, uniqueRows=TRUE)
  matcher<-match(output$name,clinsig$refsnp_id)
  output<- cbind(output,clinsig[matcher,])
  #------------------------------------------------
  write.csv(output, file = paste(dir,samp_set[i],"_SNP",".csv",sep = ""), row.names=FALSE)
}

