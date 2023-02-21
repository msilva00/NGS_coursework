setwd("E:/NGS_UCSC/Week1")
to_use<- "Galaxy260-[Normalized_counts_file_on_data_245,_data_243,_and_others].tabular"

toget <- read.table(to_use, sep = "\t", header= TRUE)
toget
countdata<-toget[,2:13]
colmnnames<-c("SRR4448112", "SRR4448113","SRR4448114","SRR4448109","SRR4448110","SRR4448111","SRR4448115","SRR4448116","SRR4448117","SRR4448121","SRR4448122","SRR4448123")

rownames(countdata)<-toget$X
colnames(countdata)<-colmnnames
head(countdata)
  
normfunct<-function(countdata,tester){
  for(a in colmnnames){
    toprint<-countdata[tester,a]
    resfile[a,1]<-toprint
    resfile[a,2]<-sum(countdata[,a])
    toprint<-countdata[tester,a]
    transcriptav<-mean(resfile$transcript)
    allav<-mean(resfile$allreads)
    resfile[,3]<-resfile$transcript/transcriptav
    resfile[,4]<-resfile$allreads/allav
  }
  print(resfile)
  return(resfile)
}

resfile<-data.frame(transcript = numeric(), allreads = numeric(), transcriptnorm = numeric(),allnorm = numeric())
outfile<-data.frame(transcript = numeric(), allreads = numeric())
tester<-c("ENSG00000111640.15",
          "ENSG00000143799.13",
          "ENSG00000196262.14",
          "ENSG00000006059.4",
          "ENSG00000107159.13",
          "ENSG00000171094.18",
          "ENSG00000134058.12")  
testgenes<-c("GAPDH","PARP1","Cyclophilin a", "KRT33a","CA9", "ALK", "CDK7")

for(b in 1:length(tester)){
  # browser()
  restest<-normfunct(countdata,tester[b]) 
  transcriptcv<-sd(restest$transcript)/mean(restest$transcript)*100
  allcv<-sd(restest$allreads)/mean(restest$allreads)*100
  outfile[b,1]<-transcriptcv;print(transcriptcv)
  transcriptcv_norm<-mean(restest$transcript)
  
  if (b==1){
    outfile[b,2]<-allcv;print(allcv)
  }
  # outfile[b,3]<-transcriptcv_norm;print(transcriptcv_norm)
}
rownames(outfile)<-testgenes
outfile
plot(density(outfile$V3))
plot((outfile$V3-mean(outfile$V3))/sd(outfile$V3))
write.csv(outfile, "NormFactTest.csv", row.names = FALSE)
