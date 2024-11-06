library(data.table)
library(Biostrings)

genome <- readDNAStringSet("/Users/mar/BIO/PROJECTS/MPOX/NCBIgenome/GCF_014621545.1_ASM1462154v1_genomic.fna")

for(sample in c("ERR10513574","ERR10963128")){
  for(caller in c("clair3rna","lowfreq","bcftools")){

    vcf <- read.csv(paste0("/Users/mar/BIO/PROJECTS/MPOX/BAM/",caller,"/",sample,"_hybrid_",caller,".vcf.gz"), comment.char = "#", sep='\t', header = FALSE)
    vcf <- data.table(vcf)
    
    sbs <- vcf[nchar(V4) == nchar(V5) & nchar(V4) == 1]
    
    # check ref nucleotide
    for(i in 1:nrow(sbs)){
      refnt <- subseq(genome[[1]],start=sbs[i,V2],end=sbs[i,V2])
      if(as.character(refnt) != sbs[i,V4]){
        stop(sbs[i,V2])
      }
    }

    # refcnt,altcnt
    if(caller == "clair3rna"){
      sbs[,tmp := tstrsplit(V10,':')[4]]
      sbs[,refcnt := tstrsplit(tmp,',')[1]]
      sbs[,refcnt := as.numeric(refcnt)]
      sbs[,altcnt := tstrsplit(tmp,',')[2]]
      sbs[,altcnt := as.numeric(altcnt)]
    } else if(caller == "lowfreq") {
      sbs[,tmp := tstrsplit(V8,';')[4]]
      sbs[,tmp := tstrsplit(tmp,"=")[2]]
      sbs[,reffwd := tstrsplit(tmp,",")[1]]
      sbs[,reffwd := as.numeric(reffwd)]
      sbs[,refrev := tstrsplit(tmp,",")[2]]
      sbs[,refrev := as.numeric(refrev)]
      sbs[,altfwd := tstrsplit(tmp,",")[3]]
      sbs[,altfwd := as.numeric(altfwd)]
      sbs[,altrev := tstrsplit(tmp,",")[4]]
      sbs[,altrev := as.numeric(altrev)]
      sbs[,refcnt := reffwd + refrev]
      sbs[,altcnt := altfwd + altrev]
    } else if (caller == "bcftools"){
      sbs[,tmp := tstrsplit(V10,':')[2]]
      sbs[,refcnt := tstrsplit(tmp,',')[1]]
      sbs[,refcnt := as.numeric(refcnt)]
      sbs[,altcnt := tstrsplit(tmp,',')[2]]
      sbs[,altcnt := as.numeric(altcnt)]
    }
    
    # depth
    if(caller == "clair3rna"){
      sbs[,depth := tstrsplit(V10,':')[3]]
      sbs[,depth := as.numeric(depth)]
    } else if(caller == "lowfreq") {
      sbs[,tmp := tstrsplit(V8,';')[1]]
      sbs[,depth := tstrsplit(tmp,"=")[2]]
      sbs[,depth := as.numeric(depth)]
    } else if(caller == "bcftools") {
      sbs[, depth := refcnt + altcnt]
    }
    
    # vaf
    if(caller == "clair3rna"){
      sbs[,vaf := tstrsplit(V10,':')[5]]
      sbs[,vaf := as.numeric(vaf)]
    } else if(caller == "lowfreq") {
      sbs[,tmp := tstrsplit(V8,';')[2]]
      sbs[,vaf := tstrsplit(tmp,"=")[2]]
      sbs[,vaf := as.numeric(vaf)]
    } else if(caller == "bcftools") {
      sbs[,tmp := tstrsplit(V10,":")[3]]
      sbs[,vaf := as.numeric(tmp)]
    }
    
    sbs[,tmp := NULL]
    
    sbs[,motif3:= ""]
    for(i in 1:nrow(sbs)){
      sbs[i, motif3 := as.character(subseq(genome[[1]],start=sbs[i,V2]-1,end=sbs[i,V2]+1))]
    }
    
    sbs[, isAPOBEC := ifelse((substr(motif3,1,2) == "TC" & V4 == "C" & V5 == "T") | (substr(motif3,2,3) == "GA" & V4 == "G" & V5 == "A"), 1, 0)]
    

    if(caller == "bcftools"){
     sbs <- sbs[depth >= 10]
    }
    
    sbs <- sbs[order(-vaf)]
    
    write.table(sbs,paste0("/Users/mar/BIO/PROJECTS/MPOX/paper/data/",sample,"_hybrid_",caller,".txt"),row.names = F, quote = F, sep='\t')
    
  }
}
