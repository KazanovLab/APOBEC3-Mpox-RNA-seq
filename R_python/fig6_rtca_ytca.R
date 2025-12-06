library(data.table)
library(stringi)
library(ggplot2)

crDNA <- function(dna)
{
  return(stri_reverse(chartr("acgtACGT","tgcaTGCA",dna)))
}

genome <- readDNAStringSet("NCBIgenome/GCF_014621545.1_ASM1462154v1_genomic.fna")

gcyt <- 0
grtca <- 0
gytca <- 0
for(i in 1:length(genome[[1]])){
  print(i)
  
  curnt <- as.character(subseq(genome[[1]],start=i,end=i))
  
  if(curnt == "C" || curnt == "G")
    gcyt <- gcyt + 1
  else
    next
  
  if(i < 3 || i > (length(genome[[1]])-3) )
    next
  
  motif5 <- as.character(subseq(genome[[1]],start=i-2,end=i+2))
  if(substr(motif5,3,3) == "C"){
    motif4 <- substr(motif5,1,4)
  }
  else if(substr(motif5,3,3) == "G"){
    motif4 <- crDNA(substr(motif5,2,5))
  }
  else
    next

  if(substr(motif4,2,4) != "TCA")
    next
  
  if(motif4 %in% c("TTCA","CTCA"))
    gytca <- gytca + 1
  else if(motif4 %in% c("GTCA","ATCA"))
    grtca <- grtca + 1
  
}

print(gcyt)
print(gytca)
print(grtca)

res <- data.table()
for(f in c("human1","human2","monkey1","monkey2","genome1","genome2","isidro","otoole")){

  print(f)
  
positions <- read.csv(paste0("positions/short_format/",f,".csv"),sep='\t')
positions <- data.table(positions)

ytca <- 0
rtca <- 0
cyt <- 0
for(i in 1:nrow(positions)){
  
  motif5 <- as.character(subseq(genome[[1]],start=positions[i,pos]-2,end=positions[i,pos]+2))
  if(substr(motif5,3,3) == "C"){
    motif4 <- substr(motif5,1,4)
    cyt <- cyt + 1
  }
  else if(substr(motif5,3,3) == "G"){
    motif4 <- crDNA(substr(motif5,2,5))
    cyt <- cyt + 1
  }
  else 
    next
  
  if(substr(motif4,2,4) != "TCA")
    next
  
  if(motif4 %in% c("TTCA","CTCA"))
    ytca <- ytca + 1
  else if(motif4 %in% c("GTCA","ATCA"))
    rtca <- rtca + 1
  
  print(motif4)
}

res <- rbind(res, data.table("posset"=f,"mut_cytosines"=cyt,"mut_ytca"=ytca,"mut_rtca"=rtca))

}

res$posset <- factor(res$posset,levels=c("human1","human2","monkey1","monkey2","genome1","genome2","isidro","otoole"))
levels(res$posset) <- c("RNAseq:human1","RNAseq:human2","RNAseq:monkey1","RNAseq:monkey2","DNAseq:human1","DNAseq:human2","Isidro et al.","O'Toole et al.")


res[,genome_cyt:=gcyt]
res[,genome_ytca:=gytca]
res[,genome_rtca:=grtca]

res[, ytca_enrich := (mut_ytca/genome_ytca)/(mut_cytosines/genome_cyt)]
res[, rtca_enrich := (mut_rtca/genome_rtca)/(mut_cytosines/genome_cyt)]

ggplot(res, aes(x = rtca_enrich, y = ytca_enrich, label = posset)) +
  geom_point(aes(fill = posset), shape = 21, size = 4, color = "blue") +         
  #geom_text(vjust = 1.5) +  
  theme_classic() +
  coord_fixed() +
  scale_x_continuous(limits = c(0, 5)) + 
  scale_y_continuous(limits = c(0, 5)) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "RTCA enrichment",
    y = "YTCA enrichment",
    fill = "Substitution set"
  )

ggsave("PICS/fig6/fig6_tca_enrich.tiff",dpi=600,units="mm",width=100,height=100,compression = "lzw")

