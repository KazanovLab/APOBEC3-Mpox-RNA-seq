library(circlize)
library(data.table)

gff <- read.csv("NCBIgenome/GCF_014621545.1_ASM1462154v1_genomic.gff", sep='\t', header=F, comment.char = "#")
gff <- data.table(gff)

genes <- gff[V3=="gene"]
genes$V3 <- as.factor(genes$V3)
genes[,strand := ifelse(V7=="+",1,0)]
genes[,strandColor := ifelse(strand==1,"red","green")]
genes[,tmp := tstrsplit(V9,"Name=")[2]]
genes[,geneName := tstrsplit(tmp,';')[1]]

genesUCSC <- read.csv("NCBIgenome/ucsc_early_late.txt",sep='\t')
genesUCSC <- data.table(genesUCSC)
genesUCSC[, chromStartFix := chromStart + 1]
genesUCSCmerge <- merge(genesUCSC,genes,by.x=c("name","chromStartFix"),by.y=c("geneName","V4"),all=T)
genesUCSCmerge[, chromEndFix := ifelse(is.na(V5),chromEnd,V5)]
genesUCSCmerge[,strandFixStr := ifelse(is.na(V7),strand.x,V7)]
genesUCSCmerge[,strandFix := ifelse(strandFixStr=="+",1,0)]
genesUCSCmerge[,stageFix := ifelse(is.na(stage),"unknown",stage)]
gfinal <- genesUCSCmerge[order(chromStartFix),.(name,chromStartFix,chromEndFix,strandFix,stageFix)]
gfinal[stageFix == "early",stageColor := "#E47B81"]
gfinal[stageFix == "intermediate",stageColor := "#60B67B"]
gfinal[stageFix == "late",stageColor := "#8EC8E2"]
gfinal[stageFix == "unknown",stageColor := "grey"]

gbed <- gfinal[,.("genome",chromStartFix,chromEndFix,name)]

#clair3rna <- read.csv("/Users/mar/BIO/PROJECTS/MPOX/paper/data_v1/old2/clair3rna.txt",sep='\t')
#clair3rna <- data.table(clair3rna)
#clair3rna[,clair3rna := 1]
#lowfreq <- read.csv("/Users/mar/BIO/PROJECTS/MPOX/paper/data_v1/old2/lowfreq.txt",sep='\t')
#lowfreq <- data.table(lowfreq)
#lowfreq[,lowfreq := 1]
#bcftools <- read.csv("/Users/mar/BIO/PROJECTS/MPOX/paper/data_v1/old2/bcftools.txt",sep='\t')
#bcftools <- data.table(bcftools)
#bcftools[,bcftools := 1]

#pos <- merge(clair3rna,lowfreq,by=c("V2","V4","V5","motif3","isAPOBEC","isAPOBECextended","vaf","altcnt"))
#pos <- merge(pos,bcftools,by=c("V2","V4","V5","motif3","isAPOBEC","isAPOBECextended","vaf","altcnt"))

for(organism in c("human","monkey")){

#positions <- read.csv("/Users/mar/BIO/PROJECTS/MPOX/paper/positions/human.txt",sep='\t')
positions <- read.csv(paste0("positions/",organism,".txt"),sep='\t')
positions <- data.table(positions)

positions[isAPOBECextended==1, isAPOBECextended_color := "#E33943"]
positions[isAPOBECextended==2, isAPOBECextended_color := "#3AB661"]
positions[isAPOBECextended==0, isAPOBECextended_color := "#65A9DC"]

positions[,antisense := 0]
for(i in 1:nrow(positions)){
  if(positions[i,isAPOBECextended] == 0)
    next
  gt <- genes[positions[i,pos] >= V4 & positions[i,pos] <= V5]
  if(nrow(gt) == 0 | length(unique(gt[,V7])) > 1){
   next
  }
  gstrand <- unique(gt[,V7])
  if((gstrand == "+" & positions[i,ref] == "G") | (gstrand == "-" & positions[i,ref] == "C")){
    positions[i,antisense := 1]
  }
}

positions[, antisense_pch := 16]
positions[antisense==1, antisense_pch := 18]

genome_length <- 197201

# Load coverage
cov <- read.csv("coverage/coverage.txt",sep='\t')
cov <- data.table(cov)
cov <- cov[pos <= genome_length]
cov[,human_cov_norm := human_cov/max(human_cov)]
cov[,monkey_cov_norm := monkey_cov/max(monkey_cov)]

# Initialize circular plot
tiff(paste0("PICS/fig2/fig2_",organism,".tiff"), width = 10000, height = 10000, res = 600, compression = "lzw")
#png("/Users/mar/BIO/PROJECTS/MPOX/paper/PICS/fig2/fig2.png", width = 2000, height = 2000, res = 100, bg = "transparent")
#png("/Users/mar/BIO/PROJECTS/MPOX/paper/PICS/fig2/fig2_monkey.png", width = 2000, height = 2000, res = 100, bg = "transparent")
circos.clear()
circos.par(start.degree = 90)
circos.initialize(sectors="genome",xlim = c(0, genome_length))

circos.track(ylim=c(0,1),bg.border=NA,track.height=0.1)
circos.genomicAxis(
  major.by = 10000,  # Tick interval (1 Mb in this example)
  labels.cex = 1.5     # Adjust the size of coordinate labels
)
circos.rect(xleft=gfinal$chromStartFix,xright=gfinal$chromEndFix,ybottom=gfinal$strandFix*0.5,ytop=0.5+gfinal$strandFix*0.5,sector.index = "genome", border=NA,col=gfinal$stageColor)
circos.genomicLabels(gbed,labels.column=4,cex=1)

circos.track(ylim=c(0,10),bg.border=NA)
circos.trackPoints(sectors=rep("genome",nrow(positions)),x=positions$pos,y=positions$mean_vaf_sample*10,col=positions$isAPOBECextended_color,pch=positions$antisense_pch,cex=2.5)
circos.trackLines(sectors=rep("genome",nrow(positions)),x=positions$pos,y=positions$mean_vaf_sample*10,type='h',col=positions$isAPOBECextended_color,lwd=2)
circos.yaxis()

circos.text(
  x = -0.5,  # X-coordinate outside the circular track
  y = 3,   # Centered vertically
  labels = "SF", 
  facing = "reverse.clockwise", 
  adj = c(1, 0),  # Adjust alignment
  cex = 1.5  # Adjust text size
)

if(organism == "human"){
 cover <- cov$human_cov_norm
} else {   
 cover <- cov$monkey_cov_norm
}

circos.track(ylim=c(0,1),bg.border=NA,track.height=0.1)
#circos.trackLines(sectors=rep("genome",nrow(cov)),x=cov$pos,y=cov$monkey_cov_norm,type='l',area=TRUE,col = "#8EC8E2",border="#65A9DC")
circos.trackLines(sectors=rep("genome",nrow(cov)),x=cov$pos,y=cover,type='l',area=TRUE,col = "#8EC8E2",border="#65A9DC")

dev.off()

}
