library(data.table)
library(ggplot2)
library(reshape2)
library(extraDistr)
library(VGAM)

refmt12 <- data.table("muttype"=c("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G"))

sample_list <- c("ERR10513574","ERR10963128",
                 "ERR11026612","ERR11026613","ERR11026615","ERR11026616","ERR11026617","ERR11026625","ERR11026630","ERR11026631","ERR11026635","ERR11026636","ERR11026637","ERR11026638",
                 "ERR11030164","ERR11030165","ERR11030167","ERR11030168","ERR11030169","ERR11030177","ERR11030182","ERR11030183","ERR11030187","ERR11030188","ERR11030189","ERR11030190",
                 #"SRR22450503","SRR22450504","SRR22450505",
                 "SRR22450506","SRR22450507","SRR22450508",
                 #"SRR22450509","SRR22450510","SRR22450511",
                 #"SRR22450515","SRR22450516","SRR22450517",
                 "SRR22450518","SRR22450519","SRR22450520",
                 #"SRR22450521","SRR22450522","SRR22450523",
                 "SRR31266107","SRR31266108","SRR31266109" #,
                 #"SRR31266110","SRR31266111","SRR31266112","SRR31266113",
                 #"SRR31266114","SRR31266115","SRR31266116","SRR31266117","SRR31266118","SRR31266119","SRR31266120",
                 #"SRR19536726","SRR19536727","SRR19536728","SRR19536729",
                 #"SRR24877167","SRR24877168"
)

project_list <- c("PRJEB56841","PRJEB56841",
                  "PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728",
                  "PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728",
                  #"PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA906618","PRJNA906618","PRJNA906618",
                  #"PRJNA906618","PRJNA906618","PRJNA906618",
                  #"PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA906618","PRJNA906618","PRJNA906618",
                  #"PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA1183318","PRJNA1183318","PRJNA1183318" #,
                  #"PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318",
                  #"PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318",
                  #"PRJNA845087","PRJNA845087","PRJNA845087","PRJNA845087",
                  #"PRJNA981509","PRJNA981509"
)

alldata <- data.table()

for(s in 1:length(sample_list)){
  
  sample <- sample_list[s]
  project <- project_list[s]
  
  bcftools <- read.csv(paste0("/Users/mar/BIO/PROJECTS/MPOX/paper/data/",sample,"_bcftools.txt"),sep='\t')
  bcftools <- data.table(bcftools)
  bcftools[, isTbefore := as.integer(substr(motif3,1,1) == "T")]
  bcftools[, isAafter := as.integer(substr(motif3,3,3) == "A")]
  
  if(project == "PRJEB60728"){
    bcftools <- bcftools[vaf >= 0.1]
  } else {
    bcftools <- bcftools[vaf >= 0.01]
  }
  
  bcftools[,muttype := paste0(ref,">",alt)]
  
  bcftools_grp <- merge(refmt12,bcftools[,.(cnt=sum(altcnt),apocnt=sum(isAPOBEC*altcnt),isTbefore=sum(isTbefore*altcnt),isAafter=sum(isAafter*altcnt)),by=muttype],by="muttype",all.x=T)
  bcftools_grp[is.na(cnt),cnt:=0]
  bcftools_grp[is.na(apocnt),apocnt:=0]
  bcftools_grp[is.na(isTbefore),isTbefore:=0]
  bcftools_grp[is.na(isAafter),isAafter:=0]
  setnames(bcftools_grp,"cnt",paste0("bcftools_cnt"))
  setnames(bcftools_grp,"apocnt",paste0("bcftools_apocnt"))
  setnames(bcftools_grp,"isTbefore",paste0("bcftools_isTbefore"))
  setnames(bcftools_grp,"isAafter",paste0("bcftools_isAafter"))
  
  refA <- nrow(bcftools[ref == "A" & alt == "G"])
  refT <- nrow(bcftools[ref == "T" & alt == "C"])   
  normA <- nrow(bcftools[ref_norm == "A" & alt_norm == "G"]) 
  normT <- nrow(bcftools[ref_norm == "T" & alt_norm == "C"])   
  
  dt <- data.table("sample"=sample,"refA"=refA,"refT"=refT,"normA"=normA,"normT"=normT)
  dt[, total1 := refA + refT]
  dt[, total2 := normA + normT]
  
  stopifnot(nrow(dt[total1 != total2]) == 0)
  
  dt[, total := total1]
  dt[, total1 := NULL]
  dt[, total2 := NULL]
  
  alldata <- rbind(alldata, dt)
  
}
  

# Estimate alpha and beta of beta-distribution, with unknown p


K <- alldata$refA
N <- alldata$total


nll <- function(par) {
  alpha <- exp(par[1])
  beta  <- exp(par[2])
  
  ll <- lchoose(N, K) + lbeta(K + alpha, N - K + beta) - lbeta(alpha, beta)
  return(-sum(ll))
}

p_init <- (sum(K) + 0.5) / (sum(N) + 1.0)  # stabilized mean
kappa_init <- 10
alpha_init <- p_init * kappa_init
beta_init  <- (1 - p_init) * kappa_init

opt <- optim(
  par = log(c(alpha_init, beta_init)),
  fn = nll,
  method = "BFGS",
  hessian = TRUE
)

alpha_hat <- exp(opt$par[1])
beta_hat  <- exp(opt$par[2])
kappa_hat <- alpha_hat + beta_hat
p_hat     <- alpha_hat / kappa_hat

# Apply beta-binomial model

alldata[, pval := (1- pbetabinom.ab(q = normA - 1, size = total, shape1 = alpha_hat, shape2 = beta_hat))]

write.table(alldata,"/Users/mar/BIO/PROJECTS/MPOX/ADAR/results/adar_genedirection_norm.txt",sep='\t',quote = F,row.names = F)

