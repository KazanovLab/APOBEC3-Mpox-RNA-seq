library(Biostrings)
library(data.table)
library(extraDistr)

reverse_complement <- function(seq) {
  # Convert to uppercase
  seq <- toupper(seq)
  
  # Define complement mapping
  complement_map <- c(
    A = "T", T = "A", G = "C", C = "G",
    R = "Y", Y = "R", S = "S", W = "W",
    K = "M", M = "K", B = "V", D = "H",
    H = "D", V = "B", N = "N"
  )
  
  # Split sequence into characters
  chars <- unlist(strsplit(seq, split = ""))
  
  # Replace each base with its complement
  comp <- complement_map[chars]
  
  # Reverse the complemented sequence
  rev_comp <- paste(rev(comp), collapse = "")
  
  return(rev_comp)
}


genome <- readDNAStringSet("NCBIgenome/GCF_014621545.1_ASM1462154v1_genomic.fna")

otherMut <- 0
apobecMut <- 0

substypes <- c("AC"=0, "CC"=0, "GC"=0, "TC"=0, "AT"=0, "CT"=0, "GT"=0, "TT"=0)

for(i in 2:(length(genome$NC_063383.1)-1))
{
 motif3 <- as.character(subseq(genome[[1]],start=i-1,end=i+1))
 nt <- substr(motif3,2,2)  
 if(nt %in% c("A","T"))
 {
   otherMut <- otherMut + 3
 }
 else
 {
   if(substr(motif3,1,2) == "TC" || substr(motif3,2,3) == "GA")
   {
     otherMut <- otherMut + 2
     apobecMut <- apobecMut + 1
   }
 }
 
 if(nt %in% c("A","G")){
   motif2 <- reverse_complement(substr(motif3,2,3))
 } else {
   motif2 <- substr(motif3,1,2)
 }
 
 substypes[motif2] <- substypes[motif2] + 1
 
}

sample_list <- c("ERR10513574","ERR10963128",
                 
                 "ERR11026612","ERR11026613","ERR11026615","ERR11026616","ERR11026617","ERR11026625","ERR11026630","ERR11026631","ERR11026635","ERR11026636","ERR11026637","ERR11026638",
                 "ERR11030164","ERR11030165","ERR11030167","ERR11030168","ERR11030169","ERR11030177","ERR11030182","ERR11030183","ERR11030187","ERR11030188","ERR11030189","ERR11030190",
                 
                 "SRR22450503","SRR22450504","SRR22450505",
                 "SRR22450506","SRR22450507","SRR22450508",
                 "SRR22450509","SRR22450510","SRR22450511","SRR22450515","SRR22450516","SRR22450517",
                 "SRR22450518","SRR22450519","SRR22450520",
                 "SRR22450521","SRR22450522","SRR22450523",
                 
                 "SRR31266107","SRR31266108","SRR31266109" #,
                 #"SRR31266110","SRR31266111","SRR31266112","SRR31266113",
                 #"SRR31266114","SRR31266115","SRR31266116","SRR31266117","SRR31266118","SRR31266119","SRR31266120",
                 
                 ##"SRR19536726","SRR19536727","SRR19536728","SRR19536729",
                 
                 ##"SRR24877167","SRR24877168"
                 )

project_list <- c("PRJEB56841","PRJEB56841",
                  
                  "PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728",
                  "PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728",
                  
                  "PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA906618","PRJNA906618","PRJNA906618",
                  
                  "PRJNA1183318","PRJNA1183318","PRJNA1183318" #,
                  #"PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318",
                  #"PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318",
                  
                  ##"PRJNA845087","PRJNA845087","PRJNA845087","PRJNA845087",
                  
                  ##"PRJNA981509","PRJNA981509"
                  )


projsampldt <- data.table(project=project_list,sample=sample_list)

alldata <- data.table()

phidttypes1 <- data.table()
phidttypes2 <- data.table()
phidtnt1 <- data.table()
phidtnt2 <- data.table()

for(s in 1:length(sample_list)){
  
  sample <- sample_list[s]
  project <- project_list[s]
 
  bcftools <- read.csv(paste0("data/",sample,"_bcftools.txt"),sep='\t')
  bcftools <- data.table(bcftools)
  if(project == "PRJEB60728"){
    bcftools <- bcftools[vaf >= 0.1]
  } else {
    bcftools <- bcftools[vaf >= 0.01]
  }

  dtless05 <- bcftools[vaf < 0.5]
  dtmore05 <- bcftools[vaf >= 0.5]
  
  dttypes1 <- data.table(subtype=names(substypes), cnt=unname(substypes))
  dttypes1[,"A":=0]
  dttypes1[,"C":=0]
  dttypes1[,"G":=0]
  dttypes1[,"T":=0]
  for(i in seq_len(nrow(dtmore05))){
    ref <- dtmore05[i,ref]
    if(ref %in% c("A","G")){
      motif2 <- reverse_complement(substr(dtmore05[i,motif3],2,3))
      alt <- reverse_complement(dtmore05[i,alt])
    } else {
      motif2 <- substr(dtmore05[i,motif3],1,2)
      alt <- dtmore05[i,alt]
    }
    dttypes1[subtype == motif2, (alt) := get(alt) + 1] 
  }
  
  dttypes1_org <- copy(dttypes1)
  dttypes1 <- dttypes1[subtype != "TC"]
  
  dttypes1[, ACGTobserved := A+C+G+T]
  N1 <- sum(dttypes1$ACGTobserved)
  N1total <- sum(dttypes1$cnt)
  dttypes1[, ACGTexpected := (cnt/N1total)*N1]
  dttypes1[, ACGTexpprop := ACGTexpected/N1]
  
  dttypes1[, prjct := project]
  dttypes1[, smpl := sample]
  
  dttypes1_org[, ACGTobserved := A+C+G+T]
  N1 <- sum(dttypes1_org$ACGTobserved)
  N1total <- sum(dttypes1_org$cnt)
  dttypes1_org[, ACGTexpected := (cnt/N1total)*N1]
  dttypes1_org[, ACGTexpprop := ACGTexpected/N1]
  dttypes1_org[, N := N1]
  
  dtnt1 <- data.table()
  for(i in 1:nrow(dttypes1_org))
   for(nt in c("A","C","G","T")){
     if(substr(dttypes1_org[i,subtype],2,2) == nt)
       next
    dtnt1 <- rbind(dtnt1, data.table("subtype"=dttypes1_org[i,subtype],"nt"=nt,"obs"=dttypes1_org[i,get(nt)],"exp"=dttypes1_org[i,ACGTexpprop]/3))
  }
  dtnt1 <- dtnt1[subtype != "TC" | nt != "T"]
  dtnt1[, N := sum(obs)]
  dtnt1[, prjct := project]
  dtnt1[, smpl := sample]
  
  phidttypes1 <- rbind(phidttypes1,dttypes1)
  phidtnt1 <- rbind(phidtnt1,dtnt1)
  
  dttypes2 <- data.table(subtype=names(substypes), cnt=unname(substypes))
  dttypes2[,"A":=0]
  dttypes2[,"C":=0]
  dttypes2[,"G":=0]
  dttypes2[,"T":=0]
  for(i in seq_len(nrow(dtless05))){
    ref <- dtless05[i,ref]
    if(ref %in% c("A","G")){
      motif2 <- reverse_complement(substr(dtless05[i,motif3],2,3))
      alt <- reverse_complement(dtless05[i,alt])
    } else {
      motif2 <- substr(dtless05[i,motif3],1,2)
      alt <- dtless05[i,alt]
    }
    dttypes2[subtype == motif2, (alt) := get(alt) + 1] 
  }
  
  dttypes2_org <- copy(dttypes2)
  dttypes2 <- dttypes2[subtype != "TC"]
  
  dttypes2[, ACGTobserved := A+C+G+T]
  N2 <- sum(dttypes2$ACGTobserved)
  N2total <- sum(dttypes2$cnt)
  dttypes2[, ACGTexpected := (cnt/N2total)*N2]
  dttypes2[, ACGTexpprop := ACGTexpected/N2]
  
  dttypes2[, prjct := project]
  dttypes2[, smpl := sample]
  
  dttypes2_org[, ACGTobserved := A+C+G+T]
  N1 <- sum(dttypes2_org$ACGTobserved)
  N1total <- sum(dttypes2_org$cnt)
  dttypes2_org[, ACGTexpected := (cnt/N1total)*N1]
  dttypes2_org[, ACGTexpprop := ACGTexpected/N1]
  dttypes2_org[, N := N1]
  
  dtnt2 <- data.table()
  for(i in 1:nrow(dttypes2_org))
    for(nt in c("A","C","G","T")){
      if(substr(dttypes2_org[i,subtype],2,2) == nt)
        next
      dtnt2 <- rbind(dtnt2, data.table("subtype"=dttypes2_org[i,subtype],"nt"=nt,"obs"=dttypes2_org[i,get(nt)],"exp"=dttypes2_org[i,ACGTexpprop]/3))
    }
  dtnt2 <- dtnt2[subtype != "TC" | nt != "T"]
  dtnt2[, N := sum(obs)]
  dtnt2[, prjct := project]
  dtnt2[, smpl := sample]
  
  phidttypes2 <- rbind(phidttypes2,dttypes2)
  phidtnt2 <- rbind(phidtnt2,dtnt2)
  
}

phidttypes1[,N:=sum(ACGTobserved),by=smpl]
phidttypes2[,N:=sum(ACGTobserved),by=smpl]


# Beta-binomial log PMF
log_beta_binom_pmf <- function(k, N, alpha, beta) {
  lchoose(N, k) + lbeta(k + alpha, N - k + beta) - lbeta(alpha, beta)
}

# Fit phi (shared dispersion) using background types

fit_phi <- function(k_vec, q_vec, N_vec) {

  # log-likelihood for a given phi
  loglik_phi <- function(phi) {
    # phi must be > 0
    if (phi <= 0) return(-Inf)
    alpha_j <- q_vec * phi
    beta_j  <- (1 - q_vec) * phi
    
    ll_j <- mapply(
      function(kj, N_total, aj, bj) {
        log_beta_binom_pmf(k = kj, N = N_total, alpha = aj, beta = bj)
      },
      kj = k_vec,
      N_total = N_vec,
      aj = alpha_j,
      bj = beta_j
    )
    sum(ll_j)
  }

  neg_loglik_over_logphi <- function(logphi) {
    phi <- exp(logphi)
    return( -loglik_phi(phi) )
  }
  
  opt <- optimize(neg_loglik_over_logphi,
                  interval = c(log(1e-4), log(1e6)))
  
  phi_hat <- exp(opt$minimum)
  phi_hat
}


phidt1 <- data.table()
phidt2 <- data.table()
phidtbynt1 <- data.table()
phidtbynt2 <- data.table()

for(s in 1:length(sample_list)){
  
  sample <- sample_list[s]
  project <- project_list[s]
  
  dts <- phidttypes1[smpl == sample]
  phi <- fit_phi(dts$ACGTobserved,dts$ACGTexpprop,dts$N)
  phidt1 <- rbind(phidt1, data.table("smpl"=sample,"prjct"= project, "phi"=phi))

  dts <- phidttypes2[smpl == sample]
  phi <- fit_phi(dts$ACGTobserved,dts$ACGTexpprop,dts$N)
  phidt2 <- rbind(phidt2, data.table("smpl"=sample,"prjct"= project, "phi"=phi))
  
  dts <- phidtnt1[smpl == sample]
  phi <- fit_phi(dts$obs,dts$exp,dts$N)
  phidtbynt1 <- rbind(phidtbynt1, data.table("smpl"=sample,"prjct"= project, "phi"=phi))
  
  dts <- phidtnt2[smpl == sample]
  phi <- fit_phi(dts$obs,dts$exp,dts$N)
  phidtbynt2 <- rbind(phidtbynt2, data.table("smpl"=sample,"prjct"= project, "phi"=phi))
  
}

#phimed1 <- median(phidt1$phi)
#phimed2 <- median(phidt2$phi)
#phimed1dt <- phidt1[,.("phimed"=median(phi)),by=prjct]
#phimed2dt <- phidt2[,.("phimed"=median(phi)),by=prjct]
phimed1 <- median(phidtbynt1$phi)
phimed2 <- median(phidtbynt2$phi)
phimed1dt <- phidtbynt1[,.("phimed"=median(phi)),by=prjct]
phimed2dt <- phidtbynt2[,.("phimed"=median(phi)),by=prjct]


alldttypes1 <- data.table()
alldttypes2 <- data.table()


for(s in 1:length(sample_list)){
  
  sample <- sample_list[s]
  project <- project_list[s]
 
  bcftools <- read.csv(paste0("data/",sample,"_bcftools.txt"),sep='\t')
  bcftools <- data.table(bcftools)
  bcftools <- bcftools[vaf >= 0.01]
  
  dtless05 <- bcftools[vaf < 0.5]
  dtmore05 <- bcftools[vaf >= 0.5]
 
  print(paste0("dtless05 size:",nrow(dtless05),", apobec: ",sum(dtless05$isAPOBEC)))
  
  # random
  genomeMax <- length(genome$NC_063383.1)-1
  nts <- c("A","G","C","T")
  apoRand <- 0
  otherRand <- 0
  for(i in 1:nrow(dtless05)){
    pos <- sample(2:genomeMax,1)
    motif3 <- as.character(subseq(genome[[1]],start=pos-1,end=pos+1))
    nt <- substr(motif3,2,2)      
    nts2 <- setdiff(nts,nt)
    alt <-sample(nts2,1)
    
   if((substr(motif3,1,2) == "TC" && alt == "T") || (substr(motif3,2,3) == "GA" && alt == "A")){
     apoRand <- apoRand + 1
   } else {
     otherRand <- otherRand + 1
   }
  }
  
  dttypes1 <- data.table(subtype=names(substypes), cnt=unname(substypes))
  dttypes1[,"A":=0]
  dttypes1[,"C":=0]
  dttypes1[,"G":=0]
  dttypes1[,"T":=0]
  for(i in seq_len(nrow(dtmore05))){
   ref <- dtmore05[i,ref]
   if(ref %in% c("A","G")){
     motif2 <- reverse_complement(substr(dtmore05[i,motif3],2,3))
     alt <- reverse_complement(dtmore05[i,alt])
   } else {
     motif2 <- substr(dtmore05[i,motif3],1,2)
     alt <- dtmore05[i,alt]
   }
   dttypes1[subtype == motif2, (alt) := get(alt) + 1] 
  }
  dttypes1[, prjct := project]
  dttypes1[, smpl := sample]
  
  
  # Chi-square 1
  
  # all 3 types of mutations together
  N1 <- nrow(dtmore05)
  if(N1 != 0){
   N1total <- sum(dttypes1$cnt)
   dttypes1[, ACGTobserved := A+C+G+T]
   dttypes1[, ACGTexpected := (cnt/N1total)*N1]
   dttypes1[, ACGTexpprop := ACGTexpected/N1]
   chipval1 <- chisq.test(dttypes1$ACGTobserved, p = dttypes1$ACGTexpprop)$p.value
   
   # each 3 types of mutation separately
   cobs1 <- c(dttypes1$A,dttypes1$G,dttypes1[substr(subtype,2,2) == "C"]$T,dttypes1[substr(subtype,2,2) == "T"]$C)
   cexp1 <- c(dttypes1$ACGTexpprop/3,dttypes1$ACGTexpprop/3,dttypes1$ACGTexpprop/3)
   chipval11 <- chisq.test(cobs1, p = cexp1)$p.value
  } else {
    chipval1 <- -1
    chipval11 <- -1
    dttypes1[, ACGTobserved := NA]
    dttypes1[, ACGTexpected := NA]
    dttypes1[, ACGTexpprop := NA]
  }
  
  alldttypes1 <- rbind(alldttypes1,dttypes1)
  
  # beta-binomial test
  phimed1 <- phimed1dt[prjct == project, phimed] #18.38 
  k <- dttypes1[subtype == "TC",T] # ACGTobserved]
  q <- dttypes1[subtype == "TC",ACGTexpprop/3] # ACGTexpprop]
  alpha <- q * phimed1
  beta <- (1 - q) * phimed1
  bbpval1 <- pbbinom(k, size = N1, alpha = alpha, beta = beta, lower.tail = FALSE)
  
  tmp <- ifelse(!is.na(bbpval1) & bbpval1 != 0 & abs(bbpval1) < 1e-3,
                     formatC(bbpval1, format = "e", digits = 2),   # scientific
                     formatC(bbpval1, format = "f", digits = 2))
  bblab1 <- paste0("p = ",tmp) 
  
  
  
  
  dttypes2 <- data.table(subtype=names(substypes), cnt=unname(substypes))
  dttypes2[,"A":=0]
  dttypes2[,"C":=0]
  dttypes2[,"G":=0]
  dttypes2[,"T":=0]
  for(i in seq_len(nrow(dtless05))){
    ref <- dtless05[i,ref]
    if(ref %in% c("A","G")){
      motif2 <- reverse_complement(substr(dtless05[i,motif3],2,3))
      alt <- reverse_complement(dtless05[i,alt])
    } else {
      motif2 <- substr(dtless05[i,motif3],1,2)
      alt <- dtless05[i,alt]
    }
    dttypes2[subtype == motif2, (alt) := get(alt) + 1] 
  }
  dttypes2[, prjct := project]
  dttypes2[, smpl := sample]
  
  # Chi-square 2
  
  N2 <- nrow(dtless05)
  if(N2 != 0){
   N2total <- sum(dttypes2$cnt)
   dttypes2[, ACGTobserved := A+C+G+T]
   dttypes2[, ACGTexpected := (cnt/N2total)*N2]
   dttypes2[, ACGTexpprop := ACGTexpected/N2]
   chipval2 <- chisq.test(dttypes2$ACGTobserved, p = dttypes2$ACGTexpprop)$p.value

   cobs2 <- c(dttypes2$A,dttypes2$G,dttypes2[substr(subtype,2,2) == "C"]$T,dttypes2[substr(subtype,2,2) == "T"]$C)
   cexp2 <- c(dttypes2$ACGTexpprop/3,dttypes2$ACGTexpprop/3,dttypes2$ACGTexpprop/3)
   chipval22 <- chisq.test(cobs2, p = cexp2)$p.value
  } else {
   chipval2 <- -1
   chipval22 <- -1
   dttypes2[, ACGTobserved := NA]
   dttypes2[, ACGTexpected := NA]
   dttypes2[, ACGTexpprop := NA]
  }
  
  alldttypes2 <- rbind(alldttypes2,dttypes2)
  
  # beta-binomial test
  phimed2 <- phimed2dt[prjct == project, phimed]
  #phimed2 <- 18.38
  k <- dttypes2[subtype == "TC",T] #ACGTobserved]
  q <- dttypes2[subtype == "TC",ACGTexpprop/3] #ACGTexpprop]
  alpha <- q * phimed2
  beta <- (1 - q) * phimed2
  bbpval2 <- pbbinom(k, size = N2, alpha = alpha, beta = beta, lower.tail = FALSE)
  print(paste0("project=",project,", sample=",sample,", k=",k,", N=",N2,", alpha=",alpha,", beta=",beta,", phi=",phi,", q=",q))

  tmp <- ifelse(!is.na(bbpval2) & bbpval2 != 0 & abs(bbpval2) < 1e-3,
                formatC(bbpval2, format = "e", digits = 2),   # scientific
                formatC(bbpval2, format = "f", digits = 2))
  bblab2 <- paste0("p = ",tmp) 
  
  print(sum(dtless05$isAPOBEC)/nrow(dtless05))
  print(apoRand/(apoRand+otherRand))      
  print(apobecMut/(otherMut+apobecMut))
 
  apo <- c(sum(dtless05$isAPOBEC),apobecMut)
  other <- c(nrow(dtless05)-sum(dtless05$isAPOBEC), otherMut)
   
  table <- matrix(c(apo, other), nrow = 2, byrow = TRUE)
  colnames(table) <- c("Success", "Failure")
  rownames(table) <- c("Group 1", "Group 2")
  
  test_result <- chisq.test(table)
  print(paste0("Project=",project,", sample=",sample,", p-value=",test_result$p.value))
  
  apocntmore05 <- sum(dtmore05$isAPOBEC)
  apocntless05 <- sum(dtless05$isAPOBEC)
  pvalmore <- 1.0 - pbinom(apocntmore05-1, size=nrow(dtmore05), prob=(apobecMut/(otherMut+apobecMut)))
  pvalless <- 1.0 - pbinom(apocntless05-1, size=nrow(dtless05), prob=(apobecMut/(otherMut+apobecMut)))
  
  plabmore <- ifelse(!is.na(pvalmore) & pvalmore != 0 & abs(pvalmore) < 1e-3,
                     formatC(pvalmore, format = "e", digits = 2),   # scientific
                     formatC(pvalmore, format = "f", digits = 2))
  plabmore <- paste0("p = ",plabmore) 
  
  plabless <- ifelse(!is.na(pvalless) & pvalmore != 0 & abs(pvalless) < 1e-3,
                     formatC(pvalless, format = "e", digits = 2),   # scientific
                     formatC(pvalless, format = "f", digits = 2))
  plabless <- paste0("p = ",plabless) 
  
  print(paste0("Nmore=",nrow(dtmore05),", apoccnt=",apocntmore05))
  print(paste0("Nless=",nrow(dtless05),", apoccnt=",apocntless05))
  print(pvalmore)
  print(pvalless)
  print(plabmore)
  print(plabless)
    
  dtplot <- data.table("type"=c("SF >= 0.5","SF >= 0.5","SF < 0.5","SF < 0.5"), #"random","random",
                       "condition"=c("apobec","other","apobec","other"), # "apobec","other",
                       "value"=c(sum(dtmore05$isAPOBEC),nrow(dtmore05)-sum(dtmore05$isAPOBEC),sum(dtless05$isAPOBEC),nrow(dtless05)-sum(dtless05$isAPOBEC))) #apoRand,otherRand,
  dtplot$type <- factor(dtplot$type,levels=c("SF >= 0.5","SF < 0.5")) # "random",  
  levels(dtplot$type) <- c("SF >= 0.5","SF < 0.5")
  
  x_more  <- which(levels(dtplot$type) == "SF >= 0.5")
  x_less  <- which(levels(dtplot$type) == "SF < 0.5")
  
  ggplot(dtplot,aes(fill=condition, y=value, x=type)) + scale_fill_manual(values = c("#E33943", "#65A9DC")) +
    geom_bar(position="fill", stat="identity", width=0.5)  + theme_light() + theme(legend.position = "none") + xlab("Substitutions") + ylab("Fraction of substitution types") +
     geom_segment(
      data = subset(dtplot, type == "SF >= 0.5"),
      aes(x = as.numeric(type) - 0.25, xend = as.numeric(type) + 0.25,
          y = 1 - apobecMut/(otherMut+apobecMut), yend = 1 - apobecMut/(otherMut+apobecMut)),
      inherit.aes = FALSE,
      linetype = "dotted", linewidth = 0.8, linetype = "11", color = "gold", lineend = "round"
    ) +
    geom_segment(
      data = subset(dtplot, type == "SF < 0.5"),
      aes(x = as.numeric(type) - 0.25, xend = as.numeric(type) + 0.25,
          y = 1 - apobecMut/(otherMut+apobecMut), yend = 1 - apobecMut/(otherMut+apobecMut)),
      inherit.aes = FALSE,
      linetype = "dotted", linewidth = 0.8, linetype = "11", color = "gold", lineend = "round"
    ) +
    annotate("text",
             x = 1, y = 1.1,
             #label = plabmore,
             #label = paste0(chipval1),#,"|",chipval11),
             label = bblab1,
             hjust = 0.5, vjust = 0,   # right & top align
             size = 3) +
    annotate("text",
             x = 2, y = 1.1,
             #label = plabless,
             #label = paste0(chipval2),#,"|",chipval22),
             label = bblab2,
             hjust = 0.5, vjust = 0,   # right & top align
             size = 3) +
    ggtitle(sample) +
    theme(plot.title = element_text(size = 10,hjust = 0.5))
    
  
#  ggsave(paste0("PICS/fig1/e/",project,"/fig1_enrichment_",project,"_",sample,".tiff"),dpi=600,compression = "lzw",width=90, height=80,units="mm")
  ggsave(paste0("PICS/fig1/e/",project,"/fig1_enrichment_",project,"_",sample,".png"),dpi=300,width=63, height=80,units="mm")
  
  dtplot[, "project" := project]
  dtplot[, "sample" := sample]
  
  alldata <- rbind(alldata, dtplot)
}

write.table(alldata, "PICS/fig1/e/alldata.txt",sep='\t',quote = F, row.names = F)



md1 <- alldttypes1[,.("cnt"=sum(cnt),"A"=sum(A),"C"=sum(C),"G"=sum(G),"T"=sum(T)),by=subtype]
md2 <- alldttypes2[,.("cnt"=sum(cnt),"A"=sum(A),"C"=sum(C),"G"=sum(G),"T"=sum(T)),by=subtype]

md1[, Af := A/cnt]
md1[, Cf := C/cnt]
md1[, Gf := G/cnt]
md1[, Tf := T/cnt]

