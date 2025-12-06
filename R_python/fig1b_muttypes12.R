library(data.table)
library(ggplot2)
library(reshape2)
library(extraDistr)

refmt12 <- data.table("muttype"=c("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G"))

sample_list <- c("ERR10513574","ERR10963128",
                 "ERR11026612","ERR11026613","ERR11026615","ERR11026616","ERR11026617","ERR11026625","ERR11026630","ERR11026631","ERR11026635","ERR11026636","ERR11026637","ERR11026638",
                 "ERR11030164","ERR11030165","ERR11030167","ERR11030168","ERR11030169","ERR11030177","ERR11030182","ERR11030183","ERR11030187","ERR11030188","ERR11030189","ERR11030190",
                 "SRR22450503","SRR22450504","SRR22450505","SRR22450506","SRR22450507","SRR22450508","SRR22450509","SRR22450510","SRR22450511",
                 "SRR22450515","SRR22450516","SRR22450517","SRR22450518","SRR22450519","SRR22450520","SRR22450521","SRR22450522","SRR22450523",
                 "SRR31266107","SRR31266108","SRR31266109","SRR31266110","SRR31266111","SRR31266112","SRR31266113",
                 "SRR31266114","SRR31266115","SRR31266116","SRR31266117","SRR31266118","SRR31266119","SRR31266120",
                 "SRR19536726","SRR19536727","SRR19536728","SRR19536729",
                 "SRR24877167","SRR24877168")

project_list <- c("PRJEB56841","PRJEB56841",
                  "PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728",
                  "PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728",
                  "PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318",
                  "PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318",
                  "PRJNA845087","PRJNA845087","PRJNA845087","PRJNA845087",
                  "PRJNA981509","PRJNA981509")

alldata <- data.table()

for(s in 1:length(sample_list)){

  sample <- sample_list[s]
  project <- project_list[s]
  
  sbs <- copy(refmt12)
  
clair3rna <- read.csv(paste0("data/",sample,"_clair3rna.txt"),sep='\t')
clair3rna <- data.table(clair3rna)
clair3rna[, isTbefore := as.integer(substr(motif3,1,1) == "T")]
clair3rna[, isAafter := as.integer(substr(motif3,3,3) == "A")]
clair3rna <- clair3rna[vaf >= 0.01]
lowfreq <- read.csv(paste0("data/",sample,"_lofreq.txt"),sep='\t')
lowfreq <- data.table(lowfreq)
lowfreq[, isTbefore := as.integer(substr(motif3,1,1) == "T")]
lowfreq[, isAafter := as.integer(substr(motif3,3,3) == "A")]
lowfreq <- lowfreq[vaf >= 0.01]
bcftools <- read.csv(paste0("data/",sample,"_bcftools.txt"),sep='\t')
bcftools <- data.table(bcftools)
bcftools[, isTbefore := as.integer(substr(motif3,1,1) == "T")]
bcftools[, isAafter := as.integer(substr(motif3,3,3) == "A")]
if(project == "PRJEB60728"){
 bcftools <- bcftools[vaf >= 0.1]
} else {
 bcftools <- bcftools[vaf >= 0.01]
}

clair3rna[,muttype := paste0(ref,">",alt)]
lowfreq[,muttype := paste0(ref,">",alt)]
bcftools[,muttype := paste0(ref,">",alt)]

clair3rna_grp <- merge(refmt12,clair3rna[,.(cnt=sum(altcnt),apocnt=sum(isAPOBEC*altcnt),isTbefore=sum(isTbefore*altcnt),isAafter=sum(isAafter*altcnt)),by=muttype],by="muttype",all.x=T)
clair3rna_grp[is.na(cnt),cnt:=0]
clair3rna_grp[is.na(apocnt),apocnt:=0]
clair3rna_grp[is.na(isTbefore),isTbefore:=0]
clair3rna_grp[is.na(isAafter),isAafter:=0]
setnames(clair3rna_grp,"cnt",paste0("clair3rna_cnt"))
setnames(clair3rna_grp,"apocnt",paste0("clair3rna_apocnt"))
setnames(clair3rna_grp,"isTbefore",paste0("clair3rna_isTbefore"))
setnames(clair3rna_grp,"isAafter",paste0("clair3rna_isAafter"))
clair3rna_grp[,muttype:=NULL]

sbs <- cbind(sbs,clair3rna_grp)

lowfreq_grp <- merge(refmt12,lowfreq[,.(cnt=sum(altcnt),apocnt=sum(isAPOBEC*altcnt),isTbefore=sum(isTbefore*altcnt),isAafter=sum(isAafter*altcnt)),by=muttype],by="muttype",all.x=T)
lowfreq_grp[is.na(cnt),cnt:=0]
lowfreq_grp[is.na(apocnt),apocnt:=0]
lowfreq_grp[is.na(isTbefore),isTbefore:=0]
lowfreq_grp[is.na(isAafter),isAafter:=0]
setnames(lowfreq_grp,"cnt",paste0("lowfreq_cnt"))
setnames(lowfreq_grp,"apocnt",paste0("lowfreq_apocnt"))
setnames(lowfreq_grp,"isTbefore",paste0("lowfreq_isTbefore"))
setnames(lowfreq_grp,"isAafter",paste0("lowfreq_isAafter"))
lowfreq_grp[,muttype:=NULL]

sbs <- cbind(sbs,lowfreq_grp)

bcftools_grp <- merge(refmt12,bcftools[,.(cnt=sum(altcnt),apocnt=sum(isAPOBEC*altcnt),isTbefore=sum(isTbefore*altcnt),isAafter=sum(isAafter*altcnt)),by=muttype],by="muttype",all.x=T)
bcftools_grp[is.na(cnt),cnt:=0]
bcftools_grp[is.na(apocnt),apocnt:=0]
bcftools_grp[is.na(isTbefore),isTbefore:=0]
bcftools_grp[is.na(isAafter),isAafter:=0]
setnames(bcftools_grp,"cnt",paste0("bcftools_cnt"))
setnames(bcftools_grp,"apocnt",paste0("bcftools_apocnt"))
setnames(bcftools_grp,"isTbefore",paste0("bcftools_isTbefore"))
setnames(bcftools_grp,"isAafter",paste0("bcftools_isAafter"))
bcftools_grp[,muttype:=NULL]

sbs <- cbind(sbs,bcftools_grp)

#}

write.table(sbs,paste0("data/",sample,"_sbs12.txt"),sep='\t',quote = F, row.names = F)

sbs[, "project" := project]
sbs[, "sample" := sample]

alldata <- rbind(alldata, sbs)

## plots

sbsplots <- read.csv(paste0("data/",sample,"_sbs12.txt"),sep='\t')
sbsplots <- data.table(sbsplots)
sbsplots$muttype <- as.factor(sbsplots$muttype)


  for(caller in c("clair3rna","lowfreq","bcftools")){
    
   if (caller %in% c("clair3rna","lowfreq"))
       next
    
   dt <- sbsplots[,.(muttype,"cnt"=get(paste0(caller,"_cnt")),
                     "apocnt"=get(paste0(caller,"_apocnt")),
                     "isTbefore"=get(paste0(caller,"_isTbefore")),
                     "isAafter"=get(paste0(caller,"_isAafter")))]  
   dt[, noapocnt := cnt - apocnt]

   N <- dt[muttype == "C>T" | muttype == "G>A", sum(cnt)]
   Nct <- dt[muttype == "C>T", cnt]
   expectedTC <- Nct / 4
   Nga <- dt[muttype == "G>A", cnt]
   expectedAG <- Nga / 4
   apoN <- dt[muttype == "C>T" | muttype == "G>A", sum(apocnt)]
   apoNct <- dt[muttype == "C>T", apocnt]
   apoNga <- dt[muttype == "G>A", apocnt]
   
   # Binomial test
   pval <- 1.0 - pbinom(apoN - 1, round(N), 1/4)
   if(pval < 0.001){
     pval_label <- "p-value < 0.001"
   } else {
     pval_label <- ""
   }
   
   # Chi-squared test
   #obs <- c(apoNct/N,apoNga/N,(Nct-apoNct)/N,(Nga-apoNga)/N)
   #p_exp <- c(1/8,1/8,3/8,3/8)
   #pval_chi <- chisq.test(x = obs, p = p_exp)$p.value
   #pval_chi <- ""
   #if(pval_chi < 0.001)
  #   pval_label_chi <- "p-chi < 0.001"
  # else
   #  pval_label_chi <- ""
   #pval_label_chi <- pval_chi
   
   dt <- dt[cnt != 0]
   if(nrow(dt) != 0){
    dt[, p_hat_before := isTbefore/cnt]
    dt[, p_hat_after := isAafter/cnt]
    dt[muttype != "C>T" & muttype != "G>A", sumcnt := sum(cnt)]
    dt[muttype != "C>T" & muttype != "G>A", w := cnt/sumcnt]
   
    dt_bg_phat <- dt[muttype != "C>T" & muttype != "G>A", .("p_hat"=p_hat_before,w)]
    dt_bg_phat <- rbind(dt_bg_phat, dt[muttype != "C>T" & muttype != "G>A", .("p_hat"=p_hat_after,w)])
    
    qL <- quantile(dt_bg_phat$p_hat, 0.05, na.rm=TRUE)
    qU <- quantile(dt_bg_phat$p_hat, 0.95, na.rm=TRUE)
    p_w <- dt_bg_phat$p_hat
    p_w[p_w < qL] <- qL
    p_w[p_w > qU] <- qU

    mean_p <- mean(p_w)
    #mean_p <- sum(dt_bg_phat$p_hat*dt_bg_phat$w)
    var_p <- var(p_w)
    #var_p <- sum(dt_bg_phat$w * (dt_bg_phat$p_hat - mean_p)^2)
    vmax <- mean_p * (1 - mean_p)
    if(!is.finite(var_p) || var_p <= 0)
      var_p <- min(vmax * 0.5, 0.25)
    if(var_p >= vmax)
      var_p <- (1 - 1e-8) * vmax
   
    k <- mean_p * (1 - mean_p) / var_p - 1
    alpha <- mean_p * k
    beta <- (1 - mean_p) * k
    
    ab <- alpha + beta
    if(ab < 20){
      s <- 20 / ab
      alpha <- alpha * s
      beta <- beta * s
    }
    
    pctga <- pbbinom(round(apoNct+apoNga), size = round(Nct+Nga), alpha = alpha, beta = beta, lower.tail = FALSE)
    if (pctga < 0)
      pctga <- 0.00
  
    pctgalab <- ifelse(!is.na(pctga) & pctga != 0 & abs(pctga) < 1e-3,
                           formatC(pctga, format = "e", digits = 2),   # scientific
                           formatC(pctga, format = "f", digits = 2))
    plab <- paste0("p = ",pctgalab) 
   }
    
   dt_orig <- copy(dt)
   
   dt[, cnt:=NULL]
   dt[, isTbefore:=NULL]
   dt[, isAafter:=NULL]
   dt[, p_hat_before:=NULL]
   dt[, p_hat_after:=NULL]
   dt[, sumcnt:=NULL]
   dt[, w:=NULL]
  
   dtmelt <- melt(dt,id.vars = "muttype")
   dtmelt <- data.table(dtmelt)
   dtmelt$variable <- factor(dtmelt$variable,levels=c("apocnt","noapocnt"))
   
   ggplot(dtmelt,aes(x=muttype,y=value,fill=variable)) + geom_bar(stat='identity') + scale_fill_manual(values = c("#E33943", "#65A9DC")) +
     theme_light() + coord_flip() + ylab("Number of single base substitutions ") + xlab("Substitution type") + theme(legend.position = "none") +
     annotate("text",
              x = Inf, y = Inf,
              label = plab,
              hjust = 1.2, vjust = 3,   # right & top align
              size = 3) +
      geom_segment(
       data = subset(dtmelt, muttype == "C>T"),
       aes(x = as.numeric(muttype) - 0.4, xend = as.numeric(muttype) + 0.4,
           y = Nct - expectedTC, yend = Nct - expectedTC),
       inherit.aes = FALSE,
       linetype = "dotted", linewidth = 0.8, linetype = "11", color = "gold", lineend = "round"
     ) +
     geom_segment(
       data = subset(dtmelt, muttype == "G>A"),
       aes(x = as.numeric(muttype) - 0.4, xend = as.numeric(muttype) + 0.4,
           y = Nga - expectedAG, yend = Nga - expectedAG),
       inherit.aes = FALSE,
       linetype = "dotted", linewidth = 0.8, linetype = "11", color = "gold", lineend = "round"
     ) +
     ggtitle(sample) +
     theme(plot.title = element_text(size = 10,hjust = 0.5))
   
   #ggsave(paste0("PICS/fig1/b/",project,"/fig1_muttypes_",caller,"_",project,"_",sample,".tiff"),dpi=600,compression = "lzw",width=80, height=100,units="mm") 
   ggsave(paste0("PICS/fig1/b/",project,"/fig1_muttypes_",caller,"_",project,"_",sample,".png"), dpi=300, width=80, height=100, units="mm") 
  }
}

write.table(alldata, "PICS/fig1/b/alldata.txt",sep='\t',quote = F, row.names = F)

