library(data.table)

sample_list <- c("ERR10513574","ERR10963128",
                 "ERR11026612","ERR11026613","ERR11026615","ERR11026616","ERR11026617","ERR11026625","ERR11026630","ERR11026631","ERR11026635","ERR11026636","ERR11026637","ERR11026638",
                 "ERR11030164","ERR11030165","ERR11030167","ERR11030168","ERR11030169","ERR11030177","ERR11030182","ERR11030183","ERR11030187","ERR11030188","ERR11030189","ERR11030190",
                 #"SRR22450503","SRR22450504","SRR22450505",
                 "SRR22450506","SRR22450507","SRR22450508",
                 #"SRR22450509","SRR22450510","SRR22450511",
                 #"SRR22450515","SRR22450516","SRR22450517",
                 "SRR22450518","SRR22450519","SRR22450520",
                 #"SRR22450521","SRR22450522","SRR22450523",
                 "SRR31266107","SRR31266108","SRR31266109",#"SRR31266110","SRR31266111","SRR31266112","SRR31266113",
                 #"SRR31266114","SRR31266115","SRR31266116","SRR31266117","SRR31266118","SRR31266119","SRR31266120",
                 "SRR19536726","SRR19536727","SRR19536728","SRR19536729",
                 "SRR24877167","SRR24877168")

project_list <- c("PRJEB56841","PRJEB56841",
                  "PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728",
                  "PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728","PRJEB60728",
                  #"PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA906618","PRJNA906618","PRJNA906618",
                  #"PRJNA906618","PRJNA906618","PRJNA906618",
                  #"PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA906618","PRJNA906618","PRJNA906618",
                  #"PRJNA906618","PRJNA906618","PRJNA906618",
                  "PRJNA1183318","PRJNA1183318","PRJNA1183318",#"PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318",
                  #"PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318",
                  "PRJNA845087","PRJNA845087","PRJNA845087","PRJNA845087",
                  "PRJNA981509","PRJNA981509")

projsampl <- data.table("project"=project_list,"sample"=sample_list)

for(p in 1:length(unique(project_list))){
  
  proj = unique(project_list)[p]
  data <- data.table()
  
  sampleDT <- projsampl[project == proj]
  for(s in 1:nrow(sampleDT)){
    
    sample <- sampleDT[s,sample]
    
    dt <- read.csv(paste0("data/",sample,"_bcftools.txt"),sep='\t')
    print(paste0("Project: ",proj,", sample: ",sample,", mutnum=",nrow(dt)))
    dt <- data.table(dt)
    dt <- dt[vaf >= 0.5]
    dt[,"sample":=sample]
    data <- rbind(data,dt)
    
  }
  
  samplesTotal <- length(unique(data$sample))
  dtgrp <- data[,.(mean_vaf=mean(vaf),samples=paste(sample,collapse = ','),samples_cnt=.N,samples_total=samplesTotal,vafs=paste(vaf,collapse = ',')),by=.(pos,ref,alt,motif3,motif21,isAPOBEC,isAPOBECextended,ref_norm,alt_norm,motif3_norm,ref_codon,var_codon,ref_aa,var_aa,aa_position,mutation_category)]
  dtgrp <- dtgrp[order(-mean_vaf)]
  
  write.table(dtgrp,paste0("positions/",proj,".txt"),sep='\t',quote = F,row.names = F)
}