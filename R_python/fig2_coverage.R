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
                 "SRR31266107","SRR31266108","SRR31266109"#,"SRR31266110","SRR31266111","SRR31266112","SRR31266113",
                 #"SRR31266114","SRR31266115","SRR31266116","SRR31266117","SRR31266118","SRR31266119","SRR31266120"
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
                  "PRJNA1183318","PRJNA1183318","PRJNA1183318"#,"PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318",
                  #"PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318","PRJNA1183318"
                  )

projsampl <- data.table("project"=project_list,"sample"=sample_list)


finalCov <- data.table(pos = 1:197209)
for(p in 1:length(unique(project_list))){
  
  proj = unique(project_list)[p]
  data <- data.table()
  
  sampleDT <- projsampl[project == proj]
  
  projCov <- data.table(pos = 1:197209)
  for(s in 1:nrow(sampleDT)){
    
    sample <- sampleDT[s,sample]
    
    for(map in c("host","hybrid","virus")) {
      
      dt <- read.csv(paste0("ALLBAMS/coverage/",proj,"/",sample,"_VirusFrom",map,"RefRealignedToVirusRef.per-base.bed"),sep='\t',header=F)
      dt <- data.table(dt)
      projCov <- merge(projCov,dt[,.(V3,V4)],by.x="pos",by.y="V3",all.x=T)  
      setnames(projCov,"V4",paste0(sample,"_",map))
    }
  }
  projCov[is.na(projCov)] <- 0
  projCov[, row_mean := rowMeans(.SD), .SDcols = 2:ncol(projCov)]
  #projCov[, log10mean := log10(row_mean+1)]
  projCov[, cov := row_mean / max(row_mean)]
  #projCov[, cov := log10mean / max(log10mean)]
  
  finalCov <- merge(finalCov,projCov[,.(pos,cov)],,by.x="pos",by.y="pos")  
  setnames(finalCov,"cov",proj)
}
    
finalCov[, human_cov := (PRJEB60728+PRJNA906618)/2]
finalCov[, monkey_cov := (PRJEB56841+PRJNA1183318)/2]

write.table(finalCov, "coverage/coverage.txt",quote = F, row.names = F, sep = '\t')



    
    
    
    