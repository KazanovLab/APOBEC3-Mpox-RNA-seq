library(data.table)
library(ggplot2)

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

for(s in 1:length(sample_list)){
  
  sample <- sample_list[s]
  project <- project_list[s]
  
  for(caller in c("clair3rna","lofreq","bcftools")){
    
    dt <- read.csv(paste0("data/",sample,"_",caller,".txt"),sep='\t')
    dt <- data.table(dt)
    #dt <- dt[vaf >= 0.1]
    if(project == "PRJEB60728"){
      dt <- dt[vaf >= 0.1]
    } else {
      dt <- dt[vaf >= 0.01]
    }
    
        
    dt <- dt[order(vaf)]
    dt[,xx:=.I]
    dt$isAPOBECextended <- as.factor(dt$isAPOBECextended)
    
    if(caller == "bcftools"){
     print(paste0("Project: ", project, ", sample: ",sample,", nrow: ",nrow(dt)))
    }
    
    ggplot(dt, aes(x=xx,y=vaf,fill=isAPOBECextended)) + geom_bar(stat="identity") + theme_light() + ylab("Variant Allele Frequency") + 
      theme(legend.position = "none") + xlab("Selected genomic positions (ranked by SF)") + ylab("Substitution frequency (SF)") + # xlab("Selected positions #") +
    scale_fill_manual(values = c("#8EC8E2", "#E47B81", "#60B67B"))
    
#    ggsave(paste0("PICS/fig1/d/",project,"/fig1_vaf_",caller,"_",project,"_",sample,".tiff"),dpi=600,units="mm",width=100,height=40,compression = "lzw") 
     ggsave(paste0("PICS/fig1/d/",project,"/fig1_vaf_",caller,"_",project,"_",sample,".png"),dpi=300,units="mm",width=200,height=80) 
    
  }
}


