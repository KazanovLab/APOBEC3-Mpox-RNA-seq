library(data.table)


for(pset in c("human","monkey","genome")){
  
 dt <- read.csv(paste0("positions/",pset,".txt"),sep='\t')
 dt <- data.table(dt)

 dtgrp <- dt[,.(cnt=.N),by=mutation_category]
 dtgrp[, total := sum(cnt)]
 dtgrp[, ratio := cnt/total]
 dtgrp[, type := "APOBEC3-like mutations"]
 if(nrow(dtgrp[mutation_category == "nonsense"]) == 0){
   dtgrp <- rbind(dtgrp,data.table("mutation_category"="nonsense","cnt"=0,"total"=unique(dtgrp$total),"ratio"=0,"type"="APOBEC3-like mutations"))
 }
 
 trg <- read.csv("data/apobec_targets_annotation.txt",sep='\t')
 trg <- data.table(trg)
 trg[, total := sum(cnt)]
 trg[, ratio := cnt/total]
 trg[, type := "TC or GA target site"]
 
 dtplot <- rbind(dtgrp,trg)
 
 p <- ggplot(dtplot, aes(x=factor(mutation_category, level = c('nonsynonymous', 'synonymous', 'nonsense', 'intergenic')), y=ratio, fill=type)) + 
   geom_bar(stat = "identity", position=position_dodge()) +
   scale_fill_manual(values=c("#8EC8E2", "#E47B81"))+
   geom_text(aes(label=cnt), vjust=-0.3, color="black",
             position = position_dodge(0.9))+ ylim(0,0.7) +
   theme_light() +
   labs(x = 'Substitution type',
        y = 'Proportion') + theme(legend.position = "none")
    # +
   #theme(text=element_text(size=24, family="Roboto Condensed"),
   #       plot.title = element_text(hjust = 0.5, size = 24),
   #       axis.text=element_text(size=24))
 
 ggsave(paste0("PICS/fig3/fig3b_",pset,".tiff"), plot = p, units="mm",height=90, width=110, dpi=300)
 
} 
