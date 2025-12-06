library(data.table)


all <- read.csv("SecondaryStructure/ss_all.txt",sep='\t')
all <- data.table(all)

allgrp <- all[,.(cnt=.N),by=.(ss_type)]
allgrp[,total := sum(cnt)]
allgrp[,ratio := cnt/total]

for(posset in c("human","monkey","genome")){

 dt <- read.csv(paste0("SecondaryStructure/ss_",posset,".txt"),sep='\t')
 dt <- data.table(dt)

 dtgrp <- dt[,.(cnt=.N),by=.(ss_type)]
 if(nrow(dtgrp[ss_type == "bulge_loop"]) == 0){
   dtgrp <- rbind(dtgrp,data.table("ss_type"="bulge_loop","cnt"=0))
 }
 
 
 dtgrp[,total := sum(cnt)]
 dtgrp[,ratio := cnt/total]

 data <- rbind(dtgrp[,.(ratio,ss_type,"pos_type"="REPs",cnt)],allgrp[,.(ratio,ss_type,"pos_type"="other",cnt)])
 data$ss_type <- as.factor(data$ss_type)
 data$pos_type <- as.factor(data$pos_type)
 data$ss_type <- factor(data$ss_type, levels=c("hairpin_loop","stem","internal_loop","bulge_loop","multifurcation_loop","external_loop"))
 levels(data$ss_type) <- c("Hairpin loop","Stem","Internal loop","Bulge","Multifurcation","External loop")
 
 p <- ggplot(data, aes(x=ss_type, y=ratio, fill=pos_type)) + geom_bar(stat = "identity", position=position_dodge(width = 0.9)) + theme_light() +
   theme(axis.text.x = element_text(angle = 45,vjust=1,hjust=1)) + geom_text(aes(label = cnt), vjust = -0.5, position = position_dodge(width = 0.9), size=2.5) +
   labs(y="Proportion",x="Secondary structure type") + ylim(0,0.7)  + theme(legend.position = "none") +
   scale_fill_manual(values=c("#8EC8E2", "#E47B81"))
 
 ggsave(paste0("PICS/fig4/fig4a_",posset,".tiff"), plot = p, units="mm",height=100, width=150, dpi=300)
 
}
