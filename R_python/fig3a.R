library(data.table)
library(ggplot2)

for(pset in c("human","monkey","genome")){

dt <- read.csv(paste0("positions/",pset,".txt"),sep='\t')
dt <- data.table(dt)
dt <- dt[isAPOBECextended == 1]
dt$isAPOBECextended <- as.factor(dt$isAPOBECextended)

dt1 <- copy(dt)
dt2 <- copy(dt)

dt1[, cat:="Original"]
dt2[, cat:="Normalized"]

dt1[, mtype := paste0(ref,"->",alt)]
dt2[, mtype := paste0(ref_norm,"->",alt_norm)]

dt <- rbind(dt1[,.(mtype,cat,isAPOBECextended)],dt2[,.(mtype,cat,isAPOBECextended)])
dt$cat <- factor(dt$cat,levels=c("Original","Normalized"))
dtgrp <- dt[,.N,by=.(cat,mtype,isAPOBECextended)]

ggplot(dtgrp, aes(x=mtype,y=N, fill=isAPOBECextended)) + geom_bar(position="stack", stat="identity") + facet_wrap(~ cat) + theme_light() + 
  scale_fill_manual(values = c("#E47B81")) +
  theme(strip.text = element_text(size = 10),legend.position="none",axis.title.x = element_blank(),strip.background=element_rect(fill="#8EC8E2")) + ylab("Number of substitutions")

ggsave(paste0("PICS/fig3/fig3a_",pset,".tiff"),dpi=300,units="mm",height=70,width=60)

}
