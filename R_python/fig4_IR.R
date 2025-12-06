library(data.table)
library(reshape2)
library(ggplot2)

data <- read.csv("DENEK/IR_resutls1.csv")
data <- data.table(data)
setnames(data,"Stem","Repeat")
setnames(data,"Loop","Spacer")
data$PosSet <- factor(data$PosSet,levels=c("total","otoole","isidro","genome2","genome1","monkey2","monkey1","human2","human1"))
levels(data$PosSet) <- c("All MPXV sites","O'Toole et al.","Isidro et al.","DNAseq:human2","DNAseq:human1","RNAseq:monkey2","RNAseq:monkey1","RNAseq:human2","RNAseq:human1")


dt <- melt(data)

ggplot(dt,aes(fill=variable, y=value, x=PosSet)) + 
  geom_bar(position="fill", stat="identity", width=0.8)  + theme_light() + xlab("Genome position set") + ylab("Fraction of substitutions") +
  coord_flip() + labs(fill = "Genome position") + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("#8EC8E2", "#60B67B", "#E47B81"))

ggsave("PICS/fig4/fig4_IR.tiff",dpi=600,units="mm",width=200,height=100,compression = "lzw")


new_colnames <- as.character(data[[1]]) 
dt <- as.data.table(t(data[, -1]))
setnames(dt, new_colnames)


chisq.test(x = dt$`RNAseq:human1`, p = dt$`All MPXV sites` / sum(dt$`All MPXV sites`))
chisq.test(x = dt$`RNAseq:human2`, p = dt$`All MPXV sites` / sum(dt$`All MPXV sites`))
chisq.test(x = dt$`RNAseq:monkey1`, p = dt$`All MPXV sites` / sum(dt$`All MPXV sites`))
chisq.test(x = dt$`RNAseq:monkey2`, p = dt$`All MPXV sites` / sum(dt$`All MPXV sites`))
chisq.test(x = dt$`DNAseq:human1`, p = dt$`All MPXV sites` / sum(dt$`All MPXV sites`))
chisq.test(x = dt$`DNAseq:human2`, p = dt$`All MPXV sites` / sum(dt$`All MPXV sites`))
chisq.test(x = dt$`Isidro et al.`, p = dt$`All MPXV sites` / sum(dt$`All MPXV sites`))
chisq.test(x = dt$`O'Toole et al.`, p = dt$`All MPXV sites` / sum(dt$`All MPXV sites`))

