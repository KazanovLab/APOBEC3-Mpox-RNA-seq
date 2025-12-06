library(UpSetR)
library(data.table)

human1 <- read.csv("/Users/mar/BIO/PROJECTS/MPOX/paper/positions/PRJEB60728.txt",sep='\t')
human1 <- data.table(human1)
human1 <- human1[samples_cnt > 4]
human2 <- read.csv("/Users/mar/BIO/PROJECTS/MPOX/paper/positions/PRJNA906618.txt",sep='\t')
human2 <- data.table(human2)
monkey1 <- read.csv("/Users/mar/BIO/PROJECTS/MPOX/paper/positions/PRJEB56841.txt",sep='\t')
monkey1 <- data.table(monkey1)
monkey2 <- read.csv("/Users/mar/BIO/PROJECTS/MPOX/paper/positions/PRJNA1183318.txt",sep='\t')
monkey2 <- data.table(monkey2)
monkey2 <- monkey2[samples_cnt > 1]
genome1 <- read.csv("/Users/mar/BIO/PROJECTS/MPOX/paper/positions/PRJNA845087.txt",sep='\t')
genome1 <- data.table(genome1)
genome2 <- read.csv("/Users/mar/BIO/PROJECTS/MPOX/paper/positions/PRJNA981509.txt",sep='\t')
genome2 <- data.table(genome2)

human1[,setname:= "human1"]
human2[,setname:= "human2"]
monkey1[,setname:= "monkey1"]
monkey2[,setname:= "monkey2"]
genome1[,setname:= "genome1"]
genome2[,setname:= "genome2"]

data_human <- rbind(human1,human2)
data_monkey <- rbind(monkey1,monkey2)
data_genomes <- rbind(genome1,genome2)
data <- rbind(human1,human2,monkey1,monkey2,genome1,genome2)

data[,freq := paste0(samples_cnt,"/",samples_total)]
data <- data[order(pos,setname)]

data_human[,freq := paste0(samples_cnt,"/",samples_total)]
data_human <- data_human[order(pos,setname)]

data_monkey[,freq := paste0(samples_cnt,"/",samples_total)]
data_monkey <- data_monkey[order(pos,setname)]

data_genomes[,freq := paste0(samples_cnt,"/",samples_total)]
data_genomes <- data_genomes[order(pos,setname)]

dt <- data[,.(pos,setname)]
dt2 <- data[,.(mean_vaf_sample=mean(mean_vaf),sets=paste(setname,collapse=","),freqs=paste(freq,collapse=",")),by=.(pos,ref,alt,motif3,motif21,isAPOBEC,isAPOBECextended,ref_norm,alt_norm,motif3_norm,ref_codon,var_codon,ref_aa,var_aa,aa_position,mutation_category)]
dt2_human <- data_human[,.(mean_vaf_sample=mean(mean_vaf),sets=paste(setname,collapse=","),freqs=paste(freq,collapse=",")),by=.(pos,ref,alt,motif3,motif21,isAPOBEC,isAPOBECextended,ref_norm,alt_norm,motif3_norm,ref_codon,var_codon,ref_aa,var_aa,aa_position,mutation_category)]
dt2_monkey <- data_monkey[,.(mean_vaf_sample=mean(mean_vaf),sets=paste(setname,collapse=","),freqs=paste(freq,collapse=",")),by=.(pos,ref,alt,motif3,motif21,isAPOBEC,isAPOBECextended,ref_norm,alt_norm,motif3_norm,ref_codon,var_codon,ref_aa,var_aa,aa_position,mutation_category)]
dt2_genomes <- data_genomes[,.(mean_vaf_sample=mean(mean_vaf),sets=paste(setname,collapse=","),freqs=paste(freq,collapse=",")),by=.(pos,ref,alt,motif3,motif21,isAPOBEC,isAPOBECextended,ref_norm,alt_norm,motif3_norm,ref_codon,var_codon,ref_aa,var_aa,aa_position,mutation_category)]

dt_wide <- dcast(
  dt, pos ~ setname,
  fun.aggregate = length, value.var = "setname"
)

#upset(dt_wide,nsets = 6, order.by = "freq")
upset(dt_wide,nsets = 2, order.by = "freq")

write.table(dt2_human,"positions/human.txt",quote=F,row.names = F,sep='\t')
write.table(dt2_monkey,"positions/monkey.txt",quote=F,row.names = F,sep='\t')
write.table(dt2_genomes,"positions/genome.txt",quote=F,row.names = F,sep='\t')

