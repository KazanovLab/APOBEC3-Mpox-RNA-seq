library(data.table)
library(plotly)

refmt12 <- data.table("muttype"=c("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G"))

  human <- read.csv(paste0("positions/human.txt"),sep='\t')
  human <- data.table(human)
  monkey <- read.csv(paste0("positions/monkey.txt"),sep='\t')
  monkey <- data.table(monkey)
  genome <- read.csv(paste0("positions/genome.txt"),sep='\t')
  genome <- data.table(genome)
  
  human[,muttype := paste0(ref,">",alt)]
  human <- human[muttype %in% c("C>T","G>A")]
  monkey[,muttype := paste0(ref,">",alt)]
  monkey <- monkey[muttype %in% c("C>T","G>A")]
  genome[,muttype := paste0(ref,">",alt)]
  genome <- genome[muttype %in% c("C>T","G>A")]
  
  human_grp <- human[,.(cnt=.N,apocnt=sum(isAPOBEC*.N)),by=.(motif3,isAPOBEC)]
  monkey_grp <- monkey[,.(cnt=sum(.N),apocnt=sum(isAPOBEC*.N)),by=.(motif3,isAPOBEC)]
  genome_grp <- genome[,.(cnt=sum(.N),apocnt=sum(isAPOBEC*.N)),by=.(motif3,isAPOBEC)]
  
  donutref <- data.table(motif3 = c("TCA","TCC","TCG","TCT",
                                    "AGA","CGA","GGA","TGA",
                                    "ACA","ACC","ACG","ACT",
                                    "CCA","CCC","CCG","CCT",
                                    "GCA","GCC","GCG","GCT",
                                    "AGC","CGC","GGC","TGC",
                                    "AGG","CGG","GGG","TGG",
                                    "AGT","CGT","GGT","TGT"
 # ), color = c('deeppink4', 'steelblue4','darkgreen','cyan4','orange4','dodgerblue4','red4','yellow4',
#                'aquamarine1','deepskyblue','green','tomato','yellow','purple','orchid','turquoise',
#                'gold','darkorange','lightsalmon','lightblue','violet','cyan','coral','hotpink',
#                'chartreuse','darkseagreen1','firebrick1','khaki','sienna1','maroon1','peachpuff','plum'))

 ), color = c('#E47B81', '#E47B81','#E47B81','#E47B81','#E47B81','#E47B81','#E47B81','#E47B81',
                '#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2',
                '#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2',
                '#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2','#8EC8E2'))

  donutref$motif3 <- as.factor(donutref$motif3)
  donutref$color <- as.factor(donutref$color)
  
  #plotdt <- merge(donutref,monkey_grp,by="motif3",all.x = T)
  #plotdt <- merge(donutref,genome_grp,by="motif3",all.x = T)
  plotdt <- merge(donutref,human_grp,by="motif3",all.x = T)
  plotdt[is.na(cnt),cnt:=0]  
  
  plot_ly(plotdt, labels=~motif3, values=~cnt, type='pie', hole=0.5, textinfo='label+percent', texttemplate = "%{label}: %{percent:.2%}",
          marker=list(colors=plotdt$color,line = list(color = 'white', width = 1) ), textposition = 'inside', showlegend = FALSE) %>%
    config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "myplot",
        width = 400,
        height = 400
      )
    )
  
  
  