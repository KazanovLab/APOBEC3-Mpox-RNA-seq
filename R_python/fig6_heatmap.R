library(data.table)
library(pheatmap)
library(reshape2)


data <- read.csv("Specificity/specificity.csv")
data <- data.table(data)

data$PosSet <- factor(data$PosSet,levels=c("human1","human2","monkey1","monkey2","genome1","genome2","isidro","otoole"))
levels(data$PosSet) <- c("RNAseq:human1","RNAseq:human2","RNAseq:monkey1","RNAseq:monkey2","DNAseq:human1","DNAseq:human2","Isidro et al.","O'Toole et al.")


# If needed, use the first column as row names
# rownames(data) <- data[[1]]
# data <- data[, -1]

df_long <- melt(data)
df_long$variable <- factor(df_long$variable, levels=c("A3G","A3F","A3B","A3A"))


# 5. Draw heatmap
ggplot(df_long, aes(x = PosSet, y = variable, fill = value)) +
  geom_tile(color = "gray90") +    # cell borders (light gray gaps)
  geom_text(aes(label = round(value, 2)), size = 3) + 
  scale_fill_gradient(low = "#deebf7", high = "#65A9DC") +  # nice blue gradient
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank()
  ) +
  labs(fill = "Pearson corr.")  # Legend title

ggsave("PICS/fig6/fig6_heatmap.tiff",dpi=600,units="mm",width=130,height=65,compression = "lzw")

