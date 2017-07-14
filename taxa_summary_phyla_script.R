require(ggplot2)
require(reshape2)

dat <- read.table("/Users/AbbyARivera/Desktop/MEGL/Corales/Data_Jul_12_17/phyla_figure_R/taxonomy_L2_R.txt", sep = "\t", header = T)

colors <- c("#99CCFF",
            "#990000",
            "#FFCC99",
            "#330066",
            "#0080FF",
            "#FFCCCC",
            "#00CC66",
            "#00CCCC")

dat$Facet <- factor(dat$Facet, c("All Taxa", "Non-Ricketsialles Taxa"))

dat.m <- melt(dat, id = c("Taxonomy", "Facet"))

plot <- ggplot(dat.m, aes(variable, fill = Taxonomy, weight = value))
plot + geom_bar(position = "fill", width = 0.98) + scale_fill_manual(values = colors[1:length(dat$Taxonomy)] ) + facet_grid(. ~ Facet) +
  scale_y_continuous(expand = c(0, 0.005)) + scale_x_discrete(expand = c(0.005, 0)) +
  theme_bw() + theme(axis.text.x = element_text(size =15, angle = 0), 
                     axis.title.x = element_text(size=15), 
                     axis.text.y = element_text(size =15, angle = 0, hjust = 1, vjust = 0.35), 
                     axis.title.y = element_text(size = 15), 
                     strip.text = element_text(size = 15), 
                     legend.title = element_text(size=15), 
                     legend.text = element_text( size = 15)) +
  xlab("Samples") + ylab("Proportion")

ggsave("/Users/AbbyARivera/Desktop/MEGL/Corales/Data_Jul_12_17/phyla_figure_R/figure_L2.tiff", height = 10, width = 20, units = "in")
