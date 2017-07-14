require(ggplot2)
require(reshape2)

dat <- read.table("/Users/AbbyARivera/Desktop/MEGL/Corales/Data_Jul_12_17/genus_figure_R/taxonomy_L6_R.txt", sep = "\t", header = T)

colors <- c("#660000",
            "#990000",
            "#FF0000",
            "#FF6666",
            "#FFCCCC",
            "#663300",
            "#994C00",
            "#FF8000",
            "#FFB266",
            "#FFE5CC",
            "#666600",
            "#999900",
            "#FFFF00",
            "#FFFF66",
            "#FFFFCC",
            "#336600",
            "#4C9900",
            "#80FF00",
            "#B2FF66",
            "#E5FFCC",
            "#006633",
            "#00CC66",
            "#33FF99",
            "#99FFCC",
            "#CCFFE5",
            "#006666",
            "#009999",
            "#00FFFF",
            "#99FFFF",
            "#CCFFFF",
            "#003366",
            "#004C99",
            "#0080FF",
            "#66B2FF",
            "#CCE5FF",
            "#000099",
            "#0000FF",
            "#6666FF",
            "#9999FF",
            "#CCCCFF",
            "#330066",
            "#6600CC",
            "#9933FF",
            "#CC99FF",
            "#E5CCFF",
            "#660066",
            "#990099",
            "#FF00FF",
            "#FF66FF",
            "#FFCCFF",
            "#660033",
            "#99004C",
            "#FF007F",
            "#FF66B2",
            "#FFCCE5",
            "#000000",
            "#404040",
            "#808080",
            "#C0C0C0",
            "#FFFFFF",
            "#AFEEEE")

dat$Facet <- factor(dat$Facet, c(""All Taxa", "Non-Ricketsialles Taxa""))

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

ggsave("/Users/AbbyARivera/Desktop/MEGL/Corales/Data_Jul_12_17/genus_figure_R/figure_L6.tiff", height = 10, width = 20, units = "in")
