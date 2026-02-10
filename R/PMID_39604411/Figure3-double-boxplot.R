#RStudio(v4.3.0) is used to generate the figures

#This Rscript is used to plot double boxplot for ONT long-read sequenced (metagenomics) read and assembly quality

#This script is inspired from the stack exchange - https://stackoverflow.com/questions/46068074/double-box-plots-in-ggplot2

setwd("Your/path/to/the/tab/separated/quality/metrics/file")

#this is to set the font style
library(extrafont) 
library(showtext)
font_add(family = "Arial",regular = "Arial.ttf")
loadfonts(device = "postscript")

#Load the necessary packages
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(plotly)
library(gridExtra)
library(patchwork)
library(svglite)

theme_set(theme_pubr())

#Load your data for read quality from NanoPlot - Example for ZymoMock community
zymo_qc=read.csv("Zymo_QC_stats.txt",sep = "\t",header = T)

#Scale conversion
zymo_qc$mean_read_length_kbp=zymo_qc$Mean.read.length/1000
zymo_qc$Total_bases_Gbp=zymo_qc$Total.bases/1000000000
zymo_qc$readn50_kbp=zymo_qc$read.n50/1000

# To plot variance for the boxplot
plot.x=ggplot(zymo_qc)+
  geom_boxplot(aes(x=Kit,y=Total_bases_Gbp))

plot.y=ggplot(zymo_qc)+
  geom_boxplot(aes(x=Kit,y=readn50_kbp))


grid.arrange(plot.x, plot.y, ncol=2)

plot.x <- layer_data(plot.x)[,1:6]
plot.y <- layer_data(plot.y)[,1:6]

colnames(plot.x) <- paste0("x.", gsub("y", "", colnames(plot.x)))
colnames(plot.y) <- paste0("y.", gsub("y", "", colnames(plot.y)))
df <- cbind(plot.x, plot.y); rm(plot.x, plot.y)
df$category <- sort(unique(zymo_qc$Kit))

df.outliers <- df %>%
  select(category, x.middle, x.outliers, y.middle, y.outliers) %>%
  data.table::data.table()
df.outliers <- df.outliers[, list(x.outliers = unlist(x.outliers), y.outliers = unlist(y.outliers)), 
                           by = list(category, x.middle, y.middle)]

#Add the calculated variance layers to the plot
zymo=ggplot(df, aes(fill = category, color = category)) +
  geom_rect(aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper), alpha = 0.3) +
  geom_rect(aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper), 
            color = "black", fill = NA) +
  geom_segment(aes(x = x.min, y = y.middle, xend = x.max, yend = y.middle)) + 
  geom_segment(aes(x = x.min, y = y.lower, xend = x.min, yend = y.upper)) + 
  geom_segment(aes(x = x.max, y = y.lower, xend = x.max, yend = y.upper))+
  geom_segment(aes(x = x.middle, y = y.min, xend = x.middle, yend = y.max)) + 
  geom_segment(aes(x = x.lower, y = y.min, xend = x.upper, yend = y.min)) + 
  geom_segment(aes(x = x.lower, y = y.max, xend = x.upper, yend = y.max))+
  labs(x="",y="",subtitle ="Zymo Mock Community (n=12)",tag=expression(bold("a")))+
  theme_bw()+
  theme(axis.text = element_text(size = 12,colour = "black",family = "Arial"))+
  theme(text=element_text(size = 12,colour = "black",family = "Arial"))+
  theme(legend.position = "none")+
  theme(plot.subtitle = element_text(hjust = 0.5,vjust = 0.5,size = 12,colour = "black",family = "Arial"))+
  geom_point(data = zymo_qc,aes(Total_bases_Gbp,readn50_kbp,
                                group=Kit,colour=Kit),
             inherit.aes = F,alpha=0.5)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 18))

#repeat the above steps for other two sample sets (ESKAPE and Swab samples)

#Arrange the plots by the sample sets

qcrlb=ggarrange(zymo+theme(axis.ticks.y = element_blank(),
                           plot.margin = margin(r = 1))
                ,eskape+theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank(),
                              plot.margin = margin(r = 1,l=1))
                ,swab+theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            plot.margin = margin(r=2,l = 1)),nrow=1,align = "h")

#Add annotations to the plot for publication
qcrlb_anot=annotate_figure(qcrlb,left=text_grob("Read n50 (Kbp)",rot = 90, vjust = 0.5,hjust = 0.5,size = 12),
                           top=text_grob("Quality Control",face="bold",color = "black",size=14,margin(0,0,-1.0,0),
                                         hjust = 0.5,vjust=0.5),
                           bottom = text_grob("Total bases generated (Gbp)",vjust=-1.0,size = 12))+
  geom_text(family="Arial",inherit.aes = T)+ggtitle("Quality control")


# Similarly load your assembly quality metrics calculated from metaQUAST

#Load the data
zassem=read.csv("Zymo_mock_assembly_stats.txt",sep = "\t",header = T)

#Scale conversion
zassem$N50_mbp=zassem$N50/1000000
zassem$Total_length_mbp=zassem$Total.length/1000000

#Calcualte variance
plot.x=ggplot(zassem)+
  geom_boxplot(aes(x=Kit,y=Total_length_mbp))

plot.y=ggplot(zassem)+
  geom_boxplot(aes(x=Kit,y=Number.of.contigs))

grid.arrange(plot.x,plot.y,ncol=2)


plot.x <- layer_data(plot.x)[,1:6]
plot.y <- layer_data(plot.y)[,1:6]

colnames(plot.x) <- paste0("x.", gsub("y", "", colnames(plot.x)))
colnames(plot.y) <- paste0("y.", gsub("y", "", colnames(plot.y)))
df <- cbind(plot.x, plot.y); rm(plot.x, plot.y)
df$category <- sort(unique(zassem$Kit))

df.outliers <- df %>%
  select(category, x.middle, x.outliers, y.middle, y.outliers) %>%
  data.table::data.table()
df.outliers <- df.outliers[, list(x.outliers = unlist(x.outliers), y.outliers = unlist(y.outliers)), 
                           by = list(category, x.middle, y.middle)]
#Add layers to the plot
zas=ggplot(df, aes(fill = category, color = category)) +
  geom_rect(aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper), alpha = 0.3) +
  geom_rect(aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper), 
            color = "black", fill = NA) +
  geom_segment(aes(x = x.min, y = y.middle, xend = x.max, yend = y.middle)) + 
  geom_segment(aes(x = x.min, y = y.lower, xend = x.min, yend = y.upper)) + 
  geom_segment(aes(x = x.max, y = y.lower, xend = x.max, yend = y.upper))+
  geom_segment(aes(x = x.middle, y = y.min, xend = x.middle, yend = y.max)) + 
  geom_segment(aes(x = x.lower, y = y.min, xend = x.upper, yend = y.min)) + 
  geom_segment(aes(x = x.lower, y = y.max, xend = x.upper, yend = y.max))+
  labs(x="",y="",subtitle = "Zymo Mock Community",tag=expression(bold("d")))+
  theme_bw()+
  theme(axis.text = element_text(size = 12,colour="black",family = "Arial"))+
  theme(text=element_text(size = 12,colour = "black",family = "Arial"))+
  theme(legend.position = "none")+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 40))+
  theme(plot.subtitle=element_text(hjust = 0.5,vjust = 0.5,size = 12,colour = "black",family = "Arial"))+
  geom_point(data = zassem,aes(Total_length_mbp,Number.of.contigs,
                               group=Kit,colour=Kit),
             inherit.aes = F,alpha=0.5)

#Arrange by sampples
as=ggarrange(zas+theme(axis.ticks.y = element_blank(),
                       plot.margin = margin(r = 12,l=-10)),
             eas+theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.title.y = element_blank(),
                       plot.margin = margin(r = 14,l=-2)),
             sas+theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.title.y = element_blank(),
                       plot.margin = margin(r = 8,l=-4)),nrow=1,align = "h")

#Annotate the plot for publication
as_anot= annotate_figure(as,left=text_grob("Number of assembled contigs",rot = 90, vjust = 0.5,hjust = 0.5,size=12),
                         top=text_grob("Assembly",face="bold",color = "black",size=14,margin(0,0,-1.0,0),
                                       hjust = 0.5,vjust = 0.5),
                         bottom = text_grob("Total length of assembled contigs (Mbp)",vjust=-1.0,size=12))+
  geom_text(family="Arial",inherit.aes = T)


#Final plot arranging both qc and assmebly together
image=ggarrange(qcrlb_anot,as_anot,nrow = 2,align = "h")
ggsave(file = "Figure3_read_and_assembly_quality.svg",width=10,height = 8)