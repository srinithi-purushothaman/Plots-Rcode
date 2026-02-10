#RStudio(v4.3.0) is used to generate the figures
#This script is to make the heat map with Mean Kraken counts for Taxa (Reads) and 
#AMR heat map with gene length coverage from CARD (Assembly contigs and Reads)

setwd("/Your/path/to/the/input/")
options(scipen=99999)

library(grid)
library(reshape2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(ggplotify)
library(patchwork)

theme_set(theme_pubr())


#Pheatmap
zymo_kraken=read.csv("heatmap_kraken2_input_Zymo.txt",sep = "\t",header=F)

#Matrix conversion
zymo_matrix=as.matrix(zymo_kraken[-1])
colnames(zymo_matrix)=c("BS(12)","CC(12)","DM(12)","PF(12)")
rownames(zymo_matrix)=c("B.subtilis",
                        "E.faecalis",
                        "L.fermentum",
                        "L.monocytogenes",
                        "S.aureus",
                        "P.aeruginosa",
                        "S.enterica",
                        "E.coli")

#This step is make the Species names in italics - inpsired from stack overflow blog
newnames1 = lapply(
  rownames(zymo_matrix),
  function(x) bquote(italic(.(x))))

#Annotate the rows
my_kit=data.frame(Kit = rep(c("Promega", "Qiagen"), c(2,2)),
                  Lysis=rep(c("Enzymatic","Mechanical"),c(3,1)))
row.names(my_kit) = colnames(zymo_matrix)

my_gram1=data.frame(Gram=rep(c("Gram_positive","Gram_negative"),c(5,3)))
row.names(my_gram1)=rownames(zymo_matrix)

#Colour scheme
annot_colors=list(Kit=c(Promega="navyblue",Qiagen="olivedrab3"))
gram_colour=list(Gram=c(Gram_positive = "#725BA6",Gram_negative = "#FF82C1"))
lysis_colour=list(Lysis=c(Enzymatic="#AF8B87",Mechanical="#FFB50A"))
brks1 = seq(3.1,5.5,length.out=6)
myplate=c(colorRampPalette(colors = c("#FFDFBF","#BD2440"))(n = length(brks1)-1))

#Pheatmap plot
pheatmap(log10(zymo_matrix),cluster_rows = F,
         cluster_cols = F,annotation_col=my_kit,
         main="Zymo Mock Community",annotation_row=my_gram1,
         color=myplate,
         annotation_colors=c(annot_colors,gram_colour,lysis_colour),
         labels_row = as.expression(newnames1),
         angle_col = 45,fontsize=12,fontsize_row = 12,fontsize_col = 12,
         border_color="black",annotation_legend=F,breaks = brks1,
         cellwidth = 20,cellheight = 20,
         filename="Pheatmap_kraken2_zymo_test1.pdf",width=5,height=5,units="in")


eskape_kraken=read.csv("eskape_mean_kraken2.txt",sep="\t",header = F)
eskape_m=as.matrix(eskape_kraken[-1])
colnames(eskape_m)=c("BS(3)","CC(3)","DM(3)","PF(3)")
rownames(eskape_m)=c("E.faecium","S.aureus",
                     "A.baumannii","K.pneumoniae",
                     "P.aeruginosa","E.coli")

newnames2 = lapply(
  rownames(eskape_m),
  function(x) bquote(italic(.(x))))


my_kit2=data.frame(Kit = rep(c("Promega", "Qiagen"), c(2,2)),
                   Lysis=rep(c("Enzymatic","Mechanical"),c(3,1)))
row.names(my_kit2) = colnames(eskape_m)

my_gram2=data.frame(Gram=rep(c("Gram_positive","Gram_negative"),c(2,4)))
row.names(my_gram2)=rownames(eskape_m)

kit_colors=list(Kit=c(Promega="navyblue",Qiagen="olivedrab3"))
gram_colour=list(Gram=c(Gram_positive = "#725BA6",Gram_negative = "#FF82C1"))
lysis_colour=list(Lysis=c(Enzymatic="#AF8B87",Mechanical="#FFB50A"))
brks2=seq(1,5.05,length.out=6)
myplate2=c(colorRampPalette(colors = c("#FFDFBF","#BD2440"))(n = length(brks2)-1))

pheatmap(log10(eskape_m),cluster_rows = F,cluster_cols = F,
         annotation_col=my_kit2,main="ESKAPE Mock",
         annotation_row = my_gram2,
         color=myplate2,
         annotation_colors=c(kit_colors,gram_colour,lysis_colour),
         labels_row = as.expression(newnames2),
         angle_col = 45,fontsize=12,fontsize_row = 12,fontsize_col = 12,
         cellwidth = 20,cellheight = 20,
         border_color="black",annotation_legend=T,breaks = brks2,
         filename="Pheatmap_kraken2_ESKAPE_test1.pdf",width=5,height=5,units = "in")

#zymo_taxa+eskape_taxa

#AMR heatmap

amr_cont=read.csv("amr_contigs_avg_coverage.txt",sep="\t",header = F)
amr_cont_m=as.matrix(amr_cont[-1])
colnames(amr_cont_m)=c("BS(3)","CC(3)","DM(3)","PF(3)")
rownames(amr_cont_m)=c("vanA","mecA","oxa-66","other oxa","kpc-2","ctx-m-65")

my_kit3=data.frame(Kit = rep(c("Promega", "Qiagen"), c(2,2)),Lysis=rep(c("Enzymatic","Mechanical"),c(3,1)))
row.names(my_kit3) = colnames(amr_cont_m)
kit_colors1=list(Kit=c(Promega="navyblue",Qiagen="olivedrab3"))
lysis_colour=list(Lysis=c(Enzymatic="#AF8B87",Mechanical="#FFB50A"))


pheatmap(amr_cont_m,cluster_rows = F,cluster_cols = F,
         annotation_col=my_kit2,main="AMR - Contigs",
         color=colorRampPalette(c("#F7EDA3","#FFAC63", "#FF7D63"))(10), #Please change the colours of your choice
         annotation_colors=c(kit_colors,lysis_colour),
         angle_col = 45,fontsize=8,fontsize_row = 8,
         fontsize_col = 8,annotation_legend=F,
         border_color="black",cellwidth = 20,cellheight = 20,
         filename="Pheatmap_AMR_contigs.pdf",width=5,height=5,units="in")


amr_reads=read.csv("amreads_1.txt",sep="\t",header = F)
amr_reads_m=as.matrix(amr_reads[-1])
colnames(amr_reads_m)=c("BS(3)","CC(3)","DM(3)","PF(3)")
rownames(amr_reads_m)=c("vanA","mecA","oxa-66","other oxa","kpc-2","ctx-m-65","oxa-486","PDC-167")

my_kit3=data.frame(Kit = rep(c("Promega", "Qiagen"), c(2,2)),Lysis=rep(c("Enzymatic","Mechanical"),c(3,1)))
row.names(my_kit3) = colnames(amr_reads_m)
kit_colors1=list(Kit=c(Promega="navyblue",Qiagen="olivedrab3"))
lysis_colour=list(Lysis=c(Enzymatic="#AF8B87",Mechanical="#FFB50A"))



pheatmap(amr_reads_m,cluster_rows = F,cluster_cols = F,
         annotation_col=my_kit3,main="AMR - Reads - CARD",
         color=colorRampPalette(c("#EFF9FD","#A2CEF5", "#5DB1D1"))(10),
         annotation_colors=c(kit_colors1,lysis_colour),
         annotation_legend=F,angle_col = 45,fontsize=8,
         fontsize_row = 8,
         fontsize_col = 8,
         border_color="black",cellwidth = 20,cellheight = 20,
         filename="Pheatmap_AMR_reads_CARD.jpeg",res=300,width=5,height=5,units = "in")

