setwd("/Your/path/to/QC-metrics/tab/separated/file/")

#Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(ggh4x)
library(patchwork)

#Load data
stats=read.csv("assembly_qualit_complete.txt",sep = "\t",header = T)

#Scale conversioon
stats$Len_mbp=as.numeric(stats$Length)/1000000
stats$n50_mbp=as.numeric(stats$N50)/1000000

#Rearrange the models, flow cells, and Basecallers
stats <- read.csv("assembly_qualit_complete.txt", sep = "\t", header = TRUE) %>%
  mutate(
    Len_mbp   = as.numeric(Length)/1e6,
    n50_mbp=as.numeric(N50)/1000000,
    Basecaller = factor(Basecaller,
                        levels = c("Guppy","Dorado","Rerio","Dorado-5kHz","Illumina","Reference")),
    Flowcell   = factor(Flowcell, levels = c("R9","R10","Illumina","Reference"))
  ) %>%
  separate(Type, into = c("Model2","Assembler"), sep = "-", extra = "merge", fill = "right") %>%
  mutate(
    Assembler = ifelse(is.na(Assembler), Model2, Assembler),
    Assembler = factor(Assembler,
                       levels = c("flye","med1.7","med2.0","srp","lrp","Illumina","Reference")),
    Model     = factor(Model2,
                       levels = c("FAST","HAC","SUP","SUP5kHz","Illumina","Reference"))
  ) %>%
  filter(!is.na(Species) & !is.na(Len_mbp))  # drop any NA rows

# your 20-colour palette
my_palette <- c(
  "#ab5852","#cb9979","brown","#7A95C4",
  "#e6194B","#3cb44b","#bfef45","#aaffc3",
  "#911eb4","#dcbeff","#fabed4","yellow",
  "#ffe119","#000075","#42d4f4","#4363d8",
  "#f58231","orange","#ffd8b1","grey"
)
names(my_palette) <- sort(unique(stats$Species))

length= ggplot(stats,
            aes(x      = Assembler,
                y      = Len_mbp,
                colour = Species,
                shape  = Model)) +
  geom_beeswarm(size = 2) +     # draws the points
  facet_nested(~ Basecaller + Flowcell,
               scales = "free_x",
               space  = "fixed") +
  scale_colour_manual(values = my_palette) +
  labs(x      = NULL,
       y      = "Assembly length (Mbp)",
       colour = "Species",
       shape  = "Model") +
  theme_bw(base_size = 14) +
  theme(
    legend.position     = "none",
    strip.text          = element_text(size = 12, family = "Arial"),
    panel.spacing       = unit(0.5, "lines"),
    panel.grid.major.x  = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank()
  )

#print(length)

ggsave("length_mbp.svg", width = 25, height = 5, units = "in", device = "svg")


############# PSEUDOGENES ###################################################

pg=ggplot(stats,
       aes(x      = Assembler,
           y      = stats$Pseudogenes,
           colour = Species,
           shape  = Model)) +
  geom_beeswarm(size = 2) +     # draws the points
  facet_nested(~ Basecaller + Flowcell,
               scales = "free_x",
               space  = "fixed") +
  scale_colour_manual(values = my_palette) +
  labs(x      = NULL,
       y      = "No. of pseudogenes",
       colour = "Species",
       shape  = "Model") +
  theme_bw(base_size = 14) +
  theme(
    legend.position     = "none",
    #axis.text.x         = element_text(angle = 45, hjust = 1, family = "Arial"),
    strip.text          = element_text(size = 12, family = "Arial"),
    panel.spacing       = unit(0.5, "lines"),
    panel.grid.major.x  = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank()
  )
pg
ggsave("pseudogenes.svg", width = 25, height = 5, units = "in", device = "svg")

###################### Assembly N50 #########################################

n50=ggplot(stats,
       aes(x      = Assembler,
           y      = stats$n50_mbp,
           colour = Species,
           shape  = Model)) +
  geom_beeswarm(size = 2) +     # draws the points
  facet_nested(~ Basecaller + Flowcell,
               scales = "free_x",
               space  = "fixed") +
  scale_colour_manual(values = my_palette) +
  labs(x      = NULL,
       y      = "Assembly N50(Mbp)",
       colour = "Species",
       shape  = "Model") +
  theme_bw(base_size = 14) +
  theme(
    legend.position     = "none",
   # axis.text.x         = element_text(angle = 45, hjust = 1, family = "Arial"),
    strip.text          = element_text(size = 12, family = "Arial"),
    panel.spacing       = unit(0.5, "lines"),
    panel.grid.major.x  = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank()
  )
n50
ggsave("n50.svg", width = 25, height = 5, units = "in", device = "svg")

############# Total rRNA ############################################################
rna=ggplot(stats,
           aes(x      = Assembler,
               y      = stats$Total_rRNA,
               colour = Species,
               shape  = Model)) +
  geom_beeswarm(size = 2) +     # draws the points
  facet_nested(~ Basecaller + Flowcell,
               scales = "free_x",
               space  = "fixed") +
  scale_colour_manual(values = my_palette) +
  labs(x      = NULL,
       y      = "Total_rRNA",
       colour = "Species",
       shape  = "Model") +
  theme_bw(base_size = 14) +
  theme(
    #legend.position     = "none",
    axis.text.x         = element_text(angle = 45, hjust = 1, family = "Arial"),
    strip.text          = element_text(size = 12, family = "Arial"),
    panel.spacing       = unit(0.5, "lines"),
    panel.grid.major.x  = element_blank(),
  )

rna
ggsave("rRNA.svg", width = 25, height = 5, units = "in", device = "svg")

############# CDS ############################################################
cds=ggplot(stats,
           aes(x      = Assembler,
               y      = stats$CDS,
               colour = Species,
               shape  = Model)) +
  geom_beeswarm(size = 2) +     # draws the points
  facet_nested(~ Basecaller + Flowcell,
               scales = "free_x",
               space  = "fixed") +
  scale_colour_manual(values = my_palette) +
  labs(x      = NULL,
       y      = "Number of CDS",
       colour = "Species",
       shape  = "Model") +
  theme_bw(base_size = 14) +
  theme(
    legend.position     = "none",
    #axis.text.x         = element_text(angle = 45, hjust = 1, family = "Arial"),
    strip.text          = element_text(size = 12, family = "Arial"),
    panel.spacing       = unit(0.5, "lines"),
    panel.grid.major.x  = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank()
  )
cds
ggsave("Total_CDS.svg", width = 25, height = 5, units = "in", device = "svg")

############ AVgCDSLength ##########################################################

avgcds=ggplot(stats,
       aes(x      = Assembler,
           y      = stats$Mean_CDS_length,
           colour = Species,
           shape  = Model)) +
  geom_beeswarm(size = 2) +     # draws the points
  facet_nested(~ Basecaller + Flowcell,
               scales = "free_x",
               space  = "fixed") +
  scale_colour_manual(values = my_palette) +
  labs(x      = NULL,
       y      = "Mean CDS length",
       colour = "Species",
       shape  = "Model") +
  theme_bw(base_size = 14) +
  theme(
    legend.position     = "none",
    strip.text          = element_text(size = 12, family = "Arial"),
    panel.spacing       = unit(0.5, "lines"),
    panel.grid.major.x  = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank()
  )
avgcds
ggsave("Avgcds.svg", width = 15, height = 10, units = "in", device = "svg")


#Arrange all the individual plots together
(n50+avgcds+pg+rna)+
  plot_layout(nrow=4,guides = "collect")
 
  
ggsave("Figure3-Assembly_quality.svg", width = 11, height = 12, units = "in", device = "svg")


