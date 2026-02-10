setwd("/your/path/to/input/file/tab/separated/with/ridom/distance/good-target/percentage")

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(ggh4x)

# Read in your data
stats <- read.delim(
  "updated-figur5.txt",
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)

# Tidy & re‐factor
stats <- stats %>%
  # 1) Force Basecaller into your custom order
  mutate(
    Basecaller = factor(Basecaller,
                        levels = c("Guppy","Dorado","Rerio","Dorado-5kHz",
                                   "Illumina","Reference")),
    Flowcell   = factor(Flowcell,
                        levels = c("R9", "R10","Illumina","Reference")),
    Model      = factor(Model,
                        levels = c("FAST","HAC","SUP","SUP5kHz","Illumina","Reference"))
  ) %>%
  # 2) Split your combined “Type” into Model2 + Assembler
  separate(
    Type,
    into  = c("Model2", "Assembler"),
    sep   = "-",
    extra = "merge",
    fill  = "right"
  ) %>%
  # 3) For Illumina & Reference rows, restore Assembler, then re‐factor
  mutate(
    Assembler = ifelse(is.na(Assembler), Model2, Assembler),
    Assembler = factor(Assembler,
                       levels = c("flye","med1.7","med2.0","srp","lrp","Illumina","Reference")),
    
    Model     = factor(Model2,
                       levels = c("FAST","HAC","SUP","SUP5kHz","Illumina","Reference"))
    
  ) %>%
  select(-Model2)


p <- ggplot(
  stats,
  aes(
    x      = Assembler,
    y      = Ridom_distance,
    colour = Species,    # outline colour = species
    shape  = Model       # shape = base-caller model
  )
) +
  
  
  # 1) draw only the <90% points faintly
  geom_quasirandom(
    data        = stats %>% filter(Good_targets_percentage < 90),
    method      = "tukey",
    dodge.width = 0.7,
    size        = 3.2,
    alpha       = 0.4,
    stroke      = 0.6,
    show.legend = FALSE
  ) +
  
  # 2) draw the >90% points on top
  geom_quasirandom(
    data        = stats %>% filter(Good_targets_percentage >= 90),
    method      = "tukey",
    dodge.width = 0.7,
    size        = 4,
    alpha       = 1,
    stroke      = 1.2
  ) +
  
  # facet by Basecaller × Flowcell in your custom order
  facet_nested(~ Basecaller + Flowcell,
               scales = "free_x",
               space = "fixed") +
  
  # x-axis discrete (levels already set)
  scale_x_discrete() +
  
  # apply your custom colour palette to Species
  scale_colour_manual(values = c(
    "#ab5852","#cb9979",
    "#e6194B","#3cb44b","#bfef45","#aaffc3",
    "#911eb4","#dcbeff","#fabed4",
    "yellow","#ffe119",
    "#000075","#42d4f4","#4363d8",
    "#f58231","orange","#ffd8b1"
  )) +
  
  # labels
  labs(
    x      = NULL,
    y      = "Allelic distance from reference",
    colour = "Species",
    shape  = "Model"
  ) +
  # grey “good” band
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 9, ymax = 24),
            fill       = "grey80",
            alpha      = 0.3,
            inherit.aes = FALSE) +
  # dashed threshold lines
  geom_hline(yintercept = c(9, 24),
             linetype   = "dashed",
             colour     = "blue") +
  
  # theme tweaks
  theme_bw(base_size = 14) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1),
    strip.placement    = "outside",
    strip.background   = element_rect(fill = "grey90"),
    panel.grid.major.x = element_blank()
  )

# Render and save
#print(p)
ggsave(
  "Figure5-withlegned.svg",
  p,
  width  = 11,
  height = 8,
  units  = "in"
)

