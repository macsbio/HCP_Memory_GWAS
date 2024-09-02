# CREATE QQ PLOTS AND GET LAMBA VALUES
library(tidyverse)
library(normentR)
library(reporter)
library(dplyr)

# Set directories
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(DATA.DIR)

iwrd <- read.table("FUMA_input_iwrd.txt", header = T, sep = "\t")
pics <- read.table("FUMA_input_pics.txt.txt", header = T, sep = "\t")
list <- read.table("FUMA_input_list.txt.txt", header = T, sep = "\t")

iwrd$Analysis <- "Penn-Word (EM)"
pics$Analysis <- "Picture-Sequence (EM)"
list$Analysis <- "List-Sorting (WM)"

input <- pics
input <- input[order(input$CHR),]

# Create cumulative positions for ordering x-axis
input$BPcum <- NA
input$BPcum <- as.numeric(input$BPcum)
s <- 0
nbp <- c()
chrnum <- unique(input$CHR)
for (i in 1:length(chrnum)) {
  nbp[i] <- max(input[input$CHR == chrnum[i],]$POS)
  input[input$CHR == chrnum[i],"BPcum"] <- input[input$CHR == chrnum[i],"POS"] + s
  s <- s + nbp[i]
}

#Create centering for each chromosome x-axis location
axis.set <- input %>%
  group_by(CHR) %>%
  summarise(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- (-log10(min(input$P)) + 1)
sig = 5E-8

axis.set$CHR[axis.set$CHR==26] <- "MT"
axis.set$CHR[axis.set$CHR==23] <- "XY"

input$CHR[input$CHR==26] <- "MT"
input$CHR[input$CHR==23] <- "XY"

input$CHR <- factor(input$CHR, levels = c(unique(input$CHR)))

# ggplot
manhattanplot <- ggplot(input, aes(x = BPcum, y = -log10(P), 
                                   colour = CHR)) +
  geom_point(size = 1) +
  geom_point(data = input, aes(x = BPcum, y = -log10(P)), size = 1) +
  scale_x_continuous(expand = c(0.05,0.05), 
                     breaks = axis.set$center, labels = axis.set$CHR,
                     guide = guide_axis(check.overlap = TRUE, n.dodge = 2)) +
  scale_color_manual(values = rep(c("#276EBF", "#183059"), 24)) + 
  geom_point(data = input[input$P<1e-5,], color="#ff9933", size = 1) +
  labs(x = "Chromosome", y = expression(paste("-log"[10],"(", plain(p),")")), title = paste0(unique(input$Analysis))) +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
    axis.text.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
    axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
    axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
    axis.ticks = element_blank()
  )
tiff(filename = paste0("pics.tif"), res = 300, width = 175, height = 85, units = "mm")
print(manhattanplot)  
dev.off()

input <- list
ci <- 0.95
nSNPs <- nrow(input)
plotdata <- data.frame(
  observed = -log10(sort(input$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))
)

plotdata_sub <- plotdata %>%
  filter(expected <= 2) %>%
  sample_frac(0.01)

plotdata_sup <- plotdata %>%
  filter(expected > 2)

plotdata_small <- rbind(plotdata_sub, plotdata_sup)

qqplot <- ggplot(plotdata_small, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey10", alpha = 0.1) +
  geom_point(shape = 1) +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 0.5, alpha = 0.5, color = "red", lineend = "round", linetype = "longdash") +
  labs(x = expression(paste("Expected -log"[10],"(", plain(p),")")),
       y = expression(paste("Observed -log"[10],"(", plain(p),")")),
       title = paste0(unique(input$Analysis))) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line = element_line(colour = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
    axis.text.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
    axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
    axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
    axis.ticks = element_blank(), 
  )
tiff(filename = paste0("list_qq.tif"), res = 300, width = 175, height = 85, units = "mm")
print(qqplot)  
dev.off()
