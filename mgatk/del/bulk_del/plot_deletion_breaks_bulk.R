#!/usr/bin/env Rscript
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(ggrepel)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))

options(warn=-1)

if(FALSE){
  clip_file <- "/Users/clareau/dat/Research/AryeeResearch/lareau_dev/mgatk/tests/barcode/pearson_bci_del6073-13095.clip.tsv"
}

#-----------------
# Command line i/o
#-----------------
args <- commandArgs(trailingOnly = TRUE)
clip_file <- args[1]


# positions that are frequently clipped
blacklist <- c(1,296,299,300,301,3108,3564,3571,10934,13753,13760,16569)

# read clip file, message to include position and sort by most frequently clipped
df <- fread(clip_file, col.names=c("pos", "coverage", "clipped", "SA"), header = TRUE) %>%
  mutate(blacklist = factor(ifelse(pos %in% blacklist, TRUE,FALSE))) %>%
  mutate(clipped = ifelse(pos == 1 | pos == 16569, 0, clipped)) %>%
  arrange(-clipped) %>%
  mutate(bin = as.numeric(cut(pos, breaks=34)))

# summarise per bin for segment mean plot on top of coverage
df2 <- df %>%
  group_by(bin) %>%
  summarise(mean = mean(coverage),
            start = min(pos),
            end = max(pos))

# plot barplot of number of clipped reads per position
# label top 10 most freq. clip positions
# color in red positions which are frequently clipped in all samples
p1 <- ggplot(df, aes(x = pos, y = clipped)) +
  geom_col(aes(fill = blacklist), width = 50,
           show.legend = FALSE) +
  scale_fill_manual(values=c("black", "lightgrey")) +
  geom_label_repel(
    data = df[1:10,],
    aes(label = pos, color = blacklist),
    size = 2,
    show.legend = FALSE) +
  scale_color_manual(values=c("black", "lightgrey")) +
  theme_bw() +
  xlab("position") +
  ylab("# of clipped reads")

# plot coverage as a line graph without smoothing
# plot segments of ~500bps with mean coverage over region
p2 <- ggplot() +
  geom_line(df, mapping = aes(x = pos, y = coverage, alpha = 0.2),
            show.legend = F) +
  xlab("position") +
  ylab("Read depth") +
  geom_segment(df2, mapping=aes(x=start, y=mean, xend=end, yend=mean)) +
  theme_bw()

# find ±2 of all positions (between 0:16569) with SA > 0 

pos <- function(position){
  pos = position + -2:2
  return(data.frame(pos=pos))
}

SA_positions=lapply((df  %>% filter(SA!=0))$pos, pos) %>% bind_rows() %>% filter(pos < 16569 & pos > 1) %>% unique()

# color: 1 if in blacklist, 3 if in SA_positions but not in blacklist, else 2
df3 = df %>% mutate(SA_sup = case_when(pos %in% SA_positions$pos ~ T, TRUE~ F))%>%
  mutate(color= case_when(SA_sup==TRUE & (blacklist==F) ~3, blacklist==T ~1, TRUE ~ 2)) %>% 
  arrange(-SA) %>% arrange(-clipped) %>% ungroup() %>%
  mutate(clipped = ifelse(pos == 1 | pos == 16569, 0, clipped)) 

p3 = ggplot(df3, aes(x = pos, y = clipped)) +
  geom_col(aes(fill = as.factor(color)), width = 50,
           show.legend = FALSE) +
  scale_fill_manual(values=c("lightgrey", "black", "red")) +
  geom_label_repel(
    data = (df3 %>% arrange(-SA))[1:4,],
    aes(label = pos, color = as.factor(color)),
    size = 2,
    show.legend = FALSE) +
  scale_color_manual(values=c("lightgrey", "black", "red")) +
  theme_bw() +
  xlab("position") +
  ylab("# of clipped reads")

# save combined plot 
out_base <- gsub(".clip.tsv","",clip_file)
ggsave(p1, file = paste0(out_base,".clipped_viz.pdf"), height = 4, width = 7)
ggsave(p2, file = paste0(out_base,".coverage_viz.pdf"), height = 4, width = 7)
ggsave(p3, file = paste0(out_base,".SA_coverage_viz.pdf"), height = 4, width = 7)