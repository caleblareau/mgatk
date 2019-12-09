#!/usr/bin/env Rscript

library(data.table)
library(ggrepel)
library(cowplot)
library(Hmisc)
library(dplyr)
library(ggplot2)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), action="store",
              help="clip.tsv"))

opt = parse_args(OptionParser(option_list=option_list))
clip <- opt$i
out_base <- basename(clip)

# positions that are frequently clipped
blacklist <- c(1,296,299,300,301,3108,3564,3571,10934,13753,13760,16569)

# read clip file, massage to include position and sort by most frequently clipped
df <- fread(clip,
            col.names=c("pos", "coverage", "clipped"), header = TRUE) %>%
  mutate(blacklist = factor(ifelse(pos %in% blacklist, TRUE,FALSE))) %>%
  mutate(clipped = ifelse(pos == 1 | pos == 16569, 0, clipped)) %>%
  arrange(-clipped) %>%
  mutate(bin = as.numeric(cut2(pos, g=34)))

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
  scale_fill_manual(values=c("black", "red")) +
  geom_label_repel(
    data = df[1:10,],
    aes(label = pos, color = blacklist),
    size = 2,
    show.legend = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  xlab("position") +
  ylab("number of clipped reads")

# plot coverage as a line graph without smoothing
# plot segments of ~500bps with mean coverage over region
p2 <- ggplot() +
  geom_line(df, mapping = aes(x = pos, y = coverage, alpha = 0.2),
            show.legend = F) +
  xlab("position") +
  ylab("read depth") +
  geom_segment(df2, mapping=aes(x=start, y=mean, xend=end, yend=mean)) +
  theme_bw()

#combine and plot
p3 <- plot_grid(p1,p2,nrow = 2)

# save combined plot and txt fike of clip position sorted so that most
# frequently clipped bases are at top
save_plot(paste0(out_base,".pdf"), p3, base_asp = 1.1)
write_tsv(df %>% dplyr::select(pos, coverage, clipped),
          paste0(out_base, ".clip.st.txt"))

