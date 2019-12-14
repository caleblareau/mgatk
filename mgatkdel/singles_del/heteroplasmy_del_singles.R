#!/usr/bin/env Rscript
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))

compute_heteroplasmy_deletion <- function(file, breakpoint_l, breakpoint_r){
  
  # Define windows based on breakpoints
  left_window <- (breakpoint_l-(66)):(breakpoint_l-28)
  right_window <- (breakpoint_r+28):(breakpoint_r+(66))
  
  # make cell name variable
  cell_id <- gsub(".stats.tsv", "", basename(file))
  
  # read in stats data and rename variables
  renames <- c("start", "end", "lc", "rc", "clip_pos", "read_name")
  read_stats <- fread(input = file) %>% rename_all(~renames)
  read_stats[is.na(read_stats)] <- 0
  
  # keep track of depth
  depth <- nrow(read_stats)
  
  # clipped reads
  read_stats <- read_stats %>%
    filter(start %in% left_window | end %in% right_window) %>%
    filter(clip_pos %in% c((breakpoint_l), (breakpoint_r), 0))
  
  # calculate heteroplasmy; ifelse so each sequenced molecule only counts once
  out <- read_stats %>%
    group_by(read_name) %>%
    summarise(total = sum(lc) + sum(rc)) %>%
    mutate(total = ifelse(total > 0, 1, 0)) %>%
    summarise(heteroplasmy = sum(100*total)/n()) %>%
    mutate(cell = cell_id, depth = depth)
  
  out
}

heteroplasmy_multiple <- function(file, bps1, bps2){
  mapply(function(x,y) compute_heteroplasmy_deletion(file, x, y), bps1, bps2, SIMPLIFY = FALSE) %>%
    bind_rows(.)
}
left <- c(13157, 9232, 8482)
right <- c(15477, 13413, 13445)
heteroplasmy_multiple(file, left, right) 