#!/usr/bin/env Rscript
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))

args <- commandArgs(trailingOnly = TRUE)

in_file = args[1]
out_file = args[2]
read_length = args[3]
window_far = args[4]
window_near = args[5]
left_coordinates = args[6]
right_coordinates = args[7]


if(FALSE){
  in_file = "/Users/clareau/dat/Research/AryeeResearch/lareau_dev/mgatk/tests/mgatk_out/temp/del/CACCACTAGGAGGCGA-1.qc.readStats.tsv"
  out_file = "/Users/clareau/dat/Research/AryeeResearch/lareau_dev/mgatk/tests/mgatk_out/temp/del/CACCACTAGGAGGCGA-1.qc.readStats.tsv"
  read_length = 72
  window_far = 6
  window_near = 28
  left_coordinates = "6073"
  right_coordinates = "13095"
  left_coordinates = "6073,4000"
  right_coordinates = "13095,8000"
}


process_coordinate <- function(string){
  as.numeric(strsplit(string, ",")[[1]])
}
read_length <- as.numeric(read_length)
window_far <- as.numeric(window_far)
window_near <- as.numeric(window_near)
left <- process_coordinate(left_coordinates)
right <- process_coordinate(right_coordinates)


compute_heteroplasmy_deletion <- function(file, breakpoint_l, breakpoint_r){
  
  # Define windows based on breakpoints
  left_window <- (breakpoint_l-(read_length-window_far)):(breakpoint_l-window_near)
  right_window <- (breakpoint_r+window_near):(breakpoint_r+(read_length-window_far))
  
  # make cell name variable
  cell_id <- gsub(".readStats.tsv", "", basename(file))
  cell_id <- gsub(".qc$", "", cell_id)
  
  # read in stats data and rename variables
  renames <- c("start", "end", "lc", "rc", "clip_pos", "read_name", "barcode", "n_clipped")
  read_stats <- fread(input = file) %>% rename_all(~renames)
  read_stats[is.na(read_stats)] <- 0
  
  # clipped reads
  read_stats <- read_stats %>%
    filter(start %in% left_window | end %in% right_window) %>%
    filter(clip_pos %in% c((breakpoint_l), (breakpoint_r), 0))
  
  # calculate heteroplasmy; ifelse so each sequenced molecule only counts once
  out <- read_stats %>%
    group_by(read_name) %>%
    summarise(total = sum(lc) + sum(rc)) %>%
    mutate(total = ifelse(total > 0, 1, 0)) %>%
    summarise(cell = cell_id, heteroplasmy = round(sum(100*total)/n(),2),
              reads_del = sum(total), reads_wt = n() - sum(total), reads_all = n(),
              total_clipped_bases = NA) %>%
    mutate(deletion = paste0("del", as.character(breakpoint_l), "-", as.character(breakpoint_r)), what = "naive")  

  # reads with barcode
  bar_reads <- read_stats %>% filter(barcode) %>%
    pull(read_name)
  
  # calculate heteroplasmy with slight modifications
  out_quant <- read_stats %>%
    filter(!(read_name %in% bar_reads)) %>%
    group_by(read_name) %>%
    mutate(count = n(),
           total = sum(lc) + sum(rc)) %>%
    ungroup() %>%
    filter((count == total | total == 0)) %>%
    dplyr::select(-count, -total) %>%
    group_by(read_name) %>%
    summarise(total = sum(lc) + sum(rc),
              sum_lc = sum(lc), sum_rc = sum(rc),
              n_clipped = max(n_clipped)) %>%
    mutate(total = ifelse(total > 0, 1, 0)) %>%
    summarise(cell = cell_id, heteroplasmy = round(sum(100*total)/n(),2),
              reads_del = sum(total), reads_wt = n() - sum(total), reads_all = n(), 
              total_clipped_bases = sum(n_clipped)) %>%
    mutate(deletion = paste0("del", as.character(breakpoint_l), "-", as.character(breakpoint_r)),  what = "improved")
  
  data.frame(rbind(out_quant, out))

}

heteroplasmy_multiple <- function(in_file, bps1, bps2){
  lapply(1:min(length(bps1), length(bps2)), function(i){
    compute_heteroplasmy_deletion(in_file, bps1[i], bps2[i])
  }) %>% rbindlist() %>% data.frame()
}

het_mult_out <- heteroplasmy_multiple(in_file, left, right) 

# Process and export
write.table(
  het_mult_out, file = out_file,
  row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
)
