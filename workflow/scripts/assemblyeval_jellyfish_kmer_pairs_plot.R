#!/bin/Rscript
library(tidyverse)

tsv <- fs::path_abs(snakemake@input[["tsv"]])
label <- snakemake@input[["tsv"]]

kmer_hist_tbl <- map_df(c(label = tsv),
    read_tsv,
    col_names = c("asm_cov", "read_cov", "count"),
    col_types = "ddd", .id = "type")

kmer_hist_tbl_filtered <- kmer_hist_tbl %>%
    filter(read_cov < 500, asm_cov < 5) %>%
    bind_rows(kmer_hist_tbl %>% filter(read_cov < 500, asm_cov >= 5) %>% group_by(type, read_cov) %>% summarise(asm_cov = 5, count = sum(count), .groups = "drop")) %>%
    mutate(asm_cov = factor(asm_cov, levels = 0:5, labels = c(seq(0, 4), "5+")))

kmer_hist_y_max <- kmer_hist_tbl_filtered %>%
    filter(read_cov > 1) %>%
    with(max(count))
kmer_hist_y_max <- 1e9
png(fs::path_abs(snakemake@output[["png"]]))
print(ggplot(kmer_hist_tbl_filtered, aes(read_cov, count, colour = asm_cov)) +
    facet_grid(cols = vars(type)) +
    geom_line() +
    scale_colour_discrete(name = "Occurences\nin assembly") +
    scale_x_continuous(breaks = seq(0, 100, by = 20)) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, kmer_hist_y_max)) +
    labs(x = "k-mer multiplicity", y = "k-mer count"))
dev.off()


kmer_hist_y_max2 <- kmer_hist_tbl %>%
    filter(asm_cov <= 5, read_cov < 500, read_cov > 1) %>%
    group_by(read_cov, type) %>%
    summarise(count = sum(count)) %>%
    with(max(count))

png(fs::path_abs(snakemake@output[["pngfill"]]))
print(ggplot(kmer_hist_tbl_filtered, aes(read_cov, count, fill = asm_cov)) +
    facet_grid(cols = vars(type)) +
    geom_area(alpha = 0.5, colour = "black", lwd = 0.1) +
    scale_fill_discrete(name = "Occurences\nin assembly") +
    scale_x_continuous(breaks = seq(0, 100, by = 20)) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, kmer_hist_y_max2)) +
    labs(x = "k-mer multiplicity", y = "k-mer count"))
dev.off()
