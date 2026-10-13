#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: export_global_union.R <target_gene_search_result.rds> <output.txt>")
results <- readRDS(args[[1]])
genes <- sort(unique(unlist(lapply(results, function(value) value$target_gene_list), use.names = FALSE)))
if (length(genes) == 0 || anyNA(genes)) stop("Invalid response-gene union")
writeLines(genes, args[[2]])
cat(length(genes), "response genes\n")

