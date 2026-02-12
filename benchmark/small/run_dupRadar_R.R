#!/usr/bin/env Rscript
# Benchmark R dupRadar on the chr6 test BAM
library(dupRadar)

bam <- "benchmark/small/test.bam"
gtf <- "benchmark/small/chr6.gtf"
outdir <- "benchmark/small/dupRadar"
stranded <- 0   # unstranded
paired <- TRUE
threads <- 1

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("=== Running R dupRadar ===\n")
cat("BAM:", bam, "\n")
cat("GTF:", gtf, "\n")
cat("Paired:", paired, "\n")
cat("Stranded:", stranded, "\n\n")

start_time <- proc.time()

# Run analyzeDuprates (this calls featureCounts 4 times internally)
dm <- analyzeDuprates(bam, gtf, stranded, paired, threads)

elapsed <- proc.time() - start_time
cat("\n=== analyzeDuprates completed in", elapsed["elapsed"], "seconds ===\n\n")

# Save duplication matrix
write.table(dm, file = file.path(outdir, "dupMatrix.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote dupMatrix.txt with", nrow(dm), "genes\n")

# Fit the model
fit <- duprateExpFit(DupMat = dm)
cat("Fit: intercept =", fit$intercept, ", slope =", fit$slope, "\n")

# Save intercept/slope
writeLines(c(
    paste0("intercept\t", fit$intercept),
    paste0("slope\t", fit$slope)
), file.path(outdir, "intercept_slope.txt"))

# Generate plots
start_plot <- proc.time()

pdf(file.path(outdir, "duprateExpDens.pdf"), width = 8, height = 6)
duprateExpDensPlot(DupMat = dm)
dev.off()
png(file.path(outdir, "duprateExpDens.png"), width = 960, height = 960, res = 144, bg = "white")
duprateExpDensPlot(DupMat = dm)
dev.off()

pdf(file.path(outdir, "duprateExpBoxplot.pdf"), width = 9, height = 6)
duprateExpBoxplot(DupMat = dm)
dev.off()
png(file.path(outdir, "duprateExpBoxplot.png"), width = 960, height = 960, res = 144, bg = "white")
duprateExpBoxplot(DupMat = dm)
dev.off()

pdf(file.path(outdir, "expressionHist.pdf"), width = 8, height = 6)
expressionHist(DupMat = dm)
dev.off()
png(file.path(outdir, "expressionHist.png"), width = 960, height = 960, res = 144, bg = "white")
expressionHist(DupMat = dm)
dev.off()

plot_elapsed <- proc.time() - start_plot
cat("Plots generated in", plot_elapsed["elapsed"], "seconds\n")

# Print summary stats
total_elapsed <- proc.time() - start_time
cat("\n=== Summary ===\n")
cat("Total genes:", nrow(dm), "\n")
cat("Genes with reads:", sum(dm$allCounts > 0), "\n")
cat("Genes with dups:", sum(dm$dupRate > 0, na.rm = TRUE), "\n")
cat("Intercept:", fit$intercept, "\n")
cat("Slope:", fit$slope, "\n")
cat("Total time:", total_elapsed["elapsed"], "seconds\n")
