library(dupRadar)

bam <- "benchmark/large/GM12878_REP1.markdup.sorted.bam"
gtf <- "benchmark/large/genes.gtf"
outdir <- "benchmark/large/dupRadar"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Starting R dupRadar benchmark...\n")
cat(paste("BAM:", bam, "\n"))
cat(paste("GTF:", gtf, "\n"))

# Time the analyzeDuprates step (the main computation)
t0 <- proc.time()
dm <- analyzeDuprates(bam, gtf, stranded = 0, paired = TRUE, threads = 1)
t1 <- proc.time()
cat(sprintf("analyzeDuprates time: %.1f seconds\n", (t1 - t0)["elapsed"]))

# Write the dup matrix
write.table(dm, file = file.path(outdir, "dupMatrix.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Fit the model
fit <- duprateExpFit(DupMatrix = dm)
cat(sprintf("Intercept: %.4f\n", fit$intercept))
cat(sprintf("Slope: %.4f\n", fit$slope))

# Generate plots (PDF and PNG)
pdf(file.path(outdir, "duprateExpDens.pdf"))
duprateExpDensPlot(DupMatrix = dm, fit = fit)
dev.off()
png(file.path(outdir, "duprateExpDens.png"), width = 960, height = 960, res = 144, bg = "white")
duprateExpDensPlot(DupMatrix = dm, fit = fit)
dev.off()

pdf(file.path(outdir, "duprateExpBoxplot.pdf"))
duprateExpBoxplot(DupMatrix = dm)
dev.off()
png(file.path(outdir, "duprateExpBoxplot.png"), width = 960, height = 960, res = 144, bg = "white")
duprateExpBoxplot(DupMatrix = dm)
dev.off()

pdf(file.path(outdir, "expressionHist.pdf"))
expressionHist(DupMatrix = dm)
dev.off()
png(file.path(outdir, "expressionHist.png"), width = 960, height = 960, res = 144, bg = "white")
expressionHist(DupMatrix = dm)
dev.off()

t2 <- proc.time()
cat(sprintf("Total time (including plots): %.1f seconds\n", (t2 - t0)["elapsed"]))
cat("Done.\n")
