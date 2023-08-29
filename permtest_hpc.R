# Set working directory
setwd("/mnt/beegfs/dcorujo/REGIONER/")

# Load packages
library(regioneReloaded)

# Load region sets
load("./regionSets/tf_sel.RData")

# Run permtest

perm_res <- crosswisePermTest(Alist = tf_sel,
                              genome = "hg38",
                              evFUN = "numOverlaps",
                              ranFUN = "resampleGenome",
                              count.once = TRUE,
                              ntimes = 1000,
                              mc.cores = 20)


perm_res <- makeCrosswiseMatrix(perm_res)

# Save results
dir.create("hpcResults", showWarnings = FALSE)
save(perm_res, file = "./hpcResults/remap_perm.RData")


pdf("./hpcResults/plotcw.pdf")
plotCrosswiseMatrix(perm_res, lineColor = "black")
dev.off()
