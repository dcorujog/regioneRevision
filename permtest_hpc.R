# Set working directory
setwd("/mnt/beegfs/dcorujo/REGIONER/regioneRevision")

# Load packages
library(regioneReloaded)

# Load region sets
load("./regionSets/ENCODE_sets.RData")

# Run permtest

permRes <- crosswisePermTest(Alist = rsList,
                             genome = "hg38",
                             evFUN = "numOverlaps",
                             ranFUN = "resampleGenome",
                             count.once = TRUE,
                             ntimes = 1000,
                             mc.cores = 20)

# Save results
dir.create("hpcResults", showWarnings = FALSE)
save(perm_res, file = "./hpcResults/permResENCODE.RData")
