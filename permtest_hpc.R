# Set working directory
setwd("/mnt/beegfs/dcorujo/REGIONER/regioneRevision")

# Load packages
library(regioneReloaded)

# Load region sets
Alist <- readRDS("./regionSets/ENCODE_sets.RData")

cores <- 60

# Run permtest

permRes <- crosswisePermTest(Alist = Alist,
                                   Blist = Alist,
                                   genome = "hg38",
                                   evFUN = "numOverlaps",
                                   ranFUN = "resampleGenome",
                                   count.once = TRUE,
                                   ntimes = 10000,
                                   mc.cores = cores)

# Save results
dir.create("hpcResults", showWarnings = FALSE)
saveRDS(permResALL, file = "./hpcResults/permResENCODE.RDS")
