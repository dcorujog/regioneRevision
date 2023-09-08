# Set working directory
setwd("/mnt/beegfs/dcorujo/REGIONER/regioneRevision")

# Load packages
library(regioneReloaded)

# Load region sets
load("./regionSets/ENCODE_sets.RData")

cores <- 60

# Run permtest

permResALL <- crosswisePermTest(Alist = allENCODE,
                                   Blist = allENCODE,
                                   genome = "hg38",
                                   evFUN = "numOverlaps",
                                   ranFUN = "resampleGenome",
                                   count.once = TRUE,
                                   ntimes = 5000,
                                   mc.cores = cores)

# Save results
dir.create("hpcResults", showWarnings = FALSE)
save(permResALL, permResTF, permResHIST, file = "./hpcResults/permResENCODE.RData")
