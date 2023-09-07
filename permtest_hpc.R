# Set working directory
setwd("/mnt/beegfs/dcorujo/REGIONER/regioneRevision")

# Load packages
library(regioneReloaded)

# Load region sets
load("./regionSets/ENCODE_sets.RData")

# Run permtest

permResALL <- crosswisePermTest(Alist = allENCODE,
                                   Blist = allENCODE,
                                   genome = "hg38",
                                   evFUN = "numOverlaps",
                                   ranFUN = "circularRandomizeRegions",
                                   count.once = TRUE,
                                   ntimes = 5000,
                                   mc.cores = 20)


permResTF <- crosswisePermTest(Alist = TF,
                             genome = "hg38",
                             evFUN = "numOverlaps",
                             ranFUN = "circularRandomizeRegions",
                             count.once = TRUE,
                             ntimes = 5000,
                             mc.cores = 20)

permResHIST <- crosswisePermTest(Alist = c(POLR2A, TF),
                               Blist = HIST,
                               genome = "hg38",
                               evFUN = "numOverlaps",
                               ranFUN = "circularRandomizeRegions",
                               count.once = TRUE,
                               ntimes = 5000,
                               mc.cores = 20)

# Save results
dir.create("hpcResults", showWarnings = FALSE)
save(permResALL, permResTF, permResHIST, file = "./hpcResults/permResENCODE_circular.RData")
