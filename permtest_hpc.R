# Set working directory
setwd("/mnt/beegfs/dcorujo/REGIONER/regioneRevision")

# Load packages
library(regioneReloaded)

# Load region sets
load("./regionSets/ENCODE_sets.RData")

# Run permtest

permResAll <- crosswisePermTest(Alist = allENCODE,
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

permResPOLR2A <- crosswisePermTest(Alist = POLR2A,
                               Blist = HIST,
                               genome = "hg38",
                               evFUN = "numOverlaps",
                               ranFUN = "circularRandomizeRegions",
                               count.once = TRUE,
                               ntimes = 5000,
                               mc.cores = 20)

# Save results
dir.create("hpcResults", showWarnings = FALSE)
save(permResTF, permResPOLR2A, file = "./hpcResults/permResENCODE_circular.RData")
