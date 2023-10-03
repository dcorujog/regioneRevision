# Set working directory
setwd("/mnt/beegfs/dcorujo/REGIONER/regioneRevision")

# Load packages
library(regioneReloaded)

# Load region sets
encode <- readRDS("regionSets/ENCODE_filtered.RDS")
universe <- readRDS("regionSets/universe_ENCODE_filtered.RDS")

cores <- 50

# Run permtest

permRes <- crosswisePermTest(Alist = Alist,
                             Blist = Alist,
                             genome = "hg38",
                             evFUN = "numOverlaps",
                             ranFUN = "resampleRegions",
                             count.once = TRUE,
                             ntimes = 5000,
                             mc.cores = cores)

# Save results
dir.create("hpcResults", showWarnings = FALSE)
saveRDS(permRes, file = "./hpcResults/permResENCODE_ReR.RDS")
