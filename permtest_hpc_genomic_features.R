# Set working directory
setwd("/mnt/beegfs/dcorujo/REGIONER/regioneRevision")

# Load packages
library(regioneReloaded)

# Load region sets
encode <- readRDS("regionSets/ENCODE_filtered.RDS")
transcripts <- readRDS("regionSets/txdbSets.RDS")
repeats <- readRDS("regionSets/lineRep.RDS")
universe <- readRDS("regionSets/universe_ENCODE_filtered.RDS")

cores <- 60

# Run permtest

Alist <- encode
Blist <- c(transcripts, "LINE1" = repeats)

permRes <- crosswisePermTest(Alist,
                         Blist,
                         ranFUN = "resampleRegions",
                         ntimes = 10,
                         universe = universe,
                         genome = "hg38",
                         count.once = TRUE,
                         ntimes = 10000,
                         mc.cores = cores)


# Save results
dir.create("hpcResults", showWarnings = FALSE)
saveRDS(permRes, file = "./hpcResults/permResENCODE_vs_features.RDS")
