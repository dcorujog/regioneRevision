# Set working directory
setwd("/mnt/beegfs/dcorujo/REGIONER/regioneRevision")

# Load packages
library(regioneReloaded)

# Load region sets
encode_filtered <- readRDS("regionSets/ENCODE_filtered.RDS")
encode_universe <- readRDS("regionSets/universe_ENCODE_filtered.RDS")

cores <- 20

# Run permtest

mlzPOL <- multiLocalZscore(encode_filtered$POLR2A_ENCFF159PYD,
                        encode_filtered, 
                        ranFUN = "resampleRegions", 
                        ntimes = 5000,
                        genome = "hg38", 
                        universe = encode_universe,
                        window = 2000,
                        step = 50,
                        mc.cores = cores)

mlzCTCF <- multiLocalZscore(encode_filtered$CTCF_ENCFF199YFA,
                           encode_filtered, 
                           ranFUN = "resampleRegions", 
                           ntimes = 5000,
                           genome = "hg38", 
                           universe = encode_universe,
                           window = 2000,
                           step = 50,
                           mc.cores = cores)

# Save results
dir.create("hpcResults", showWarnings = FALSE)
saveRDS(mlzPOL, file = "./hpcResults/POLR2A_localZS.RDS")
saveRDS(mlzCTCF, file = "./hpcResults/CTCF_localZS.RDS")
