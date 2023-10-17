# Set working directory
setwd("/mnt/beegfs/dcorujo/REGIONER/regioneRevision")

# Load packages
library(regioneReloaded)

# Load region sets
setX <- readRDS("regionSets/ENCODE_filterScore.RDS")

# Generate universe/genome
gR <- GRanges() # gR is the genome calculated for the regions of newSets
for(i in 1:length(setX)){ 
  gR <- mergeRegions(gR, setX[[i]])
}

# Function to run permTest in lapply

nzPerm <- function(A, sets, ranfun = "randomizeRegions", nt = 5000, uni = NULL) {
  
  # Vector of n of regions to use
  vecLength <- round(seq(100, length(A), length.out = 10))
  # vecLength <- round(quantile(seq(1, length(A))))

  # Generate a list of 10 subsamples of increasing size from a region set

  listA<-list(A01 = A[sample(vecLength[1])],
              A02 = A[sample(vecLength[2])],
              A03 = A[sample(vecLength[3])],
              A04 = A[sample(vecLength[4])],
              A05 = A[sample(vecLength[5])],
              A06 = A[sample(vecLength[6])],
              A07 = A[sample(vecLength[7])],
              A08 = A[sample(vecLength[8])],
              A09 = A[sample(vecLength[9])],
              A10 = A[sample(vecLength[10])])
  
  # Run permtest
  permRes <- crosswisePermTest(Alist = listA, 
                               Blist = sets, 
                               ntimes = nt, 
                               ranFUN = ranfun, 
                               evFUN = "meanDistance",
                               genome = "hg38",
                               universe = uni,
                               mc.cores = 30,
                               count.once = T)
}

# Subset a few region sets to run the test on
setNames <- c("MAFF_ENCFF005YUC", "FOXA1_ENCFF011QFM", "RAD21_ENCFF155CEQ", "POLR2A_ENCFF159PYD", 
              "CTCF_ENCFF199YFA", "H3K9me3_ENCFF372HCL", "H3K27ac_ENCFF392KDI", "H3K4me3_ENCFF982DUT")
subSetX <- setX[names(setX) %in% setNames]

# Run tests with different randomization strategies
permResList_RaR <- lapply(subSetX, FUN = nzPerm, sets = setX, nt = 5000)

# Store results
saveRDS(permResList_RaR, file = "hpcResults/permResList_RaR_meanDist.RDS")
