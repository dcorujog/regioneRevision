library(regioneReloaded)
library(ggplot2)
library(egg)

# ENCODE PEAKS MATRIX ----

mpt <- readRDS("10_vs_all.rds")
getParameters(mpt)

selectRow <- c("RAD21", 
               "CTCF",
               "FOXA1_ENCFF081USG",
               "FOXA2",
               "MAFK", 
               "MAFF",
               "JUND_ENCFF562",
               "JUN_",
               "POLR2A_",
               "POLR2G_",
               "H3K4me3",
               "H3K27ac",
               "H3K9me3")

selectCol <- c("FOXA1_ENCFF081USG",
               "JUN_",
               "POLR2A",
               "CTCF",
               "RAD21",
               "MAFF")

mpt <- makeCrosswiseMatrix(mpt, 
                           symm_matrix = F,
                           selectRow = selectRow,
                           selectCol = selectCol, 
                           pvcut = 0.01, 
                           hc.method = "average")

pmat <- plotCrosswiseMatrix(mpt,
                    lineColor = "black") +
  labs(x = "",
       y = "") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

pdf("plots/fig1_association_matrix.pdf", height = unit(14, "cm"), widt = unit(6, "cm"))
print(pmat)
dev.off()

mat <- getMatrix(mpt)
mat

ord_row <- c(1:length(rownames(mat)))
names(ord_row) <- rownames(mat)
ord_row

ord_col <- c(1:length(colnames(mat)))
names(ord_col) <- colnames(mat)
ord_col

# ENCODE PEAKS DIMRED ----

pca <- plotCrosswiseDimRed(mpt, ellipse = F, nc = 6, labSize = 1.5) +
  coord_cartesian(xlim = c(-3, 3),
                  ylim = c(-3, 3)) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

set.seed(42)
tsne <- plotCrosswiseDimRed(mpt, ellipse = F, nc = 6, type = "tSNE", 
                            perplexity = 3, labSize = 1.5) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8))
tsne

pdf("plots/ENCODE_PCA.pdf")
grid.arrange(grobs = list(set_panel_size(pca, width = unit(4, "cm"), height = unit(4, "cm"))))
dev.off()

pdf("plots/ENCODE_tSNE.pdf")
grid.arrange(grobs = list(set_panel_size(tsne, width = unit(4, "cm"), height = unit(4, "cm"))))
dev.off()

# POLR2A LOCAL Z-SCORE

mlz <- readRDS("hpcResults/POLR2A_localZS.RDS")
mlz <- makeLZMatrix(mlz, normalize = TRUE, scale = TRUE)
  
pmlz <- plotLocalZScoreMatrix(mlz, maxVal = 1) +
  labs(x = "",
       y = "") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))

pdf("plots/ENCODE_local.pdf")
grid.arrange(grobs = list(set_panel_size(pmlz, width = unit(6.5, "cm"), height = unit(12, "cm"))))
dev.off()
