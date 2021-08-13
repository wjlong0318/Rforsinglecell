BiocManager::install("idpr")
library(idpr)
pep="CQQYVNSDYTF"
meanScaledHydropathy(pep)
scaledHydropathyGlobal(pep)
