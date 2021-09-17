

require(plyranges)
require(GenomicRanges)
require(ArchR)




callPeaks <- function( ){


}

candidatePeaks <- scMACS::determine_dynamic_range(sampleFrags, cellSubsetArchR,500,doBin=FALSE)
int_matrix_1 <- scMACS::calculate_intensities(fragMat = sampleFrags[[1]],candidatePeaks=candidatePeaks,
			NBDistribution = NBDistribution_global,normalizeBins=FALSE)

int_matrix_2 <- scMACS::calculate_intensities(fragMat = sampleFrags[[2]],candidatePeaks=candidatePeaks,
			NBDistribution = NBDistribution_global,normalizeBins=FALSE)

